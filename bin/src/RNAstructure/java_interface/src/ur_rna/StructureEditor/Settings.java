package ur_rna.StructureEditor;

import ur_rna.StructureEditor.models.DrawSettings;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.Convert;
import ur_rna.Utilities.EventSource;
import ur_rna.Utilities.ObjTools;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import static ur_rna.Utilities.ObjTools.contains;

/**
 * Holds program options and settings.
 */
public class Settings {
    private static final String nullValue = "{~NULL}";
    public static final Object Missing = new Object();


    // These public fields will be automatically be loaded and saved.
    // -----------------------------------------------------------------
    public boolean ShowDashboard = true;
    public boolean DrawClockwise = false;
    public boolean ShowHints = true;
    public boolean SnapGuides = true;
    public boolean ExpandMotif = false;
    public boolean DisableRotation = false;
    public boolean DisableLoopResize = false;
    public boolean DisableBranchSlide = false;
    // -------------End of public fields--------------------------------

    private final EventSource.NoArg _changeEvent = new EventSource.NoArg();
    //private final static Field[] settingsFields = Settings.class.getFields();
    //private final static Field[] drawSettingsFields = DrawSettings.class.getFields();
    private final static Field[] fields;

//    private static class Setting {
//        public Field field;
//        public Object source;
//    }
    //private static Map<String, Setting> _settings = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);

    static {
        Field[] sf = Settings.class.getFields();
        Field[] df = DrawSettings.class.getFields();
        fields = new Field[sf.length + df.length];
        int count = 0;
        for(Field f : sf)
            if ((f.getModifiers() & Modifier.TRANSIENT) == 0)
                fields[count++] = f;
    }

    private DrawSettings drawSettings = new DrawSettings();


    public void addChangeHandler(Runnable handler) {
        _changeEvent.add(handler);
    }
    public void removeChangeHandler(Runnable handler) {
        _changeEvent.remove(handler);
    }
    public void notifyChanged() {
        _changeEvent.invoke();
    }
    public Object get(String name) {
        Field f = getSettingsField(name);
        if (f != null)
            return get(f);
        logError("Error retrieving program setting: Invalid setting name: \"%s\".", name);
        return Missing;
    }

    public Object get(Field f) {
        try {
            if (f.getDeclaringClass() == Settings.class)
                return f.get(this);
            else if (f.getDeclaringClass() == DrawSettings.class)
                return f.get(drawSettings);
            else
                throw new IllegalAccessException("Unknown Declaring Class");
        } catch (Exception ex) {
            logError("Error retrieving program setting '%s': %s", f.getName(), ex.getMessage());
            return Missing;
        }
    }
    public void set(String name, Object value) {
        Field f = getSettingsField(name);
        if (f == null)
            logError("Error updating program setting: Invalid setting name: \"%s\".", name);
        else
            set(f, value);
    }

    private boolean suspendChanges;
    private boolean settingsWereChanged;
    private void set(Field f, Object value) {
        try {
            // Note: value could be null. toType converts null to default values (e.g. 0 for numbers, '\0' for char, null for String.
            value = Convert.toType(value, f.getType());
            if (f.getDeclaringClass() == Settings.class)
                f.set(this, value);
            else if (f.getDeclaringClass() == DrawSettings.class)
                f.set(drawSettings, value);
            else
                throw new IllegalAccessException("Unknown Declaring Class");
            if (suspendChanges)
                settingsWereChanged=true;
            else
                notifyChanged();
        } catch (Exception ex) {
            logError("Error loading setting '%s':  to value '%s': %s", f.getName(), value, ex.getMessage());
        }
    }
    private void logError(String format, final Object ... args) {
        AppLog.getDefault().error(format, args);
    }

    private Field getSettingsField(String name) {
        for (Field f : fields) {
            if (f == null)
                break;
            if (f.getName().equals(name))
                return f;
        }
        return null;
    }

    public void load(final Preferences node) {
        String[] savedFields = ObjTools.EMPTY_STRING_ARRAY;
        suspendChanges = true;
        try {
            savedFields = node.keys();
        } catch (BackingStoreException ex) {
            ex.printStackTrace();
        }
        for (Field f : fields) {
            if (f == null)
                break;
            if (contains(savedFields, f.getName())) {
                String value = node.get(f.getName(), null);
                if (value == null)
                    continue;
                set(f, nullValue.equals(value) ? null : value);
            }
        }
        suspendChanges = false;
    }

    public void save(final Preferences node) {
        try {
            node.clear();
        } catch (BackingStoreException ex) {
            ex.printStackTrace();
        }
        for (Field f : fields) {
            if (f == null) break;
            try {
                Object value = get(f);
                if (value != Missing)
                    node.put(f.getName(), Convert.toString(value, nullValue));
            } catch (Exception ex) {
                System.out.println(String.format("Error saving setting '%s': %s", f.getName(), ex.getMessage()));
            }
        }
        try {
            node.flush();
        } catch (BackingStoreException ex) {
            ex.printStackTrace();
        }
    }
    public void suspendChangeNotification() {
        suspendChanges = true;
    }
    public void resumeChangeNotification() {
        suspendChanges = false;
        if (settingsWereChanged) {
            settingsWereChanged = false;
            notifyChanged();
        }
    }
    public DrawSettings getDrawSettings() {
        return drawSettings.copy();
    }
}
