package ur_rna.RNAstructureUI.utilities;

import java.util.HashMap;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

/**
 * @author Richard M. Watson
 */
public class UserOptions {
    public static final String DRAW_ALWAYS = "Y";
    public static final String DRAW_NEVER = "N";
    public static final String DRAW_PROMPT = "P";
    public static HashMap<String, String> overrides = new HashMap<>();
    public boolean useNativeLAF = true; 
    public boolean useNativeMenus = true; 
    public String drawStructures;
    private Preferences prefs;

    public void load(final Preferences prefs) {
        this.prefs = prefs;
        drawStructures = get(prefs, "draw-structures", DRAW_ALWAYS);
    }

    public void save(final Preferences prefs) {
        put(prefs, "draw-structures", drawStructures);
        try {
            prefs.flush();
        } catch (BackingStoreException ex) {
            //
        }
    }

    private String get(Preferences prefs, String name, String defaultValue) {
        String val = overrides.get(name);
        if (val != null) return val;
        return prefs.get(name, defaultValue);
    }
    private boolean put(Preferences prefs, String name, String value) {
        if (overrides.containsKey(name)) return false; // do NOT set it if overridden.
        prefs.put(name, value);
        return true;
    }

    public void save() {
        save(this.prefs);
    }
    public static void override(final String name, final String value) {
        overrides.put(name, value);
    }
}
