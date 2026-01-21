package ur_rna.StructureEditor.services;

import ur_rna.Utilities.prefs.FilePreferences;

import java.io.File;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

/**
 * Assists with storing User Preferences to file.
 */
public class UserPrefs {
    public static void register() {
        FilePreferences.Factory.register();
    }

    public static final String DEFAULT_USER = "default";
    private final FilePreferences root;
    private Preferences appNode;
    private Preferences userNode;

    public UserPrefs(String appName) {
        this(FilePreferences.Factory.getPreferencesFile().toPath().getParent().resolve(appName + ".prefs").toFile());
    }
    public UserPrefs(File settingsFile) {
        root = FilePreferences.createRoot(settingsFile.toString());
        appNode = root;
        setUser(DEFAULT_USER); // sets userNode to default.
    }

    public void setUser(String name) { userNode = appNode.node(name); }
    public String getUser() { return userNode.name(); }

    /**
     * Store application-wide settings that are not user-specific here.
     */
    public Preferences app() { return appNode; }
    /**
     * Store user-specific settings here.
     */
    public Preferences user() { return userNode; }

    /**
     * Save user preferences to file.
     * This happens automatically upon proper Java shutdown, and is only
     * really necessary in scenarios where the JVM is terminated unexpectedly.
     */
    protected void save() {
        try {
            root.flush();
        } catch (BackingStoreException ex) {
            ex.printStackTrace();
        }
    }
    public String get(String name) { return get(name, null); }
    public String get(String name, String defaultValue) {
        return userNode.get(name, defaultValue);
    }
    public void put(String name, String value) {
        userNode.put(name, value);
    }
}
