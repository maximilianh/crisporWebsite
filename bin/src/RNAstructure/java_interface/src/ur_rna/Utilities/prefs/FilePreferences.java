package ur_rna.Utilities.prefs;

import ur_rna.Utilities.AppInfo;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;
import java.util.logging.Logger;
import java.util.prefs.AbstractPreferences;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.prefs.PreferencesFactory;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * Preferences implementation that stores to a user-defined file.
 */
public class FilePreferences extends AbstractPreferences  {
    private static final Logger log = Logger.getLogger(FilePreferences.class.getName());
    private static final String[] EMPTY_NAME_ARRAY = new String[0];

    private boolean modified;
    private Map<String, String> prefs;
    private Set<String> childNames;
    protected final RootNode root;

    protected static class RootNode extends FilePreferences {
        private File file;
        private Properties cachedProps;

        public File getFile() {
            return file == null ? Factory.getPreferencesFile() : file;
        }
        public void setFile(final File file) {
            this.file = file;
        }

        public RootNode(File file) {
            super(null, "");
            this.file = file;
            Thread t = new Thread() { @Override public void run() { autoFlush(); } };
            t.setPriority(Thread.MAX_PRIORITY);
            Runtime.getRuntime().addShutdownHook(t);
        }

//        public void notifyUpdated() {
//
//        }

        protected void autoFlush() {
            try {
                flush();
            } catch (BackingStoreException ex) {
                log.severe(String.format("Error while auto-flushing %s to file \"%s\" : %s",
                        FilePreferences.class.getSimpleName(),
                        getFile(),
                        ex));
            }
        }

        public void loadProps() throws BackingStoreException {
            cachedProps = new Properties();
            final File file = getFile();
            if (!file.exists()) return;
            try (FileInputStream fs = new FileInputStream(file)){
                cachedProps.load(fs);
            } catch (IOException e) {
                throw new BackingStoreException(e);
            }
        }
        public Properties getProps() {
            return cachedProps;
        }
        public void storeProps() throws BackingStoreException {
            final File file = getFile();
            try (FileOutputStream fs = new FileOutputStream(file)){
                cachedProps.store(fs, FilePreferences.class.getSimpleName());
            } catch (IOException e) {
                throw new BackingStoreException(e);
            }
        }

        public void removeProps(final String basePath, boolean removeChildren) {
            // Make a list of all direct children of this node to be removed
            final Enumeration<?> propNames = cachedProps.propertyNames();
            List<String> toRemove = new ArrayList<String>();
            int len = basePath.length();

            if (removeChildren && basePath.equals("/")) {
                cachedProps = new Properties();
                return;
            }

            while (propNames.hasMoreElements()) {
                String propKey = (String) propNames.nextElement();
                // remove the key only if
                //   (A) we are removing this node and all child nodes or
                //   (B) we are removing keys (NOT child nodes) and and this property represents a KEY (i.e. has no slash after the basePath)
                if (propKey.startsWith(basePath) && (removeChildren || propKey.indexOf('/', len) == -1))
                    toRemove.add(propKey);
            }

            // Remove them now that the enumeration is done
            for (String propKey : toRemove)
                cachedProps.remove(propKey);
        }

        /**
         * Returns the absolute path and ensures that it ends with a slash (/)
         * For child nodes this is the same as {@code absolutePath() + '/'} but for
         * the root node it is just '/'.
         */
        @Override
        public String nodePrefix() {
            return "/";
        }
    }

    protected FilePreferences(AbstractPreferences parent, String name) {
        super(parent, name);
        log.finest("Instantiating node " + name);
        prefs = new TreeMap<>(); // LinkedHashMap -- order doesn't (and can't) matter. (Because underlying Properties map is unordered)
        childNames = new TreeSet<>(); //LinkedHash
        //children = new TreeMap<String, FilePreferences>();
        root = parent == null ? (RootNode)this : ((FilePreferences)parent).rootNode();
    }

    public static FilePreferences createRoot()  {        return createRoot(null);    }
    public static FilePreferences createRoot(String preferencesFilePath) {
        FilePreferences.RootNode prefs = new RootNode(preferencesFilePath == null ? null : new File(preferencesFilePath));
        try {
            prefs.loadProps();
            prefs.readProps();
        } catch (BackingStoreException ex) {
            log.severe("Exception initializing preferences from file. " + prefs.getPreferencesFile() + " - " + ex.toString());
        }
        return prefs;
    }

    public File getPreferencesFile() {
        return rootNode().getFile();
    }

    public void setPreferencesFile(File file) {
        rootNode().setFile(file);
    }

    public FilePreferences getRoot() {
        return root;
    }

    protected RootNode rootNode() {
        return root;
    }

    public boolean isModified() { return modified || isRemoved(); }

    protected void putSpi(String key, String value) {
        modified = true;
        prefs.put(key, value);
    }

    protected String getSpi(String key) {
        return prefs.get(key);
    }

    protected void removeSpi(String key) {
        modified = true;
        prefs.remove(key);
    }

    protected void removeNodeSpi() throws BackingStoreException {
        ((FilePreferences)parent()).childNames.remove(name());
        rootNode().removeProps(nodePrefix(), true);
        rootNode().storeProps();
    }

    protected String[] keysSpi() throws BackingStoreException {
        return prefs.keySet().toArray(new String[prefs.keySet().size()]);
    }

    protected String[] childrenNamesSpi() throws BackingStoreException {
         return childNames.size() == 0 ? EMPTY_NAME_ARRAY : childNames.toArray(new String[childNames.size()]);
    }

    protected FilePreferences childSpi(String name) {
        FilePreferences child = new FilePreferences(this, name);
        if (childNames.contains(name))
            child.readProps();
        else {
            child.newNode = true;
            child.modified = true; // If the node is new it is, by definition modified.
        }
        return child;
    }

    @Override
    public void sync() throws BackingStoreException {
        rootNode().loadProps();
        super.sync();
    }

    protected void syncSpi() throws BackingStoreException {
        readProps();
    }
    protected void readProps() {
        Properties p = root.getProps();
        String path = nodePrefix();
        int len = path.length();
        final Enumeration<?> names = p.propertyNames();
        while (names.hasMoreElements()) {
            String propKey = (String) names.nextElement();
            if (propKey.startsWith(path)) {
                String subKey = propKey.substring(len);
                int pos = subKey.indexOf('/');
                if (pos == -1)
                    prefs.put(subKey, p.getProperty(propKey));
                else if (pos == 0)
                    throw new IndexOutOfBoundsException("Property key name must not be empty: " + propKey);
                else
                    childNames.add(subKey.substring(0, pos));
            }
        }
        modified = false;
    }

    /**
     * Returns the absolute path and ensures that it ends with a slash (/)
     * For child nodes this is the same as {@code absolutePath() + '/'} but for
     * the root node it is just '/'.
     */
    public String nodePrefix() {
        return absolutePath() + '/';
    }

    @Override
    public void flush() throws BackingStoreException {
        rootNode().removeProps(nodePrefix(), false);
        super.flush();
        rootNode().storeProps();
    }

    protected void flushSpi() throws BackingStoreException { //throw new InternalError("Abstract member not used."); } // should never be called.
    //protected void flush(Properties p) throws BackingStoreException {
        Properties p = root.getProps();
        String path = nodePrefix();
        for (String s : prefs.keySet())
            p.setProperty(path + s, prefs.get(s));
    }

    /**
     * PreferencesFactory implementation that stores the preferences in an XML file.
     * (N.B. The system property <tt>java.util.prefs.PreferencesFactory</tt> must be set to
     * <tt>ur_rna.Utilities.prefs.FilePreferences.Factory</tt>
     * <p>
     * The file defaults to <tt>[user.home]/[ProgramTitle].prefs</tt>, but may be overridden with
     * the system property <tt>ur_rna.Utilities.prefs.FilePreferences.Factory</tt>
     *
     * @author Richard M. Watson
     */
    public static class Factory implements PreferencesFactory {
        private static final Logger log = Logger.getLogger(Factory.class.getName());

        FilePreferences rootPreferences;

        /**
         * Sets the system property {@code java.util.prefs.PreferencesFactory} to the {@link Factory}
         * class so that the {@link Preferences} class will use it for the global singleton PreferencesFactory.
         * This affects calls such as {@link Preferences#userRoot()} and {@link Preferences#userNodeForPackage(Class)}
         * etc.
         */
        public static void register() {
            System.setProperty("java.util.prefs.PreferencesFactory", Factory.class.getName());
        }
        public static void register(File preferencesFile) {
            register();
            setPreferencesFile(preferencesFile);
        }

        public Preferences systemRoot() {
            return userRoot();
        }

        public Preferences userRoot() {
            if (rootPreferences == null) {
                log.finer("Instantiating root preferences");
                rootPreferences = FilePreferences.createRoot();
            }
            return rootPreferences;
        }

        private static File preferencesFile;
        public static void setPreferencesFile(String systemPropertyFile) {
            setPreferencesFile(new File(systemPropertyFile));
        }
        public static void setPreferencesFile(File systemPropertyFile) {
            preferencesFile = systemPropertyFile;
        }
        public static File getPreferencesFile() {
            if (preferencesFile == null) {
                //String prefsFile = System.getProperty(SYSTEM_PROPERTY_FILE);
                //if (prefsFile == null || prefsFile.length() == 0) {
                String title = System.getProperty("program.title");
                if (isEmpty(title)) {
                    Class c = AppInfo.getMainClass();
                    title = c == null ? "java" : c.getName();
                }
                String prefsFile = System.getProperty("user.home") + File.separator + title + ".prefs";
                //}
                preferencesFile = new File(prefsFile).getAbsoluteFile();
                log.finer("Preferences file is " + preferencesFile);
            }
            return preferencesFile;
        }
    }
}
