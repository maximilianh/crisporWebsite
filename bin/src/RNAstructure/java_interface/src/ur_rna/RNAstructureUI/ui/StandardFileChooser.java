package ur_rna.RNAstructureUI.ui;

import ur_rna.RNAstructureUI.AppMainFrame;
import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.PathTools;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.io.File;
import java.util.HashMap;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * @author Richard M. Watson
 */
public class StandardFileChooser extends JFileChooser {
    private static final String DEFAULT_CONTEXT = "(default)";
    public static final String PREFS = "recent-file-locations";
    private boolean fileShouldExist;
    private boolean isCurDirSet;
    private String context;
    private static HashMap<String,File> recentDirs;
    private boolean modifyDir = true;
    private static Preferences recentPathPreferences;

    public  static void loadRecentPaths() { loadRecentPaths(null); }
    public  static void loadRecentPaths(Preferences storage) {
        try {
            if (storage == null) storage = getDefaultPreferenceStorage();
            Preferences p = recentPathPreferences = storage;
            if (recentDirs == null)
                recentDirs = new HashMap<>();
            else
                recentDirs.clear();

            for (String key : p.keys())
                setRecentDir(key, p.get(key, null));
            if (getRecentDir(null, true) == null) // if it isn't set OR if it doesn't exist, it will return null.
                setRecentDir(null, PathTools.getHomeDir());
        } catch (BackingStoreException ex) {
            AppLog.getDefault().error("Failed to load recent paths.", ex);
        }
    }

    public  static void saveRecentPaths() { saveRecentPaths(recentPathPreferences); }
    public static void saveRecentPaths(Preferences storage) {
        try {
            if (storage == null) storage = getDefaultPreferenceStorage();
            for (String key : recentDirs.keySet()) {
                File f = recentDirs.get(key);
                if (f != null) storage.put(key, f.toString());
            }
            storage.flush();
        } catch (Exception ex) {
            AppLog.getDefault().error("Failed to save recent paths.", ex);
        }
    }
    private static Preferences getDefaultPreferenceStorage() {
        return Preferences.userNodeForPackage(StandardFileChooser.class).node(PREFS);
    }

    public static File getRecentDir(String context) { return getRecentDir(context, false, null); }
    public static File getRecentDirOrDefault(String context, boolean mustExist) {
        return getRecentDir(context, mustExist, recentDirs.get(DEFAULT_CONTEXT));
    }
    public static File getRecentDir(String context, boolean mustExist) { return getRecentDir(context, mustExist, null); }
    public static File getRecentDir(String context, boolean mustExist, File defaultIfNotSet) {
        if (recentDirs == null) loadRecentPaths();
        if (isEmpty(context)) context = DEFAULT_CONTEXT;
        File f = recentDirs.get(context);
        if (f == null || mustExist && !f.exists())
            return defaultIfNotSet;
        return f;
    }

    public static void setRecentDir(String context, String value) {         setRecentDir(context, PathTools.fileFromPath(value));     }
    public static void setRecentDir(String context, File value) {
        if (recentDirs == null) loadRecentPaths();
        if (isEmpty(context)) context = DEFAULT_CONTEXT;
        if (value == null && !context.equals(DEFAULT_CONTEXT))
            recentDirs.remove(context);
        else
            recentDirs.put(context, value);
        saveRecentPaths();
    }

    public StandardFileChooser() { this(CUSTOM_DIALOG, null, null); }
    public StandardFileChooser(boolean isSaveDialog, String filters, String title) { this (isSaveDialog ? SAVE_DIALOG : OPEN_DIALOG, filters, title); }
    public StandardFileChooser(int dialogType, String filters, String title) {
        super.setDialogType(dialogType);
        super.setFileSelectionMode(JFileChooser.FILES_ONLY);
        if (filters == null || filters.trim().length() == 0)
            super.setAcceptAllFileFilterUsed(true);
        else {
            super.setAcceptAllFileFilterUsed(false);
            addFilters(filters);
        }
        if (title != null && title.length() != 0)
            super.setDialogTitle(title);
        setFileShouldExist(dialogType == OPEN_DIALOG);
    }

    public static File getDefaultDir() {
        return getRecentDir(null, true, PathTools.getHomeDir());
    }
    public static void setDefaultDir(final File dir) {        setRecentDir(null, dir);    }
    public StandardFileChooser setDefaultFile(String path) {
        return setDefaultFile(PathTools.fileFromPath(path, false));
    }
    public StandardFileChooser setDefaultFile(File f) {
          if (f == null || f.isDirectory())
            setStartDir(f);
        else {
            super.setSelectedFile(f);
            File parent = f.getParentFile();
            if (parent != null && parent.exists())
                    setStartDir(parent);
        }
        return this;
    }
    public StandardFileChooser setStartDir(String path) { return setStartDir(path == null? null : new File(path)); }
    public StandardFileChooser setStartDir(File dir) {
        isCurDirSet = true;
        super.setCurrentDirectory(dir);
        return this;
    }

    /**
     * Gets the property 'modifyDir', which affects whether the user's actions (i.e. changing the directory
     * of the file chooser) should affect the system by storing the new directory as the default for the
     * given context.
     * @return  Returns the current value of the property: True if the users actions should affect the recent
     * directory stored for the current context or False if the system's recent directories should not be affected.
     */
    public boolean modifyDir() {
        return modifyDir;
    }
    /**
     * Sets the property 'modifyDir', which affects whether the user's actions (i.e. changing the directory
     * of the file chooser) should affect the system by storing the new directory as the default for the
     * given context.
     * @param modifyDir Set to True if the users actions should affect the recent directory stored for the current
     *                  context. Set to false if the system's recent directories should not be affected.
     * @return The current StandardFileChooser for a fluent interface.
     */
    public StandardFileChooser setModifyDir(final boolean modifyDir) {
        this.modifyDir = modifyDir;
        return this;
    }
    /**
     * Show the File Chooser and return the selected file path, if successful, or null if the user
     * cancels.
     * @return The selected file path, if successful, or null if the user cancels.
     */
    public File showFileDialog(Component parent) {
        // Get the current file chooser home directory as a file.
        File cd = getRecentDir(context, true);
        if (cd == null) {
            if (getDialogType() == SAVE_DIALOG)
                cd = PathTools.getHomeDir();
            else
                cd = getDefaultDir();
        }
        if (!isCurDirSet && cd != null) setCurrentDirectory( cd );
        int success = super.showDialog(parent, this.getApproveButtonText());
        if( success != JFileChooser.APPROVE_OPTION ) { return null; }

        if (modifyDir)
            StandardFileChooser.setRecentDir(context, getCurrentDirectory());
        return getSelectedFile();
    }

    /**
     * Set the context for default startup folder.
     * @param context
     * @return Returns this StandardFileChooser for chained calls during configuration.
     */
    public StandardFileChooser setContext(String context) {
        this.context = context;
        return this;
    }

    public String getContext() {
        return context;
    }

    /**
     * Set the text that appears on the accept (OK) button.
     * @param text The text to display on the button.
     * @return Returns this StandardFileChooser for chained calls during configuration.
     */
    public StandardFileChooser setButtonText(String text) {
        this.setApproveButtonText(text);
        return this;
    }

    /**
     *  Adds a set of filter-specs to the list of file filters for the FileChooser.
     *  A "filter-spec" is a string containing a File Type Description followed by a bar ({@code |}) and one or more
     *  File Extensions, separated by semi-colons ({@code ;}). I.e.: {@code Description|ext1[;ext2...]}
     *  Both the description and at least one extension are required and cannot be blank.
     *  A "set" of filter-specs consists of one or more filter-specs separated by a bar. E.g.:
     *  {@code filter-spec_1[|filter-spec_2...]} i.e. {@code Description|ext1_1[;ext1_2...]|Description2|ext2_1[;ext2_2...] }
     *  A blank filter-spec is allowed and will simply be ignored. So the following will result in
     *  the addition of exactly two filters (despite the empty file-specs at the start, middle and end).
     *  ${ code |Description|ext1;ext2||Description2|ext1;ext2| }
     *  However the following filter-spec set would raise an error, because the list of extensions is effectively blank.
     *  (the ext1... would be interpreted as the start of the <em>second</em> filter-spec).
     *  ${ code Description||ext1;ext2 }
     *  As a special case, the filter-specs "*" or "*|*" are both interpreted as the accept-all-files filter, i.e.
     *  "All Files (*.*)"
     * @param filterSpecSet
     * @return Returns this StandardFileChooser for chained calls during configuration.
     */
    public StandardFileChooser addFilters(String filterSpecSet) {
        String[] parts = filterSpecSet.split("\\|");
        for (int i = 1; i < parts.length; i++) {
            String desc = parts[i-1].trim();
            if (desc.length() == 0)
                continue; //ignore blanks between file-specs
            if (desc.equals("*")) {
                super.setAcceptAllFileFilterUsed(true);
                if (parts[i].equals("*"))
                    i++;
                continue;
            }
            String[] ext = parts[i].split(";");

            desc += String.format(" Files (*.%s)", String.join(", *.", ext));
            this.addChoosableFileFilter(new FileNameExtensionFilter(desc, ext));
            i++; // move to the next filter-spec (two elements are consumed for each filter-spec)
        }
        return this;
    }

    /**
     * Clear any existing filters.
     * This is useful in any case where a file filter has been added, but is no longer desired.
     * For example when the default constructor is used, the all-files-filter is automatically added.
     * If that filter is not desired, then calling clearFilters will remove it along with any other existing filters.
     * If the only undesired filter is the all-files-filter, then simply calling
     * {@link #setAcceptAllFileFilterUsed(boolean)} with the value false would be more efficient
     * than clearing all filters.
     * @return Returns this StandardFileChooser for chained calls during configuration.
     */
    public StandardFileChooser clearFilters() {
        this.setAcceptAllFileFilterUsed(false);
        super.resetChoosableFileFilters();
        return this;
    }

//    public static class Save extends StandardFileChooser {
//        public Save() {super(true, null, null);}
//        public Save(String filters) { super(true, filters, null); }
//
//    }
//    public static class Open extends StandardFileChooser {
//        public Open() { super(false, null, null);  }
//        public Open(String filters) { super(false, filters, null); }
//    }
    public static String getSaveName(String filters, String defaultPath, String title, String context, Component parent) {
        return getFileName(true, parent, title, defaultPath, filters, context);
    }
    public static String getSaveName(String filters) { return getSaveName(filters, null); }
    public static String getSaveName(String filters, String defaultPath) { return getSaveName(filters, defaultPath, null, null, null); }

    public static String getOpenName(String filters, String defaultPath, String title, String context, Component parent) {
        return getFileName(false, parent, title, defaultPath, filters, context);
    }

    public static String getOpenName(String filters) { return getOpenName(filters, null); }
    public static String getOpenName(String filters, String defaultPath) { return getOpenName(filters, defaultPath, null, null, null); }

    public static String getFileName(boolean forSave, Component parent, String title, String defaultPath, String filters, String context) { return getFileName(forSave, parent, title, defaultPath, filters, context, false);    }
    public static String getFileName(boolean forSave, Component parent, String title, String defaultPath, String filters, String context, boolean disableDirChange ) {
        // See preemptNextFileNameRequest
        if (preemptTimeout != 0 &&
                (preemptFiltersRequired == null || preemptFiltersRequired.equals(filters))) {
            // FYI: preemptTimeout is milliseconds after 1970-01-01, and because it is a `long`, it will not rollover for 292 million years
            long time = System.currentTimeMillis() - preemptTimeout;
            preemptTimeout = 0; //whether we return the preempted value or not, this preemption is canceled/completed
            if (time <= 100)
                return preemptFileName;
        }

        if (AppMainFrame.getUseSimpleFileChooser()) {
            if (isEmpty(title))
                title = String.format("RNAstructure - Get %s File Name", forSave ? "Save" : "Open");
            return Dialogs.getInput(title, "Enter a file path.", defaultPath);
        }
        StandardFileChooser fc = new StandardFileChooser(forSave, filters, title);
        if (context == null) context = title;
        if (context != null)
            fc.setContext(context);
        if (defaultPath != null)
            fc.setDefaultFile(defaultPath);

        if (disableDirChange)
            fc.setModifyDir(false);

        File f = fc.showFileDialog(parent);
        if (f == null)
            return null;
        return f.getAbsolutePath();
    }
    /**
     * Causes the next call to getFileName to return with the specified file WITHOUT showing the dialog.
     * This is used by the "Examples" button to preempt the display of the dialog when a user selects an example file from the list.
     * A preemption of this type will only be valid for 100 ms, which should be plenty of time for the subsequent call to getSaveName.
     * @param file The file to return when getSaveName is called.
     * @param filters The filter string that should be used to restrict preemption.
     */
    public static void preemptNextFileNameRequest(String file, String filters) {
        preemptTimeout = System.currentTimeMillis();
        preemptFiltersRequired = filters;
        preemptFileName = file;
    }
    private static String preemptFileName, preemptFiltersRequired; private static long preemptTimeout;

    public boolean getFileShouldExist() { return fileShouldExist; }
    public void setFileShouldExist(final boolean fileShouldExist) {
        this.fileShouldExist = fileShouldExist;
    }

    public String getDefaultExtension() {
        FileFilter f = this.getFileFilter();
        if (f instanceof FileNameExtensionFilter) {
            FileNameExtensionFilter ff = (FileNameExtensionFilter)f;
            String[] exts = ff.getExtensions();
            if (exts.length != 0)
                return exts[0];
        }
        return "";
    }
    public boolean passesFilter(File f) {
        FileFilter ff = this.getFileFilter();
        return ff == null || ff.accept(f);
    }
    /**
     * Returns the selected file. This can be set either by the
     * programmer via <code>setSelectedFile</code> or by a user action, such as
     * either typing the filename into the UI or selecting the
     * file from a list in the UI.
     *
     * @return the selected file
     * @see #setSelectedFile
     */
    @Override
    public File getSelectedFile() {
        File f = super.getSelectedFile();
        if (f == null || f.exists()) return f;
        // Get the file, and make sure it ends with the proper extension.
        if (!passesFilter(f))
            f = new File(f.getPath() + "." + getDefaultExtension());
        return f;
    }
    @Override
    public void approveSelection(){
        File f = getSelectedFile();
        if (f == null) {
            Dialogs.showMessage( "Please select a file.");
            return;
        }
        if(!f.exists() && getFileShouldExist()){
            Dialogs.showMessage( "The specified file does not exist: \n" + f.getPath());
            return;
        } else if(f.exists() && getDialogType() == SAVE_DIALOG){
            int result = JOptionPane.showConfirmDialog(this,"The file exists, overwrite?\n\nDetails:\nFile: " + f.getPath(),"Existing file",JOptionPane.YES_NO_CANCEL_OPTION);
            switch(result){
                case JOptionPane.YES_OPTION:
                    super.approveSelection();
                    return;
                case JOptionPane.NO_OPTION:
                    return;
                case JOptionPane.CLOSED_OPTION:
                    return;
                case JOptionPane.CANCEL_OPTION:
                    cancelSelection();
                    return;
            }
        }
        super.approveSelection();
    }
}
