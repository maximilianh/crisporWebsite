package ur_rna.Utilities.swing;

import ur_rna.Utilities.AppLog;
import ur_rna.Utilities.PathTools;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.basic.BasicFileChooserUI;
import java.awt.*;
import java.beans.PropertyChangeEvent;
import java.io.File;
import java.util.HashMap;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * A file-open and -save dialog.
 */
public class FileDialog extends JFileChooser {
    private static final String DEFAULT_CONTEXT = "(default)";
    public static final String PREFS = "recent-file-locations";
    private static FileFilter lastChosenFilter;
    private boolean fileShouldExist;
    private boolean isCurDirSet;
    private String context;
    private static HashMap<String,File> recentDirs;
    private boolean modifyDir = true;
    private static Preferences recentPathPreferences;
    private Component owner;

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
        return Preferences.userNodeForPackage(FileDialog.class).node(PREFS);
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

    public FileDialog() { this(CUSTOM_DIALOG, null, null); }
    public FileDialog(boolean isSaveDialog, String filters, String title) { this (isSaveDialog ? SAVE_DIALOG : OPEN_DIALOG, filters, title); }
    public FileDialog(int dialogType, String filters, String title) {
        super.addPropertyChangeListener(FILE_FILTER_CHANGED_PROPERTY, evt->filterChanged(evt) );
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
    private void filterChanged(final PropertyChangeEvent evt) {
        if (getDialogType() != SAVE_DIALOG) return;
        // TODO: getSelectedFile is always NULL here. Try to find a way to get the value of the TextBox.
        if (!(this.getUI() instanceof BasicFileChooserUI)) return;
        BasicFileChooserUI ui = (BasicFileChooserUI)this.getUI();
        String name = ui.getFileName();
        if (isEmpty(name)) return;
        PathTools.FileName parsed = PathTools.parse(name);
        //FileFilter f1 = (FileFilter)evt.getOldValue();
        FileFilter f2 = (FileFilter)evt.getNewValue();
        if (f2 instanceof FileNameExtensionFilter) {
            FileNameExtensionFilter fn = (FileNameExtensionFilter)f2;
            if (!fn.accept(new File(name))) {
                String[] exts = fn.getExtensions();
                if (exts.length!=0)
                    ui.setFileName(parsed.baseName() + '.' + exts[0]);
            }
        }
    }

    public static File getDefaultDir() {
        return getRecentDir(null, true, PathTools.getHomeDir());
    }
    public static void setDefaultDir(final File dir) {        setRecentDir(null, dir);    }
    public FileDialog setDefaultFile(String path) {
        return setDefaultFile(PathTools.fileFromPath(path, false));
    }
    public FileDialog setDefaultFile(File f) {
        if (f == null || f.isDirectory())
            setStartDir(f);
        else {
            super.setSelectedFile(f);
            File parent = f.getParentFile();
            if (parent != null && parent.exists())
                setStartDir(parent);
            setFilterFromFile(f);
        }
        return this;
    }
    public FileDialog setFilterIndex(int filterIndex) {
        setFileFilter(getChoosableFileFilters()[filterIndex]); return this;
    }
    public FileDialog setFilterFromFile(File f) {
        FileFilter ff = findMatchingFilter(f);
        if (ff != null && !ff.accept(new File(""))) setFileFilter(ff);
        return this;
    }

    private FileFilter findMatchingFilter(File f) {
        for (FileFilter ff : super.getChoosableFileFilters())
            if (ff.accept(f))
                return ff;
        return  null;
    }

    public FileDialog setFilterFromExt(String findExt) {
        if (findExt == null)
            setFileFilter(null);
        else {
            for (FileFilter ff : super.getChoosableFileFilters())
                if (ff instanceof FileNameExtensionFilter)
                    for (String ext : ((FileNameExtensionFilter) ff).getExtensions())
                        if (findExt.equalsIgnoreCase(ext)) {
                            setFileFilter(ff);
                            return this;
                        }
        }
        return this;
    }

    public FileDialog setStartDir(String path) { return setStartDir(path == null? null : new File(path)); }
    public FileDialog setStartDir(File dir) {
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
     * @return The current FileDialog for a fluent interface.
     */
    public FileDialog setModifyDir(final boolean modifyDir) {
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
            FileDialog.setRecentDir(context, getCurrentDirectory());
        return getSelectedFile();
    }

    /**
     * Set the context for default startup folder.
     * @param context
     * @return Returns this FileDialog for chained calls during configuration.
     */
    public FileDialog setContext(String context) {
        this.context = context;
        return this;
    }

    public String getContext() {
        return context;
    }

    /**
     * Set the text that appears on the accept (OK) button.
     * @param text The text to display on the button.
     * @return Returns this FileDialog for chained calls during configuration.
     */
    public FileDialog setButtonText(String text) {
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
     * @return Returns this FileDialog for chained calls during configuration.
     */
    public FileDialog addFilters(String filterSpecSet) {
        String[] parts = filterSpecSet.split("\\|");
        for (int i = 0; i < parts.length; i++) {
            String desc = parts[i].trim();
            if (desc.length() == 0)
                continue; //ignore blanks between file-specs
            if (desc.equals("*")) {
                super.setAcceptAllFileFilterUsed(true);
                // A '*' is allowed to be just the description OR both the desciption and the extension.
                // e.g.:   "Text files|txt|*"   or   "Text files|txt|*|*"
                if (parts.length > i+1 && parts[i+1].equals("*"))
                    i++;
                continue;
            }
            i++; // move to the next filter-spec (two elements are consumed for each filter-spec)
            if (parts.length <= i)
                throw new IllegalArgumentException("Wrong number of filter components. Each filter description must have a corresponding list of extensions."+
                        "\nExample: \"JPEG Image|jpg;jpeg|PNG Image|png|*\""+
                        "\nAttempted: " + filterSpecSet);
            String[] ext = parts[i].split(";");
            desc += String.format(" (*.%s)", String.join(", *.", ext));
            this.addChoosableFileFilter(new FileNameExtensionFilter(desc, ext));
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
     * @return Returns this FileDialog for chained calls during configuration.
     */
    public FileDialog clearFilters() {
        this.setAcceptAllFileFilterUsed(false);
        super.resetChoosableFileFilters();
        return this;
    }

    //    public static class Save extends FileDialog {
//        public Save() {super(true, null, null);}
//        public Save(String filters) { super(true, filters, null); }
//
//    }
//    public static class Open extends FileDialog {
//        public Open() { super(false, null, null);  }
//        public Open(String filters) { super(false, filters, null); }
//    }
    public static String getSaveName(String filters, String defaultPath, String title, String context, Component parent) {
        return getSaveName(filters, defaultPath, title, context, parent, null);
    }
    public static String getSaveName(String filters, String defaultPath, String title, String context, Component parent, String defaultExtension) {
        return getFileName(true, parent, title, defaultPath, filters, context, false, defaultExtension);
    }
    public static String getSaveName(String filters) { return getSaveName(filters, null); }
    public static String getSaveName(String filters, String defaultPath) { return getSaveName(filters, defaultPath, null, null, null); }

    public static String getOpenName(String filters, String defaultPath, String title, String context, Component parent) {
        return getFileName(false, parent, title, defaultPath, filters, context);
    }

    public static String getOpenName(String filters) { return getOpenName(filters, null); }
    public static String getOpenName(String filters, String defaultPath) { return getOpenName(filters, defaultPath, null, null, null); }

    public static String getFileName(boolean forSave, Component parent, String title, String defaultPath, String filters, String context) { return getFileName(forSave, parent, title, defaultPath, filters, context, false, null);    }
    public static String getFileName(boolean forSave, Component parent, String title, String defaultPath, String filters, String context, boolean disableDirChange, String defaultExtension) {
//        // See preemptNextFileNameRequest
//        if (preemptTimeout != 0 &&
//                (preemptFiltersRequired == null || preemptFiltersRequired.equals(filters))) {
//            // FYI: preemptTimeout is milliseconds after 1970-01-01, and because it is a `long`, it will not rollover for 292 million years
//            long time = System.currentTimeMillis() - preemptTimeout;
//            preemptTimeout = 0; //whether we return the preempted value or not, this preemption is canceled/completed
//            if (time <= 100)
//                return preemptFileName;
//        }
//
//        if (AppMainFrame.getUseSimpleFileChooser()) {
//            if (isEmpty(title))
//                title = String.format("RNAstructure - Get %s File Name", forSave ? "Save" : "Open");
//            return Dialogs.getInput(title, "Enter a file path.", defaultPath);
//        }
        FileDialog fc = new FileDialog(forSave, filters, title);

        if (parent != null)
            fc.setOwner(parent);
        if (context == null) context = title;
        if (context != null)
            fc.setContext(context);

        // Set default extension before defaultPath, because if defaultPath has a known extension, it will override defaultExtension.
        if (!isEmpty(defaultExtension))
            fc.setFilterFromExt(defaultExtension);

        if (defaultPath != null)
            fc.setDefaultFile(defaultPath);

        if (disableDirChange)
            fc.setModifyDir(false);

        File f = fc.showFileDialog(parent);
        if (f == null)
            return null;
        lastChosenFilter = fc.getFileFilter();
        return f.getAbsolutePath();
    }
    public static FileFilter getLastChosenFilter() {
        return lastChosenFilter;
    }

    //    /**
//     * Causes the next call to getFileName to return with the specified file WITHOUT showing the dialog.
//     * This is used by the "Examples" button to preempt the display of the dialog when a user selects an example file from the list.
//     * A preemption of this type will only be valid for 100 ms, which should be plenty of time for the subsequent call to getSaveName.
//     * @param file The file to return when getSaveName is called.
//     * @param filters The filter string that should be used to restrict preemption.
//     */
//    public static void preemptNextFileNameRequest(String file, String filters) {
//        preemptTimeout = System.currentTimeMillis();
//        preemptFiltersRequired = filters;
//        preemptFileName = file;
//    }
//    private static String preemptFileName, preemptFiltersRequired; private static long preemptTimeout;

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
            JOptionPane.showMessageDialog(null, "Please select a file.", null, JOptionPane.INFORMATION_MESSAGE);
            //Dialogs.showMessage( "Please select a file.");
            return;
        }
        if(!f.exists() && getFileShouldExist()){
            JOptionPane.showMessageDialog(null, "The specified file does not exist: \n" + f.getPath(), null, JOptionPane.INFORMATION_MESSAGE);
            //Dialogs.showMessage( "The specified file does not exist: \n" + f.getPath());
            return;
        } else if(f.exists() && getDialogType() == SAVE_DIALOG){
            int result = JOptionPane.showConfirmDialog(this,"The file exists, overwrite?\n\nDetails:\nFile: " + f.getPath(),"Existing file", JOptionPane.YES_NO_CANCEL_OPTION);
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
    public void setOwner(final Component owner) {
        this.owner = owner;
    }
    public Component getOwner() {
        return owner;
    }
}
