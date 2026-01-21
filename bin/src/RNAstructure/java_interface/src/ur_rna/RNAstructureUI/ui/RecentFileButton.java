package ur_rna.RNAstructureUI.ui;

import ur_rna.RNAstructureUI.RNAstructure;
import ur_rna.RNAstructureUI.utilities.FileSelectedEvent;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.function.Consumer;

import static ur_rna.Utilities.ObjTools.indexOf;

/**
 * @author Richard M. Watson
 * Creates a button that allows users to select recently-used files.
 */
public class RecentFileButton extends JButton {
    public String filters;
    public Consumer<File> onFileSelected;

    public RecentFileButton() { }
    public RecentFileButton(String text, String filters, ActionListener fileSelectedListener, String fileSelectedCommand) {
        this(text, filters, null);
        onFileSelected = file->fileSelectedListener.actionPerformed(new FileSelectedEvent(this, file, fileSelectedCommand));
    }

    public RecentFileButton(String text, String filters, Consumer<File> fileSelectedListener) {
        super(text);
        this.filters = filters;
        onFileSelected = fileSelectedListener;
        addActionListener(e -> showFileSelectMenu());
    }

    // Called on click
    public void showFileSelectMenu() {
        showRecentFileMenu(this, this.filters, this.onFileSelected);
    }

    public static void showRecentFileMenu(final JComponent displayTarget, final String filters, final Consumer<File> onFileSelected) {
        FileFilter filter = parseFilters(filters);
        List<File> examples = filterFiles(getExamples(), filter);
        List<File> recent = filterFiles(getRecent(), filter);
        JPopupMenu menu = new JPopupMenu("Select a File");
        Font headerFont = menu.getFont().deriveFont(Font.BOLD);
        menu.add(createMenuItem("Recent Files:", null, headerFont));
        menu.addSeparator();
        if (recent == null || recent.size() == 0)
            menu.add(createMenuItem("  (No recent files yet.)", null));
        else {
            JPopupMenu recentItems = menu; // start by adding directly to the main drop-down menu (later switch to an overflow sub-menu)
            // recent files are listed in reverse order.
            final int MAX_ITEMS = 8;
            int added = 0;
            for (int i = recent.size() -1; i >=0 ; i--) {
                if (added++ == MAX_ITEMS) {
                    // if we exceed the maximum number of recent items, create an overflow menu.
                    JMenu extra = new JMenu("More Recent Files");
                    menu.add(extra);
                    recentItems = extra.getPopupMenu();
                }
                recentItems.add(createFileItem(recent.get(i), onFileSelected, true));
            }
        }
        menu.addSeparator();
        menu.add(createMenuItem("Example Files:", null, headerFont));
        menu.addSeparator();
        if (examples == null)
            menu.add(createMenuItem("  (Example files Missing! Click to create.)", createExampleFilesClicked));
        else if (examples.size() == 0)
            menu.add(createMenuItem("  (There are no examples of this file type.)", null));
        else {
            for(File f : examples)
                menu.add(createFileItem(f, onFileSelected, false));
        }
        menu.show(displayTarget, 0, displayTarget==null?0:displayTarget.getHeight());
    }

    public static JMenuItem createMenuItem(final String text, final ActionListener listener) { return createMenuItem(text, listener, null); }
    public static JMenuItem createMenuItem(final String text, final ActionListener listener, Font font) {
        JMenuItem mi = new JMenuItem(text);
        if (listener==null)
            mi.setEnabled(false);
        else
            mi.addActionListener(listener);
        if (font!=null) mi.setFont(font);
        return mi;
    }
    public static JMenuItem createFileItem(final File f, final Consumer<File> onFileSelected, boolean showDirectory) {
        String text = "  " + f.getName();
        if (showDirectory) text += " in " + (RNAstructure.isExampleFile(f) ? "<Examples>" : f.getParent());
        return createMenuItem(text, e->onFileSelected.accept(f));
    }

    public static class ExtensionFileFilter implements FileFilter {
        private String[] extensions, normalizedExt;
        private boolean acceptAll, mustExist, includeDirs;

        public ExtensionFileFilter(final boolean mustExist, final String... extensions) {
            this.acceptAll = indexOf(extensions, "*") != -1;
            this.mustExist = mustExist;
            setExtensions(extensions);
        }
        public ExtensionFileFilter(final String... extensions) { this(false, extensions); }

        public static ExtensionFileFilter AcceptAll(boolean mustExist) {
            ExtensionFileFilter f = new ExtensionFileFilter();
            f.setAcceptAll(true);
            f.setMustExist(mustExist);
            return f;
        }


        /**
         * Tests whether or not the specified pathname should be
         * included in a pathname list.
         * @param pathname The pathname to be tested
         * @return <code>true</code> if and only if <code>pathname</code>
         * should be included
         */
        @Override
        public boolean accept(final File pathname) {
            if (pathname == null) return false;
            if (mustExist && !pathname.exists()) return false;
            if (acceptAll) return true;
            if (pathname.isDirectory() && !includeDirs) return false;
            // NOTE: we tested implementations using Maps, binary search
            // on a sorted list and this implementation. All implementations
            // provided roughly the same speed, most likely because of
            // overhead associated with java.io.File. Therefor we've stuck
            // with the simple lightweight approach.
            String fileName = pathname.getName();
            int pos = fileName.lastIndexOf('.');
            int len = normalizedExt.length;
            if (pos > 0 && pos < fileName.length() - 1) {
                String ext = fileName.substring(pos + 1).toLowerCase(Locale.ENGLISH);
                for (int i = 0; i < len; i++)
                    if (ext.equals(normalizedExt[i]))
                        return true;
            }
            return false;
        }
        public void setExtensions(final String[] extensions) {
            this.extensions = new String[extensions.length];
            this.normalizedExt = new String[extensions.length];
            for (int i = 0; i < extensions.length; i++) {
                if (extensions[i] == null || extensions[i].length() == 0) {
                    throw new IllegalArgumentException("Each extension must be non-null and not empty");
                }
                this.extensions[i] = extensions[i];
                normalizedExt[i] = extensions[i].toLowerCase(Locale.ENGLISH);
            }
        }
         public ExtensionFileFilter() { }
         public String[] getExtensions() { return extensions; }
         public String[] getNormalizedExtensions() { return normalizedExt; }
         public boolean acceptAll() { return acceptAll; }
         public void setAcceptAll(final boolean acceptAll) { this.acceptAll = acceptAll; }
         public boolean mustExist() { return mustExist; }
         public void setMustExist(final boolean mustExist) { this.mustExist = mustExist; }
         public boolean includeFolders() { return includeDirs; }
         public void setIncludeFolders(final boolean include) { this.includeDirs = include; }
    }
    private static FileFilter filterAcceptAll = new FileFilter() {
        @Override public boolean accept(final File pathname) {
        return pathname != null && pathname.exists() && !pathname.isDirectory();
    }};

    private static FileFilter parseFilters(String filters) { //parses filters of the form "Description|ext;ext|Description2|ext1;ext2;"
        if (filters == null || filters.trim().length() == 0)
            return filterAcceptAll;
        ArrayList<String> exts = new ArrayList<>();
        String[] parts = filters.split("\\|");
        for (int i = 1; i < parts.length; i++) {
            String desc = parts[i - 1].trim();
            if (desc.length() == 0)
                continue; //ignore blanks between file-specs
            if (desc.equals("*"))
                return filterAcceptAll; // null means no filter--accept all files.
            for (String ext : parts[i].split(";"))
                exts.add(ext.trim());
            i++; // move to the next filter-spec (two elements are consumed for each filter-spec)
        }
        return new ExtensionFileFilter(true, exts.toArray(new String[exts.size()]));
    }

    public static List<File> getExamples() {
        File examples = RNAstructure.getExamplesDir();
        if (examples == null) return null;
        File[] list = examples.listFiles();
        return list == null ? null : Arrays.asList(list);
    }
    public static List<File> getRecent() {
        ArrayList<File> list = new ArrayList<>();
        for (String path : RNAstructure.MRUFiles.getList())
            list.add(new File(path));
        return list;
    }

    public static List<File> filterFiles(List<File> files, FileFilter filter) {
        if (files == null) return null;
        ArrayList<File> list = new ArrayList<>(files.size());
        for (File f : files)
            if (filter.accept(f))
                list.add(f);
        return list;
    }


    //ActionListener menuItemClicked = e -> raiseFileSelectedEvent(((FileMenuItem)e.getSource()).file);

    static ActionListener createExampleFilesClicked = e -> RNAstructure.ExtractExamples();

//    private static class FileMenuItem extends JMenuItem {
//        public final File file;
//        public FileMenuItem(File f) {            this(f, f.getName(), null);        }
//        public FileMenuItem(File f, String text, ActionListener clicked) {
//            super(text);
//            file = f;
//            if (clicked != null)
//                this.addActionListener(clicked);
//        }
//
//        public static String stripExt(String name) {
//            int pos = name.lastIndexOf(".");
//            if (pos > 0)
//                name = name.substring(0, pos);
//            return name;
//        }
//    }
//
//    // Called when a fl
//    public void raiseFileSelectedEvent(File f) {
//        if (f != null) {
//            ActionEvent ev = new FileSelectedEvent(this, f, command);
//            onFileSelected.actionPerformed(ev);
//        }
//    }

}
