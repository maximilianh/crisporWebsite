package ur_rna.StructureEditor.menus;

import ur_rna.StructureEditor.AppActions;
import ur_rna.Utilities.swing.MergeMenu;

public class FileMenu extends MainMenu {
    public static final int OPEN_FILE_POS = -100;
    public static final int SAVE_FILE_POS = -80;
    public static final int CLOSE_FILE_POS = -70;
    public static final int BOTTOM_POS = 10000;
    private final MergeMenu recent;
    public FileMenu() {
        super("&File");
        addItem(AppActions.NEW_FILE, OPEN_FILE_POS);
        addItem(AppActions.SHOW_OPEN_FILE, OPEN_FILE_POS);
        recent = new MergeMenu("Open &Recent", OPEN_FILE_POS);
        recent.addSeparator(1000);
        recent.add(AppActions.CLEANUP_RECENT_FILES, 1000);
        recent.add(AppActions.CLEAR_RECENT_FILES, 1000);
        add(recent);
        addSeparator(BOTTOM_POS);
        addItem(AppActions.EXIT_PROGRAM, BOTTOM_POS);
    }
    public MergeMenu recentFilesMenu() {
        return recent;
    }
}
