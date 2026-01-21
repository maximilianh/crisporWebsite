package ur_rna.StructureEditor.menus;

import ur_rna.StructureEditor.windows.DrawWindow;
import ur_rna.Utilities.swing.MergeMenu;

/**
 * @author Richard M. Watson
 */
public class DrawingEditMenu extends MergeMenu {
    private final DrawWindow owner;
    public DrawingEditMenu(final DrawWindow frame) {
        super("&Edit");
        this.owner = frame;
    }
}
