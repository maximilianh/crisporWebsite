package ur_rna.StructureEditor.menus;

import ur_rna.StructureEditor.windows.DrawWindow;
import ur_rna.Utilities.swing.MergeMenu;

/**
 * @author Richard M. Watson
 */
public class DrawingViewMenu extends MergeMenu {
    private final DrawWindow owner;
    public DrawingViewMenu(final DrawWindow frame) {
        super("&View");
        this.owner = frame;
    }
}
