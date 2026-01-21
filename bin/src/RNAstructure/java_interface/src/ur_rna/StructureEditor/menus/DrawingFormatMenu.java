package ur_rna.StructureEditor.menus;

import ur_rna.StructureEditor.windows.DrawWindow;
import ur_rna.Utilities.swing.MergeMenu;

/**
 * @author Richard M. Watson
 */
public class DrawingFormatMenu extends MergeMenu {
    private final DrawWindow owner;
    public DrawingFormatMenu(final DrawWindow frame) {
        super("For&mat");
        this.owner = frame;
    }
}
