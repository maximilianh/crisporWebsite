package ur_rna.StructureEditor.menus;

import ur_rna.StructureEditor.windows.DrawWindow;
import ur_rna.Utilities.swing.MergeItem;
import ur_rna.Utilities.swing.MergeMenu;
import ur_rna.Utilities.swing.UiAction;

import javax.swing.*;

/**
 *
 */
public class DrawingFileMenu extends MergeMenu {
    private final DrawWindow owner;
    public DrawingFileMenu(final DrawWindow frame)  {
        super("&File");
        this.owner = frame;
        add(new MergeItem("Bananas are Middle"));
        add(new MergeItem("Bananas are TOP", -2000));
        add(new MergeItem(new UiAction("Bananas are Bottom", e->((AbstractButton)e.getSource()).setText("ZZZZZZZZZZ")), 20000));
    }
}
