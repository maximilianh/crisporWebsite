package ur_rna.StructureEditor.ui;

import javax.swing.*;
import javax.swing.plaf.basic.BasicInternalFrameUI;

/**
 *
 */
public class ToolWindow extends JInternalFrame {
    public ToolWindow() {
        super("", false, true, false, false);
        putClientProperty("JInternalFrame.isPalette", Boolean.FALSE);
        BasicInternalFrameUI bi = (BasicInternalFrameUI)getUI();
        bi.setNorthPane(null);
    }
}
