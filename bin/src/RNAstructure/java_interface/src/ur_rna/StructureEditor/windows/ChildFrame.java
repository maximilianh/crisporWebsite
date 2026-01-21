package ur_rna.StructureEditor.windows;

import ur_rna.Utilities.swing.MdiChildFrame;

import javax.swing.*;
import java.awt.*;
import java.beans.PropertyVetoException;
import java.util.Collection;

public class ChildFrame extends MdiChildFrame {
    public ChildFrame() { super(); }
    public ChildFrame(String title) { this(title, true); }
    public ChildFrame(String title, boolean resizeable) {
        super(title, resizeable);
    }
//    public MainFrame getParent() { return (MainFrame)super.getParentFrame(); }
//    public void setParent(MainFrame frame) {super.setParentFrame(frame); }
    public Collection<? extends JMenu> getMenus() { return null; }
    public Collection<? extends Component> getToolbarButtons() { return null; }

    public void close() {
        try {
            setClosed(true);
        } catch (PropertyVetoException ex) {
            //ex.printStackTrace();
        }
    }
    /** Called when Program Settings have changed. ChildFrame descendants can override this to perform custom updates as necessary. */
    public void programSettingsUpdated() { }
}
