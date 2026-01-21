package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.function.Consumer;

/**
 * @author Richard M. Watson
 */
public class MdiChildFrame extends JInternalFrame {
    private MdiParentFrame parentFrame;
    protected void setParentFrame(MdiParentFrame parent) {
        parentFrame = parent;
        parent.addChild(this);
    }
    protected MdiParentFrame getParentFrame() { return parentFrame; }

    public MdiChildFrame() {
        this("", true);
        this.setTitle(this.getClass().getSimpleName());
    }
    public MdiChildFrame(String title) { this(title, true); }
    public MdiChildFrame(String title, boolean resizeable) {
        super(title, resizeable, true, resizeable, resizeable);
    }

    /** For use only from within MdiParentFrame */
    void setParentFrame_Internal(final MdiParentFrame frame) {
        parentFrame = frame;
    }

    protected void addAction(final String actionName, final Consumer<ActionEvent> handler) {
        ActionHelper.addAction(this, actionName, handler);
    }
    protected void addKeyBinding(final char keyChar, final String actionName) {
        ActionHelper.addKeyBinding(this, keyChar, actionName);
    }
    protected void addKeyBinding(final int keyCode, final int modifiers, final String actionName) {
        ActionHelper.addKeyBinding(this, keyCode, modifiers, actionName);
    }
}
