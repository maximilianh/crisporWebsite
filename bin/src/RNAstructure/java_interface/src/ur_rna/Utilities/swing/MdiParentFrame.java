package ur_rna.Utilities.swing;

import ur_rna.Utilities.ObjTools;

import javax.swing.*;
import javax.swing.event.InternalFrameEvent;
import javax.swing.event.InternalFrameListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * @author Richard M. Watson
 */
public abstract class MdiParentFrame extends JFrame implements InternalFrameListener {
    protected abstract JDesktopPane getDesktop();

    protected MdiParentFrame() { }
    protected MdiParentFrame(String title) { super(title);  }
    protected MdiParentFrame(String title, GraphicsConfiguration gc) { super(title, gc);  }

    public List<JInternalFrame> getChildFrames(Class ofType) {
        JInternalFrame[] all = getDesktop().getAllFrames();
        List<JInternalFrame> list = new ArrayList<>(all.length);
        for (JInternalFrame f : all)
            if (ofType.isAssignableFrom(f.getClass()))
                list.add(f);
        return list;
    }

    protected List<MdiChildFrame> getMdiChildFrames(Class ofType) {
        JInternalFrame[] all = getDesktop().getAllFrames();
        List<MdiChildFrame> list = new ArrayList<>(all.length);
        for (JInternalFrame f : all)
        if (f instanceof MdiChildFrame)
            list.add((MdiChildFrame)f);
        return list;
    }

    protected MdiChildFrame getActiveChild() {
        JInternalFrame active = getDesktop().getSelectedFrame();
        if (active instanceof MdiChildFrame)
            return (MdiChildFrame)active;
        return null;
    }
    protected void activateChild(MdiChildFrame child) {
        getDesktop().setSelectedFrame(child);
    }

    protected void addChild(MdiChildFrame child) {
        // add the listener (if not already added)
        if (!ObjTools.contains(child.getInternalFrameListeners(), this))
            child.addInternalFrameListener(this);

        // add the child window to the desktop (if this hasn't already happened)
        if (!ObjTools.contains(getDesktop().getAllFrames(), child))
            getDesktop().add(child);

        child.setParentFrame_Internal(this);
    }

    /**
     * Invoked when a internal frame has been opened.
     *
     * @param e
     * @see JInternalFrame#show
     */
    @Override
    public void internalFrameOpened(final InternalFrameEvent e) {

    }
    /**
     * Invoked when an internal frame is in the process of being closed.
     * The close operation can be overridden at this point.
     *
     * @param e
     * @see JInternalFrame#setDefaultCloseOperation
     */
    @Override
    public void internalFrameClosing(final InternalFrameEvent e) {

    }
    /**
     * Invoked when an internal frame has been closed.
     *
     * @param e
     * @see JInternalFrame#setClosed
     */
    @Override
    public void internalFrameClosed(final InternalFrameEvent e) {

    }
    /**
     * Invoked when an internal frame is iconified.
     *
     * @param e
     * @see JInternalFrame#setIcon
     */
    @Override
    public void internalFrameIconified(final InternalFrameEvent e) {

    }
    /**
     * Invoked when an internal frame is de-iconified.
     *
     * @param e
     * @see JInternalFrame#setIcon
     */
    @Override
    public void internalFrameDeiconified(final InternalFrameEvent e) {

    }
    /**
     * Invoked when an internal frame is activated.
     *
     * @param e
     * @see JInternalFrame#setSelected
     */
    @Override
    public void internalFrameActivated(final InternalFrameEvent e) {

    }
    /**
     * Invoked when an internal frame is de-activated.
     *
     * @param e
     * @see JInternalFrame#setSelected
     */
    @Override
    public void internalFrameDeactivated(final InternalFrameEvent e) {

    }
    public void close() {
        this.dispose();
    }

    protected void addAction(final String actionName, final Consumer<ActionEvent> handler) {
        ActionHelper.addAction(getDesktop(), actionName, handler);
    }
    protected void addKeyBinding(final char keyChar, final String actionName) {
        ActionHelper.addKeyBinding(getDesktop(), keyChar, actionName);
    }
    protected void addKeyBinding(final int keyCode, final int modifiers, final String actionName) {
        ActionHelper.addKeyBinding(getDesktop(), keyCode, modifiers, actionName);
    }
}
