package ur_rna.StructureEditor.menus;

import ur_rna.Utilities.swing.AcceleratorKey;
import ur_rna.Utilities.swing.MergeItem;
import ur_rna.Utilities.swing.UiAction;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.InputEvent;
import java.util.ArrayList;
import java.util.function.Predicate;


public class WindowsMenu extends MainMenu {
    UiAction nextWindow = new UiAction("Switch to Next Window", ()->gotoNext(1)).setKeyStroke("*TAB");
    UiAction prevWindow = new UiAction("Switch to Previous Window", ()->gotoNext(-1)).setKeyStroke("*+TAB");
    private final JDesktopPane desktop;
    private final JPopupMenu popup;
    private final Predicate<JInternalFrame> windowFilter;
    private final JMenuItem noWindowsItem;

    public static void enableTabNext(JComponent c) {
        AcceleratorKey.removeTabTraversalKeys(c, AWTKeyStroke.getAWTKeyStroke(9, InputEvent.CTRL_DOWN_MASK), AWTKeyStroke.getAWTKeyStroke(9, InputEvent.CTRL_DOWN_MASK|InputEvent.SHIFT_DOWN_MASK));
//        c.getActionMap().put(nextWindow.name(), nextWindow);
//        c.getActionMap().put(prevWindow.name(), prevWindow);
//        c.getInputMap(WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(nextWindow.keyStroke(), nextWindow.name());
//        c.getInputMap(WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(prevWindow.keyStroke(), prevWindow.name());
//        c.getInputMap(WHEN_FOCUSED).put(nextWindow.keyStroke(), nextWindow.name());
//        c.getInputMap(WHEN_FOCUSED).put(prevWindow.keyStroke(), prevWindow.name());
//        c.getInputMap(WHEN_IN_FOCUSED_WINDOW).put(nextWindow.keyStroke(), nextWindow.name());
//        c.getInputMap(WHEN_IN_FOCUSED_WINDOW).put(prevWindow.keyStroke(), prevWindow.name());
    }

    public WindowsMenu(JDesktopPane desktop, Predicate<JInternalFrame> windowFilter) {
        super("&Windows");
        this.desktop = desktop;
        this.windowFilter = windowFilter;
        popup = getPopupMenu();

        noWindowsItem = new MergeItem("(No windows are currently open.)");
        noWindowsItem.setEnabled(false);

        addSeparator(10);
        add(nextWindow, 10);
        add(prevWindow, 10);
//        JMenuItem mi = new JMenuItem("Banans");
//        mi.addActionListener(e-> {
////                    enableTabNext(Program.getInstance().getMainFrame().getRootPane());
////                    enableTabNext(Program.getInstance().getMainFrame().getLayeredPane());
////                    enableTabNext(Program.getInstance().getMainFrame().getToolBar());
////                    enableTabNext(Program.getInstance().getMainFrame().getJMenuBar());
////                    enableTabNext((JComponent) Program.getInstance().getMainFrame().getGlassPane());
////                    enableTabNext(Program.getInstance().getMainFrame().desktop);
//                    for (JInternalFrame f : desktop.getAllFrames())
//                        enableTabNext(f);
//                });
//        add(mi);
        popup.addPopupMenuListener(new PopupMenuListener() {
            @Override
            public void popupMenuWillBecomeVisible(final PopupMenuEvent e) {
                //System.out.println("popupMenuWillBecomeVisible");
                updateWindowList();
            }
            @Override
            public void popupMenuWillBecomeInvisible(final PopupMenuEvent e) {
                //System.out.println("popupMenuWillBecomeInvisible");
            }
            @Override
            public void popupMenuCanceled(final PopupMenuEvent e) {
                //System.out.println("popupMenuCanceled");
            }
        });
    }
    private void removeWindowItems() {
        for (int i = popup.getComponentCount()-1; i >=0; i--) {
            if (popup.getComponent(i) instanceof  WindowItem)
                popup.remove(i);
        }
    }
    private void updateWindowList() {
        removeWindowItems();
        int added = 0;
        for(JInternalFrame f : desktop.getAllFrames())
            if (windowFilter.test(f)) {
                popup.add(createWindowItem(f), added++);
            }
        if (added == 0)
            popup.add(noWindowsItem, 0);
        else
            popup.remove(noWindowsItem);

        nextWindow.setEnabled(added>1);
        prevWindow.setEnabled(added>1);
    }
    private static class WindowItem extends MergeItem {
        public final JInternalFrame window;
        public WindowItem(final JInternalFrame window) {
            this.window = window;
            setText(window.getTitle());
        }
    }
    public WindowItem createWindowItem(JInternalFrame window) {
        WindowItem w = new WindowItem(window);
        w.addActionListener(this::windowItem_click);
        return w;
    }
    private void windowItem_click(final ActionEvent event) {
        WindowItem w = (WindowItem)event.getSource();
        w.window.toFront();
    }

    void gotoNext(int direction) {
        JInternalFrame[] list = desktop.getAllFrames();
        ArrayList<JInternalFrame> frames = new ArrayList<>(list.length);
        for(JInternalFrame f : list)
            if (windowFilter.test(f))
                frames.add(f);
        if (frames.size()<2) return;
        int pos = frames.indexOf(desktop.getSelectedFrame());
        if (pos==-1)
            pos = 0;
        else
            pos = (pos + direction) % frames.size();
        if (pos < 0) pos += frames.size();
        final JInternalFrame f = frames.get(pos);
        f.toFront();
        desktop.setSelectedFrame(f);
        SwingUtilities.invokeLater(()-> {
            if (desktop.getSelectedFrame() == f)
                try {
                    f.revalidate();
                    desktop.revalidate();
                    f.repaint();
                    desktop.repaint();
                    f.updateUI();
                    desktop.updateUI();
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
        });
    }
}
