package ur_rna.StructureEditor.services;

import ur_rna.StructureEditor.AppActions;
import ur_rna.Utilities.swing.MergeItem;
import ur_rna.Utilities.swing.MergeMenu;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.util.HashMap;

/**
 * Handles recent file menu.
 */
public class RecentFileMenuManager {
    private final MergeMenu parent;
    public Container getComponent() { return parent; }
    public MergeMenu getMenu() { return parent; }

    public HashMap<MergeItem, RecentFileList.RecentFile> menuMap = new HashMap<>();
    public void setRecent(RecentFileList list) {
        for(MergeItem m : menuMap.keySet())
            parent.remove(m);
        int i = 100;
        for(RecentFileList.RecentFile f : list.getFiles()) {
            MergeItem item = new MergeItem(getMenuText(f), i++);
            menuMap.put(item, f);
            item.addActionListener(this::menuClickListener);
            item.setVisible(true);
            parent.add(item);
            //item.setToolTipText("Open " + toFriendlyName(p.item2) + ": " + p.item1);
        }
    }
    private void menuClickListener(final ActionEvent event) {
        @SuppressWarnings("SuspiciousMethodCalls") // menuMap.get will return null if getSource is null or not an MergeItem. This is OK, because the return value is checked.
        RecentFileList.RecentFile f =  menuMap.get(event.getSource());
        if (f != null)
            AppActions.openFile(f.path, f.type);
    }
    private String getMenuText(final RecentFileList.RecentFile f) {
        File file = new File(f.path);
        return file.getName() + "  in  " + file.getParent() + (f.type==null?"":"   [" + f.type.name() + "]");
    }

    public RecentFileMenuManager(final MergeMenu parent) {
        if (parent == null)
            throw new IllegalArgumentException("RecentFileMenuManager parent cannot be null.");
        this.parent = parent;
    }
}
