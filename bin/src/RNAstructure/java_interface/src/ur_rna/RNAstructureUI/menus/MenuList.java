package ur_rna.RNAstructureUI.menus;

import ur_rna.Utilities.swing.*;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *  A list of MainMenu items that should be incorporated into the application MenuBar when a Window receives focus.
 */
public class MenuList extends ArrayList<MergeMenu> implements IMenuItem {
    int subItemMergePos = MergeItem.MergePosInvalid;
    public void unsetSubItemMergePos() { subItemMergePos = MergeItem.MergePosInvalid; }
    public void setSubItemMergePos(int pos) {
        subItemMergePos = pos;
    }
    public void add(MergeMenu m, int mergePos) {
        add(m);
        IMergeItem.setMergePos(m, mergePos);
    }
    public void add(MergeMenu... menus) {
        for (MergeMenu m : menus)
            IMergeItem.setInitialMergePos(m, subItemMergePos);
        super.addAll(Arrays.asList(menus));
    }
    @Override
    public boolean add(final MergeMenu menu) {
        IMergeItem.setInitialMergePos(menu, subItemMergePos);
        return super.add(menu);
    }
    @Override
    public void add(final int index, final MergeMenu element) {
        super.add(index, element);
        IMergeItem.setInitialMergePos(element, subItemMergePos);
    }
    @Override
    public Component getComponent() {
        return null;
    }
    @Override
    public boolean hasSubItems() {
        return size()!=0;
    }
    @Override
    public List<MergeMenu> getMenus() {
        return this;
    }
    @Override
    public List<MergeMenu> getMenuItems() {
        return this;
    }
    public void enableMenus() {
        Menus.setEnabled(this, true);
    }
}
