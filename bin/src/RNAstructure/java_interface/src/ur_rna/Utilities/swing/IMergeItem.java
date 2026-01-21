package ur_rna.Utilities.swing;

import ur_rna.Utilities.annotation.NotNull;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * Provides an interface for a UI component that be included in a hierarchical list. In addition,
 * items can be ordered, grouped and merged (i.e. items from two parents with the same name are combined).
*/
public interface IMergeItem extends Cloneable {
    /** The default merge position for NON-merging items or merge items that have never been positioned. */
    int MergePosDefault = 0;
    /** The minimum allowed merge position. Values lower than this are reserved. */
    int MergePosMin = Integer.MIN_VALUE+16;
    /** The maximum allowed merge position. Values above than this are reserved. */
    int MergePosMax = Integer.MAX_VALUE;
    /** Indicates an invalid merge position.
     *  If a MergeMenu's subItemIndex is set to this value, it will NOT modify the mergePositions of
     *  child items added to it.  */
    int MergePosInvalid = Integer.MIN_VALUE+1;

    /** The merge position for newly created items. This distinguishes items that have never been positioned from those set to e.g. 0. */
    int MergePosInitial = Integer.MIN_VALUE+1;

//    /** Indicates that a child item should keep its current mergePosition when being added to a parent */
//    int MergePosPreserve = Integer.MIN_VALUE+1;
//    /** Indicates that a child item should obtain its parent's current subItemIndex when being added to a parent. */
//    int MergePosSubItem = Integer.MIN_VALUE+2;
//    /** Indicates that a child item have its merge position set to the default value when being added to a parent. */
//    int MergePosUseDefault = Integer.MIN_VALUE+3;

    int getMergePos();
    // Returns fluent reference to self
    IMergeItem setMergePos(int value);
    Component getComponent();
    default String getMergeName() { return Components.getNameOrText(getComponent()); }

//    String getMergeName();
//    default String getMergeGroup( ) { return null; }
//    default int getMergeOrientation( ) { return null; }
    default void mergeStart() {}
    default void mergeEnd() {}
    //default int compareTo(IMergeItem other) { return Integer.compare(getMergePos(), other.getMergePos()); }

    @NotNull
    static String getMergeName(Component c) {
        String s;
        if (c instanceof IMergeItem && !isEmpty(s = ((IMergeItem) c).getMergeName()))
            return s;
        return Components.getNameOrText(c);
    }

    static int getValidMergePos(Component c) {
        return fixInvalid(getMergePos(c));
    }
    static int getMergePos(Component c) {
        return c instanceof IMergeItem ? ((IMergeItem) c).getMergePos() : MergePosDefault;
    }
    static int compareByMergePos(Component c1, Component c2) {   return getValidMergePos(c1) - getValidMergePos(c2);  }
    static int compareByMergePos(IMergeItem c1, IMergeItem c2) { return fixInvalid(c1.getMergePos()) - fixInvalid(c2.getMergePos());    }
    static int fixInvalid(int pos) { return pos==MergePosInitial?MergePosDefault:pos; }

    static int setMergePos(IMergeItem mi, int newPos) {
        if (newPos!=MergePosInvalid)
            mi.setMergePos(newPos);
//        switch (newPos) {
//            case MergePosPreserve: break;
//            case MergePosSubItem: c.setMergePos(parentSubItemPos); break;
//            case MergePosUseDefault: c.setMergePos(MergePosDefault); break;
//            default: c.setMergePos(newPos);
//        }
        return mi.getMergePos();
    }
    static int setInitialMergePos(IMergeItem mi, int newPos) {
        if (mi.getMergePos()==MergePosInitial)
            mi.setMergePos(newPos);
        return mi.getMergePos();
    }

    class Proxy implements IMergeItem {
        private final Component component;
        private int mergePos;
        public Proxy(Component c) { this(c, MergePosInitial); }
        public Proxy(Component c, int mergePos) { this.component = c; this.mergePos = mergePos; }
        @Override
        public int getMergePos() {
            return mergePos;
        }
        @Override
        public Proxy setMergePos(final int value) { mergePos = value; return this; }
        @Override
        public Component getComponent() {
            return component;
        }
//        @Override
//        public String getMergeName() {
//            return getMergeName(component);
//        }
    }

    static IMergeItem asMergeItem(Component c) {
        if (c instanceof IMergeItem) return (IMergeItem)c;
        return new IMergeItem() {
            @Override
            public int getMergePos() {
                return MergePosDefault;
            }
            @Override
            public IMergeItem setMergePos(final int value) {
                throw new UnsupportedOperationException("The merge position of a proxy IMergeItem cannot be changed.");
            }
            @Override
            public Component getComponent() {
                return c;
            }
//            @Override
//            public String getMergeName() {
//                return getMergeName(c);
//            }
        };
    }

    static List<IMergeItem> asMergeItems(Collection<? extends Component> items) {
        List<IMergeItem> list = new ArrayList<>(items.size());
        for(Component c : items)
            list.add(asMergeItem(c));
        return list;
    }
}
