package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;

/**
 * A container for a merge-able menu item.
 */
public class MergeItem extends JMenuItem implements IMergeItem {
//    public static final List<MergeItem> EMPTY_LIST = Collections.emptyList();
//    public static final JMenuItem SEPARATOR = new JMenuItem();

    private int pos = MergePosInitial;

    /** allows storage of user-defined additional information */
    public Object tag;

    public JMenuItem getMenuItem() { return this; }

    public MergeItem()  { }
    public MergeItem(int mergePos)  { setMergePos(mergePos); }
    public MergeItem(Action action)  { this(action, MergePosInitial); }
    public MergeItem(Action action, int mergePos)  { super(action); setMergePos(mergePos); }
    public MergeItem(Action action, String text, int mergePos)  { super(action); setMergePos(mergePos); setText(text); }
    public MergeItem(String text)  { setText(text); }
    public MergeItem(String text, int mergePos)  { this(mergePos); setText(text); }
    public MergeItem(String text, ActionListener action, String actionCommand,  int mnemonic, KeyStroke key)  {
        super(new UiAction(text, action).setMnemonic(mnemonic).setKeyStroke(key).setCommand(actionCommand));
    }

    @Override
    public void setText(String text) {
        Menus.setTextWithMnemonic(this, text, super::setText);
    }

//    @Override
//    @NotNull
//    public String getName() {
//        String s;
//        if (!isEmpty(s = super.getName())) return s;
//        if (null != (s = getText())) return s;
//        return "";
//    }

    @Override
    public int getMergePos() { return pos; }
    @Override
    public MergeItem setMergePos(final int value) {
        if (value != MergePosInvalid)
            pos = value;
        return this;
    }
    public MergeItem setKeyStroke(final char keyStroke) {
        setKeyStroke(AcceleratorKey.getKey(keyStroke));
        return  this;
    }
    public MergeItem setKeyStroke(final String keyStroke) {
        setKeyStroke(AcceleratorKey.getKey(keyStroke));
        return  this;
    }
    public MergeItem setKeyStroke(final KeyStroke keyStroke) { this.getMenuItem().setAccelerator(keyStroke);
        return  this;
    }
    public KeyStroke getKeyStroke() {
        return this.getMenuItem().getAccelerator();
    }

//    /**
//     * This function always returns an empty list.
//     * It is overridden in MergeMenu to return a true list of sub-items.
//     */
//    public List<MergeItem> getSubItems() { return EMPTY_LIST; }
//
//    /**
//     * Always throws {@link IndexOutOfBoundsException}
//     * This is overridden in MergeMenu to return the sub-item at the specified position.
//     */
//    public MergeItem getSubItem(int pos) { throw new IndexOutOfBoundsException(); }
//
//    /**
//     * This function always returns null.
//     * It is overridden in MergeMenu to search for a sub-item with the
//     * specified name and return it if found.
//     */
//    public MergeItem findByName(final String name) { return null; }
//    /**
//     * This function always returns null.
//     * It is overridden in MergeMenu to search for a sub-item with the
//     * specified name and return it if found.
//     */
//    public MergeItem findByTreePath(final String name) { return null; }

    // (this instanceof JSeparator) ? "(SEPARATOR)"
    @Override
    public String toString() {
        String type = getClass().getSimpleName();
        return String.format("%s [ Name: \"%s\"  Text: \"%s\" Pos: %d ]", type, getName(), getText(), this.getMergePos());
    }

    public static class Separator extends JPopupMenu.Separator implements IMergeItem {
        private int pos = MergePosInitial;
        public Separator() {}
        public Separator(int mergePos) { pos = mergePos; }
        @Override
        public int getMergePos() {
            return pos;
        }
        @Override
        public Separator setMergePos(final int value) {
            pos = value;
            return this;
        }
        @Override
        public Component getComponent() { return this; }
    }
}
