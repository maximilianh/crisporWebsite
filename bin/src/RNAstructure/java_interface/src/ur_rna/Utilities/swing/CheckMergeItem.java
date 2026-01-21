package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.event.ActionListener;

/**
 * A container for a merge-able menu item.
 */
public class CheckMergeItem extends JCheckBoxMenuItem implements IMergeItem {
    private int pos = MergePosInitial;

    /** allows storage of user-defined additional information */
    public Object tag;

    public JCheckBoxMenuItem getMenuItem() { return this; }

    public CheckMergeItem()  { }
    public CheckMergeItem(int mergePos)  { pos = mergePos; }
    public CheckMergeItem(Action action)  { this(action, MergePosDefault); }
    public CheckMergeItem(Action action, int mergePos)  { super(action); this.pos = mergePos; }
    public CheckMergeItem(Action action, String text, int mergePos)  { super(action); this.pos = mergePos; setText(text); }
    public CheckMergeItem(String text)  { setText(text); }
    public CheckMergeItem(String text, int mergePos)  { this(mergePos); setText(text); }
    public CheckMergeItem(String text, ActionListener action, String actionCommand)  { this(text, action, actionCommand, 0, null); }
    public CheckMergeItem(String text, ActionListener action, String actionCommand,  int mnemonic, KeyStroke key)  {
        super(new UiAction(text, action).setMnemonic(mnemonic).setKeyStroke(key).setCommand(actionCommand));
    }

    @Override
    public void setText(String text) {
        Menus.setTextWithMnemonic(this, text, super::setText);
    }

    @Override
    public int getMergePos() { return pos; }
    @Override
    public CheckMergeItem setMergePos(final int value) { pos = value; return this; }

    public CheckMergeItem setKeyStroke(final char keyStroke) {
        setKeyStroke(AcceleratorKey.getKey(keyStroke));
        return  this;
    }

    public CheckMergeItem setKeyStroke(final String keyStroke) {
        setKeyStroke(AcceleratorKey.getKey(keyStroke));
        return  this;
    }

    public CheckMergeItem setKeyStroke(final KeyStroke keyStroke) { this.getMenuItem().setAccelerator(keyStroke);
        return  this;
    }

    public KeyStroke getKeyStroke() {
        return this.getMenuItem().getAccelerator();
    }

    @Override
    public String toString() {
        String type = getClass().getSimpleName();
        return String.format("%s [ Name: \"%s\"  Text: \"%s\" Pos: %d ]", type, getName(), getText(), this.getMergePos());
    }

}
