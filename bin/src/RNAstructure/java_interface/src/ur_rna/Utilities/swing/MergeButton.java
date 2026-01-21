package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;

/**
 *
 */
public class MergeButton extends JButton implements IMergeItem {
    int mergePos = MergePosInitial;

    /**
     * Creates a button with no set text or icon.
     */
    public MergeButton() {
    }
    /**
     * Creates a button with text.
     *
     * @param text the text of the button
     */
    public MergeButton(final String text) {
        super(text);
    }
    /**
     * Creates a button where properties are taken from the
     * <code>Action</code> supplied.
     *
     * @param a the <code>Action</code> used to specify the new button
     * @since 1.3
     */
    public MergeButton(final Action a) {
        super(a);
    }
    /**
     * Creates a button with initial text and an icon.
     *
     * @param text the text of the button
     * @param icon the Icon image to display on the button
     */
    public MergeButton(final String text, final Icon icon) {
        super(text, icon);
    }

    /**
     * Creates a button with text.
     *
     * @param text the text of the button
     */
    public MergeButton(final String text, final int mergePos) {
        super(text);
        this.mergePos = mergePos;
    }
    /**
     * Creates a button where properties are taken from the
     * <code>Action</code> supplied.
     *
     * @param a the <code>Action</code> used to specify the new button
     * @since 1.3
     */
    public MergeButton(final Action a, final int mergePos) {
        super(a);
        this.mergePos = mergePos;
    }
    @Override
    public int getMergePos() {
        return mergePos;
    }
    @Override
    public MergeButton setMergePos(final int value) {
        mergePos = value;
        return this;
    }
    @Override
    public Component getComponent() {
        return this;
    }

    public static class Separator extends JToolBar.Separator implements IMergeItem {
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
