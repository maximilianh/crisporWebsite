package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;

/**
 *
 */
public class MergeToggleButton extends JToggleButton implements IMergeItem {
    int mergePos = MergePosInitial;

    /**
     * Creates a button with no set text or icon.
     */
    public MergeToggleButton() {
    }
    /**
     * Creates a button with text.
     *
     * @param text the text of the button
     */
    public MergeToggleButton(final String text) {
        super(text);
    }
    /**
     * Creates a button where properties are taken from the
     * <code>Action</code> supplied.
     *
     * @param a the <code>Action</code> used to specify the new button
     * @since 1.3
     */
    public MergeToggleButton(final Action a) {
        super(a);
    }
    /**
     * Creates a button with initial text and an icon.
     *
     * @param text the text of the button
     * @param icon the Icon image to display on the button
     */
    public MergeToggleButton(final String text, final Icon icon) {
        super(text, icon);
    }

    /**
     * Creates a button with text.
     *
     * @param text the text of the button
     */
    public MergeToggleButton(final String text, final int mergePos) {
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
    public MergeToggleButton(final Action a, final int mergePos) {
        super(a);
        this.mergePos = mergePos;
    }
    @Override
    public int getMergePos() {
        return mergePos;
    }
    @Override
    public MergeToggleButton setMergePos(final int value) {
        mergePos = value;
        return this;
    }
    @Override
    public Component getComponent() {
        return this;
    }
}
