package ur_rna.Utilities.swing;

import javax.swing.*;
import java.util.Collections;
import java.util.List;

/**
 * Extends {@link javax.swing.ButtonGroup} to allow selecting items by index.
 *
 */
public class ButtonGroupList extends ButtonGroup {
    public boolean isSelected(int index) {
        return isSelected(buttons.get(index).getModel());
    }
    public void select(int index) {
        setSelected(index, true);
    }
    public void setSelected(int index, boolean selected) {
        setSelected(buttons.get(index).getModel(), selected);
    }
    public List<AbstractButton> getButtons() {
        return Collections.unmodifiableList(buttons);
    }
    public void removeAll() {
        AbstractButton[] arr = buttons.toArray(new AbstractButton[buttons.size()]);
        for (int i = arr.length - 1; i > -1; i--)
            remove(arr[i]);
    }
}
