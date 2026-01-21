package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;
import java.util.Collections;
import java.util.List;

import static ur_rna.Utilities.Strings.isEmpty;

/**
 * Utility to help with Swing Components
 */
public abstract class Components {
    public static Component[] EMPTY_COMPONENT_ARRAY = new Component[0];
    public static List<Component> EMPTY_COMPONENT_LIST = Collections.emptyList();

    public static void addAll(Container parent, Component[] list) { for(Component c : list) parent.add(c);  }
    public static void addAll(Container parent, Iterable<? extends Component> list) { for(Component c : list) parent.add(c);  }
    /**
     * Returns true if the component c is equal to possibleAncestor or is a child or descendant of possibleAncestor.
     */
    public static boolean isAncestor(Component c, Container possibleAncestor) {
        if (c == possibleAncestor) return true;
        Container parent = c.getParent();
        while(parent!=null) {
            if (parent==possibleAncestor) return true;
            parent = parent.getParent();
        }
        return false;
    }

    /**
     * Return the first non-null non-empty String of the following: the item's name, text, or actionCommand.
     * @param c
     * @return
     */
    public static String getNameOrText(Component c) {
        String s;
        if (!isEmpty(s = c.getName())) return s;
        if (c instanceof AbstractButton) {
            if (!isEmpty(s = ((AbstractButton)c).getText())) return s;
            if (!isEmpty(s = ((AbstractButton)c).getActionCommand())) return s;
        }
        return null;
    }

}

