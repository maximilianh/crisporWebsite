package ur_rna.Utilities.swing;

import ur_rna.Utilities.OSInfo;
import ur_rna.Utilities.Strings;

import javax.swing.*;
import java.awt.*;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

/**
 * Utility for getting KeyStrokes from Strings or chars
 */
public class AcceleratorKey {
    /**
     * The integer key code that represents the standard modifier key for menu shortcuts.
     * On Windows and Linux, this is the Control Key ({@link Event#CTRL_MASK}), while on
     * Mac, this is the Command Key ({@link Event#META_MASK}).
     */
    public static int SystemMenuKeyModifier = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(); //on Windows and Linux, this is Event.CTRL_MASK. On Mac, this is Event.META_MASK

    public static KeyStroke getKey(char key) { return KeyStroke.getKeyStroke(key, SystemMenuKeyModifier); }
    public static KeyStroke getKey(String key) {
        if (Strings.isEmpty(key))
            return null;
        int modifier = 0;
        for (int i = 0; i < key.length(); i++) {
            char c = key.charAt(i);
            switch (c) {
                case '*':
                    modifier |= SystemMenuKeyModifier;
                    break;
                case '^':
                    modifier |= InputEvent.CTRL_MASK;
                    break;
                case '#':
                    modifier |= InputEvent.META_MASK;
                    break;
                case '+':
                    modifier |= InputEvent.SHIFT_MASK;
                    break;
                case '%':
                    modifier |= InputEvent.ALT_MASK;
                    break;
                default:
                    key = key.substring(i); // removes modifiers before this position in the string.
                    KeyStroke ks;
                    if (key.length() == 1)
                        ks = KeyStroke.getKeyStroke(Character.toUpperCase(key.charAt(0)), modifier);
                    else {
                        ks = KeyStroke.getKeyStroke(key.toUpperCase());
                        if (ks != null && modifier != 0)
                            ks = KeyStroke.getKeyStroke(ks.getKeyCode(), modifier);
                    }
                    if (ks == null)
                        throw new RuntimeException(String.format("No KeyStroke can be obtained from the string '%s'", key));
                    return ks;
            }
        }
        // We shouldn't get here unless the only characters in the string were modifier characters.
        return KeyStroke.getKeyStroke(modifier, 0);
    }
    /**
     * Return a nice string representation of the given keystroke
     * (the toString method of KeyStroke isn't nice.)
     */
    public static String toString(final KeyStroke key) {
        int modifiers = key.getModifiers();
        StringBuilder buffer = new StringBuilder();

        if (0!=(modifiers & Event.CTRL_MASK))
            buffer.append("Ctrl-");
        if (0!=(modifiers & Event.SHIFT_MASK))
            buffer.append("Shift-");
        if (0!=(modifiers & Event.META_MASK)) {
            if (OSInfo.isMac())
                buffer.append("Command-");
            else if (OSInfo.isWin())
                buffer.append("Win-");
            else
                buffer.append("Meta-");
        }
        if (0!=(modifiers & Event.ALT_MASK))
            buffer.append(OSInfo.isMac()?"Option-":"Alt-");
        buffer.append(KeyEvent.getKeyText(key.getKeyCode()));
        return buffer.toString();
    }

    public static String getShortcutModifierKeyName() { return getShortcutModifierKeyName(false); }
    public static String getShortcutModifierKeyName(boolean abbreviated) {
        return getModifierKeyName(SystemMenuKeyModifier, abbreviated);
    }

    public static String getAltKeyName() { return getAltKeyName(false); }
    public static String getAltKeyName(boolean abbreviated) {
        return OSInfo.isMac()?"Option":"Alt";
    }

    public static String getMetaKeyName() { return getMetaKeyName(false); }
    public static String getMetaKeyName(boolean abbreviated) {
        if (OSInfo.isMac())
            return "Command";
        if (OSInfo.isWin())
            return abbreviated?"Win":"Windows";
        if (OSInfo.isNix())
            return "Super";
        return "Meta";
    }

    public static String getModifierKeyName(int modifierKeyCodeOrMask) { return getModifierKeyName(modifierKeyCodeOrMask, false); }
    public static String getModifierKeyName(int modifierKeyCodeOrMask, boolean abbreviated) {
        switch (modifierKeyCodeOrMask) {
            case KeyEvent.VK_CONTROL:
            case KeyEvent.CTRL_DOWN_MASK:
            case KeyEvent.CTRL_MASK: return abbreviated?"Ctrl":"Control";

            case KeyEvent.VK_SHIFT:
            case KeyEvent.SHIFT_DOWN_MASK:
            case KeyEvent.SHIFT_MASK: return "Shift";

            case KeyEvent.VK_ALT:
            case KeyEvent.ALT_DOWN_MASK:
            case KeyEvent.ALT_MASK: return getAltKeyName(abbreviated);

            case KeyEvent.VK_META:
            case KeyEvent.META_DOWN_MASK:
            case KeyEvent.META_MASK: return getMetaKeyName(abbreviated);

            default: {
                String s = KeyStroke.getKeyStroke(Character.valueOf('A'), modifierKeyCodeOrMask).toString();
                // s should be something like "ctrl shift typed A", where the text before " typed" is the modifier keys in string form.
                int pos = s.indexOf(" typed");
                return pos == -1 ? s : s.substring(0, pos);
            }
        }
    }

    /**
     * Some Components (e.g. JTextArea) disable the default Tab and Shift+Tab behavior,
     * which is to move to the next/previous component respectively.
     * This function resets those behaviors.
     * @param components A list of components (e.g. JTextArea etc) to fix.
     */
    public static void resetTabTraversalKeys(Component ... components) {
        for(Component c : components) {
            c.setFocusTraversalKeys(KeyboardFocusManager.FORWARD_TRAVERSAL_KEYS, null);
            c.setFocusTraversalKeys(KeyboardFocusManager.BACKWARD_TRAVERSAL_KEYS, null);
        }
    }


    /**
     * Disable the default behavior of TAB, Shift+TAB, Ctrl+TAB, Ctrl+Shift+TAB etc,
     * which is to move to the next/previous component respectively.
     * @param components A list of components (e.g. JTextArea etc) to modify.
     */
    public static void removeAllTabTraversalKeys(Component ... components) {
        for(Component c : components) {
            c.setFocusTraversalKeys(KeyboardFocusManager.FORWARD_TRAVERSAL_KEYS, Collections.emptySet());
            c.setFocusTraversalKeys(KeyboardFocusManager.BACKWARD_TRAVERSAL_KEYS, Collections.emptySet());
        }
    }


    /**
     * Disable the default behavior of TAB, Shift+TAB, Ctrl+TAB, Ctrl+Shift+TAB etc,
     * which is to move to the next/previous component respectively.
     * @param component A component (e.g. JTextArea etc) to modify.
     */
    public static void removeTabTraversalKeys(Component component, AWTKeyStroke...removeKeys) {
        HashSet<AWTKeyStroke> set = new HashSet<>(component.getFocusTraversalKeys(KeyboardFocusManager.FORWARD_TRAVERSAL_KEYS));
        set.removeAll(Arrays.asList(removeKeys));
        component.setFocusTraversalKeys(KeyboardFocusManager.FORWARD_TRAVERSAL_KEYS, set);

        set = new HashSet<>(component.getFocusTraversalKeys(KeyboardFocusManager.BACKWARD_TRAVERSAL_KEYS));
        set.removeAll(Arrays.asList(removeKeys));
        component.setFocusTraversalKeys(KeyboardFocusManager.BACKWARD_TRAVERSAL_KEYS, set);
    }
}
