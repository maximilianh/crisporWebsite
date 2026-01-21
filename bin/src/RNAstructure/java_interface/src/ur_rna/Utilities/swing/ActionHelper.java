package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.function.Consumer;

/**
 * Helps add actions to a Frame/Window etc.
 */
public abstract class ActionHelper {
    public static final int COMMAND_MODIFIER_MASK = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();
    public enum KeyBindFocus {
        /** The component itself must have the focus */
        Self(JComponent.WHEN_FOCUSED),
        /** The component's top-level window must have the focus */
        Window(JComponent.WHEN_IN_FOCUSED_WINDOW),
        /** The component or any of its child components must have the focus */
        Child(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
        public final int value;
        KeyBindFocus(int numericValue) {
            value = numericValue;
        }
    }
    public static void addAction(JComponent c, final String actionName, final Consumer<ActionEvent> handler) {
        c.getActionMap().put(actionName, new AbstractAction() {
                @Override
                public void actionPerformed(final ActionEvent e) {
                handler.accept(e);
            }
        });
    }
    public static void addKeyBinding(JComponent c, final char keyChar, final String actionName) {
        addKeyBinding(c, KeyStroke.getKeyStroke(keyChar, COMMAND_MODIFIER_MASK), actionName, KeyBindFocus.Child);
    }
    public static void addKeyBinding(JComponent c, final int keyCode, final int modifiers, final String actionName) {
        addKeyBinding(c, KeyStroke.getKeyStroke(keyCode, modifiers), actionName, KeyBindFocus.Child);
    }
    public static void addKeyBinding(JComponent c, final KeyStroke key, final String actionName) {
        addKeyBinding(c, key, actionName, KeyBindFocus.Child);
    }
    public static void addKeyBinding(JComponent c, final KeyStroke key, final String actionName, final KeyBindFocus focusCondition) {
        c.getInputMap(focusCondition.value).put(key, actionName);
    }
    /** Combine addKeyBinding with addAction */
    public static void addKeyAction(JComponent c, final KeyStroke key, final String actionName, final Consumer<ActionEvent> handler) {
        addKeyAction(c, key, actionName, handler, KeyBindFocus.Child);
    }
    /** Combine addKeyBinding with addAction */
    public static void addKeyAction(JComponent c, final KeyStroke key, final String actionName, final Consumer<ActionEvent> handler, final KeyBindFocus focusCondition) {
        addAction(c, actionName, handler);
        addKeyBinding(c, key, actionName, focusCondition);
    }
}
