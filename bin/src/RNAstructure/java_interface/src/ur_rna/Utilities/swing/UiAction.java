package ur_rna.Utilities.swing;

import ur_rna.Utilities.Convert;
import ur_rna.Utilities.ObjTools;
import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.annotation.Nullable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashSet;
import java.util.Set;

/**
 * Extension of the swing {@link AbstractAction} class that accepts various methods of handling the action.
 */
 public class UiAction extends AbstractAction {
    public static final String TEXT_KEY = "HandlerTextKey";

    private static ActionListener[] EMPTY_LISTENERS = new ActionListener[0];
    private ActionListener listener;
    private Set<ActionListener> listeners = null;
    public UiAction(@Nullable String name, @NotNull final Runnable runnable) { this(name, runnable, null); }
    public UiAction(@Nullable String name, @NotNull final Runnable runnable, @Nullable Icon icon) { this(name, e->runnable.run(), icon); }
    public UiAction(@Nullable String name, @Nullable final ActionListener listener) { this(name, listener, null); }
    public UiAction(@Nullable String name, @Nullable ActionListener listener, @Nullable Icon icon) {
        this.listener = listener;
        if (name != null) {
            setName(name);
            setDesc(KeyMnemonic.stripMnemonics(name));
        }
        if (icon != null) setIcon(icon);
    }
    /**
     * Get the name of this Action. This is the same as {@code getValue(Action.NAME) }
     * if that property has been set to a String.
     */
    public String name() { return getName(this); }
    //public String text() { return getText(this); }
    public int mnemonic() { return  getMnemonic(this); }
    public Icon icon() { return  getIcon(this); }
    public String desc() { return  getDesc(this); }
    public KeyStroke keyStroke() { return getKeyStroke(this); }
    public String command() { return getCommand(this); }

    public UiAction setName(String text) {
        int m = KeyMnemonic.getMnemonic(text);
        if (m != 0) {
            setMnemonic(m);
            text = KeyMnemonic.stripMnemonics(text);
        }
        putValue(UiAction.NAME, text);
        return this;
    }

    public UiAction setMnemonic(int value) { putValue(Action.MNEMONIC_KEY, value); return this; }
    public UiAction setCommand(String actionCommand) { putValue(Action.ACTION_COMMAND_KEY, actionCommand); return this; }
    public UiAction setIcon(Icon value) { putValue(Action.SMALL_ICON, value); return this; }
    public UiAction setIcon(Image value) { setIcon(new ImageIcon(value)); return this; }
    public UiAction setDesc(String desc) { putValue(Action.SHORT_DESCRIPTION, desc);return this; }
    public UiAction setLongDesc(String desc) { putValue(Action.LONG_DESCRIPTION, desc); return this; }
    public UiAction setKeyStroke(KeyStroke keyStroke) { putValue(Action.ACCELERATOR_KEY, keyStroke); return this; }
    public UiAction setKeyStroke(String keyString) { putValue(Action.ACCELERATOR_KEY, AcceleratorKey.getKey(keyString)); return this; }
    public UiAction setKeyStroke(char keyChar) { putValue(Action.ACCELERATOR_KEY, AcceleratorKey.getKey(keyChar)); return this; }

    /**
     * Get the main listener (which is usually specified in the constructor).
     */
    public ActionListener getListener() {
        return listener;
    }

    /**
     * Get a list of all listeners of this action, including the main listener
     * (which is usually specified in the constructor).
     */
    public ActionListener[] getListeners() {
        if (listeners == null || listeners.size() == 0)
            return listener == null ? EMPTY_LISTENERS : new ActionListener[] { listener };
        if (listener == null) return listeners.toArray(new ActionListener[listeners.size()]);
        ActionListener[] all = new ActionListener[listeners.size() + 1];
        all[0] = listener;
        ObjTools.copyTo(listeners, all, 1);
        return all;
    }

    /**
     * Add an ActionListener to be notified when this Action is invoked.
      * @param newListener The ActionListener to add.
     */
    public void addListener(ActionListener newListener) {
        if (listener == null)
            listener = newListener;
        else {
            if (listeners == null)
                listeners = new HashSet<>();
            listeners.add(newListener);
        }
    }
    /**
     * Remove an ActionListener that should no longer be notified when this Action is invoked.
     * @param existingListener The ActionListener to remove.
     */
    public void removeListener(ActionListener existingListener) {
        if (listener == existingListener)
            listener = null;
        else if (listeners != null)
            listeners.remove(existingListener);
    }

    @Override
    public void actionPerformed(final ActionEvent e) {
        if (listener != null)
            listener.actionPerformed(e);
    }

    public void invoke(Object source) {
        ActionEvent e = new ActionEvent(source, ActionEvent.ACTION_PERFORMED, name(), System.currentTimeMillis(), 0);
        actionPerformed(e);
    }
    public void invoke(Object source, String command) {
        ActionEvent e = new ActionEvent(source, ActionEvent.ACTION_PERFORMED, command, System.currentTimeMillis(), 0);
        actionPerformed(e);
    }

//    public static String getText(final Action a) {
//        Object o = a.getValue(UiAction.TEXT_KEY);
//        if (o == null)
//            o = a.getValue(Action.NAME);
//        return Convert.toString(o);
//    }
    public static String getName(final Action a) {
        return Convert.toString(a.getValue(Action.NAME));
    }
    public static String getDesc(final Action a) {
        return Convert.toString(a.getValue(Action.SHORT_DESCRIPTION));
    }
    public static Icon getIcon(final Action a) {
        Object icon = a.getValue(Action.SMALL_ICON);
        if (icon instanceof Icon)
            return (Icon)icon;
        return null;
    }
    public static int getMnemonic(final Action a) {
        Object o = a.getValue(Action.MNEMONIC_KEY);
        if (o instanceof Number)
            return ((Number)o).intValue();
        if (o instanceof Character)
            return (Character)o;
        return 0;
    }
    public static KeyStroke getKeyStroke(final Action a) {
        Object o = a.getValue(Action.ACCELERATOR_KEY);
        if (o instanceof KeyStroke)
            return (KeyStroke)o;
        if (o instanceof String)
            return AcceleratorKey.getKey((String)o);
        if (o instanceof Character)
            return AcceleratorKey.getKey((Character) o);
        return null;
    }
    public static String getCommand(final Action a) {
        return Convert.toString(a.getValue(Action.ACTION_COMMAND_KEY));
    }
    public static boolean isSelected(final Action a) {
        return Boolean.TRUE.equals(a.getValue(Action.SELECTED_KEY));
    }
    public static void setSelected(final Action a, boolean selected) {
        a.putValue(Action.SELECTED_KEY, selected);
    }
}
