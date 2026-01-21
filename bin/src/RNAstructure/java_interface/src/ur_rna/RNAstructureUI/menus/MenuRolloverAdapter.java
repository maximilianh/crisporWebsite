package ur_rna.RNAstructureUI.menus;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.Serializable;
import java.util.HashMap;
import java.util.function.Consumer;

/**
 * An inner class that creates a rollover listener for menus in the
 * RNAstructure GUI.
 * <br><br>
 * The rollover text is placed in the message bar at the bottom of the main
 * application frame, whenever a user interacts with the menu item this
 * listener is attached to, and removed to make way for the default message
 * text when the user stops interacting with the menu item.
 *
 * @author Jessica S. Reuter
 * @author Richard M. Watson
 */
public class MenuRolloverAdapter implements MouseListener, Serializable {
    private static final long serialVersionUID = 20120802;

    private HashMap<Object, String> rollovers = new HashMap<>();
    private final Consumer<String> rolloverHandler;

    public MenuRolloverAdapter(Consumer<String> rolloverHandler) { this.rolloverHandler = rolloverHandler; }

    public void add(JMenuItem item, String rolloverText) {
        rollovers.put(item, rolloverText);
        item.addMouseListener(this);
    }

    public static MenuRolloverAdapter fromToolTips(Consumer<String> rolloverHandler, JComponent parentMenu) {
        MenuRolloverAdapter ra = new MenuRolloverAdapter(rolloverHandler);
        ra.addFromToolTips(parentMenu);
        return ra;
    }

    /**
     * Convert tooltips to rollover text.
     */
    public void addFromToolTips(JComponent parentMenu) {
        for( Component c : parentMenu.getComponents()) {
            if( c instanceof JMenuItem ) {
                JMenuItem m = (JMenuItem)c;
                String tooltip = m.getToolTipText();
                if (tooltip != null && tooltip.length() != 0) {
                    // m.setToolTipText(null);
                    add(m, tooltip);
                }
                if (m.getComponentCount() != 0)
                    addFromToolTips(m);
            }
        }
    }

    @Override
    public void mouseClicked(MouseEvent e) { reset(); }

    @Override
    public void mouseEntered(MouseEvent e) { set(e.getSource()); }

    @Override
    public void mouseExited(MouseEvent e) { reset(); }

    @Override
    public void mousePressed(MouseEvent e) { reset(); }

    @Override
    public void mouseReleased(MouseEvent e) { reset(); }

    /**
     * Reset the main frame message bar text to its default.
     */
    private void reset() {
        rolloverHandler.accept(null);
    }

    /**
     * Set the main frame message bar text to this listener's specific
     * rollover text.
     */
    private void set(Object ref) {
        rolloverHandler.accept(rollovers.get(ref));
    }
}
