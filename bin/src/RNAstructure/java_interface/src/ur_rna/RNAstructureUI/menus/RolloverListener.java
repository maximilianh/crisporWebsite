package ur_rna.RNAstructureUI.menus;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.Serializable;
import java.util.function.Consumer;

/**
 * A class that creates a rollover listener for menus in the RNAstructure GUI.
 * <br><br>
 * The rollover text is placed in the message bar at the bottom of the main
 * application frame, whenever a user interacts with the menu item this
 * listener is attached to, and removed to make way for the default message
 * text when the user stops interacting with the menu item.
 *
 * @author Jessica S. Reuter
 * @author Richard M. Watson
 */
public class RolloverListener implements MouseListener, Serializable {
    private static final long serialVersionUID = 20160316;

    //private HashMap<Object, String> rollovers = new HashMap<>();
    private final Consumer<String> handler;
    private final String text;
    public RolloverListener(Consumer<String> rolloverHandler, String rolloverText) {
        handler = rolloverHandler;
        text = rolloverText;
    }

    @Override
    public void mouseClicked(MouseEvent e) { reset(); }

    @Override
    public void mouseEntered(MouseEvent e) { set(text); }

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
        handler.accept(null);
    }

    /**
     * Set the main frame message bar text to this listener's specific
     * rollover text.
     */
    private void set(String rolloverText) {        handler.accept(rolloverText);    }
}
