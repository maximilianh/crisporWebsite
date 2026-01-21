package ur_rna.Utilities.swing;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;

/**
 * A label that reacts to clicks by calling {@link Action#actionPerformed(ActionEvent)} or {@link ActionListener#actionPerformed(ActionEvent)}
 * with a custom actionCommand.
 * This can be useful for links etc.
 */
public class JActionLabel extends JLabel {
    private String link;
    private URI linkURI;
    private Action action;
    private String actionCommand;

    private ActionListener listener;
    private List<ActionListener> listeners = null;

    public JActionLabel(final Icon icon) {
        this(null, icon, CENTER);
    }
    public JActionLabel(final String text) {
        this(text, text);
    }
    public JActionLabel(final String text, final String actionCommand) {
        this(text, null, CENTER);
        this.actionCommand = actionCommand;
    }
    public JActionLabel(final Action a) {
        this(UiAction.getName(a), UiAction.getCommand(a));
        setAction(a);
    }
    public JActionLabel(final String text, final Icon icon) {
        super(text, icon, SwingConstants.CENTER);
    }
    public JActionLabel(final String text, final Icon icon, int horizontalAlignment) {
        super(text, icon, horizontalAlignment);
        addClickListener();
    }

    public JActionLabel() { addClickListener(); }

    private void addClickListener() {
        this.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(final MouseEvent e) {
                fireActionPerformed(new ActionEvent(this, ActionEvent.ACTION_PERFORMED, actionCommand, System.currentTimeMillis(), e.getModifiersEx()));
            }
        });
    }

    public void setAction(Action a) {
        action = a;
        updateCursor();

        if (action == null) return;

        String name = UiAction.getName(a);

        if (name!=null && (getText()==null||getText().isEmpty()))
            setText(name);

        if (getName()==null||getName().isEmpty())
            setName(name);

        if (getActionCommand()==null)
            setActionCommand(UiAction.getCommand(a));

        if (getIcon()==null) {
            Icon i = UiAction.getIcon(a);
            if (i != null) setIcon(i);
        }

        if (getToolTipText()==null) {
            String s = UiAction.getDesc(a);
            if (s != null) setToolTipText(s);
        }

        if (getDisplayedMnemonic() <= 0) {
            int m = UiAction.getMnemonic(a);
            if (m != 0) setDisplayedMnemonic(m);
        }
    }

    protected void fireActionPerformed(final ActionEvent event) {
        if (linkURI!=null)
            try {
                Desktop.getDesktop().browse(linkURI);
            } catch (IOException ex) {
                ex.printStackTrace();
            }

        if (action != null)
            action.actionPerformed(event);

        if (listener != null)
            listener.actionPerformed(event);

        if (listeners != null)
            for (ActionListener l : listeners)
                l.actionPerformed(event);
    }

    public void setLink(final String link) {
        this.link = link;
        try {
            linkURI = new URI(link);
        } catch (URISyntaxException ex) {
            linkURI = null;
        }
        updateCursor();
    }

    public String getLink() { return link; }
    public URI getLinkUri() { return linkURI; }
    public void setLinkUri(final URI linkURI) {
        this.linkURI = linkURI;
        link = linkURI.toString();
        updateCursor();
    }
    public Action getAction() {
        return action;
    }
    public String getActionCommand() {
        return actionCommand;
    }
    public void setActionCommand(final String actionCommand) {
        this.actionCommand = actionCommand;
    }

    public void addActionListener(ActionListener l) {
        if (listener==null)
            listener = l;
        else {
            if (listeners == null)
                listeners = new ArrayList<>();
            listeners.add(l);
        }
        updateCursor();
    }

    private void updateCursor() {
        if (hasClickAction())
            setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        else
            setCursor(Cursor.getDefaultCursor());
    }

    /**
     * Returns true if some action is performed on click. I.e. if any ActionListeners are registered or if
     * an Action or Link has been set.
     */
    public boolean hasClickAction() {
        return listener!=null||(listeners!=null&&listeners.size()!=0)||action!=null||linkURI!=null;
    }

    public boolean removeActionListener(ActionListener l) {
        boolean found = false;
        if (listener==l) {
            listener = null;
            found = true;
        } else if (listeners != null)
            found = listeners.remove(l);
        updateCursor();
        return found;
    }

    public ActionListener[] getActionListeners() {
        if ((listeners == null || listeners.size()==0)) {
            if (listener == null)
                return new ActionListener[0];
            return new ActionListener[]{listener};
        } else {
            // listeners is NOT empty.
            if (listener == null)
                return listeners.toArray(new ActionListener[listeners.size()]);

            ActionListener[] all = new ActionListener[listeners.size() + 1];
            all[0] = listener;
            for (int i = listeners.size() - 1; i >= 0; i--)
                all[i + 1] = listeners.get(0);
            return all;
        }
    }
}
