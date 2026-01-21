package ur_rna.Utilities.swing;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

/**
 * A button that lacks the typical boarder and styling of a button.
 * This can be useful for links etc.
 */
public class JFlatButton extends JButton {
    private String link;
    private URI linkURI;

    public JFlatButton(final Icon icon) {
        super(icon);
        makeFlat();
    }
    public JFlatButton(final String text) {
        super(text);
        makeFlat();
    }
    public JFlatButton(final Action a) {
        super(a);
        makeFlat();
    }
    public JFlatButton(final String text, final Icon icon) {
        super(text, icon);
        makeFlat();
    }
    public JFlatButton() {
        makeFlat();
    }
    private void makeFlat() {
        setFocusPainted(false);
        setMargin(new Insets(0, 0, 0, 0));
        setContentAreaFilled(false);
        setBorderPainted(false);
        setOpaque(false);
        setBorder(new EmptyBorder(0,0,0,0));
        setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
    }


    @Override
    protected void fireActionPerformed(final ActionEvent event) {
        if (linkURI!=null)
            try {
                Desktop.getDesktop().browse(linkURI);
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        super.fireActionPerformed(event);
    }
    public void setLink(final String link) {
        this.link = link;
        try {
            linkURI = new URI(link);
        } catch (URISyntaxException ex) {
            linkURI = null;
        }
    }
    public String getLink() { return link; }
    public URI getLinkUri() { return linkURI; }
    public void setLinkUri(final URI linkURI) {
        this.linkURI = linkURI;
        link = linkURI.toString();
    }
}
