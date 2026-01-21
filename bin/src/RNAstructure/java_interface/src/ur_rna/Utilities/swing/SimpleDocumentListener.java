package ur_rna.Utilities.swing;

import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.JTextComponent;

/**
 * Similar to TextDocListener, but this is lighter weight and allows for lambda syntax.
 * Usage: {@code
 * field.getDocument().addDocumentListener((SimpleDocumentListener) e -> {
 *   // handler code
 * });
 * }
  */
@FunctionalInterface
public interface SimpleDocumentListener extends DocumentListener {
    void update(DocumentEvent e);

    static void listen(JTextComponent c, SimpleDocumentListener listener) {
        c.getDocument().addDocumentListener(listener);
    }

    @Override
    default void insertUpdate(DocumentEvent e) {
        update(e);
    }
    @Override
    default void removeUpdate(DocumentEvent e) {
        update(e);
    }
    @Override
    default void changedUpdate(DocumentEvent e) {
        update(e);
    }
}