package ur_rna.Utilities.swing;

import javax.swing.*;

/**
 * @author Richard M. Watson
 */
public class FormValidationException extends Exception {
    private static final long serialVersionUID = 20170810L;
    private JComponent control;

    public FormValidationException() { }
    public FormValidationException(String message) {
        this(message, null);
    }
    public FormValidationException(String message, JComponent control) {
        super(message);
        this.control = control;
    }
    public FormValidationException(String message, JComponent control, Throwable cause) {
        super(message, cause);
        this.control = control;
    }
    public JComponent getControl() { return control; }
}
