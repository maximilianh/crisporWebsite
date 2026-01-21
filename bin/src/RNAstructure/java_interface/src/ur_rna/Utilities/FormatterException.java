package ur_rna.Utilities;

/**
 * Represents an error converting an Object or Value to a String.
 */
public class FormatterException extends Exception {
    private static final long serialVersionUID = 1L;

    public FormatterException() { }
    public FormatterException(String message) {
        super(message);
    }
    public FormatterException(Throwable cause) {
        super(cause);
    }
    public FormatterException(String message, Throwable cause) {
        super(message, cause);
    }
}
