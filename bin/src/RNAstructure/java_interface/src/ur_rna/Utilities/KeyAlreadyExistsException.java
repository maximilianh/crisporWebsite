package ur_rna.Utilities;

public class KeyAlreadyExistsException extends RuntimeException {
    private static final long serialVersionUID = 1L;
    protected static String defaultMessage = "The specified key already exists in the collection.";
    private Object existingKey = this;
    private Object existingValue;

    public KeyAlreadyExistsException() { this(defaultMessage); }
    public KeyAlreadyExistsException(String message) {
        super(message);
    }
    public KeyAlreadyExistsException(Throwable cause) {
        super(defaultMessage, cause);
    }
    public KeyAlreadyExistsException(String message, Throwable cause) {
        super(message, cause);
    }
    public KeyAlreadyExistsException initDetails(Object existingKey, Object existingValue) {
        if (this.existingKey != this)
            throw new IllegalStateException("Cannot specify the key after it has already been set.", this);
        if (existingKey == this)
            throw new IllegalArgumentException("Cannot set the key to this " + this.getClass().getSimpleName() + ".", this);
        this.existingKey = existingKey;
        this.existingValue = existingValue;
        return this;
    }
    public Object getExistingKey() {
        return existingKey == this ? null : existingKey;
    }
    public Object getExistingValue() {
        return existingKey == this ? null : existingValue;
    }
}