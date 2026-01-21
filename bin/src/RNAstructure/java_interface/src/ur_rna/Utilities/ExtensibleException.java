package ur_rna.Utilities;

import java.io.Serializable;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Exception that allows additional details to be added.
 */
public class ExtensibleException extends Exception {
    private static final long serialVersionUID = 1L;

    public ExtensibleException() {}
    protected String updatedMessage;
    protected HashMap<String, Serializable> data;
    protected Map<String, Serializable> readOnlyData;

    public ExtensibleException(String message) {
        super(message);
    }
    public ExtensibleException(Throwable cause) { super(cause); }

    public ExtensibleException(String message, Throwable cause) {
        super(message);
        initCause(cause);
    }

    @SuppressWarnings("unchecked")
    public Map<String, Serializable> getData() {
        if (data == null)
            return Collections.EMPTY_MAP;
        if (readOnlyData == null)
            readOnlyData = Collections.unmodifiableMap(data);
        return readOnlyData;
    }
    public synchronized ExtensibleException putData(String key, Serializable value) {
        return putData(key, value, true);
    }
    public synchronized ExtensibleException putData(String key, Serializable value, boolean overwrite) {
        if (data == null) data = new HashMap<>();
        if (overwrite)
            data.put(key, value);
        else
            data.putIfAbsent(key, value);
        return this;
    }
    /**
     * Returns the value stored in the data table for the specified key if it exists, or null otherwise.
     */
    public synchronized Serializable getData(String key) {
        return getData(key, null);
    }
    /**
     * Returns the value stored in the data table for the specified key if it exists, or defaultValue otherwise.
     */
    public synchronized Serializable getData(String key, Serializable defaultValue) {
        if (data == null)
            return defaultValue;
        return data.getOrDefault(key, defaultValue);
    }

    /**
     * Allows the calling code to update the message stored in this exception.
     * @param value The new message.
     */
    public void setMessage(String value) {
        updatedMessage = value;
    }
    /**
     * This returns the original message passed in the constructor (to Throwable).
     * This will only differ from getMessage if {@code setMessage} has been used
     * to store an updated message.
     */
    public String getOriginalMessage() {
        return super.getMessage();
    }
    /**
     * If a new message has been stored (using {@code setMessage}) this will
     * revert the message back to the original (passed in the constructor.)
     */
    public void resetMessage() {
        setMessage(null);
    }
    /**
     * If a new message has been stored (using {@code setMessage}) this will return
     * the new message. Otherwise, the result is the same as calling {@link Throwable#getMessage()}
     */
    @Override
    public String getMessage() {
        return updatedMessage == null ? super.getMessage() : updatedMessage;
    }
}
