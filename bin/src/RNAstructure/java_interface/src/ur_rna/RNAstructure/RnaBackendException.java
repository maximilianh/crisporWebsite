package ur_rna.RNAstructure;

import java.util.Objects;

/** Encapsulates errors raised in RNABackendJNI functions, e.g. loading and saving CT files etc. */
public class RnaBackendException extends Exception {
    private static final long serialVersionUID = 1L;

    private String file;
    private int errorCode;

    public RnaBackendException() { }
    public RnaBackendException(String message) { super(message);  }
    public RnaBackendException(Throwable cause) {
        super(cause);
    }
    public RnaBackendException(String message, Throwable cause) {
        super(message, cause);
    }
    public RnaBackendException(String message, int errCode) {
        this(message, errCode, null, null);
    }
    public RnaBackendException(String message, int errCode, String filePath) {
        this(message, errCode, filePath, null);
    }
    public RnaBackendException(String message, int errCode, String filePath, Throwable cause) {
        super(message);
        if (filePath != null) this.file = filePath;
        this.errorCode = errCode;
        if (cause != null) initCause(cause);
    }

    public String getFile() {
        return file;
    }
    public void initFile(final String file) {
        if (this.file != null)
            throw new IllegalStateException("Can't overwrite file with " + Objects.toString(file, "a null") + ".", this);
        this.file = file;
    }
    public int getErrorCode() {
        return errorCode;
    }
    public void initErrorCode(final int errorCode) {
        if (errorCode != 0)
            throw new IllegalStateException("Can't overwrite error code with " + errorCode + ".", this);
        this.errorCode = errorCode;
    }

    @Override
    public String toString() {
        return getLocalizedMessage() + (errorCode==0?"":" (RNAstructure ErrorCode: " + errorCode+")");
    }
}