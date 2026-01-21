package ur_rna.GUITester;

import ur_rna.GUITester.ScriptParser.ScriptNode;
import ur_rna.Utilities.ExtensibleException;

public class ScriptRuntimeException extends ExtensibleException {
    private static final long serialVersionUID = 1L;
    private ScriptNode _node;
    /**
     * Default constructor.
     */
    public ScriptRuntimeException() {
    }

    /**
     * @param message
     */
    public ScriptRuntimeException(String message) {
        super(message);
    }

    /**
     * @param cause
     */
    public ScriptRuntimeException(Throwable cause) {
        super(cause);
    }

    /**
     * @param message
     * @param cause
     */
    public ScriptRuntimeException(String message, Throwable cause) {
        super(message, cause);
    }

    /**
     * @param message
     * @param cause
     */
    public ScriptRuntimeException(String message, Throwable cause, ScriptNode currentNode) {
        super(message);
        _node = currentNode;
        if (cause != null) initCause(cause);
    }

    /**
     * Get the current script node at the time of the error.
     * @return the current script node at the time of the error.
     */
    public ScriptNode getNode() { return _node;}

    public ScriptRuntimeException setNode(ScriptNode node) { return setNode(node, true); }
    public ScriptRuntimeException setNode(ScriptNode node, final boolean overwrite) {
        if (this._node==null || overwrite)
            this._node = node;
        return this;
    }

//    public String getMessage() {
//        if (super.getMessage() == null && argumentName != null)
//            return String.format("Invalid or malformed command-line argument '%s'.", argumentName);
//        return super.getMessage();
//    }
}
