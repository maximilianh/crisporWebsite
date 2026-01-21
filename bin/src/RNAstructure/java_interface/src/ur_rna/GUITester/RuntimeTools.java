package ur_rna.GUITester;

import static ur_rna.Utilities.ObjTools.toDisplayString;

public abstract class RuntimeTools {
    public static void expectArgs(String action, Object[] args, int minCount, int maxCount) throws
            ActionArgumentException {
        if (args.length < minCount || (maxCount > -1 && args.length > maxCount)) {
            String message = String.format(
                    "The action \"%s\" expects %s argument(s), but it received %d.",
                    action,
                    minCount == 0 && maxCount <= minCount ? "no" :
                            (minCount == maxCount ? "exactly " + minCount :
                                    (maxCount == -1 ? minCount + " or more" :
                                            "between " + minCount + " and " + maxCount)
                            ),
                    args.length
            );
            throw new ActionArgumentException(message);
        }
    }

    @SuppressWarnings("unchecked") // isAssignableFrom results in a warning about using unchecked or unsafe operations.
    public static void expectArgType(String action, Object[] args, int ordinal, String argName, boolean allowNull, Class... expected)
            throws ActionArgumentException {
        Object value = args[ordinal -1];
        Class foundClass = value == null ? null : value.getClass();
        if (foundClass == null && allowNull) return;
        if (foundClass != null) {
            for (Class c : expected)
                if (c.isAssignableFrom(foundClass))
                    return;
        }

        StringBuilder sbTypes = new StringBuilder();
        for (int i = 0; i < expected.length; i++) {
            if (i != 0)
                sbTypes.append(" or ");
            sbTypes.append(expected[i].getSimpleName());
        }
        String message = String.format(
                "The action \"%s\" expects an argument of type %s for parameter %d%s, but %s was found instead.",
                action,
                sbTypes.toString(),
                ordinal,
                ((argName == null || argName.equals("")) ? "" : " (" + argName + ")"),
                foundClass == null ? "null" : "the type " + foundClass.getSimpleName()
        );
        throw new ActionArgumentException(message);
    }

    public static void getInvalidArgMessage(String action, Object[] args, int ordinal, String argName)
            throws ActionArgumentException {
        Object value = args[ordinal -1];
        String message = String.format(
                "Invalid argument value for action \"%s\" in parameter %d%s. Value: %s.",
                action,
                ordinal,
                ((argName == null || argName.equals("")) ? "" : " (" + argName + ")"),
                toDisplayString(value)
        );
        throw new ActionArgumentException(message);
    }

    /**
     * Base class for exceptions thrown from  {@code doAction}
     */
    public static class ActionException extends ScriptRuntimeException {
        public ActionException() {}
        public ActionException(String message) {super(message); }
        public ActionException(Throwable cause) { super(cause); }
        public ActionException(String message, Throwable cause) {
            super(message);
            if (cause != null) this.initCause(cause);
        }
    }

    /**
     * Thrown if the action is unknown or not implemented.
     */
    public static class UnsupportedActionException extends ActionException {
        public UnsupportedActionException() { }
        public UnsupportedActionException(String action) {
            super("Unsupported action: \"" + action + "\".");
        }
        public UnsupportedActionException(Throwable cause) {
            super(cause);
        }
        public UnsupportedActionException(String fullMessage, Throwable cause) { super(fullMessage, cause); }
    }

    /**
     * Thrown if the number or type of arguments is invalid for the given action.
     */
    public static class ActionArgumentException extends ActionException {
        public ActionArgumentException() { }
        public ActionArgumentException(String message) {
            super(message);
        }
        public ActionArgumentException(Throwable cause) {
            super(cause);
        }
        public ActionArgumentException(String message, Throwable cause) {
            super(message, cause);
        }
    }
}
