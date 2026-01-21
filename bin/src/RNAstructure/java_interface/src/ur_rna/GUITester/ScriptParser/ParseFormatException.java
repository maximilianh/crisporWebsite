package ur_rna.GUITester.ScriptParser;

import static ur_rna.Utilities.Strings.escapeStringLiteral;

/**
 * This exception is thrown when high-level parse errors, such as 
 * evaluation of a literal are encountered.
 */
public class ParseFormatException extends ParseException {

  /**
   * The version identifier for this Serializable class.
   * Increment only if the <i>serialized</i> form of the
   * class changes.
   */
  private static final long serialVersionUID = 1L;

  /**
   * This constructor is used by the method "generateParseException"
   * in the generated parser.  Calling this constructor generates
   * a new object of this type with the fields "currentToken",
   * "expectedTokenSequences", and "tokenImage" set.
   */
  public ParseFormatException(String message, Token relatedToken, Throwable cause)
  {
    this(message + "  " + getInfoFromToken(relatedToken));
	if (cause != null) initCause(cause);
    currentToken = relatedToken;
  }
  
  public ParseFormatException(String message, Token relatedToken) { this(message, relatedToken, null); }

  public ParseFormatException() {
    super();
  }

  /** Constructor with message. */
  public ParseFormatException(String message) {
    super(message);
  }

  /**
   * Uses "currentToken" to generate additional information to help
   * the user locate the source of the error.
   */
  public static String getInfoFromToken(Token t) {
	String tokenInfo = GuiTestScriptParserConstants.tokenImage[t.kind];
	if (!tokenInfo.startsWith("\""))
		tokenInfo += ": \"" + escapeStringLiteral(t.image) + "\"";

    return String.format("Offending token: %s at line %d, column %d.", tokenInfo, t.beginLine, t.beginColumn);
  }
}