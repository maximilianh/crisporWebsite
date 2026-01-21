/**
 * 
 */
package ur_rna.Utilities;

import ur_rna.Utilities.annotation.Nullable;

/**
 * Indicates that there is a format or syntax error in the string argument passed to the method.
 * <p>
 * This should be thrown only from methods that parse a string (e.g. a command-line parser, etc) whenever the
 * format of the string argument is wrong. This should NOT be used for internal parsing errors, but rather when
 * string itself has an identifiable syntax error.
 * These derive from {@link Exception} (rather than the unchecked {@link RuntimeException} because the caller
 * should be aware of the potential for errors in user-provided input.
 * </p>
 */
public class SyntaxErrorException extends Exception {
	private static final long serialVersionUID = 1L;

	/**
	 * {@link SyntaxErrorException#leftCodeQuote } and {@link SyntaxErrorException#rightCodeQuote } are the quote
	 * characters used to distinguish the offending segment of code in the error message.
	 */
	public static String leftCodeQuote = ">>";
	/**
	 * {@link SyntaxErrorException#leftCodeQuote } and {@link SyntaxErrorException#rightCodeQuote } are the quote
	 * characters used to distinguish the offending segment of code in the error message.
	 */
	public static String rightCodeQuote = "<<";

	/** A brief description or name of the error, e.g. "Invalid number format" or "Unterminated string constant". */
	protected String errorDescription;

	/** The part of the subject string most closely associated with the error. */
	protected String offendingText;

	/**
	 * Information about the location of the offending segment in the subject, so the user can find it easily.
	 * If set to any value other than null, this overrides the default return value of {@link #getLocation}.
	 */
	protected String location;

	/** The estimated line where the syntax error occurred. */
	protected int lineOffset = -1;

	/** The estimated column number or character offset of the start of the syntax error. */
	protected int columnOffset = -1;

	/**
	 * Default constructor. Provides no information about the details of the error.
	 */
	public SyntaxErrorException() {	}

	/**
	 * Basic exception constructor, delegated to super class. {@link Throwable#Throwable(String) }
	 * Does not fill in details of the error.
	 */
	public SyntaxErrorException(String message) { super(message); }
	public SyntaxErrorException(Throwable cause) {
		super(cause);
	}

	/**
	 * Basic exception constructor, delegated to super class. {@link Throwable#Throwable(String,Throwable) }
	 * Does not fill in details of the error.
	 */
	public SyntaxErrorException(String message, Throwable cause) {
		super(message, cause);
	}

	public SyntaxErrorException(@Nullable String message, @Nullable final String briefErrorDescription,
								@Nullable final String offendingText, @Nullable final String errorLocationDescription,
								final int line, final int column, @Nullable final Throwable cause) {
		super(
				message == null || message.isEmpty() ?
						buildMessage(briefErrorDescription, offendingText, line, column, null) :
						message
		);
		if (cause != null) initCause(cause);
		this.errorDescription = briefErrorDescription;
		this.offendingText = offendingText;
		this.location = errorLocationDescription;
		this.lineOffset = line;
		this.columnOffset = column;
	}

	public SyntaxErrorException(@Nullable String errorDescription, @Nullable String offendingText, int line, int column) {
		this(null, errorDescription, offendingText, null, line, column, null);
	}
	public SyntaxErrorException(@Nullable String errorDescription, @Nullable String offendingText, String locationDescription) {
		this(null, errorDescription, offendingText, locationDescription, -1, -1, null);
	}
	private static String buildMessage(final String error, final String text, final int line, final int column, final String location) {
		/*   Example:  Unterminated String Constant (line 1, col 2) >>"hello there<<
		*/
		StringBuilderEx sb = new StringBuilderEx();
		if (error != null)
			sb.append(error).trimEnd(' ').appendSeparator("."); //add period at end (if there isn't already one there)

		buildLocationInfo(sb, line, column, location);

		if (text != null) {
			sb.appendSeparator(" ")
					.append(leftCodeQuote)
					.append(Strings.escapeStringLiteral(text, CharType.NON_PRINTABLE))
					.append(rightCodeQuote);
		}
		return sb.toString();
	}

	protected static void buildLocationInfo(final StringBuilderEx sb, final int line, final int column, final String location) {
		if (location != null)
			sb.appendSeparator(" ").append(location);

		if (line != -1 || column != -1) {
			sb.appendSeparator(" ").append("(");
			if (line != -1)
				sb.append("line ").append(line);

			if (column != -1) {
				if (line != -1)
					sb.append(", ");
				sb.append("col ").append(column);
			}
			sb.append(")");
		}
	}

	public String getErrorDescription() {
		return errorDescription;
	}
	public String getOffendingText() {
		return offendingText;
	}
	public String getLocationDescription() { return location; }
	public String getLocation() {
		StringBuilderEx sb = new StringBuilderEx();
		buildLocationInfo(sb, lineOffset, columnOffset, location);
		return sb.toString();
	}
	public boolean hasLocationInfo() {
		return lineOffset != -1 || columnOffset != -1 || location != null;
	}
	public int getLine() {
		return lineOffset;
	}
	public int getColumn() {
		return columnOffset;
	}
}
