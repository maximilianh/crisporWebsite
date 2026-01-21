package ur_rna.Utilities;

import ur_rna.Utilities.annotation.MagicConstant;
import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.annotation.Nullable;

import java.io.*;
import java.nio.charset.Charset;
import java.security.NoSuchAlgorithmException;
import java.text.ParseException;
import java.util.Arrays;
import java.util.function.Function;

/**
 * A collection of String utilities.
 */
public abstract class Strings {

    public static String fmt(String format, Object... args) {
        return String.format(format, args);
    }
    /**
     * Used to convert raw characters to their escaped version
     * when these raw version cannot be used as part of a Java
     * string literal.
     * <p>
     * Escapes NULL, BKSP, TAB, CR, LF, FF, backslash, double (") and single (') quotes
     * to their java escape codes. Also escapes all other control characters (1 to 31) and
     * all characters above 127 to the corresponding Java unicode escape code (\u0001, \u26BD etc. )
     * </p>
     *
     * @param subject The string to escape.
     * @return The escaped string.
     * @see #escapeStringLiteral(String,int)
     */
    public static String escapeStringLiteral(String subject) {
        return escapeStringLiteral(subject, CharType.ALL);
    }
    /**
     * Used to convert raw characters to their escaped version
     * when these raw version cannot be used as part of a Java
     * string literal.
     * <p>
     * Escapes NULL, BKSP, TAB, CR, LF, FF, backslash, double (") and single (') quotes
     * to their java escape codes. Also escapes all other control characters (1 to 31) and
     * all characters above 127 to the corresponding Java unicode escape code (\u0001, \u26BD etc. )
     * </p>
     *
     * @param subject The string to escape.
     *
     * @param escapeOptions
     *   A bitfield of {@link CharType} values indicating which types
     *   of characters to escape.
     *   <p>
     *   However note {@link CharType#CONTROL} characters are ALWAYS escaped, while
     *   {@link CharType#ALPHA}, {@link CharType#DIGITS}, and most {@link CharType#SYMBOLS}
     *   are NEVER escaped.
     *   The only character types for which escaping can be enabled or disabled are the following:
     *   {@link CharType#TAB}, {@link CharType#NEWLINE}, {@link CharType#DOUBLE_QUOTES},
     *   {@link CharType#SINGLE_QUOTE}, {@link CharType#BACKSLASH}, and {@link CharType#EXT_ASCII}.
     *   </p>
     *
     * @return The escaped string.
     * @see #escapeStringLiteral(String, int)
     */
    public static String escapeStringLiteral(String subject,
            @MagicConstant(flagsFromClass = CharType.class) int escapeOptions) {
        boolean bDQ = ObjTools.isSet(CharType.DOUBLE_QUOTES, escapeOptions);
        boolean bSQ = ObjTools.isSet(CharType.SINGLE_QUOTE, escapeOptions);
        boolean bEA = ObjTools.isSet(CharType.EXT_ASCII, escapeOptions);
        boolean bSL = ObjTools.isSet(CharType.BACKSLASH, escapeOptions);

        StringBuilder sb = new StringBuilder(subject.length());
        char c;
        for (int i = 0; i < subject.length(); i++) {
            switch (c = subject.charAt(i)) {
                case 0:
                    sb.append("\\0");
                    continue;
                case CharConstants.BKSP:
                    sb.append("\\b");
                    continue;
                case CharConstants.TAB:
                    sb.append("\\t");
                    continue;
                case CharConstants.LF:
                    sb.append("\\n");
                    continue;
                case CharConstants.FF:
                    sb.append("\\f");
                    continue;
                case CharConstants.CR:
                    sb.append("\\r");
                    continue;
                case CharConstants.DBL_QUOT:
                    sb.append(bDQ ? "\\\"" : c);
                    continue;
                case CharConstants.SNG_QUOT:
                    sb.append(bSQ ? "\\\'" : c);
                    continue;
                case CharConstants.BACK_SLASH:
                    sb.append(bSL ? "\\\\" : c);
                    continue;
                case CharConstants.DEL:
                    sb.append("\\u007f");
                default:
                    if (c < CharConstants.SP || bEA && c >= CharConstants.UPPER_ASCII || c > CharConstants.LAST_ASCII) {
                        String s = "000" + Integer.toString(c, Constants.HEX_RADIX);
                        sb.append("\\u");
                        sb.append(s, s.length() - 4, s.length());
                    } else {
                        sb.append(c);
                    }
            }
        }
        return sb.toString();
    }

    /**
     * UN-escapes a string by replacing backslash-escaped sequences with
     * the character literals they represent.
     * Handles the following cases:
     *   {@code \\  \"  \' } are replaced by backslash (\), double-quote ("), and apostrophe (') respectively.
     *   {@code \r \n \t \f \b \a } are replaced by CR, LF, TAB, FormFeed, and BackSpace respectively.
     *   {@code \e  \a } are replaced by ESC and ALERT respectively
     *   {@code \cX } is replaced by the control character represented by X  (e.g. \cM is the carriage return (CR, 0x0D, #13)
     *   {@code \O \OO \OOO } are replaced by their octal character values (where 'O' can be any digit 0-7)
     *   {@code \xHH } is replaced by its hexadecimal character value (where 'H' is a hex digit: 0-9, A-F, or a-f)
     *   {@code \\uHHHH } is replaced by its unicode character value ('H' as above).
     *   {@code \UHHHHHHHH } is replaced by its unicode character value ('H' as above).
     *
     *   Unknown backslash escape sequences are passed through as-is. They do not result in errors.
     *   E.g. "c:\Windows" yields "c:\Windows"   whereas "c:\Users" yields ERROR (invalid \U escape)
     *   Similarly, a trailing backslash is not considered an error and is simply passed through. "hi\" yields "hi\"
     *
     * @param subject
     * @return The un-escaped string
     * @throws ParseException if any of the {@code \c \O \x \\u or \U } escapes is NOT followed by a valid
     *         sequence of valid characters of the appropriate length.
     */
    @NotNull
    public static String unescapeStringLiteral(@NotNull String subject) throws ParseException {
        StringBuilder sb = new StringBuilder(subject.length());
        boolean prev_backslash = false;

        for (int i = 0; i < subject.length(); i++) {
            int cp = subject.codePointAt(i);
            if (subject.codePointAt(i) > Character.MAX_VALUE)
                i++; /* UTF-16 */

            if (!prev_backslash) {
                if (cp == '\\')
                    prev_backslash = true;
                else
                    sb.append(Character.toChars(cp));
                continue;
            }

            switch (cp) {
                case '\\': sb.append('\\'); break;
                case '\"':  sb.append('\"'); break;
                case '\'':  sb.append('\''); break;
                case 'r':  sb.append('\r'); break;
                case 'n':  sb.append('\n'); break;
                case 't':  sb.append('\t'); break;
                case 'f':  sb.append(CharConstants.FF); break;
                case 'b':  sb.append(CharConstants.BKSP); break;
                case 'a':  sb.append(CharConstants.ALERT); break;
                case 'e':  sb.append(CharConstants.ESC); break;

            /*
             * A "control" character: \cX
             * X = cp xor 64 (aka '@')
             * e.g. \cM is the carriage return (CR, 0x0D, #13) because (('M'=77) ^ 64) yields 13
             * Note that cp need not actually be a control character, e.g.  '\c{' == ';'
             */
                case 'c':   {
                    if (++i == subject.length()) { throw parseError("trailing \\c", i-1); }
                    cp = subject.codePointAt(i);
                     if (cp > 0x7f) { throw parseError("expected ASCII character after \\c", i); }
                    sb.append(Character.toChars(cp ^ 64));
                    break;
                } /* end of control chars */

//                case '8':
//                case '9':
//                    throw parseError("illegal octal digit", i);
                 /*
                 * Octal: 1 to 3 octal digits, e.g \0 \33 \775
                 */
                case '0': case '1': case '2': case '3': case '4': case '5': case '6':
                case '7': {
                    int digits = 1;
                    while (digits < 3) {
                        if (i+digits == subject.length())
                            break;
                        int ch = subject.charAt(i+digits);
                        if (ch < '0' || ch > '7')
                            break;
                        digits++;
                    }
                    try {
                        sb.append(Character.toChars(Integer.parseInt(subject.substring(i, i+digits), Constants.OCTAL_RADIX)));
                        i += digits-1;
                    } catch (NumberFormatException nfe) {
                        throw parseError("invalid octal value for \\0 escape", i);
                    }
                    break;
                } /* end of octals */

                case 'x':  {
                    if (i+2 >= subject.length())
                        throw parseError("string too short for \\x escape", i);
                    try {
                        i++;
                        sb.append(Character.toChars(Integer.parseInt(subject.substring(i, i+2), Constants.HEX_RADIX)));
                        i+=1;
                    } catch (NumberFormatException nfe) {
                        throw parseError("invalid hex value for \\x escape", i);
                    }
                    break;
                }

                case 'u': {
                    if (i+4 >= subject.length())
                        throw parseError("string too short for \\u escape", i);
                    try {
                        i++;
                        sb.append(Character.toChars(Integer.parseInt( subject.substring(i, i+4), Constants.HEX_RADIX)));
                        i+=3;
                    } catch (NumberFormatException nfe) {
                        throw parseError("invalid hex value for \\u escape", i);
                    }
                    break;
                }

                case 'U': {
                    if (i+8 >= subject.length())
                        throw parseError("string too short for \\U escape", i);
                    try {
                        i++;
                        sb.append(Character.toChars(Integer.parseInt( subject.substring(i, i+8), Constants.HEX_RADIX)));
                        i+=7;
                    } catch (NumberFormatException nfe) {
                        throw parseError("invalid hex value for \\U escape", i);
                    }
                    break;
                }

                default:
                    sb.append('\\').append(Character.toChars(cp));
                    break;
            }
            prev_backslash = false;
        }

        /* final character was a backslash */
        if (prev_backslash)
            sb.append('\\');

        return sb.toString();
    }

//    /*
//     * Return a string "U+XX.XXX.XXXX" etc, where each XX set is the
//     * digits of the logical Unicode code point.
//     */
//    public static String uniplus(String s) {
//        if (s.length() == 0)
//            return "";
//        /* This is just the minimum; sb will grow as needed. */
//        StringBuilder sb = new StringBuilder(2 + 3 * s.length());
//        sb.append("U+");
//        for (int i = 0; i < s.length(); i++) {
//            sb.append(String.format("%X", s.codePointAt(i)));
//            if (s.codePointAt(i) > Character.MAX_VALUE)
//                i++; /* UTF-16 */
//            if (i+1 < s.length()) {
//                sb.append(".");
//            }
//        }
//        return sb.toString();
//    }

    private static ParseException parseError(String message, int pos) {
        return new ParseException("Error parsing string: " + message + " at char " + pos + ".", pos);
    }

    /**
     * Returns true if the string is null or empty ("").
     * @param str
     * @return
     */
    public static boolean isEmpty(final String str) {
        return str == null || str.isEmpty();
    }

    /**
     * Returns true if the string is null, empty or is composed entirely of
     * whitespace (as determined by Character.is
     * (e.g. space, tab, newlines, etc).
     * @param str
     * @return
     */
    public static boolean isWhiteSpace(@Nullable final String str) {
        if (str == null) return true;
        int length = str.length();
        for (int i = 0; i < length; i++) {
            if (!Character.isWhitespace(str.charAt(i)))
                return false;
        }
        return true;
    }
    /**
    * Attempts to convert a string representing a conceptual "true/false" into a boolean.
    * @return  {@code defaultIfEmpty } if {@code boolString} is null, empty (""), or whitespace.
    *           Returns FALSE if {@code boolString} == "false", "no" or "0" (case insensitive, trimmed).
    *           True otherwise.
    */
    public static boolean asBool(String boolString, boolean defaultIfEmpty) {
        if (boolString == null) return defaultIfEmpty;
        boolString = boolString.trim().toLowerCase();
        if ("".equals(boolString)) return defaultIfEmpty; //check length again AFTER trim()
        switch (boolString) {
            case "false":
            case "0":
            case "no":
                return false;
        }
        return true; //anything else, return true.
    }
    public static boolean asBool(String boolString) {
        return asBool(boolString, false);
    }
    /**
     * Returns the first string in the list of parameters that is not null or empty ("")
     * Returns empty ("") if all strings are either null or empty.
     */
    @NotNull
    public static String firstNonEmpty(final String... names) {
        for (String name : names) {
            if (!isEmpty(name))
                return name;
        }
        return "";
    }
    public static String escapeHTML(String s) {
        StringBuilder out = new StringBuilder(Math.max(16, s.length()));
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                //case '"': out.append("&quot;"); break;
                case '<': out.append("&lt;"); break;
                case '>': out.append("&gt;"); break;
                case '&': out.append("&amp;"); break;
                default:
//                    if (c >= 0x7F) {
//                        out.append("&#");
//                        out.append((int) c);
//                        out.append(';');
//                    } else
                        out.append(c);
            }
        }
        return out.toString();
    }
    /**
     * Converts a camel-case name into space separated words.
     *
     *        firstName yields "First Name"
     *        getID yields "Get ID"
     *        SSNNumber yields "SSN Number"
     *        iAmNice   yields "I Am Nice"
     *        _myName yields "My Name"
     *        name1 yields "Name 1"
     *        first_name yields First Name
     */
    public static String toFriendlyName(final String name) {
        if (isEmpty(name)) return "";
        StringBuilder sb = new StringBuilder();
        char[] chars = name.toCharArray();
        boolean prevUpper = false, prevDigit = false, prevSpace = true, doubleUpper = false;

        for (int i = 0; i < chars.length; i++) {
            char c = chars[i];

            if ((c == '_' || c == ' ')) {
                if (!prevSpace) {
                    prevSpace = true;
                    doubleUpper = prevUpper = prevDigit = false;
                    sb.append(' ');
                }
                continue;
            }

            boolean isDigit = Character.isDigit(c);
            boolean isUpper = Character.isUpperCase(c);
            if (prevSpace) {
                c = Character.toUpperCase(c); //add the first letter after a space as upper-case.
                prevSpace = false;
            } else {
                if (isUpper && !prevUpper || isDigit && !prevDigit) { // the current character is upper case or digit and the previous character was not. So insert a space before this character.
                    sb.append(' ');
                } else if (prevDigit && !isDigit) {
                    sb.append(' ');
                    c = Character.toUpperCase(c); //add the first letter after a number as upper-case.
                } else if (doubleUpper && !isUpper) {
                    // Insert a space before the previous character.
                    // e.g.   SSNNumber yields SSN Number
                    int prev = sb.length() - 1;
                    sb.append(sb.charAt(prev));
                    sb.setCharAt(prev, ' ');
                }
            }
            sb.append(c);
            doubleUpper = prevUpper && isUpper; // Both this character and the previous were upper-case, so do NOT insert a space (unless a lower-case character follows this)
            prevUpper = isUpper;
            prevDigit = isDigit;
        }
        int last = sb.length() - 1;
        if (last != -1 && sb.charAt(last) == ' ') // remove trailing space if it exists.
            sb.setLength(last);
        return sb.toString();
    }
    public static String replaceNonPrintable(final String s, final char replaceControlChars, final char replaceWhiteSpace) {
        int len = s.length();
        char c;
        boolean doReplace = false;
        for (int i = 0; i < len; i++) {
            c = s.charAt(i);
            if (Character.isHighSurrogate(c) && i < len - 1)
                // This char is part of a surrogate pair (two subsequent chars that together form a single unicode character.)
                i++; // skip this and the next character
            else if (c < 33 || c == 127) {
                doReplace = true;
                break;
            }
        }
        if (!doReplace) return s;
        char[] ss = s.toCharArray();
        for (int i = 0; i < len; i++) {
            c = ss[i];
            if (Character.isHighSurrogate(c) && i < len - 1)
                // This char is part of a surrogate pair (two subsequent chars that together form a single unicode character.)
                i++; // skip this and the next character
            else {
                switch (c) {
                    case '\t': case '\r': case '\n': case ' ':
                        ss[i] = replaceWhiteSpace;
                        break;
                    default:
                        if (c < 32 || c == 127)
                            ss[i] = replaceControlChars;
                        break;
                }
            }
        }
        return new String(ss);
    }

    /** Same as {@link #replaceNonPrintable } but includes support for unicode characters, including Surrogate pairs. */
    public static String replaceNonPrintableU(final String s, final char replaceControlChars, final char replaceWhiteSpace) {
        int len = s.length();
        char c;
        int code;
        boolean doReplace = false;
        for (int i = 0; i < len; i++) {
            c = s.charAt(i);
            if (Character.isHighSurrogate(c) && i < len - 1)
                // This char is part of a surrogate pair (two subsequent chars that together form a single unicode character.)
                code = Character.toCodePoint(c, s.charAt(i++));
            else
                code = c;
            if (Character.isWhitespace(code) || Character.isISOControl(code)) {
                doReplace = true;
                break;
            }
        }
        if (!doReplace) return s;
        char[] ss = s.toCharArray();
        for (int i = 0; i < len; i++) {
            c = ss[i];
            if (Character.isHighSurrogate(c) && i < len - 1)
                // This char is part of a surrogate pair (two subsequent chars that together form a single unicode character.)
                code = Character.toCodePoint(c, ss[i++]);
            else
                code = c;
            if (Character.isWhitespace(code))
                ss[i] = replaceWhiteSpace;
            else if (Character.isISOControl(code))
                ss[i] = replaceControlChars;
        }
        return new String(ss);
    }

//    public static String padr(final int targetLength, final String s) { return pad(true, targetLength, s, ' ', 0); }
//    public static String padr(final int targetLength, final String s, final char padWith, final int truncate) { return pad(true, targetLength, s, padWith, truncate);    }
//    public static String padl(final int targetLength, final String s) { return pad(false, targetLength, s, ' ', 0); }
//    public static String padl(final int targetLength, final String s, final char padWith, final int truncate) { return pad(false, targetLength, s, padWith, truncate); }
    /**
     * Pad a string with the specified character.
     * The padding can be placed on the right or left side of the string, depending on the sign of the width parameter:
     * A negative width places the padding on the left, while a positive width pads to the right side.
     * @param targetLength Determines both the side to pad on and the target number of characters to return.
     *              A negative value causes padding to be placed on the left side, while a positive one pads on the right.
     *              Padding is added until the string length is at least abs(width) long. (Where abs is the absolute value function).
     * @param subject The string to pad.
     * @param padWith The character to pad with.
     * @param truncate An int that specifies whether the string should be truncated if its current length is larger than the desired width.
     *                 0 indicates that no truncating should occur. If the string is longer than width, it is returned unmodified.
     *                 1 indicates that the string should be truncated by removing excess characters from the RIGHT side of the string.
     *                 -1 indicates that the string should be truncated by removing excess characters from the LEFT side of the string.
     * @return If truncate is 1 or -1, the length of the returned string will be exactly targetLength. If truncate = 0,
     * the string will have a minimum length of targetLength.
     */
    public static String pad(int targetLength, final String subject, final char padWith, final int truncate) {
        int len = subject.length();
        boolean padLeft = targetLength < 0;
        if (padLeft) targetLength = -targetLength;
        if (len > targetLength && truncate!=0)
            return truncate > 0 ? subject.substring(0, targetLength) :  subject.substring(subject.length()-targetLength); // <0: Truncate left, >0: Truncate right.
        else if (len < targetLength) {
            String pad = fromChar(padWith, targetLength - len);
            return padLeft ? pad+subject : subject+pad;
        }
        return subject;
    }
    public static String pad(int targetLength, final String subject) { return pad(targetLength, subject, ' ', 0); }
    public static String pad0s(int targetLength, final String subject) { return pad(targetLength, subject, '0', 0); }
    /** Pad-Right. Synonym for pad, but negative lengths are not allowed. */
    public static String padr(int targetLength, final String subject) {
        if (targetLength<0)
            throw new IllegalArgumentException("Length cannot be negative.");
        return pad(targetLength, subject, ' ', 0);
    }
    /** Pad-Left. Synonym for pad(-targetLength), but negative lengths are not allowed. */
    public static String padl(int targetLength, final String subject) {         if (targetLength<0)
        throw new IllegalArgumentException("Length cannot be negative.");
        return pad(-targetLength, subject, ' ', 0);
    }
    /**
     * Creates a string of the specified length, composed by repeating the specified character.
     * @param fillWith The character used to fill the string.
     * @param length The desired length of the string.
     * @return a string of the specified length, composed by repeating the specified character.
     */
    public static String fromChar(final char fillWith, final int length) {
        if (length==0) return "";
        char[] buf = new char[length];
        for (int i = 0; i < length; i++) buf[i] = fillWith;
        return new String(buf);
    }

    public static <T> String join(final String delimiter, final Iterable<T> values, @Nullable Function<? super T, String> formatter) {
        StringBuilder sb = new StringBuilder();
        for(T obj : values)
            sb.append(formatter==null?ObjTools.toStr(obj):formatter.apply(obj)).append(delimiter);
        if (sb.length() >= delimiter.length())
            sb.setLength(sb.length() - delimiter.length());
        return sb.toString();
    }
    public static <T> String join(final String delim, final Iterable<T> values) {
        return join(delim, values, null);
    }
    public static <T> String join(final String delim, final T[] values) {
        return join(delim, Arrays.asList(values), null);
    }

    /** Returns a 32-character hexadecimal-formatted MD5 hash of the subject string. */
    public static String md5(final String input) {
        try {
            java.security.MessageDigest md5 = java.security.MessageDigest.getInstance("MD5");
            byte[] hash = md5.digest(input.getBytes("UTF-8"));
            StringBuilder sb = new StringBuilder(32);
            for (byte b : hash) sb.append(String.format("%02x", b));
            return sb.toString();
        } catch (NoSuchAlgorithmException | UnsupportedEncodingException e) {
            e.printStackTrace();
            return null;
        }
    }
    public static final int DefaultBufferSize = 4 * 1024;
    public static String readAll(final InputStream input, Charset charset) throws IOException{ return readAll(input, charset, DefaultBufferSize); }
    public static String readAll(final InputStream input, String charset, int bufferSize) throws IOException{ return readAll(input, Charset.forName(charset), bufferSize); }
    public static String readAll(final InputStream input, String charset) throws IOException{ return readAll(input, charset, DefaultBufferSize); }
    public static String readAll(final InputStream input) throws IOException{ return readAll(input, Charset.defaultCharset(), DefaultBufferSize); }
    public static String readAll(final InputStream input, Charset charset, int bufferSize) throws IOException {
        ByteArrayOutputStream result = new ByteArrayOutputStream();
        byte[] buffer = new byte[bufferSize];
        int length;
        while ((length = input.read(buffer)) != -1) {
            result.write(buffer, 0, length);
        }
        return result.toString(charset.name());
    }
    public static void writeAll(final String content, final OutputStream output, Charset charset) throws IOException{ writeAll(content, output, charset, DefaultBufferSize); }
    public static void writeAll(final String content, final OutputStream output, String charset, int bufferSize) throws IOException{ writeAll(content, output, Charset.forName(charset), bufferSize); }
    public static void writeAll(final String content, final OutputStream output, String charset) throws IOException{ writeAll(content, output, charset, DefaultBufferSize); }
    public static void writeAll(final String content, final OutputStream output) throws IOException{ writeAll(content, output, Charset.defaultCharset(), DefaultBufferSize); }
    public static void writeAll(final String content, final OutputStream output, Charset charset, int bufferSize) throws IOException {
        try(BufferedOutputStream stream = new BufferedOutputStream(output, bufferSize)) {
            stream.write(content.getBytes(charset));
            stream.flush();
        }
    }
}
