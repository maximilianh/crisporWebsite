package ur_rna.Utilities;

import ur_rna.Utilities.annotation.Nullable;

import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;
import java.text.ParseException;
import java.util.*;
import java.util.regex.MatchResult;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * SimpleDataSerializer
 * Serializes Lists, Maps and primitive types (Strings, Numbers, Booleans) to
 * a simple text format.
 * Capable of reading JSON as well as SSON, a simplified form of JSON defined here.
 *
 * SSON (Simple Structured Object Notation) is less strict than JSON. SSON allows the following:
 *   - In Maps, quotes are not required for property names (as long as they conform to the rules of unquoted strings described later).
 *   - Commas are optional between elements of Lists and Maps
 *   - C-style Single-line and multi-line comments are allowed.
 *   - Strings do not require quotes if they start with a letter and contain only letters, numbers, and the symbols underscore (_) dot (.) and hyphen (-)
 *   - Arrays can have an optional primitive type constraints:
 *           e.g. [%i 1 2 3 4 ]  [%d 1.1 0.2 3.7 4 ]  [%s "a" "b" "c" ]
 *           i=int L=long f=float d=double k=byte b=bool s=string (default)
 * ------ Not yet implemented: -----------
 *   - Addition of Tables, which represent an array of same-structured objects in a compact form.
 *           [%% "A"   "B"   "C"  |
 *               O1.A  O1.B  O1.C |
 *               O2.A  O2.B  O2.C |
 *               O3.A  O3.B  O3.C ]
 *     Instead of:
 *           [ { A:O1.A, B:O1.B, C:O1.C}
 *             { A:O2.A, B:O2.B, C:O2.C}
 *             { A:O3.A, B:O3.B, C:O3.C} ]
 */
public class SimpleDataSerializer {
    /**
     * Parses text into an SONMap
     * characterSourceIndex is an array that provides information about the
     * original line and column of a particular character index in the final input code.
     * The format is [ index1, line1, col1, index2, line2, col2, ... ] such that
     * if an error is found at index=531, the characterSourceIndex can be searched until
     * one index_a &lt; 531 &lt; index_b. Then the original line is line_a and the original
     * column is (531 - index_a + column_a).
     *
     */
    public Object parse(CharSequence input, @Nullable int[] characterSourceIndex) throws SyntaxErrorException {
        if (characterSourceIndex==null)
            characterSourceIndex = buildSourceIndex(input);
        Token[] tokens = tokenize(input, characterSourceIndex, true, true);
        ParseContext ctx = new ParseContext(tokens, 0, characterSourceIndex);
        Object value = parseSONValue(ctx);
        if (ctx.hasNext()) {
            Token next = ctx.next();
            throw syntaxErr("unexpected " + next.type().name() + " (expected END-OF-DOCUMENT)", next.start(), ctx.sourceIndex, next.text());
        }
        return value;
    }
    private int[] buildSourceIndex(final CharSequence input) {
        int len = input.length();
        int[] index = new int[30];
        int lineCount = 0;
        for (int i = 0; i <len; i++) {
            if (input.charAt(i)=='\n') {
                lineCount++;
                if (index.length < (lineCount+1) * 3)
                    index = Arrays.copyOf(index, index.length * 2);
                index[lineCount*3] = i+1;
                index[lineCount*3+1] = lineCount;
            }
        }
        return index;
    }

    public static class OutputFormatSettings  implements Cloneable {
        public final static OutputFormatSettings SSON;
        public final static OutputFormatSettings PrettySSON;
        public final static OutputFormatSettings JSON;
        public final static OutputFormatSettings PrettyJSON;
         //public static OutputFormatSettings getDefault() {            return Default.clone();        }

        static {
            SSON = new OutputFormatSettings(0, false, false, true, true, false, false, false);
            JSON = new OutputFormatSettings(0, true, true, false, false, true, true, true);
            PrettySSON = SSON.clone(); PrettySSON.indent = 2;
            PrettyJSON = JSON.clone(); PrettyJSON.indent = 2;
        }
        public OutputFormatSettings(){}
        public OutputFormatSettings(int indent, boolean strict){
            this(indent, strict,strict, !strict, !strict, strict, strict, strict);
        }
        public OutputFormatSettings(final int indent, final boolean quoteAlways, final boolean quoteProps, final boolean inlineLists, final boolean inlineMaps, final boolean commasInLists, final boolean commasInMaps, final boolean spaceAfterProp) {
            this.indent = indent;
            this.quoteAlways = quoteAlways;
            this.quoteProps = quoteProps;
            this.inlineLists = inlineLists;
            this.inlineMaps = inlineMaps;
            this.commasInLists = commasInLists;
            this.commasInMaps = commasInMaps;
            this.spaceAfterProp = spaceAfterProp;
        }
        public int indent;
        public boolean quoteAlways;
        public boolean quoteProps;
        public boolean inlineLists;
        public boolean inlineMaps;
        //public Func numbers;
        public boolean commasInLists;
        public boolean commasInMaps;
        public boolean spaceAfterProp;

        public OutputFormatSettings clone() {
            try {
                return (OutputFormatSettings)super.clone();
            } catch (CloneNotSupportedException ex) {
                throw new  InternalError(ex);
            }
        }
        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            OutputFormatSettings settings = (OutputFormatSettings) o;

            if (indent != settings.indent) return false;
            if (quoteAlways != settings.quoteAlways) return false;
            if (quoteProps != settings.quoteProps) return false;
            if (inlineLists != settings.inlineLists) return false;
            if (inlineMaps != settings.inlineMaps) return false;
            if (commasInLists != settings.commasInLists) return false;
            if (commasInMaps != settings.commasInMaps) return false;
            return spaceAfterProp == settings.spaceAfterProp;
        }
        @Override
        public int hashCode() {
            int result = indent;
            result = 31 * result + (quoteAlways ? 1 : 0);
            result = 31 * result + (quoteProps ? 1 : 0);
            result = 31 * result + (inlineLists ? 1 : 0);
            result = 31 * result + (inlineMaps ? 1 : 0);
            result = 31 * result + (commasInLists ? 1 : 0);
            result = 31 * result + (commasInMaps ? 1 : 0);
            result = 31 * result + (spaceAfterProp ? 1 : 0);
            return result;
        }
    }

    public String toString(Object value, OutputFormatSettings settings) throws FormatterException {
        try {
            StringWriter w = new StringWriter();
            this.write(value, w, settings);
            return w.toString();
        } catch (IOException ex) {
            throw new FormatterException("Unexpected IO Error", ex);
        }
    }

    public void write(Object value, Writer out, OutputFormatSettings settings) throws IOException, FormatterException {
        Queue<String> trace = new ArrayDeque<>();
        if (settings == null) settings = OutputFormatSettings.SSON;
        write(value, out, settings, trace);
    }
    private void write(Object value, Writer out, OutputFormatSettings settings, Queue<String> trace) throws IOException, FormatterException {
        if (value == null)
            out.write("null");
        else if (value instanceof String)
            out.write(writeString((String)value, settings.quoteAlways));
        else if (value instanceof Boolean)
            out.write(((Boolean)value) ? "true" : "false");
        else if (value instanceof Number)
            out.write(writeNumber((Number)value, settings));
        else if (value instanceof Character)
            out.write(writeString(value.toString(), settings.quoteAlways));
        else if (value instanceof List<?>)
            writeList((List<?>) value, out, settings, trace);
        else if (value instanceof Map<?,?>)
            writeMap((Map<?,?>) value, out, settings, trace);
        else
            throw new FormatterException(String.format("Cannot serialize a value of type %s.\nAt: %s", value.getClass().getSimpleName(), formatTrace(trace)));

    }
    private void writeMap(final Map<?, ?> value, final Writer out, final OutputFormatSettings settings, final Queue<String> trace)  throws IOException, FormatterException {
        int len = value.size();

        if (len == 0) {
            out.write("{ }");
            return;
        }

        boolean indent = false;
        int indentLevel = trace.size() + 1;
        if (settings.indent != 0) {
            if (settings.inlineMaps) {
                for (Object val : value.values()) {
                    if (isListOrMap(val)) {
                        indent = true;
                        break;
                    }
                }
            } else
                indent = true;
        }
        out.write('{');
        int remaining = len;
        for (Map.Entry<?,?> kv : value.entrySet()) {
            if (indent) {
                out.write(System.lineSeparator());
                writeIndent(out, settings, indentLevel);
            } else
                out.write(' ');

            String key = kv.getKey().toString();
            out.write(writeString(key, settings.quoteProps));
            out.write(':');
            if (settings.spaceAfterProp) out.write(' ');
            trace.add("map:"+key);
            write(kv.getValue(), out, settings, trace);
            trace.remove();
            if (settings.commasInMaps && --remaining > 0)
                out.write(',');
        }
        if (indent) {
            out.write(System.lineSeparator());
            writeIndent(out, settings, indentLevel - 1);
        } else
            out.write(' ');
        out.write('}');
    }
    private void writeList(final List<?> value, final Writer out, final OutputFormatSettings settings, final Queue<String> trace) throws IOException, FormatterException  {
        int len = value.size();

        if (len == 0) {
            out.write("[ ]");
            return;
        }

        boolean indent = false;
        int indentLevel = trace.size() + 1;
        if (settings.indent != 0) {
            if (settings.inlineLists) {
                for (int i = 0; i < len; i++) {
                    if (isListOrMap(value.get(i))) {
                        indent = true;
                        break;
                    }
                }
            } else
                indent = true;
        }

        out.write('[');
        for (int i = 0; i < len; i++) {
            if (indent) {
                out.write(System.lineSeparator());
                writeIndent(out, settings, indentLevel);
            } else
                out.write(' ');

            trace.add("list#"+i);
            write(value.get(i), out, settings, trace);
            trace.remove();
            if (settings.commasInLists && i < len - 1)
                out.write(',');
        }
        if (indent) {
            out.write(System.lineSeparator());
            writeIndent(out, settings, indentLevel - 1);
        } else
            out.write(' ');
        out.write(']');
    }

    private void writeIndent(Writer out, OutputFormatSettings settings, int level) throws IOException {
        int count = settings.indent * level;
        for (int j = 0; j < count; j++)
            out.write(' ');
    }

    private boolean isListOrMap(final Object o) {
        return o instanceof List<?> || o instanceof Map<?,?>;
    }
    private String formatTrace(final Queue<String> trace) {
        if (trace.size() == 0) return "/";
        if (trace.size() == 1) return trace.poll();
        StringBuilder sb = new StringBuilder(trace.size() * 5); // rough estimate of size. each item is "list", "map", an index (0 to 99 in most cases), or a property-name (which can be any length, but is limited to 10 characters). plus "/"  for each item.
        for (String s : trace)
            sb.append('/').append(Strings.replaceNonPrintable(s, '?', '_'));
        return sb.toString();
    }

    private String writeString(String s, boolean forceQuotes) {
        if (forceQuotes)
            return '"' +  Strings.escapeStringLiteral(s, STRING_ESCAPE_OPTIONS) + '"';
        return quoteIfNeeded(s);
    }
    private String writeNumber(Number n, OutputFormatSettings settings) {
        //return (settings.numbers == null ? OutputFormatSettings.Default : settings).numbers.format(n);
        return n.toString();
    }


    Token[] tokenize(final CharSequence input, final int[] index, boolean dropWhitespace, boolean dropComments) throws SyntaxErrorException {
        ArrayList<Token> tokens = new ArrayList<>();
        TokenPattern[] types = TokenPattern.values();
        int pos = 0, end = input.length();
        outer:
        while (pos < end) {
            for (TokenPattern type : types) {
                Token t = type.match(input, pos);
                if (t != null) {
                    pos = t.end();
                    boolean skip = dropWhitespace && type == TokenPattern.WHITESPACE
                            || dropComments && type == TokenPattern.COMMENT;
                    if (!skip)
                        tokens.add(t);
                    continue outer;
                }
            }
            // the character data did not match anything
            throw syntaxErr("syntax error", pos, index, input.subSequence(pos, Math.min(pos + 8, end)) + "...");
        }
        return tokens.toArray(new Token[tokens.size()]);
    }

    SyntaxErrorException syntaxErr(String message, int pos, int[] sourceIndex, @Nullable String offendingText) {
        int line = 0, col = pos;
        if (sourceIndex != null && sourceIndex.length >= 3) {
            int index = sourceIndex.length - 3;
            for (int i = 3; i < sourceIndex.length; i += 3) {
                if (sourceIndex[i] > pos) {
                    index = i - 3;
                    break;
                }
            }
            line = sourceIndex[index + 1];
            col = pos - sourceIndex[index] + sourceIndex[index + 2];
        }
        return new SyntaxErrorException("Error parsing object tree: " + message, offendingText, line+1, col+1);
    }
    private transient static Matcher identifierMatcher;
    private static final int STRING_ESCAPE_OPTIONS = CharType.BACKSLASH | CharType.DOUBLE_QUOTES;
    public static String quoteIfNeeded(final String s) {
        if (identifierMatcher == null)
            identifierMatcher = TokenPattern.IDENTIFIER.pattern().matcher(s);
        else
            identifierMatcher.reset(s);
        if (identifierMatcher.matches())
            return s;
        return '"' +  Strings.escapeStringLiteral(s, STRING_ESCAPE_OPTIONS) + '"';
    }


    private static final String END_TOKEN = "(?=\\s|[:{}()\\[\\]\"]|$)";
    //private static Pattern COMMAND = Pattern.compile("\\s*(?<pfx>[!~]+)(?<name>[\\w-]+)");
    enum TokenPattern {
        // Parsing a number -- The regex "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?" is a good, narrow regex.
        // But then expressions like 5+6 or 6frog are tokenized into [5, +6] and [ 6, frog ] respectively,
        // when they should probably actually result in an error.
        // So the negative lookahead "(?![-+\w])" was added to prevent this. (note that \b doesn't work, because '+' would match)
        NUMBER(Pattern.compile("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"+END_TOKEN)),
        // A double-quote followed by one or more of:
        //    1) any character that is neither \ nor "
        //    2) \ followed by any of: b f n r t " \ or ' (e.g. \n, \0, \\)
        //    3) \xHH or \\uHHHH (where H is a hexadecimal number)
        //    4) \O \OO \OOO where the 'O's represent an octal number (O = 0-7)
        STRING(Pattern.compile("\"(?:[^\\\\\"]|\\\\(?:[abfnrt\"\\\\']|c[@-~]|x[A-Za-z0-9]{2}|u[A-Za-z0-9]{4}|[0-7]{1,3}))*+\"")),  // uses backslash as escape character. Valid escapes include \\ \r \n \t \f \b \a \"  \'  \cX \xHH \\uHHHH and \O \OO \OOO ('O'=0-7)
        BOOL(Pattern.compile("true|false", Pattern.CASE_INSENSITIVE)), // i.e. "(?:[^"\]|\"|\)+"
        NULL(Pattern.compile("null", Pattern.CASE_INSENSITIVE)), // i.e. "(?:[^"\]|\"|\)+"
        IDENTIFIER(Pattern.compile("\\p{L}[\\p{L}\\p{N}\\p{Pc}_.-]*"+END_TOKEN)), // i.e. "(?:[^"\]|\"|\)+"   // must start with a letter and include only letters (and unicode-connecting-punctuation), numbers, underscore, dot, and hyphen
        WHITESPACE(Pattern.compile("\\s+")),
        COMMENT(Pattern.compile("//[^\\n]*+|/\\*(?:[^*]|\\*(?!/))*+\\*/")), // both single-line and multi-line C-style comments are supported
        OPEN_BRACKET(Pattern.compile("\\[")), // begins an array or table
        CLOSE_BRACKET(Pattern.compile("\\]")), // ends an array or table
        OPEN_BRACE(Pattern.compile("\\{")), // begins an object
        CLOSE_BRACE(Pattern.compile("\\}")), // ends an object
        // START_P(Pattern.compile("\\(")),
        // END_P(Pattern.compile("\\)")),
        COMMA(Pattern.compile(",")),
        COLON(Pattern.compile(":")),     // separates property name from value in objects
        SYM_PERCENT(Pattern.compile("%")),   // indicates a type constraint
        SYM_BAR(Pattern.compile("\\|"));     // Used for separation of "rows" of values in tables.

        private final Pattern p;
        TokenPattern(Pattern p) {
            this.p = p;
        }
        public Pattern pattern() { return p;}
        public Token match(CharSequence s, int pos) { return match(s, pos, s.length()); }
        public Token match(CharSequence s, int start, int end) {
            Matcher m = p.matcher(s).region(start, end);
            if (m.lookingAt())
                return new Token(m, this);
            return null;
        }
    }

    static class Token {
        public TokenPattern type() { return p; }
        public int start() { return m.start(); }
        public int end() { return m.end(); }
        public int length() { return m.end() - m.start(); }
        public String text() { return m.group(0); }
        public String group(int index) { return m.group(index); }
        public String group(String name) { return m.group(name); }
        public MatchResult result() { return m; }
        private final Matcher m;
        private final TokenPattern p;
        public Token(final Matcher m, final TokenPattern p) {
            this.m = m;
            this.p = p;
        }
        @Override
        public String toString() {
            return "{" + this.getClass().getSimpleName() + ": " + type().name() + " \"" + Strings.escapeStringLiteral(text()) + "\"}";
        }
    }
//    private static class Node extends TreeNode<Node> {
//        public final int startPos;
//        public final int endPos;
//        public final Expr expr;
//        public final Token token;
//        public Node(final int startPos, final int endPos, final Expr expr, final Matcher matcher) {
//            this.startPos = startPos;
//            this.endPos = endPos;
//            this.expr = expr;
//            this.matcher = matcher;
//        }
//    }

    private class ParseContext {
        public int pos;
        public int[] sourceIndex;
        public Token[] tokens;
//        public ParseContext copy() {
//            return new ParseContext(tokens, pos, sourceIndex);
//        }
        public ParseContext(final Token[] tokens, final int pos, final int[] sourceIndex) {
            this.pos = pos;
            this.sourceIndex = sourceIndex;
            this.tokens = tokens;
        }
        public Token next() throws SyntaxErrorException {
            if (this.pos >= tokens.length)
                throw syntaxErr("unexpected end of document", last() == null ? 0 : last().end(), sourceIndex, null);
            return tokens[this.pos++];
        }
        public boolean hasNext() {
            return this.pos < tokens.length;
        }
        public Token peek() {
            return tokens[this.pos];
        }
        public Token last() {
            return tokens.length == 0 ? null : tokens[tokens.length - 1];
        }
    }
    private Object parseSONValue(ParseContext ctx) throws SyntaxErrorException {
        Token t = ctx.next();
        switch (t.type()) {
            case BOOL:
                return Boolean.valueOf(t.text());
            case NUMBER:
                return parseNumber(t, ctx);
            case IDENTIFIER:
            case STRING:
                return parseString(t, ctx);
            case NULL:
                return null;
            case OPEN_BRACKET:
                return parseArray(ctx);
            case OPEN_BRACE:
                return parseObject(ctx);
            default:
                throw syntaxErr("unexpected " + t.type().name(), t.start(), ctx.sourceIndex, t.text());
        }
    }

    private Object parseArray(ParseContext ctx) throws SyntaxErrorException {
        List<Object> arr = new SSONList();
        Token openBracket = ctx.tokens[ctx.pos - 1];
        while (ctx.hasNext()) {
            Token t = ctx.peek();
            if (t.type() == TokenPattern.CLOSE_BRACKET) {
                // we found the end-bracket.
                ctx.next(); //consume the token
                return arr;
            } else
                arr.add(parseSONValue(ctx));
            if (ctx.hasNext() && ctx.peek().type() == TokenPattern.COMMA)
                ctx.next(); //consume the separator if one exists.
        }
        throw syntaxErr("array is missing its end-bracket", openBracket.start(), ctx.sourceIndex, null);
    }

    private Object parseObject(ParseContext ctx) throws SyntaxErrorException {
        Map<String,Object> obj = new SSONMap();
        Token openBrace = ctx.tokens[ctx.pos - 1];
        boolean started = false;
        while (ctx.hasNext()) {
            Token t = ctx.next();
            if (t.type() == TokenPattern.CLOSE_BRACE) {
                // we found the end-brace.
                return obj;
            } else {
                if (started) {
                    if (t.type() == TokenPattern.COMMA)
                        t = ctx.next();
                    //else
                    //    throw syntaxErr("missing comma in property list", t.start(), ctx.sourceIndex, t.text());
                }
                String name = parseString(t, ctx);
                if (!ctx.hasNext() || ctx.next().type() != TokenPattern.COLON)
                    throw syntaxErr("expected a colon (:) after property " + t.text(), t.end(), ctx.sourceIndex, t.text());
                else if (!ctx.hasNext())
                    throw syntaxErr("unexpected end of document -- expected a value for property \"" + t.text(), t.start(), ctx.sourceIndex, t.text());
                obj.put(name, parseSONValue(ctx));
                started = true;
            }
        }
        throw syntaxErr("object is missing its end-brace", openBrace.start(), ctx.sourceIndex, null);
    }

    //private ParsePosition pp = new ParsePosition(0);
    private Number parseNumber(Token t, ParseContext ctx) throws SyntaxErrorException {
        try {
            return Double.parseDouble(t.text());
//            pp.setIndex(0); pp.setErrorIndex(-1);
//            Number n = NumberFormat.getInstance().parse(t.text(), pp);
//            if (pp.getErrorIndex() != -1)
//                throw new ParseException("", pp.getErrorIndex());
//            if (pp.getIndex() != t.length())
//                throw new ParseException("", pp.getIndex());
//            return n;
//        } catch (ParseException ex) {
//            throw syntaxErr("error parsing number", t.start() + ex.getErrorOffset(), ctx.sourceIndex, t.text());
//        }
        } catch (NumberFormatException ex) {
            throw syntaxErr("error parsing number (" + ex.getMessage() + ")", t.start(), ctx.sourceIndex, t.text());
        }
    }
    private String parseString(Token t, ParseContext ctx) throws SyntaxErrorException {
        if (t.type() == TokenPattern.STRING) {
            String s = t.text();
            try {
                return Strings.unescapeStringLiteral(s.substring(1, s.length() - 1));
            } catch (ParseException ex) {
                throw syntaxErr("error parsing escape sequence in string", t.start() + ex.getErrorOffset() + 1, ctx.sourceIndex, t.text());
            }
        }
        if (t.type() == TokenPattern.IDENTIFIER)
            return t.text();
        throw syntaxErr("unexpected " + t.type().name(), t.start(), ctx.sourceIndex, t.text());
    }

    public static class SSONList extends ArrayList<Object> {
        public SSONMap addMap() {
            SSONMap m = new SSONMap();
            add(m);
            return m;
        }
        public SSONList addList() {
            SSONList list = new SSONList();
            add(list);
            return list;
        }
    }
    public static class SSONMap extends LinkedHashMap<String,Object> {
        public int get(String name, int defaultValue) {  return Convert.toInt(get(name), defaultValue); }
        public double get(String name, double defaultValue) {  return Convert.toDouble(get(name), defaultValue); }
        public String get(String name, String defaultValue) { return Convert.toString(get(name), defaultValue); }
        public float get(String name, float defaultValue) {  return Convert.toFloat(get(name), defaultValue); }
        public boolean get(String name, boolean defaultValue) {  return Convert.toBool(get(name), defaultValue); }
        public SSONMap getMap(String name) {
            Object val = get(name);
            if (val instanceof SSONMap)
                return (SSONMap)val;
            throw new IllegalArgumentException(String.format("The value of property '%s' is of type %s, not SSONMap.", name, val == null ? "Null" : val.getClass().getSimpleName()));
        }
        public SSONList getList(String name) {
            Object val = get(name);
            if (val instanceof SSONList)
                return (SSONList)val;
            throw new IllegalArgumentException(String.format("The value of property '%s' is of type %s, not SSONList.", name, val == null ? "Null" : val.getClass().getSimpleName()));
        }

        public SSONList putList(final String key) {
            SSONList list = new SSONList();
            put(key, list);
            return list;
        }
        public SSONMap putMap(final String key) {
            SSONMap m = new SSONMap();
            put(key, m);
            return m;
        }
        public boolean containsAnyKey(final String...keys) {
            for (String s : keys)
                if (this.containsKey(s)) return true;
            return false;
        }
    }
}
