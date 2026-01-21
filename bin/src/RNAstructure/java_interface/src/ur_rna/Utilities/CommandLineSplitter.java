package ur_rna.Utilities;

import java.util.ArrayList;
import java.util.List;

/**
 * Split a string into an array in a manner similar to how a console command line is split.
 *
 */
public class CommandLineSplitter {
    public static String[] split(String text) {
        StringBuilder word = new StringBuilder();
        List<String> args = new ArrayList<>();
        char[] chars = text.toCharArray();
        // State vars
        boolean escapeNext = false;
        boolean inQuote = false;
        boolean hardQuotes = false; // if true, escapes will be ignored.
        char endQuote = 0;
        boolean forceAddEmpty = false; // set to true if the word buffer is empty,
                                       // but an argument should be added anyway because
                                       // a quoted string was present.

        // Loop over each char in the string
        for (int i = 0; i < chars.length; i++) {
            char c = chars[i];
            if (escapeNext) {
                // append un-escaped character to the current word.
                escapeNext = false;
                switch (c) {
                    case 'n': c = '\n'; break;
                    case 'r': c = '\r'; break;
                    case 't': c = '\t'; break;
                    case '0': c = 0; break;
                    case 'f': c = '\f'; break;
                    case 'v': c = 11; break;
                }
                word.append(c);
            } else if (inQuote) {
                // We are currently in a quoted string.
                if (c == endQuote) { // Test for the end-quote
                    inQuote = false;
                    forceAddEmpty = word.length() == 0;
                }
                else if (c == '\\' && !hardQuotes) // respect the escape char unless hardQuotes is on.
                    escapeNext = true;
                else
                    word.append(c);
            } else {
                switch(c) {
                    case '\\':
                        escapeNext = true;
                        break;
                    case '\'':
                    case '"':
                        inQuote = true;
                        endQuote = c;
                        hardQuotes = c == '\'';
                        break;
                    default:
                        if (Character.isWhitespace(c)) {
                            // Whitespace ends the current token
                            if (word.length() != 0 || forceAddEmpty)
                                args.add(word.toString());
                            word.setLength(0);
                            forceAddEmpty = false;
                        } else
                            word.append(c);
                }
            }
        }

        // If the last char was a backslash, add it.
        if (escapeNext)
            word.append('\\');

        // Add the last argument (unless already added by presence of whitespace)
        if (word.length() != 0 || forceAddEmpty)
            args.add(word.toString());

        return args.toArray(new String[args.size()]);
    }
}
