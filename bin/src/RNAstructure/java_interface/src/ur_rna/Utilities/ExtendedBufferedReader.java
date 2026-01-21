package ur_rna.Utilities;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

/**
 * A special buffered reader which tracks the current character position and line index.
 * <p>
 * The reader also contains a peek method, which allows you to see the next char returned by
 * </p>
 */
final class ExtendedBufferedReader extends BufferedReader {
    public static final int CR = 10;
    public static final int LF = 13;
    public static final int EOF = -1;
    public static final int UNKNOWN = -2;
    private static final int defaultGuessLineLength = 256;

    /**
     * The last char returned from read() (or an overload)
     * Usage: If lastReadChar is CR and a LF is encountered next, the currentLine should NOT be incremented.
     */
    private int lastReadChar = UNKNOWN;

    /** The current line index. This is equal to the number of EOLs (CR/LF/CRLF) read thus far */
    private long currentLine;

    /** The position, which is number of characters read so far */
    private long position;
    /** Whether or not this reader has been closed. */
    private boolean closed;

    /**
     * Created extended buffered reader using default buffer-size
     */
    ExtendedBufferedReader(final Reader reader) {
        super(reader);
    }

    @Override
    public int read() throws IOException {
        final int ichr = super.read();
        if (ichr == CR || (ichr == LF && lastReadChar != CR)) {
            currentLine++;
        }
        lastReadChar = ichr;
        this.position++;
        return lastReadChar;
    }

    /**
     * Returns the last character that was read as an integer (0 to 65535). This will be the last character returned by
     * any of the read methods. This will not include a character read using the {@link #peek()} method. If no
     * character has been read then this will return {@link #UNKNOWN}. If the end of the stream was reached
     * on the last read then this will return {@link #EOF}.
     *
     * @return the last character that was read
     */
    int getLastChar() {
        return lastReadChar;
    }

    @Override
    public int read(final char[] buf, final int offset, final int length) throws IOException {
        if (length == 0)
            return 0;

        final int len = super.read(buf, offset, length);

        if (len > 0) {
            for (int i = offset; i < offset + len; i++) {
                if (buf[i] == LF) {
                    // if i == 0, we are at the start of the buffer, and the character before this was lastReadChar.
                    // conversely, if i > 0, the preceding character is simply buf[i - 1]
                    if (CR != (i == 0 ? lastReadChar : buf[i - 1]))
                        currentLine++;
                } else if (buf[i] == CR)
                    currentLine++;
            }
            lastReadChar = buf[offset + len - 1]; //update to the last character in the buffer.
        } else if (len == -1) {
            lastReadChar = EOF;
        }

        position += len;
        return len;
    }

    /**
     * Does not call {@link BufferedReader#readLine()} because that method drops the line terminator(s).
     * <p>
     * Increments {@link #currentLine}
     * <p>
     * Sets {@link #lastReadChar} to {@link #EOF} at EOF, otherwise to LF
     *
     * @return the line that was read, or null if reached EOF.
     */
    @Override
    public String readLine() throws IOException {
        return readLine(true, defaultGuessLineLength);
    }
    public String readLine(boolean dropEolChars) throws IOException {
        return readLine(dropEolChars, defaultGuessLineLength);
    }
    public String readLine(final boolean dropEolChars, final int guessLineLength) throws IOException {
        //char[] buf = new char[guessLineLength == 0 ? defaultGuessLineLength : guessLineLength];
        //int len;
//        while((len = super.read(buf)) > 0) {
//            for
//        }

        StringBuilder sb = new StringBuilder();
        boolean testNL = false;
        int i;
        while (0 <= (i = read())) {
            sb.append((char) i);
            if (i == LF)
                break;

            if (testNL) {
                sb.setLength(sb.length() - 1); //remove last character because it wasn't a newline.
                reset(); //move the reader back one character.
                break;
            }

            if (i == CR) {
                // the next character might be a LF.
                // let's let the loop continue one more time to see.
                // Save the current location in the buffered reader so we can return to it
                mark(1);
                testNL = true;
            }
        }
        return sb.length() == 0 ? null : sb.toString();

        //from extended:
//        if (line != null) {
//            lastReadChar = LF; // needed for detecting start of line
//            currentLine++;
//        } else {
//            lastReadChar = EOF;
//        }
//
//        return line;
    }

    /**
     * Returns the next character in the current reader without consuming it. So the next call to {@link #read()} will
     * still return this value. Does not affect line number or last character.
     *
     * @return the next character
     * @throws IOException if there is an error in reading
     */
    int peek() throws IOException {
        super.mark(1);
        final int c = super.read();
        super.reset();

        return c;
    }

    /**
     * Returns the current line number
     *
     * @return the current line number
     */
    long getCurrentLineNumber() {
        // Check if we are at EOL or EOF or just starting
        if (lastReadChar == CR || lastReadChar == LF || lastReadChar == UNKNOWN || lastReadChar == EOF) {
            return currentLine; // counter is accurate
        }
        return currentLine + 1; // Allow for counter being incremented only at EOL
    }

    /**
     * Gets the character position in the reader.
     *
     * @return the current position in the reader (counting characters, not bytes since this is a Reader)
     */
    long getPosition() {
        return this.position;
    }

    public boolean isClosed() {
        return closed;
    }

    /**
     * Closes the stream.
     *
     * @throws IOException If an I/O error occurs
     */
    @Override
    public void close() throws IOException {
        // Set ivars before calling super close() in case close() throws an IOException.
        closed = true;
        lastReadChar = EOF;
        super.close();
    }
}