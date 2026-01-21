package ur_rna.Utilities;

public final class CharType {
    private CharType() {} //disallow instantiation.

    public final static int NONE = 0;

    /** All character types. */
    public final static int ALL = ~NONE;

    /** ASCII Control codes (#0 to #31 and #127, but NOT including the white space charcters TAB=#9, LF=#10, and CR=#13). */
    public final static int CONTROL = 1; // bit-0
    //reserved bits 1 - 3

    /** The tab character (#9) */
    public final static int TAB = 1<<4;
    /** CR and LF */
    public final static int NEWLINE = 1<<5;

    /** The space character (#32) */
    //reserved bits 6,7
    public final static int SPACE = 1<<8;
    //reserved bit 9

    /** A-Z */
    public final static int UPPER_ALPHA = 1<<10;
    /** '0'-'9' */
    public final static int DIGITS = 1<<11;
    /** a-z */
    public final static int LOWER_ALPHA = 1<<12;
    //reserved bits 13 - 15

    /** Single quotes (') */
    public final static int SINGLE_QUOTE = 1<<16;
    /** Double quotes (") */
    public final static int DOUBLE_QUOTES = 1<<17;
    public final static int SLASH = 1<<18;
    public final static int BACKSLASH = 1<<19;
    public final static int DOT = 1<<20;
    public final static int COMMA = 1<<21;
    public final static int UNDERSCORE = 1<<22;
    //reserved bits 23 - 26
    /** Lower-ASCII Symbols not accounted for by any other flag. */
    public final static int MISC_SYMBOLS = 1<<27;
    //reserved bit 28

    /** Extended ASCII characters (129-255) */
    public final static int EXT_ASCII = 1<<29;
    //reserved bits 30, 31

    /** Whitespace characters: TAB, CR, LF, and SPACE */
    public final static int WHITESPACE = TAB | NEWLINE | SPACE;
    public final static int ALPHA = LOWER_ALPHA | UPPER_ALPHA;
    public final static int ALPHA_NUM = LOWER_ALPHA | UPPER_ALPHA | DIGITS;
    public final static int DECIMAL_NUMBER = DIGITS | DOT | COMMA;
    public final static int WORD_CHARS = LOWER_ALPHA | UPPER_ALPHA | UNDERSCORE;
    public final static int PUNCTUATION = DOT | COMMA;

    /** Single (') and double (") quotes. */
    public final static int QUOTES = SINGLE_QUOTE | DOUBLE_QUOTES;

    /** All Lower-ASCII symbols. */
    public final static int SYMBOLS =  0x0FFF<<16;  // bit 31--> 0000 1111 1111 1111 0000 0000 0000 0000 <--bit 0

    /** Non-printable characters, including control characters, TAB, CR, and LF. */
    public final static int NON_PRINTABLE = CONTROL | TAB | NEWLINE;
}
