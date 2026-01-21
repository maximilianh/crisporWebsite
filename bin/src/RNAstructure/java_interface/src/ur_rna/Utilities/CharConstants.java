package ur_rna.Utilities;

public final class CharConstants {
    private CharConstants(){} //prevent instantiation

    public static final char NUL = '\0';
    public static final char BKSP = '\b';
    public static final char TAB = '\t';
    public static final char LF =  '\n';
    public static final char FF =  '\f';
    public static final char CR =  '\r';

    public static final char SP = ' ';
    public static final char DEL = 127; // '\u007f'
    public static final char UPPER_ASCII = 128; // '\u00FF'
    public static final char LAST_ASCII = 255; // '\u00FF'

    public static final char SNG_QUOT =  '\"';
    public static final char DBL_QUOT =  '\'';
    public static final char BACK_SLASH =  '\\';
    public static final char SLASH =  '/';
    public static final char ESC =  27; // '\033'
    public static final char ALERT = '\007';
    public static final String CRLF = "\r\n";

    public static final int HEX_RADIX=16;
}
