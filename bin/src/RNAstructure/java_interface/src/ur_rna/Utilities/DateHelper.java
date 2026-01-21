package ur_rna.Utilities;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
/**
 * @author Richard M. Watson
 */
public class DateHelper {
    public final static String ISO8601Format = "yyyy-MM-dd'T'HH:mm:ss";
    public final static String FileDateFormat = "yyyy-MM-dd";
    public final static String FileDateTimeFormat = "yyyy-MM-dd.HHmmss";
    public final static String EN_US_Format = "MM/dd/yyyy HH:mm:ss";

    public final static int ISO8601 = 1;
    public final static int FileDate = 2;
    public final static int FileDateTime = 3;
    public final static int EN_US = 4;

    public static Date getDate() {
        return new Date();
    }
    public static String getFormattedDate(String format) {
        return format(new Date(), format);
    }
    public static String getFormattedDate(int knownStyle) {
        return format(new Date(), knownStyle);
    }
    public static String format(Date d, int knownStyle) {
        String fmt;
        switch (knownStyle) {
            case ISO8601: fmt = ISO8601Format; break;
            case FileDate: fmt = FileDateFormat; break;
            case FileDateTime: fmt = FileDateTimeFormat; break;
            case EN_US: fmt = EN_US_Format; break;
            default:
                return "";
        }
        return format(d, fmt);
    }
    public static String format(Date d, String format) {
        DateFormat dateFormat = new SimpleDateFormat(format);
        return dateFormat.format(d);
    }
}
