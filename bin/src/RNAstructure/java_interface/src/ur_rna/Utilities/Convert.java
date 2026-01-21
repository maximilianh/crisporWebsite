package ur_rna.Utilities;

import ur_rna.Utilities.annotation.NotNull;
import ur_rna.Utilities.swing.FontUtil;

import java.awt.*;

/**
 * @author Richard M. Watson
 */
public class Convert {
    private Convert() {}
    public static boolean toBool(Object obj, boolean defaultIfNullOrEmpty) {
        if (obj == null || obj.equals("")) return defaultIfNullOrEmpty;
        return toBool(obj);
    }
    public static boolean toBool(Object obj) {
        if (obj == null) return false;
        if (obj instanceof Boolean)
            return (Boolean) obj;
        if (obj instanceof Number)
            return ((Number)obj).doubleValue() == 0d;
        if (obj instanceof String)
            return toBool((String)obj);
        if (obj instanceof Void)
            return false;
        if (obj instanceof java.util.Optional<?>)
            return ((java.util.Optional<?>)obj).isPresent();
        return true;
    }
    public static boolean toBool(String s) {
        if (s == null || s.length() == 0) return false;
        s = s.toUpperCase();
        return !(
                s.equals("F")
                        || s.equals("FALSE")
                        || s.equals("NO")
                        || s.equals("0")
        );
    }
//    public static int toInt(String s, int defaultValue) {
//        if (s == null) return defaultValue;
//        try {
//            return Integer.parseInt(s);
//        }catch (NumberFormatException e) {
//            return defaultValue;
//        }
//    }
//    public static double toDouble(String s, double defaultValue) {
//        if (s == null) return defaultValue;
//        try {
//            return Double.parseDouble(s);
//        }catch (NumberFormatException e) {
//            return defaultValue;
//        }
//    }
    public static int toInt(Object value) { return toInt(value, 0, false); }
    public static int toInt(Object value, int defaultValue) { return toInt(value, defaultValue, false); }
    public static int toInt(Object value, int defaultValue, boolean strict) {
        if (value == null) return defaultValue;
        if (value instanceof Number)
            return ((Number) value).intValue();
        if (value instanceof String)
            try {
                return value.equals("") ? defaultValue : Integer.parseInt((String) value);
            } catch (NumberFormatException ex) {
                if (strict) throw ex;
                return defaultValue;
            }

        if (value instanceof Color)
            return ((Color)value).getRGB();

        if (strict)
            throw createConvertException(value, Integer.class);
        return defaultValue;
    }
    public static float toFloat(Object value) { return toFloat(value, 0, false); }
    public static float toFloat(Object value, float defaultValue) { return toFloat(value, defaultValue, false); }
    public static float toFloat(Object value, float defaultValue, boolean strict) {
        if (value == null) return defaultValue;
        if (value instanceof Number)
            return ((Number)value).floatValue();
        if (value instanceof String)
            try {
                return value.equals("") ? defaultValue :Float.parseFloat((String) value);
            } catch (NumberFormatException ex) {
                if (strict) throw ex;
                return defaultValue;
            }
        if (strict)
            throw createConvertException(value, Float.class);
        return defaultValue;
    }
    public static double toDouble(Object value) { return toDouble(value, 0, false); }
    public static double toDouble(Object value, double defaultValue) { return toDouble(value, defaultValue, false); }
    public static double toDouble(Object value, double defaultValue, boolean strict) {
        if (value == null) return defaultValue;
        if (value instanceof Number)
            return ((Number)value).doubleValue();
        if (value instanceof String)
            try {
                return value.equals("") ? defaultValue : Double.parseDouble((String) value);
            } catch (NumberFormatException ex) {
                if (strict) throw ex;
                return defaultValue;
            }
        if (strict)
            throw createConvertException(value, Double.class);
        return defaultValue;
    }
    public static char toChar(Object value) { return toChar(value, '\0', false); }
    public static char toChar(Object value, char defaultValue) { return toChar(value, defaultValue, false); }
    public static char toChar(Object value, char defaultValue, boolean strict) {
        if (value == null) return defaultValue;
        if (value instanceof Character)
            return (Character)value;
        if (value instanceof Number)
            return (char)((Number)value).intValue();
        if (value instanceof String)
            return value.equals("") ? defaultValue : ((String) value).charAt(0);
        if (strict)
            throw createConvertException(value, Character.class);
        return defaultValue;
    }
    public static String toString(Object value) { return toString(value, null); }
    public static String toString(Object value, String valueIfNull) {
        if (value == null) return valueIfNull;
        if (value instanceof String)
            return (String)value;

        if (value instanceof Color)
            return Colors.getName((Color)value);

        if (value instanceof Font)
            return FontUtil.encode((Font)value);

        return value.toString();
    }
    public static Object toType(Object value, Class targetType) {
        if (value!=null&&value.getClass()==targetType) return value;

        if (targetType == Integer.class || targetType == int.class)
            return toInt(value);
        if (targetType == Float.class || targetType == float.class)
            return toFloat(value);
        if (targetType == Double.class || targetType == double.class)
            return toDouble(value);
        if (targetType == Character.class || targetType == char.class)
            return toChar(value);
        if (targetType == Boolean.class || targetType == boolean.class)
            return toBool(value);
        if (targetType == String.class)
            return toString(value);
        if (targetType == Font.class)
            return toFont(value);
        if (targetType == Color.class)
            return toColor(value);
        throw createConvertException(value, targetType);
    }
    @SuppressWarnings("unchecked")
    public static <T> T toType(Object value, @NotNull T defaultValue) {
        // these are all final classes, so if defaultValue is an instance of one of them,
        // then T is exactly equal to that class (Not a sub-class)
        // Therefore a cast error will never occur.
        if (defaultValue instanceof Integer)
            return (T)(Integer)toInt(value, (Integer) defaultValue);
        if (defaultValue instanceof Float)
            return (T)(Float)toFloat(value);
        if (defaultValue instanceof Double)
            return (T)(Double)toDouble(value);
        if (defaultValue instanceof Character)
            return (T)(Character)toChar(value);
        if (defaultValue instanceof Boolean)
            return (T)(Boolean)toBool(value);
        if (defaultValue instanceof String)
            return (T)toString(value);
        throw createConvertException(value, defaultValue.getClass());
    }

    protected static RuntimeException createConvertException(Object val, Class target) {
        return new RuntimeException(String.format("Unable to convert from %s to %s.",
                val == null ? "NULL" : val.getClass().getSimpleName(),
                target.getSimpleName()));
    }
    /**
     * Converts a standard POSIX Shell globbing pattern into a regular expression
     * pattern. The result can be used with the standard {@link java.util.regex} API to
     * recognize strings which match the glob pattern.
     * <p>
     * See also, the POSIX Shell language:
     * http://pubs.opengroup.org/onlinepubs/009695399/utilities/xcu_chap02.html#tag_02_13_01
     * By: Neil Traft - http://stackoverflow.com/questions/1247772/is-there-an-equivalent-of-java-util-regex-for-glob-type-patterns
     * </p>
     * @param pattern A glob pattern.
     * @return A regex pattern to recognize the given glob pattern.
     *
     */
    public static String globToRegex(String pattern) {
        StringBuilder sb = new StringBuilder(pattern.length());
        int inGroup = 0;
        int inClass = 0;
        int firstIndexInClass = -1;
        char[] arr = pattern.toCharArray();
        for (int i = 0; i < arr.length; i++) {
            char ch = arr[i];
            switch (ch) {
                case '\\':
                    if (++i >= arr.length) {
                        sb.append('\\');
                    } else {
                        char next = arr[i];
                        switch (next) {
                            case ',':
                                // escape not needed
                                break;
                            case 'Q':
                            case 'E':
                                // extra escape needed
                                sb.append('\\');
                            default:
                                sb.append('\\');
                        }
                        sb.append(next);
                    }
                    break;
                case '*':
                    if (inClass == 0)
                        sb.append(".*");
                    else
                        sb.append('*');
                    break;
                case '?':
                    if (inClass == 0)
                        sb.append('.');
                    else
                        sb.append('?');
                    break;
                case '[':
                    inClass++;
                    firstIndexInClass = i+1;
                    sb.append('[');
                    break;
                case ']':
                    inClass--;
                    sb.append(']');
                    break;
                case '.':
                case '(':
                case ')':
                case '+':
                case '|':
                case '^':
                case '$':
                case '@':
                case '%':
                    if (inClass == 0 || (firstIndexInClass == i && ch == '^'))
                        sb.append('\\');
                    sb.append(ch);
                    break;
                case '!':
                    if (firstIndexInClass == i)
                        sb.append('^');
                    else
                        sb.append('!');
                    break;
                case '{':
                    inGroup++;
                    sb.append('(');
                    break;
                case '}':
                    inGroup--;
                    sb.append(')');
                    break;
                case ',':
                    if (inGroup > 0)
                        sb.append('|');
                    else
                        sb.append(',');
                    break;
                default:
                    sb.append(ch);
            }
        }
        return sb.toString();
    }


    public static Font toFont(Object value) { return toFont(value, null, false); }
    public static Font toFont(Object value, Font defaultValue) { return toFont(value, defaultValue, false); }
    public static Font toFont(Object value, Font defaultValue, boolean strict) {
        if (value == null) return defaultValue;
        if (value instanceof Font)
            return (Font)value;
        if (value instanceof String)
            return Font.decode((String) value);
        if (value instanceof Number)
            return FontUtil.defaultFont.deriveFont(((Number)value).floatValue());
        if (strict)
            throw createConvertException(value, Font.class);
        return defaultValue;
    }
    public static Color toColor(Object value) { return toColor(value, null, false); }
    public static Color toColor(Object value, Color defaultValue) { return toColor(value, defaultValue, false); }
    public static Color toColor(Object value, Color defaultValue, boolean strict) {
        if (value == null) return defaultValue;
        if (value instanceof Color)
            return (Color)value;
        if (value instanceof String) {
            Color c = Colors.getColor((String) value);
            return c == null ? defaultValue : c;
        }
        if (value instanceof Number) {
            int num = ((Number) value).intValue();
            return new Color(num, num > 0xFFFFFFL || num < 0);
        }
        if (strict)
            throw createConvertException(value, Color.class);
        return defaultValue;
    }
}
