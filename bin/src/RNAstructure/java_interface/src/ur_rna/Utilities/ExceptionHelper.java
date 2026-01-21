package ur_rna.Utilities;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.regex.Pattern;

/**
 * Tools for helping with exceptions.
 */
public class ExceptionHelper {
    /** Gets the stack trace as a String (i.e. the text that would be printed by printStackTrace) */
    public static String getStackTrace(Throwable ex) {
        StringWriter errWriter = new StringWriter();
        ex.printStackTrace(new PrintWriter(errWriter));
        return errWriter.toString();
    }

    /** Java exceptions often contain no message or a few details *without* describing the actual error.
     * For example {@link IndexOutOfBoundsException} has the message "Index ... Size ..." but does not
     * include "Index Out Of Bounds" in its message.
     * {@link java.io.FileNotFoundException} includes the file path, but does not include "File not found".
     *
     * So when displaying the message to a user, one is forced to either show {@link Exception#getMessage()}  which shows
     * meaningless information or {@link Exception#toString()} which includes the full type name,
     * e.g. "java.lang.nio.NoSuchFileException"
     */
    public static String getMessage(Throwable ex) {
        String name = ex.getClass().getName();
        String msg = ex.getLocalizedMessage();
        String s = msg==null ? name : name + ": " + msg;
        String actual = ex.toString();
        if (s.equals(actual)) {
            name = ex.getClass().getSimpleName();
            if (name.endsWith("Exception"))
                name = name.substring(0, name.length() - 9);
            name = Strings.toFriendlyName(name);
            return  msg==null ? name : name + ": " + msg;
        } else {
            String[] classStart = "java javax".split(" ");
            // Attempt to remove classpath
            for (String pkg : classStart)
                if (actual.startsWith(pkg + ".")) {
                    Pattern rex = Pattern.compile("^(\\w+\\.)+");
                    actual = rex.matcher(actual).replaceFirst("");
                    break;
                }
            return actual;
        }
    }
}
