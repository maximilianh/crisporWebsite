package ur_rna.Utilities;

import java.io.*;
import java.util.HashMap;

import static ur_rna.Utilities.Strings.asBool;

/**
 * A utility class to streamline and standardize logging across applications.
 */
public class AppLog {
    protected static PrintStream nullStream = new NullPrintStream();
    public static Verbosity DefaultLogImportance = Verbosity.Info;
    private static AppLog defaultAppLog;
    public static AppLog getDefault() {
        if (defaultAppLog == null)
            defaultAppLog = new AppLog();
        return defaultAppLog;
    }
    public static void setDefault(AppLog defLog) {
        defaultAppLog = defLog;
    }

    //** Use of the following stream references allows the application to easily redirect. */
    private PrintStream[] streams = new PrintStream[] {
            AppLog.nullStream, // Verbosity.Silent
            System.err, // Verbosity.Error
            System.err, // Verbosity.Warn
            System.out, // Verbosity.Info
            System.out, // Verbosity.Debug
            System.out  // Verbosity.Trace
    };
    // In general, messages are logged if their importance is <= the set verbosity.
    // However, by setting forceEnabled[i] to true, a message can be shown if its importance is
    // exactly i even if i > verbosity.
    private boolean[] forceEnabled = new boolean[ObjTools.max(Verbosity.values(), v->v.importance)+1];

    public AppLog() {
        //readSystemProperties();
    }

    /** Checks for the existence of a system property in multiple places.
     * 1) System.getProperty(key)
     * 2) System.getenv(key)
     * 3) System.getenv(APPNAME_key)
     *
     * Upper/lower case is also checked for getenv.
     * If key contains hyphens, they are converted to underscores for getenv.
     */
    public static String getSysProp(String... possibleKeys) {
        if (possibleKeys.length == 1)
            return getSysProp(possibleKeys[0]);
        for(String key : possibleKeys) {
            String value = getSysProp(key);
            if (key != null) return value;
        }
//        String name = System.getProperty("APPNAME");
//        if (name == null) return null;
        return null;
    }

    public static String getSysProp(String key) {
        String value = System.getProperty(key); if (value != null) return value;

        // If not found as a system property, check environment.
        value = System.getenv(key); if (value != null) return value;

        // If getting from the environment, also check <APP-NAME>_<PROPERTY>
        String appPrefix = AppInfo.getAppName().toLowerCase() + "_";
        value = System.getenv(appPrefix+key); if (value != null) return value;

        // hyphen may be unavailable for use in environment variables. convert to underscore
        if (key.indexOf('-')!=-1) {
            key = key.replace('-', '_');
            value = System.getenv(key); if (value != null) return value;
            value = System.getenv(appPrefix+key); if (value != null) return value;
        }

        // Try UCase and LCase from environment.
        value = System.getenv(key.toLowerCase()); if (value != null) return value;
        value = System.getenv(key.toUpperCase()); if (value != null) return value;
        value = System.getenv((appPrefix + key).toLowerCase()); if (value != null) return value;
        value = System.getenv((appPrefix + key).toUpperCase());

        return value;
    }

    public void readSystemProperties() {
        String val = getSysProp("log-verbose", "verbose", "log-verbosity", "verbosity");
        if (val!=null) setVerbosity(Convert.toInt(val, DefaultLogImportance.importance));

        PrintStream out = (null!=(val=getSysProp("log-out-file")))?getCustomStream(val):null;
        PrintStream err = (null!=(val=getSysProp("log-err-file")))?getCustomStream(val):null;

        for(Verbosity v : Verbosity.values()) {
            // Check for e.g. log-debug  or  <app_name>_log-debug
            val = getSysProp("log-" + v.name().toLowerCase());
            if (val!=null) setEnabledSpecial(v, asBool(val));

            // Check for e.g. log-debug-file  or  <app_name>_log-debug-file
            val = getSysProp("log-" + v.name().toLowerCase() + "-file");
            if (val!=null) setStream(v, getCustomStream(val));

            if (out!=null&&getStream(v)==System.out) setStream(v, out);
            if (err!=null&&getStream(v)==System.err) setStream(v, err);
        }
    }

    private HashMap<File,PrintStream> _customStreams = new HashMap<File, PrintStream>();
    private PrintStream getCustomStream(String path) {
        switch (path) {
            case "!null": case "!close": case "!none": case "!silent": return nullStream;
            case "!out": case "stdout": return System.out;
            case "!err": case "stderr": return System.err;
            default:
                try {
                    // Store the stream, so that if two verbosity levels use the same file, only a single stream will be created.
                    // (otherwise there might be sharing violations) or overwriting.
                    File f = new File(path).getCanonicalFile();
                    PrintStream ps = _customStreams.get(f);
                    if (ps == null) {
                        ps = new PrintStream(new FileOutputStream(f, true), true,"UTF-8");
                        _customStreams.put(f, ps);
                    }
                    return ps;
                } catch (IOException e) {
                    e.printStackTrace();
                    return System.out;
                }
        }
    }

    /** Returns true if Trace log messages are enabled, either via the current Verbosity setting or by force-enabling
     * them specifically (e.g. via {@link #setEnabledSpecial(Verbosity, boolean)}.
     * @return True if Trace is enabled (by any means).
     */
    public boolean isTraceEnabled() { return isEnabled(Verbosity.Trace); }

    /** Returns true if Debug log messages are enabled, either via the current Verbosity setting or by force-enabling
     * them specifically (e.g. via {@link #setEnabledSpecial(Verbosity, boolean)}.
     * @return True if Debug is enabled (by any means).
     */
    public boolean isDebugEnabled() { return isEnabled(Verbosity.Debug); }

    /** Returns true if the specified log-level is enabled, either via the current Verbosity setting or by
     * force-enabling the specific log-level (e.g. via {@link #setEnabledSpecial(Verbosity, boolean)}.
     * @param verbosity The log-level to check.
     * @return True if the specified log-level is enabled (by any means).
     */
    public boolean isEnabled(final Verbosity verbosity) {
        return verbosity.importance <= _requiredVerbosity.importance || forceEnabled[verbosity.importance];
    }

    /** Returns true if the specified log-level is enabled, either via the current Verbosity setting or by
     * force-enabling the specific log-level (e.g. via {@link #setEnabledSpecial(Verbosity, boolean)}.
     * @param verbosity The log-level to check.
     * @return True if the specified log-level is enabled (by any means).
     */
    public boolean isEnabled(final int verbosity) {
        return verbosity <= _requiredVerbosity.importance || forceEnabled[verbosity];
    }

    /** Returns true if the specified log-level has been force-enabled
     * (e.g. via {@link #setEnabledSpecial(Verbosity, boolean)}.
     * Note that this log-level may be enabled by the Verbosity setting, regardless of whether it has been
     * force-enabled.
     * @param verbosity The log-level to check.
     * @return True if the specified log-level is force-enabled.
     */
    public boolean isEnabledSpecial(final Verbosity verbosity) {
        return forceEnabled[verbosity.importance];
    }
    public void setEnabledSpecial(final Verbosity verbosity, boolean enabled) {
        forceEnabled[verbosity.importance] = enabled;
    }

    protected Verbosity _requiredVerbosity = Verbosity.Info;

//    public boolean getDebug() { return isSpecificallyEnabled(Verbosity.Debug); }
//    public void setDebug(boolean value) { setSpecificallyEnabled(Verbosity.Debug, value); }
//    public boolean getTrace() { return isSpecificallyEnabled(Verbosity.Trace); }
//    public void setTrace(boolean value) { setSpecificallyEnabled(Verbosity.Trace, value); }

    /** Returns the log-level (aka importance or verbosity) required for log entries to be processed.
     * Log entries with importance lower than this will be ignored (unless they are force-enabled).
     */
    public Verbosity getVerbosity() { return _requiredVerbosity; }
    public void setVerbosity(int value) throws IllegalArgumentException {
        _requiredVerbosity = Verbosity.fromImportance(value);
    }
    public void setVerbosity(Verbosity value) throws IllegalArgumentException {
        getStream(value); //throws error if invalid.
        _requiredVerbosity = value;
    }

    public void log(String s) {
        log(s, DefaultLogImportance);
    }
    public void log(String format, Object ... args) {
        logFmt(format, DefaultLogImportance, args);
    }
    private void logFmt(String format, Verbosity importance, Object ... args) {
        if (args == null || args.length == 0)
            log(format, importance);
        else
            log(String.format(format, args), importance);
    }
    public void log(String s, Verbosity messageVerbosity) {
        try {
            getStream(messageVerbosity).println(s);
        } catch (Throwable e) {
            e.printStackTrace(System.err);
            System.err.println(s);
        }
    }

    public PrintStream getStream(Verbosity streamVerbosity) {
        int i = streamVerbosity.importance;
        if (i < 0 || i >= streams.length)
            throw new InternalError("Invalid log type: " + i); //should never happen because it is not one of the defined values.
        return isEnabled(i) ? streams[i] : nullStream;
    }

    public void setStream(Verbosity streamVerbosity, PrintStream stream) {
        streams[streamVerbosity.importance] = stream;
    }

    public PrintStream getDbgStream() { return getStream(Verbosity.Debug); }
    public PrintStream getErrStream() { return getStream(Verbosity.Error); }
    public PrintStream getTrStream() { return getStream(Verbosity.Trace); }

    public void trace(String s) { log(s, Verbosity.Trace); }
    public void debug(String s) { log(s, Verbosity.Debug); }
    public void info(String s) { log(s, Verbosity.Info); }
    public void warn(String s) { log(s, Verbosity.Warn); }
    public void error(String s) { log(s, Verbosity.Error); }

    public void trace(String format, Object... args) { logFmt(format, Verbosity.Trace, args); }
    public void debug(String format, Object... args) { logFmt(format, Verbosity.Debug, args); }
    public void info(String format, Object... args) { logFmt(format, Verbosity.Info, args); }
    public void warn(String format, Object... args) { logFmt(format, Verbosity.Warn, args); }
    public void error(String format, Object... args) { logFmt(format, Verbosity.Error, args); }
    public void error(String s, Throwable e) { log(s, Verbosity.Error); log(getErrorInfo(e), Verbosity.Error); }

    public static String getErrorInfo(Throwable t) {
        java.io.StringWriter sw = new java.io.StringWriter();
        t.printStackTrace(new java.io.PrintWriter(sw));
        return sw.toString();
    }
}
