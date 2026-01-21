package ur_rna.Utilities;

/**
 * @author Richard M. Watson
 */
public enum Verbosity {
    Silent(0),
    Error(1),
    Warn(2),
    Info(3),
    Debug(4),
    Trace(5);

    public final int importance;
    Verbosity(int importance) { //, PrintStream ps) {
        this.importance = importance;
        //this.ps=ps;
    }
    public static Verbosity fromImportance(int importance) {
        Verbosity[] values = Verbosity.values();
        if (importance > -1 && importance < values.length)
            return values[importance];
        throw new IllegalArgumentException("Invalid log type: " + importance);
    }

//        private PrintStream ps;
//        public void setStream(PrintStream redirectToStream) {
//            this.ps = redirectToStream==null?nullStream:redirectToStream;
//        }
    //public PrintStream getStream() { return ps; }
}
