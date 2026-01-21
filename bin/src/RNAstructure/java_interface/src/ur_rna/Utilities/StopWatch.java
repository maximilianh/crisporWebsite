package ur_rna.Utilities;

import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Class used for debugging and diagnostics to display the ellapsed time for sections of code at runtime.
 * Use the start, stop, reset, and restart methods to determine when the timer is running.
 * Use the elapsed, println and prinf functions to determine the ellapsed time and/or output debug messages.
 * @author Richard M. Watson
 */
public class StopWatch {
    // Global stopwatches that are not local to any function
    public static PrintWriter defaultOutput = new PrintWriter(System.out, true);
    public static final Map<String, StopWatch> watches = new LinkedHashMap<>();

    /** Get named StopWatch if it exists. Otherwise create it and optionally start it if createStarted is true. */
    public static StopWatch get(String name, boolean createStarted) {
        StopWatch w = watches.get(name);
        if (w == null) w = create(name, createStarted);
        return w;
    }
    /** Get named StopWatch if it exists. Otherwise return null. */
    public static StopWatch get(String name) { return watches.get(name); }
    /** Create a new named StopWatch and optionally start it if `started` is true.
     * A KeyAlreadyExistsException is thrown if the named StopWatch already exists.  */
    public static StopWatch create(String name, boolean started) {
        StopWatch w = new StopWatch(started);
        if (watches.containsKey(name)) throw new KeyAlreadyExistsException("StopWatch "+name+" already exists.");
        watches.put(name, w);
        return w;
    }
    /** Write a message to the named StopWatch. */
    public static StopWatch println(String watchName, String message) { return watches.get(watchName).println(message); }
    /** Write a message to the named StopWatch. */
    public static StopWatch printf(String watchName, String message, Object... args) { return watches.get(watchName).printf(message, args); }
    /** Get the elapsed time of the named StopWatch. */
    public static long elapsed(String watchName) { return watches.get(watchName).elapsed(); }
    /** Gets the named StopWatch and starts it. If it doesn't already exist, it will be created. */
    public static StopWatch enter(String watchName) { return get(watchName, true).start(); }
    /** Gets the named StopWatch and stops it. If it doesn't already exist, it will be created. */
    public static StopWatch exit(String watchName) {
        StopWatch w = get(watchName, true);
        if (!w.isStarted()) throw new IllegalStateException("StopWatch "+watchName+"is not yet started.");
        return w.stop();
    }

    /**
     * Prints the elapsed time and current state of all named StopWatches.
     * Also calculates the percentage of each StopWatch relative to the totalElapsed time.
     * If totalElapsed is less than or equal to 0, percentages are not calculated.
     */
    public void printAll() { printAll(0); }
    /** Prints the elapsed time and current state of all named StopWatches */
    public void printAll(long totalElapsed) { printAll(totalElapsed==0?"\t%1$5d\t%4$5s%%\t%2$3s\t%3$s%n":"\t%1$5d\t%2$3s\t%3$s%n", totalElapsed); }
    /**
     * Prints the elapsed time and current state of all named StopWatches using the specified format.
     * Use positional formatting 1=Elapsed Time, 2=State (on/off), 3=Name
     * e.g. {@code "\t%1$5d\t%2$3s\t%2$s%n"}
     */
    public void printAll(String format, long totalElapsed) {
        for (Map.Entry<String, StopWatch> entry : watches.entrySet()) {
            StopWatch w = entry.getValue();
            defaultOutput.printf(format, w.elapsed(), w.isStarted() ? "on" : "off", entry.getKey(), totalElapsed == 0 ? 0 : 100 & w.elapsed() / totalElapsed);
        }
    }
    /** Remove all named StopWatches */
    public void clearNamed() { watches.clear(); }


    // Stores the time at which the StopWatch was started.
    // If the watch is stopped, it stores the elapsed time as a negative number.
    // When restarted, the time (as a negative number) is added to the CURRENT time so that the total elapsed will be correct.
    // So:
    //    (time == 0) --  Not yet started. Elapsed time is 0.
    //    (time > 0)  --  Already started and currently running. Elapsed time is `time() - time`
    //    (time < 0)  --  Started previously, but currently stopped. Elapsed time is `-time`
    private long time;
    private PrintWriter out;

    public StopWatch(){ this(false, defaultOutput);}
    public StopWatch(boolean start){ this(start, defaultOutput); }
    public StopWatch(boolean start, PrintWriter output){ out = output; if (start) start(); }
    public StopWatch(boolean start, PrintStream output){ this(start, new PrintWriter(output, true)); }

    /**
     * Start the timer. If it is already running, this has no effect.
     */
    public StopWatch start() {
        if (time == 0)
            time = time();
        else if (time < 0)
            time = time + System.currentTimeMillis();
        return this;
        // otherwise it is already started.
    }
    /**
     * Stop the timer, but store the elapsed time so that when {@link #start()} is called,
     * the elapsed time will continue to increase.
     * If the timer is already stopped, this has no effect.
     */
    public StopWatch stop() {
        if (time > 0)
            time = time - System.currentTimeMillis(); // get amount of time run, but store it as negative to indicate that it is a duration.
        return this;
    }
    /**
     * The total elapsed time in milliseconds.
     */
    public long elapsed() {
        if (time == 0) return 0;
        return (time < 0 ? 0 : time()) - time;
    }
    /**
     * Reset elapsed time to 0 and stop the timer.
     */
    public StopWatch reset() {
        time = 0;
        return this;
    }
    /**
     * Reset elapsed time to 0 and start the timer.
     */
    public StopWatch restart() {
        time = time();
        return this;
    }
    public boolean isStarted() { return time > 0; }

    public static long time() {
        return System.currentTimeMillis();
    }

    public StopWatch println(String message) {
        out.printf("@\t%5d\t%s%n", elapsed(), message);
        return this;
    }
    public StopWatch printf(String message, Object... args) {
        out.printf("@\t%5d\t%s", elapsed(), String.format(message, args));
        return this;
    }
}
