package ur_rna.Utilities;

/**
 * This class operates as a proxy for a {@link Runnable }.
 * When {@link #run()} is called initially, it is delegated to {@link Runnable#run()}, and
 * subsequent calls to {@link #run()} are be ignored (i.e. they do not result in a call to {@link Runnable#run()}).
 * The method {@link #hasRun()} can be used to determine whether or not {@link #run()} has been called.
 * The {@link #reset()} method can be called to re-enable {@link #run()} to make another one-time call to
 * {@link Runnable#run()}).
 */
public class RunOnce implements RunnableOnce {
    private boolean has_run = false;
    private final Runnable inner;
    public RunOnce() { inner = null; }
    public RunOnce(Runnable innerRunnable) {
        inner = innerRunnable;
    }

    /**
     * Indicates whether {@link #run()} has been called.
     *
     * @return True if {@link #run()} has already been called. This can be reset by a call to {@link #reset()}.
     */
    public boolean hasRun() { return has_run; }

    /**
     * Resets the internal state of the RunnableOnce instance, so {@link #run()} can be called again,
     * even if it has already been called once.
     */
    public void reset() {
        has_run = false;
    }

    /**
     * When an object implementing interface <code>Runnable</code> is used
     * to create a thread, starting the thread causes the object's
     * <code>run</code> method to be called in that separately executing
     * thread.
     * <p>
     * The general contract of the method <code>run</code> is that it may
     * take any action whatsoever.
     *
     * @see Thread#run()
     */
    @Override
    public void run() {
        has_run = true;
        if (inner != null) inner.run();
    }

    /**
     * Fakes a call to run. After calling this, {@link #hasRun()} will return true,
     * but {@link #run()} is never actually called.
     */
    public void mockRun() {
        has_run = true;
    }
}
