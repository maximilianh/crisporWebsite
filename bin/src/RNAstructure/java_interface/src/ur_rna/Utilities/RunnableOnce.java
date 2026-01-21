package ur_rna.Utilities;

/**
 * Indicates that a call to {@link #run()} will be effective at most once.
 */
public interface RunnableOnce extends Runnable {
    /**
     * Indicates whether {@link #run()} has been called.
     *
     * @return True if {@link #run()} has already been called. This can be reset by a call to {@link #reset()}.
     */
    boolean hasRun();

    /**
     * Resets the internal state of the RunnableOnce instance, so {@link #run()} can be called again,
     * even if it has already been called once.
     */
    void reset();
}
