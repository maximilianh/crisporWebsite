package ur_rna.Utilities;

/**
 * Represents an event listener.
 */
@FunctionalInterface
public interface Listener<T> {
    void handleEvent(T eventInfo);
}
