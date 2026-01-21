package ur_rna.Utilities;

import java.util.Map;
import java.util.Objects;

/** General purpose key-value pair. */
public class Pair<K,V> implements Map.Entry {
    public final K key;
    public final V value;
    public Pair(K key, V value) { this.key = key; this.value = value; }

    /**
     * Returns the key corresponding to this entry.
     *
     * @return the key corresponding to this entry
     * @throws IllegalStateException implementations may, but are not
     *                               required to, throw this exception if the entry has been
     *                               removed from the backing map.
     */
    @Override
    public K getKey() {
        return key;
    }
    /**
     * Returns the value corresponding to this entry.  If the mapping
     * has been removed from the backing map (by the iterator's
     * <tt>remove</tt> operation), the results of this call are undefined.
     *
     * @return the value corresponding to this entry
     * @throws IllegalStateException implementations may, but are not
     *                               required to, throw this exception if the entry has been
     *                               removed from the backing map.
     */
    @Override
    public V getValue() {
        return value;
    }
    /**
     * Replaces the value corresponding to this entry with the specified
     * value (optional operation).  (Writes through to the map.)  The
     * behavior of this call is undefined if the mapping has already been
     * removed from the map (by the iterator's <tt>remove</tt> operation).
     *
     * @param value new value to be stored in this entry
     * @return old value corresponding to the entry
     * @throws UnsupportedOperationException if the <tt>put</tt> operation
     *                                       is not supported by the backing map
     * @throws ClassCastException            if the class of the specified value
     *                                       prevents it from being stored in the backing map
     * @throws NullPointerException          if the backing map does not permit
     *                                       null values, and the specified value is null
     * @throws IllegalArgumentException      if some property of this value
     *                                       prevents it from being stored in the backing map
     * @throws IllegalStateException         implementations may, but are not
     *                                       required to, throw this exception if the entry has been
     *                                       removed from the backing map.
     */
    @Override
    public V setValue(final Object value) {
        throw new UnsupportedOperationException("The Pair is immutable -- changes cannot " +
                "be made to the key or the value after creation.");
    }

    @Override
    public int hashCode() {
        return (Objects.hashCode(key) * 16777619) ^ Objects.hashCode(value);
    }
    @Override
    public boolean equals(final Object other) {
        return other instanceof Map.Entry && equals((Map.Entry) other);
    }
    public boolean equals(final Map.Entry other) {
        return this == other || Objects.equals(key, other.getKey()) && Objects.equals(value, other.getValue());
    }
}
