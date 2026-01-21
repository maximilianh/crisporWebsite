package ur_rna.Utilities;

import java.util.Iterator;
import java.util.NoSuchElementException;

public abstract class Filter<T> {
    public abstract boolean passes(T object);
    public Iterator<T> filter(Iterator<T> iterator) {
        return new FilterIterator(iterator);
    }
    public Iterable<T> filter(final Iterable<T> iterable) {
        return new Iterable<T>() {
            public Iterator<T> iterator() {
                return filter(iterable.iterator());
            }
        };
    }
    private class FilterIterator implements Iterator<T> {
        private static final int ITERATING = 0;
        private static final int AWAITING_CONSUMER = 1;
        private static final int AT_END = -1;
        private int state = ITERATING;
        private Iterator<T> iterator;
        private T foundItem;
        private FilterIterator(Iterator<T> iterator) {
            this.iterator = iterator;
        }
        public boolean hasNext() {
            //moves to next only if the last value has been consumed. returns true if another element has been found.
            return state != AT_END && moveToNext();
        }
        public T next() {
            moveToNext();
            if (state == AWAITING_CONSUMER)
                state = ITERATING; //indicate we have consumed the item.
            return foundItem;
        }
        public void remove() {
            throw new UnsupportedOperationException();
        }
        /**
         * Advances the underlying iterator and searches for the next passing item.
         * <p>
         * The underlying iterator is only advanced if the previously found item has been consumed
         * (by a call to next() function). This means that calls to hasNext() will only cause
         * the iterator to advance on the first call, and all subsequent calls will NOT advance the iterator
         * until next() has been called.
         * </p>
         *
         * @return True if a passing item was found or false otherwise. The latter result also indicates that the
         * underlying iterator reached its end, and subsequent calls to hasNext() should return false.
         */
        private boolean moveToNext() {
            switch (state) {
                case ITERATING:
                    while (iterator.hasNext()) {
                        T item = iterator.next();
                        if (passes(item)) {
                            foundItem = item;
                            state = AWAITING_CONSUMER; // indicate an item was found.
                            return true;
                        }
                    }
                    state = AT_END;
                    return false;
                case AWAITING_CONSUMER:
                    //  There is still a next item pending consumption.
                    // I.e. hasNext() was called and advanced the iterator, but the found item has not
                    // yet been consumed by a call to next()
                    // Therefore, there is item remaining to be consumed, so return true.
                    return true;
                case AT_END:
                default:
                    //tried to call next() when there are no more items.
                    throw new NoSuchElementException();
            }
        }
    }
}