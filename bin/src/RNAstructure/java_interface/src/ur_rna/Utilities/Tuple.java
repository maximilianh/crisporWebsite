package ur_rna.Utilities;

import java.util.Comparator;
import java.util.Objects;

/**
 * Implementation of Tuples that have mutable items.
 * Note that equals and hashCode have NOT been overridden, so two distinct tuples with identical items are NOT considered equal.
 * (The items should be immutable for equals and hashCode to be based on them, otherwise we risk changing the hashCode of a tuple AFTER
 *  it has been stored in a hash-table, which would lead to unspecified and probably wrong results).
 */
public abstract class Tuple {
//    public static <T1, T2> Pair<T1, T2> create(T1 t1, T2 t2) {
//        return new Pair<T1, T2>(t1, t2);
//    }
//    public static <T1, T2, T3> Triplet<T1, T2, T3> create(T1 t1, T2 t2, T3 t3) {
//        return new Triplet<T1, T2, T3>(t1, t2, t3);
//    }

    public static final class Singlet<T1> {
        public T1 item1;
        public Singlet(final T1 item1) {
            this.item1 = item1;
        }
        public Singlet() { }

    }
    public static final class Pair<T1, T2> implements Cloneable {
        public T1 item1;
        public T2 item2;
        public Pair(final T1 item1, final T2 item2) {
            this.item1 = item1;
            this.item2 = item2;
        }
        public Pair() { }

        public static <T extends Comparable<? super T>> Comparator<Pair<? extends T,?>> compareFirstItem() {
            return (o1, o2) -> o1.item1.compareTo(o2.item1);
        }
        public static <T1 extends Comparable<? super T1>, T2 extends Comparable<? super T2>> Comparator<Pair<? extends T1,? extends T2>> compareBothItems() {
            return (o1, o2) -> {
                int cmp = o1.item1.compareTo(o2.item1);
                return cmp != 0 ? cmp : o1.item2.compareTo(o2.item2);
            };
        }

        public boolean itemsEqual(final Pair<T1, T2> obj) {
            return obj == this || Objects.equals(item1, obj.item1) && Objects.equals(item2, obj.item2);
        }
        @Override
        @SuppressWarnings("unchecked")
        protected Pair<T1, T2> clone()  {
            try {
                return (Pair<T1, T2>) super.clone();
            } catch (CloneNotSupportedException ex) {
                throw new InternalError(ex);
            }
        }
        @Override
        public String toString() {
            return String.format("%s { %s, %s }", Pair.class.getSimpleName(),  item1, item2);
        }
    }
    public static final class Triplet<T1, T2, T3>  implements Cloneable {
        public T1 item1;
        public T2 item2;
        public T3 item3;
        public Triplet(final T1 item1, final T2 item2, T3 item3) {
            this.item1 = item1;
            this.item2 = item2;
            this.item3 = item3;
        }
        public Triplet() {
        }
        public boolean itemsEqual(final Triplet<T1, T2, T3> obj) {
            return obj == this || Objects.equals(item1, obj.item1) && Objects.equals(item2, obj.item2) && Objects.equals(item3, obj.item3);
        }
        @Override
        @SuppressWarnings("unchecked")
        protected Triplet<T1, T2, T3> clone()  {
            try {
                return (Triplet<T1, T2, T3>) super.clone();
            } catch (CloneNotSupportedException ex) {
                throw new InternalError(ex);
            }
        }
        @Override
        public String toString() {
            return String.format("%s { %s, %s, %s }", Triplet.class.getSimpleName(),  item1, item2, item3);
        }
    }
    public static final class Quad<T1, T2, T3, T4> {
        public T1 item1;
        public T2 item2;
        public T3 item3;
        public T4 item4;
        public Quad(final T1 item1, final T2 item2, T3 item3, T4 item4) {
            this.item1 = item1;
            this.item2 = item2;
            this.item3 = item3;
            this.item4 = item4;
        }
        public Quad() {
        }
        public boolean itemsEqual(final Quad<T1, T2, T3, T4> obj) {
            return obj == this || Objects.equals(item1, obj.item1) && Objects.equals(item2, obj.item2) && Objects.equals(item3, obj.item3) && Objects.equals(item4, obj.item4);
        }
        @Override
        @SuppressWarnings("unchecked")
        protected Quad<T1, T2, T3, T4> clone()  {
            try {
                return (Quad<T1, T2, T3, T4>) super.clone();
            } catch (CloneNotSupportedException ex) {
                throw new InternalError(ex);
            }
        }
        @Override
        public String toString() {
            return String.format("%s { %s, %s, %s, %s }", Quad.class.getSimpleName(),  item1, item2, item3, item4);
        }
    }
}