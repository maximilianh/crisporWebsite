package ur_rna.Utilities;

import ur_rna.Utilities.annotation.ToFriendlyString;

import java.awt.*;
import java.lang.annotation.Annotation;
import java.lang.reflect.Array;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;
import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;

public abstract class ObjTools {
    public static final String[] EMPTY_STRING_ARRAY = new String[0];
    public static final Object[] EMPTY_OBJECT_ARRAY = EMPTY_STRING_ARRAY;
    public static final int[] EMPTY_INT_ARRAY = new int[0];
    public static final boolean[] EMPTY_BOOL_ARRAY = new boolean[0];

    /**
     * An object that can be used to represent a missing value.
     * For example, a Map's get(key) function will return null under either of two circumstances:
     * (1) the key was not found or (2) the key was associated with a null value.
     * Instead of using null as an ambiguous result, MISSING can be returned in the case of
     * the missing key in such circumstances.
     */
    public static Object MISSING = new Object();

    private ObjTools(){} //prevent instantiation

    @SuppressWarnings("unchecked") // getMethod results in a warning about using unchecked or unsafe operations.
    public static String toDisplayString(Object value) {
        if (value == null)
            return "null";
        if (value instanceof String)
            return "\"" + value + "\"";
        if (value instanceof Double || value instanceof Integer)
            return value.toString(); // + value.getClass().getSimpleName().substring(0,1); // D,F,I,B

        Class cls =value.getClass();
        String name = cls.getSimpleName();

        if (value instanceof Number)
            return value.toString() + name.substring(0, 1).toLowerCase(); // D,F,I,B
        if (value instanceof Component) {
            Component c = (Component) value;
            Rectangle r = c.getBounds();
            return "[GUI " + name + " \"" + c.getName() + "\" Bounds: (x:" + r.x + ", y:" + r.y + ", w:" + r.width + ", h:" + r.height + ") ]";
        }
        if (value instanceof Throwable) {
            return "[ERROR: " + ((Throwable) value).getMessage() + "]";
        }
        if (value == MISSING)
            return "[MISSING]";

        if (value instanceof Displayable)
            return ((Displayable) value).toDisplayString();

		/* List classes here that can be converted by using toString() directly, without prepending the class name */
        Class[] knownTypes = new Class[]{Boolean.class};
        for (Class c : knownTypes) {
            if (c.isInstance(value))
                return value.toString();
        }

        Annotation[] notes = cls.getAnnotationsByType(ToFriendlyString.class);
        if (notes.length != 0) {
            ToFriendlyString a = (ToFriendlyString)notes[0];
            String method = a.method();
            if ("toString".equals(method))
                return value.toString();
            if (!method.isEmpty()) {
                try {
                    Method m = cls.getMethod(method);
                    if (m != null)
                        return m.invoke(value).toString();
                } catch (Exception ex) {
                    // do nothing. use another method.
                }
            }
        }

        String s = value.toString();
        if (s.startsWith("[") && s.endsWith("]") || s.startsWith("{") && s.endsWith("}"))
            return s;
        return "[" + name + ": " + s + "]";
    }
    public static String toStr(Object value) {
        return toStr(value, "");
    }
    public static String toStr(Object value, String nullValue) {
        if (value == null) return nullValue;
        return value.toString();
    }

    public static boolean isSet(final int flagToFind, final int valueToSearch) {
        return (valueToSearch & flagToFind) == flagToFind;
    }
    public static boolean isAnySet(final int flagToFind, final int valueToSearch) {
        return (valueToSearch & flagToFind) != 0;
    }

    /**
     * If obj is a boolean, its value is returned.
     * If obj is null, an empty array, or a number, char or byte representing 0, the return value is False.
     * If obj is a string, then it is evaluated as described in AsBool(string)
     * Otherwise the return value is the the result of AsBool(string) applied to obj.ToString()
     */
    public static boolean asBool(final Object obj) {
        if (obj == null) return false;
        if (obj instanceof String) return Strings.asBool((String) obj, false);
        if (obj instanceof Boolean)
            return (Boolean) obj;
        if (obj instanceof Byte)
            return 0 != (byte) obj;
        if (obj instanceof Integer)
            return 0 != (int) obj;
        if (obj instanceof Float)
            return 0F != (float) obj;
        if (obj instanceof Long)
            return 0L != (long) obj;
        if (obj instanceof Character)
            return 0 == (Character)obj;
        if (obj instanceof Array)
            return 0 != Array.getLength(obj);
        return Strings.asBool(obj.toString(), false); //uses string version. "", false, no, "0" evaluate to false. All others are true, including "{TypeName}" etc.
    }

    /** Copies from the list to the array */
    public static <T> void copyTo(java.util.List<? extends T> list, T[] arr) { copyTo(list, arr, 0, list.size()); }
    /** Copies from the list to the array. The starting destination array index is specified by the {@code start} parameter. */
    public static <T> void copyTo(java.util.List<? extends T> list, T[] arr, int start) { copyTo(list, arr, start, list.size()); }
    /** Copies the specified number of elements from the list to the array. The starting destination array index is specified by the {@code start} parameter. */
    public static <T> void copyTo(java.util.List<? extends T> list, T[] arr, int start, int length) {
        for (int i = 0; i < length; i++)
            arr[start + i] = list.get(i);
    }
    /** Copies from the collection to the array */
    public static <T> void copyTo(java.util.Collection<? extends T> list, T[] arr) { copyTo(list, arr, 0, list.size()); }
    /** Copies from the collection to the array. The starting destination array index is specified by the {@code start} parameter. */
    public static <T> void copyTo(java.util.Collection<? extends T> list, T[] arr, int start) { copyTo(list, arr, start, list.size()); }
    /** Copies the specified number of elements from the collection to the array. The starting destination array index is specified by the {@code start} parameter. */
    public static <T> void copyTo(java.util.Collection<? extends T> list, T[] arr, int start, int length) {
        int i = 0;
        for (T c : list) {
            arr[start + i++] = c;
            if (i == length) break;
        }
    }

    public static <K,V,K2,V2> Map<K2,V2> mapMap(Map<K,V> src, Function<? super K, ? extends K2> keyMap, Function<? super V, ? extends V2> valueMap) {
        return mapMapEntries(src, e->keyMap.apply(e.getKey()), e->valueMap.apply(e.getValue()));
    }
    public static <K,V,K2,V2> Map<K2,V2> mapMapEntries(Map<K,V> src, Function<? super Map.Entry<K,V>, ? extends K2> keyMap, Function<? super Map.Entry<K,V>, ? extends V2> valueMap) {
        Map<K2,V2> dest;
        if (src instanceof TreeMap)
            dest = new TreeMap<K2, V2>();
        else if (src instanceof LinkedHashMap)
            dest = new LinkedHashMap<K2, V2>(src.size());
        else
            dest = new HashMap<>(src.size());
        return mapMapEntries(src, dest, keyMap, valueMap);
    }
    public static <K,V,K2,V2> Map<K2,V2> mapMapEntries(Map<K,V> src, Map<K2,V2> dest, Function<? super Map.Entry<K,V>, ? extends K2> keyMap, Function<? super Map.Entry<K,V>, ? extends V2> valueMap) {
        for(Map.Entry<K,V> e : src.entrySet())
            dest.put(keyMap.apply(e), valueMap.apply(e));
        return dest;
    }
    public static <K,V,K2,V2> Map<K2,V2> mapMap(Map<K,V> src, Map<K2,V2> dest, Function<? super K, ? extends K2> keyMap, Function<? super V, ? extends V2> valueMap) {
        for(Map.Entry<K,V> e : src.entrySet())
            dest.put(keyMap.apply(e.getKey()), valueMap.apply(e.getValue()));
        return dest;
    }
    public static <T,T2> List<T2> mapList(Collection<T> src, Function<T,T2> converter) {
        List<T2> dest;
        if (src instanceof LinkedList)
            dest = new LinkedList<>();
        else
            dest = new ArrayList<>(src.size());
        return mapList(src, dest, converter);
    }
    public static <T,T2> List<T2> mapList(Collection<T> src, List<T2> dest, Function<? super T,? extends T2> converter) {
        for(T t : src)
            dest.add(converter.apply(t));
        return dest;
    }


    @SuppressWarnings("unchecked")
    public static <T> ArrayList<T> safeClone(ArrayList<T> list) {
        return (ArrayList<T>)list.clone();
    }
    private static java.lang.reflect.Method cloneMethod = getCloneMethod();
    private static Method getCloneMethod() {
        try {
            return Object.class.getDeclaredMethod("clone");
        }catch (NoSuchMethodException ex) {
            for (java.lang.reflect.Method m : Object.class.getDeclaredMethods())
                AppLog.getDefault().warn("Found Method: " + m.getName());
            //should never happen because Object has a method called clone.
            throw new InternalError(ex);
        }
    }
    @SuppressWarnings("unchecked")
    public static <T> T safeClone(T obj) {
        try {
            return (T) cloneMethod.invoke(obj);
        } catch (IllegalAccessException | InvocationTargetException ex) {
            return  null;
        }
    }
    @SuppressWarnings("unchecked")
    public static <T> T firstOfType(final Object[] arr, final Class<T> type) {
        Objects.requireNonNull(type);
        Objects.requireNonNull(arr);
        for(Object item : arr) {
            if (item != null && type.isInstance(item))
                return (T)item;
        }
        return null;
    }

    public static <T> boolean contains(final T[] subject, final T find) {
        return indexOf(subject, find)!=-1;
    }
    public static <T> boolean contains(final T[] subject, final Predicate<T> test) {
        return indexOf(subject, test)!=-1;
    }
    public static <T> int indexOf(final T[] subject, final Predicate<T> test) {
        for (int i = 0; i < subject.length; i++)
            if (test.test(subject[i]))
                return i;
        return -1;
    }
    public static <T> int indexOf(final T[] subject, final T find) {
        if (find == null) {
            for (int i = 0; i < subject.length; i++)
                if (subject[i] == null) return i;
        } else {
            for (int i = 0; i < subject.length; i++)
                if (find.equals(subject[i]))
                    return i;
        }
        return -1;
    }

    public static int indexOf(final String[] subject, final String find) { return indexOf(subject, find, false); }
    public static int indexOf(final String[] subject, final String find, boolean ignoreCase) {
        if (find == null) {
            for (int i = 0; i < subject.length; i++)
                if (subject[i] == null) return i;
        } else if (ignoreCase) {
            for (int i = 0; i < subject.length; i++)
                if (find.equalsIgnoreCase(subject[i]))
                    return i;
        } else {
            for (int i = 0; i < subject.length; i++)
                if (find.equals(subject[i]))
                    return i;
        }
        return -1;
    }

    public static <K, V> Map<K, V> toMap(K[] keys, V[] values) {
        Map<K, V> map = new LinkedHashMap<K, V>(keys.length);
        for (int i = 0; i < keys.length; i++) {
            map.put(keys[i], values[i]);
        }
        return map;
    }
    /** Returns true if the test function returns true for any of the items in the iterable. */
    public static <T> boolean any(final Iterable<T> items, final Predicate<T> test) {
        for (T item : items)
            if (test.test(item)) return true;
        return false;
    }

    /** Returns true if the test function returns true for all of the items in the iterable. */
    public static <T> boolean all(final Iterable<T> items, final Predicate<T> test) {
        for (T item : items)
            if (!test.test(item)) return false;
        return true;
    }
    public static <T> T first(final Iterable<T> list) {
        //noinspection LoopStatementThatDoesntLoop
        for(T t : list)
            return t;
        throw new IndexOutOfBoundsException("The list has no first element.");
    }
    public static <T> T first(final Iterable<T> list, T valueIfEmpty) {
        //noinspection LoopStatementThatDoesntLoop
        for(T t : list)
            return t;
        return valueIfEmpty;
    }

    public static <T> T maxOf(T[] values, ToIntFunction<T> eval) {
        if (values.length==0) throw new IllegalArgumentException("The list of values passed to max() must not be empty.");
        T vmax = values[0];
        int max = eval.applyAsInt(vmax);
        for(int i=1; i<values.length;i++) {
            int result = eval.applyAsInt(values[i]);
            if (result > max) {
                max = result;
                vmax = values[i];
            }
        }
        return vmax;
    }
    public static <T> int max(T[] values, ToIntFunction<T> eval) {
        if (values.length==0) throw new IllegalArgumentException("The list of values passed to max() must not be empty.");
        int max = eval.applyAsInt(values[0]);
        for(int i=1; i<values.length;i++) {
            int result = eval.applyAsInt(values[i]);
            if (result > max) max = result;
        }
        return max;
    }
    public static <T> T minOf(T[] values, ToIntFunction<T> eval) {
        if (values.length==0) throw new IllegalArgumentException("The list of values passed to min() must not be empty.");
        T vmin = values[0];
        int min = eval.applyAsInt(vmin);
        for(int i=1; i<values.length;i++) {
            int result = eval.applyAsInt(values[i]);
            if (result < min) {
                min = result;
                vmin = values[i];
            }
        }
        return vmin;
    }
    public static <T> int min(T[] values, ToIntFunction<T> eval) {
        if (values.length==0) throw new IllegalArgumentException("The list of values passed to min() must not be empty.");
        int min = eval.applyAsInt(values[0]);
        for(int i=1; i<values.length;i++) {
            int result = eval.applyAsInt(values[i]);
            if (result < min) min = result;
        }
        return min;
    }

    /** Use to make up for the lack of pass-by-reference in Java */
    public static class RefInt {
        public RefInt(){}
        public RefInt(int init){ value = init; }
        public int value;
    }
    /** Use to make up for the lack of pass-by-reference in Java */
    public static class RefBool {
        public RefBool(){}
        public RefBool(boolean init){ value = init; }
        public boolean value;
    }
    /** Use to make up for the lack of pass-by-reference in Java */
    public static class RefStr {
        public RefStr(){}
        public RefStr(String init){ value = init; }
        public String value;
    }

    public static List<Double> toList(double[] arr) { return Arrays.asList(toObjArr(arr));}
    public static List<Float> toList(float[] arr) { return Arrays.asList(toObjArr(arr));}
    public static List<Integer> toList(int[] arr) { return Arrays.asList(toObjArr(arr));}

    public static Double[] toObjArr(double[] arr) {
        Double[] narr = new Double[arr.length];
        for (int i = 0; i < arr.length; i++) {
            narr[i] = arr[i];
        }
        return narr;
    }
    public static Float[] toObjArr(float[] arr) {
        Float[] narr = new Float[arr.length];
        for (int i = 0; i < arr.length; i++) {
            narr[i] = arr[i];
        }
        return narr;
    }
    public static Integer[] toObjArr(int[] arr) {
        Integer[] narr = new Integer[arr.length];
        for (int i = 0; i < arr.length; i++) {
            narr[i] = arr[i];
        }
        return narr;
    }
}
