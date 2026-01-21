package ur_rna.Utilities;

/**
 * Extended math functions.
 */
public class MathUtil {
    /**
     * Returns the smallest integer x such that 2^x is greater than or equal to the specified value.
     * e.g. {@code Pow2(37)} returns 6 (because 2^6 is the smallest power of two greater than or equal to 37)
     * @param value
     * @return An integer that represents the exponent of the smallest power of two that is greater than or equal to the specified value.
     */
    public static int log2i (double value) {
        return (int)Math.ceil(log2(value));
    }

    /** Returns the base 2 logarithm of a number. i.e. {@code log(value)/log(2)} */
    public static double log2 (double value) {
        return Math.log(value)/Math.log(2);
    }

    /**
     * Returns the smallest power of 2 greater than or equal to the specified value.
     * e.g. {@code Pow2(37)} returns 64, which is the smallest power of 2 greater than or equal to 37.
     * Note: If the return value might exceed the range of an integer,
     *       then {@link #pow2L(double)} or {@link #pow2(double)} should be used instead.
     */
    public static int pow2i (double value) {
        return 1<<log2i(value);
    }
    /**
     * Returns the smallest power of 2 greater than or equal to the specified value.
     * e.g. {@code Pow2(37)} returns 64, which is the smallest power of 2 greater than or equal to 37.
     */
    public static  long pow2L (double value) {
        return 1L<<log2i(value);
    }
    /**
     * Returns the smallest power of 2 greater than or equal to the specified value.
     * e.g. {@code Pow2(37)} returns 64  because  32 &lt; 37 &lt;= 64
     */
    public static double pow2 (double value) {
        return Math.pow(2,log2i(value));
    }
}
