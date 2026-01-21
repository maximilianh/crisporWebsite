package ur_rna.Utilities;

import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

/**
 * @author Richard M. Watson
 */
public class AppInfoTest {
    @Test
    public void testGetMainClass() throws Exception {
        assertNull(AppInfo.getMainClass());
    }

    @Test
    public void testBitStr() throws Exception {
        assertEquals("BANANA", getBitStr(1,6,16,22,33,34,35,38,48,54,65,66,67,70,80,86));
        assertEquals("X# 5~`,./\"9)Hhy_+=\\|][{':", getBitStr(3,4,6,16,17,21,37,48,50,52,53,65,66,67,68,69,70,85,86,98,99,101,113,114,115,117,128,129,130,131,133,145,149,160,163,164,165,176,179,181,195,198,211,213,214,224,227,228,229,230,240,241,242,243,244,246,256,257,259,261,272,274,275,276,277,290,291,292,294,306,307,308,309,310,320,322,323,324,326,336,337,339,340,342,352,353,355,356,357,358,368,369,370,373,385,387,388,389));
        char[] c = new char[5000/16+1];
        c[c.length-1] = 1<<(5000 - 16*(5000/16));
        assertEquals(new String(c), getBitStr(5000));
    }

    static char[] _selStrBuilder = new char[16];
    final static int BITS_PER_CHAR = 16; // 0xFFFF => [0000 0000 0000 0000]
    private static String getBitStr(int... values) {
        char[] arr = _selStrBuilder;
        int used = _selStrBuilder.length; while (used-- > 0) arr[used]=0; // clear the array

        for (int index : values) {
            int pos = index / BITS_PER_CHAR;
            int val = index - (pos * BITS_PER_CHAR);
            if (pos > used) {
                used = pos;
                if (used >= arr.length)
                    arr = _selStrBuilder = Arrays.copyOf(arr, Math.max(arr.length * 2, MathUtil.pow2i(used+1)));
            }
            arr[pos] |= 1<<val;
        }
        return new String(arr, 0, used+1);
    }
}