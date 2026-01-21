package ur_rna.RNAstructure.backend;

import org.junit.Before;
import org.junit.Test;
import ur_rna.StructureEditor.Program;
import ur_rna.Utilities.StopWatch;

public class RNATest {
    protected static Boolean nativeLibLoaded;

    public static void loadNativeLib() throws UnsatisfiedLinkError {
        if (nativeLibLoaded == null) {
            nativeLibLoaded = false;
            Program.loadNativeLib();
            nativeLibLoaded = true;
        } else if (!nativeLibLoaded)
            throw new UnsatisfiedLinkError("Failed to load native lib.");
    }

    @Before
    public void setUp() throws Exception {
       loadNativeLib();
    }

    @Test
    public void TestSwigArrays() throws Exception {
        final int ARR_SIZE = 1000;
        final int TRIALS = 10000;
        int[] ints = new int[ARR_SIZE];

        StopWatch s = new StopWatch(true);
        RNA rna = new RNA("ACGT", true);
        s.println("RNA startup").restart();


        rna.getInts(ints, ARR_SIZE); // load the library and call it once before the stopwatch
        s.println("ints startup").restart();

        for (int i = 0; i < TRIALS; i++) {
            rna.getInts(ints, ARR_SIZE);
        }
        s.stop().println("getInts").restart();

        rna.initInts();
        rna.getIntAt(0); // call it once before the stopwatch to load the array
        s.stop().println("initInts").restart();

        for (int i = 0; i < TRIALS; i++) {
            for (int j = 0; j < ARR_SIZE; j++) {
                rna.getIntAt(j);
            }
        }
        s.println("getIntAt");
    }


    @Test
    public void TestSwigArrays2() throws Exception {

    }
}