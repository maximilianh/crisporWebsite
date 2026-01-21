package ur_rna.StructureEditor.services.fileIO;

import org.junit.Before;
import org.junit.Test;
import ur_rna.RNAstructure.backend.RNA;
import ur_rna.RNAstructure.backend.RNATest;
import ur_rna.StructureEditor.models.RnaScene;
import ur_rna.StructureEditor.services.fileIO.RnaFileIO.NucIndexConverter;

import static org.junit.Assert.assertEquals;
import static ur_rna.Utilities.testing.TestUtils.assertError;
import static ur_rna.Utilities.testing.TestUtils.assertInstanceOf;

/**
 * Test the NucIndexConverter in RNAFileIO
 */
public class NucIndexConverterTest {
    @Before
    public void setUp() throws Exception {
        RNATest.loadNativeLib();
    }

    @Test
    public void TestCreate() throws Exception {
        assertInstanceOf(NucIndexConverter.SingleStrand.getClass(), NucIndexConverter.create("AUGCAUGCAUGC"));
        assertInstanceOf(NucIndexConverter.DoubleStrand.class, NucIndexConverter.create("AUGCAUIIIGCAUGC"));
        assertInstanceOf(NucIndexConverter.MultiStrand.class, NucIndexConverter.create("AUGIIICAUGCIIIAUGC"));
    }

    @Test
    public void testSingleStrand() throws Exception {
        final String seq = "AUGCAUGCAUGC";
        RNA rna = new RNA(seq, true);
        RnaScene scene = RnaFileIO.createSceneFromRNA(rna, 1, false);
        assertEquals(seq, rna.GetSequence());
        assertEquals(seq.length(), scene.getNucCount());
        assertEquals(seq.length(), scene.strands.get(0).size());
        NucIndexConverter n = NucIndexConverter.create(rna.GetSequence());
        for(int i = 0; i < seq.length(); i++) {
            assertEquals(scene.strands.get(0).get(i), n.getSceneNuc(scene, i + 1));
            assertEquals(i+1, n.getBackendIndex(scene.strands.get(0).get(i)));
        }
        assertError(()->n.getSceneNuc(scene, seq.length()+1), IndexOutOfBoundsException.class);
        assertError(()->n.getSceneNuc(scene, -1), IndexOutOfBoundsException.class);
    }
    @Test
    public void testDoubleStrand() throws Exception {
        final String seq = "AUG"+"III"+"TAUGCAUGC";
        RNA rna = new RNA(seq, true);
        RnaScene scene = RnaFileIO.createSceneFromRNA(rna, 1, false);
        assertEquals(seq, rna.GetSequence());
        assertEquals(seq.length()-3, scene.getNucCount());
        assertEquals(3, scene.strands.get(0).size());
        assertEquals(9, scene.strands.get(1).size());
        assertEquals(12, scene.getNucCount());
        NucIndexConverter n = NucIndexConverter.create(rna.GetSequence());
        for(int i = 0; i < 3; i++) {
            assertEquals(scene.strands.get(0).get(i), n.getSceneNuc(scene, i + 1));
            assertEquals(i + 1, n.getBackendIndex(scene.strands.get(0).get(i)));
        }
        assertEquals('I', rna.GetNucleotide(6));
        assertEquals('T', rna.GetNucleotide(7));
        assertEquals('T', scene.strands.get(1).get(0).symbol.charAt(0));
        for(int i = 0; i < 9; i++) {
            assertEquals(scene.strands.get(1).get(i), n.getSceneNuc(scene, i + 7));
            assertEquals(i + 7, n.getBackendIndex(scene.strands.get(1).get(i)));
        }
        assertEquals(scene.strands.get(1).get(scene.strands.get(1).size()-1), n.getSceneNuc(scene, seq.length())); //verify last
        for(int i: new int[]{4,5,6}) {
            assertEquals('I', rna.GetNucleotide(i));
            assertError(() -> n.getSceneNuc(scene, i), IndexOutOfBoundsException.class);
        }
        assertError(()->n.getSceneNuc(scene, seq.length()+1), IndexOutOfBoundsException.class);
        assertError(()->n.getSceneNuc(scene, -1), IndexOutOfBoundsException.class);
    }

    @Test
    public void testMultiStrand() throws Exception {
        final String seq = "AUG"+"III"+"TAUGCAUGC"+"III"+"uccGX";
        RNA rna = new RNA(seq, true);
        RnaScene scene = RnaFileIO.createSceneFromRNA(rna, 1, false);
        assertEquals(seq, rna.GetSequence());
        assertEquals(seq.length()-6, scene.getNucCount());
        assertEquals(3, scene.strands.get(0).size());
        assertEquals(9, scene.strands.get(1).size());
        assertEquals(5, scene.strands.get(2).size());
        assertEquals(17, scene.getNucCount());
        NucIndexConverter n = NucIndexConverter.create(rna.GetSequence());
        for(int i = 0; i < 3; i++) {
            assertEquals(scene.strands.get(0).get(i), n.getSceneNuc(scene, i + 1));
            assertEquals(i + 1, n.getBackendIndex(scene.strands.get(0).get(i)));
        }

        for(int i: new int[]{4,5,6,16,17,18}) {
            assertEquals('I', rna.GetNucleotide(i));
            assertError(()->n.getSceneNuc(scene, i), IndexOutOfBoundsException.class);
        }

        assertEquals('T', rna.GetNucleotide(7));
        assertEquals('T', scene.strands.get(1).get(0).symbol.charAt(0));
        for(int i = 0; i < 9; i++) {
            assertEquals(scene.strands.get(1).get(i), n.getSceneNuc(scene, i + 7));
            assertEquals(i + 7, n.getBackendIndex(scene.strands.get(1).get(i)));
        }

        assertEquals('u', rna.GetNucleotide(19));
        assertEquals('u', scene.strands.get(2).get(0).symbol.charAt(0));
        for(int i = 0; i < 5; i++) {
            assertEquals(scene.strands.get(2).get(i), n.getSceneNuc(scene, i + 19));
            assertEquals(i + 19, n.getBackendIndex(scene.strands.get(2).get(i)));
        }

        assertEquals(scene.strands.get(2).get(scene.strands.get(2).size()-1), n.getSceneNuc(scene, seq.length())); //verify last
        assertEquals("X", n.getSceneNuc(scene, seq.length()).symbol); //verify last
        assertError(()->n.getSceneNuc(scene, seq.length()+1), IndexOutOfBoundsException.class);
        assertError(()->n.getSceneNuc(scene, -1), IndexOutOfBoundsException.class);
    }

}