package ur_rna.StructureEditor.models;

import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

/**
 * Test StrandList
 */
public class StrandListTest {
    @Test
    public void AddToSingleStrandTest() throws Exception {
        RnaScene rna = new RnaScene();
        StrandList strands = rna.strands;

        assertEquals(1, strands.size());
        Strand s1 = strands.first();
        assertEquals(s1, strands.get(0));
        assertEquals(s1, strands.last());

        assertEquals(0, strands.allNucs().size());
        assertEquals(0, s1.getSceneStart());
        assertEquals(0, s1.size());

        s1.add();

        assertEquals(1, strands.allNucs().size());
        assertEquals(0, s1.getSceneStart());
        assertEquals(1, s1.size());

        s1.add();
        assertEquals(2, strands.allNucs().size());
        assertEquals(0, s1.getSceneStart());
        assertEquals(2, s1.size());

        verifyNucIndex(strands);
    }
    private void verifyNucIndex(final StrandList strands) {
        int ni = 0, si = 0;
        for (Strand s : strands) {
            assertEquals(si++, s.index);
            for (Nuc n : s) {
                assertEquals(ni++, n.indexInScene);
                assertEquals(s, n.strand);
            }
        }
    }

    @Test
    public void AddToTwoStrandsTest() throws Exception {
        RnaScene rna = new RnaScene();
        StrandList strands = rna.strands;

        assertEquals(1, strands.size());
        Strand s1 = strands.first();
        Strand s2 = strands.add();
        assertEquals(2, strands.size());
        assertEquals(s2, strands.second());
        assertEquals(s2, strands.last());
        assertEquals(s2, strands.get(1));
        assertEquals(s1, strands.first());
        assertEquals(s1, strands.get(0));

        for (int i = 0; i < 5; i++)
            s1.add(""+i);

        for (int i = 5; i < 8; i++)
            s2.add(""+i);

        verifyNucIndex(strands);

        assertEquals(8, strands.allNucs().size());
        assertEquals(0, s1.getSceneStart());
        assertEquals(5, s2.getSceneStart());
        assertEquals(5, s1.size());
        assertEquals(3, s2.size());
        assertEquals(5, s1.getSceneEnd());
        assertEquals(8, s2.getSceneEnd());

        for (int i = 8; i < 14; i++)
            s1.add(""+i);
        s1.add(""+14);

        verifyNucIndex(strands);


        assertEquals(15, strands.allNucs().size());
        assertEquals(0, s1.getSceneStart());
        assertEquals(12, s2.getSceneStart());
        assertEquals(12, s1.size());
        assertEquals(3, s2.size());
        assertEquals(12, s1.getSceneEnd());
        assertEquals(15, s2.getSceneEnd());
    }

    @Test
    public void DivideStrandTest() throws Exception {
        RnaScene rna = new RnaScene();
        StrandList strands = rna.strands;
        for (int i = 0; i < 15; i++)
            strands.first().add(""+i);
        verifyNucIndex(strands);

        Strand sA = strands.first();
        assertEquals(15, strands.allNucs().size());
        assertEquals(0, sA.getSceneStart());
        assertEquals(15, sA.size());
        assertEquals(15, sA.getSceneEnd());

        strands.divideStrand(sA, 5);
        verifyNucIndex(strands);

        Strand sB = strands.second();
        assertEquals(15, strands.allNucs().size());
        assertEquals(0, sA.getSceneStart());
        assertEquals(5, sA.size());
        assertEquals(5, sA.getSceneEnd());
        assertEquals(10, sB.size());
        assertEquals(5, sB.getSceneStart());
        assertEquals(15, sB.getSceneEnd());
    }

    @Test
    public void StrandInsertTest() throws Exception {
        RnaScene rna = new RnaScene();
        StrandList strands = rna.strands;
        for (int i = 0; i < 15; i++)
            strands.first().add(""+i);
        verifyNucIndex(strands);

        Strand sA = strands.first();
        strands.divideStrand(sA, 5);
        verifyNucIndex(strands);

        Strand sB = strands.second();
        Strand sC = strands.add(0);
        verifyNucIndex(strands);

        assertEquals(0, sC.index);
        assertEquals(1, sA.index);
        assertEquals(2, sB.index);
        assertEquals(sC, strands.get(0));
        assertEquals(sA, strands.get(1));
        assertEquals(sB, strands.get(2));

        assertEquals(0, sC.size());
        assertEquals(0, sC.getSceneStart());
        assertEquals(0, sA.getSceneStart());
        assertEquals(5, sB.getSceneStart());

        sC.addAll(Arrays.asList(new Nuc("C1"), new Nuc("C2"), new Nuc("C3")));
        verifyNucIndex(strands);

        assertEquals(3, sC.size());
        assertEquals(0, sC.getSceneStart());
        assertEquals(3, sA.getSceneStart());
        assertEquals(8, sB.getSceneStart());
        assertEquals(18, sB.getSceneEnd());
        assertEquals(18, strands.allNucs().size());
    }
}