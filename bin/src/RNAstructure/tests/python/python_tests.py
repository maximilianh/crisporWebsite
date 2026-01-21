import unittest as ut
import RNAstructure as r

# a tRNA sequence
RR1664 = 'GCGCCCGUAGCUCAGCUGGAUAGAGCGCUGCCCUCCGGAGGCAGAGGUCUCAGGUUCGAAUCCUGUCGGGCGCGCCA'
# a short hairpin sequence
short = "GGGGGAAACCCCC"

class testRNAMethods(ut.TestCase):
    def setUp(self):
        self.rna = r.RNA.fromString(RR1664)

    def test_FoldSingleStrand(self):
        self.rna.FoldSingleStrand()

    def test_PartitionFunction(self):
        self.rna.PartitionFunction()

    def test_MaxExpect(self):
        self.rna.MaximizeExpectedAccuracy()

class testProbScanMethods(ut.TestCase):
    def setUp(self):
        self.ps = r.ProbScan.fromString(RR1664)

    def test_Hairpin(self):
        self.ps.probability_of_all_hairpins()

    def test_InternalLoops(self):
        self.ps.probability_of_all_internal_loops()

    def test_Helices(self):
        self.ps.probability_of_all_helices(length=3)

class testErrorHandling(ut.TestCase):
    def setUp(self):
        self.rna = r.RNA.fromString(short)

    def test_BoundsChecking(self):
        with self.assertRaises(IndexError):
            self.rna.GetNucleotide(self.rna.GetSequenceLength()+1)
        with self.assertRaises(IndexError):
            self.rna.GetNucleotide(0)
        with self.assertRaises(IndexError):
            self.rna.GetPair(self.rna.GetSequenceLength()+1)

    def test_ProbabilityRequiresPartitionFunction(self):
        with self.assertRaises(ValueError):
            self.rna.GetPairProbability(1,5)

    def test_PairRequiresFolding(self):
        with self.assertRaises(ValueError):
            self.rna.GetPair(1,2)

class testIterationRNA(ut.TestCase):
    def setUp(self):
        self.rna = r.RNA.fromString(short)
    def test_RNA_iterator(self):
        self.assertEqual("".join(x for x in self.rna),short)
    def test_RNA_index_iterator(self):
        self.assertEqual(list(self.rna.iterIndices()),
                         list(range(1,self.rna.GetSequenceLength()+1)))

class testIterationProbScan(ut.TestCase):
    def setUp(self):
        self.rna = r.ProbScan.fromString(short)
    def test_RNA_iterator(self):
        self.assertEqual("".join(x for x in self.rna),short)
    def test_RNA_index_iterator(self):
        self.assertEqual(list(self.rna.iterIndices()),
                         list(range(1,self.rna.GetSequenceLength()+1)))



if __name__ == "__main__":
    ut.main(verbosity=2)
"""
suite = ut.TestLoader().loadTestsFromTestCase(testRNAMethods)
ut.TextTestRunner(verbosity=2).run(suite)
"""
