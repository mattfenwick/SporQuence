import model
import unittest




def getCodons(bases):
    '''[Base] -> [Codons]'''
    if len(bases) % 3 != 0:
        raise ValueError("number of bases must be divisible by 3")

    codons = []
    x = 0
    while x < len(bases):
        codons.append(model.Codon(bases[x:x + 3]))
        x += 3

    assert len(bases) / 3 == len(codons)

    return codons


STARTS = ["GTG", "ATG", "TTG"]
STOPS = ["TAA", "TAG", "TGA"]

def getOrfs(codons):
    '''[Codons] -> [ORFs]'''
    inORF = False
    orfs = []
    ix = 0
    for c in codons:
        if not inORF and c.bases in STARTS:
            orfs.append(model.ORF(ix, [], None))
            inORF = True
        if inORF and c.bases in STOPS:
            inORF = False
            orfs[-1].stop = c
        if inORF:
            orfs[-1].codons.append(c)
        ix += 1
    return orfs


############################

class AnalyzeTest(unittest.TestCase):

    def setUp(self):
        pass

    @unittest.expectedFailure
    def testGetCodonsBadLength(self):
        bases = 'ACGTCCTG'
        codons = getCodons(bases)

    def testGetCodons(self):
        bases = 'ACGTCCTGA'
        codons = getCodons(bases)
        self.assertEqual(codons[0].bases, "ACG")
        self.assertEqual(codons[2].bases, "TGA")

    def testGetOrfs(self):
        bases = 'ACTGTGACCTCATATTAGGGGTTT'
        codons = getCodons(bases)
        orfs = getOrfs(codons)
        self.assertEqual(4, len(orfs[0].codons))
        self.assertEqual("TAG", orfs[0].stop.bases)
        self.assertEqual(1, orfs[0].index, "ORF indices are 0-indexed!")

    def testGetOrfsTwoStarts(self):
        codons = getCodons('ACTGTGACCTCATATCAGGGGTTT' + 'ACTGTGACCTCATATTAGGGGTTT' + 'ACTTTGACCTCATATTGGTGAATT')
        orfs = getOrfs(codons)
        self.assertEqual(2, len(orfs))
        self.assertEqual(1, orfs[0].index)
        self.assertEqual(12, len(orfs[0].codons))
        self.assertEqual(17, orfs[1].index)
        self.assertEqual(5, len(orfs[1].codons))

    def testGetOrfsTwoStops(self):
        # 1, start, 7, start, 4, stop, 2, stop, 3, start, stop, 1 -> 2 ORFs
        bases = ['ACT', 'GTG', 'ACCTCATATCAGGGGTTTACT', 'GTG', 'ACCTCATATCCC', 
                 'TGA', 'CCCACT', 'TAA', 'ACCTCATAT', 'TTG', 'TGA', 'ATT']
        codons = getCodons(''.join(bases))
        orfs = getOrfs(codons)
        self.assertEqual(2, len(orfs))
        self.assertEqual(1, orfs[0].index)
        self.assertEqual(13, len(orfs[0].codons))
        self.assertEqual(21, orfs[1].index)
        self.assertEqual(1, len(orfs[1].codons))

    def testGetOrfsNoEndEnd(self):
        bases = ['ACC', 'GTG', 'CAC', 'TTT']
        orfs = getOrfs(getCodons(''.join(bases)))
        self.assertEqual(3, len(orfs[0].codons))
        self.assertEqual(None, orfs[0].stop)
    


def getSuite():
    suite1 = unittest.TestLoader().loadTestsFromTestCase(AnalyzeTest)
    return unittest.TestSuite([suite1])
