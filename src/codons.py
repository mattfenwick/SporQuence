import model
import unittest




def makeCodons(bases):
    '''[Base] -> [Codon]'''
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

def getOrfEnds(codons):
    '''[Codon] -> [(Int, Int)]'''
    inORF, orfEnds = False, []
    ix, start = 0, None

    for c in codons:
        if not inORF and c.bases in STARTS:
            start = ix
            inORF = True
        if inORF and c.bases in STOPS:
            inORF = False
            orfEnds.append((start, ix))
        ix += 1

    if inORF:
        orfEnds.append((start, None))

    return orfEnds


############################

class CodonsTest(unittest.TestCase):

    def setUp(self):
        pass

    @unittest.expectedFailure
    def testGetCodonsBadLength(self):
        bases = 'ACGTCCTG'
        codons = makeCodons(bases)

    def testGetCodons(self):
        bases = 'ACGTCCTGA'
        codons = makeCodons(bases)
        self.assertEqual(codons[0].bases, "ACG")
        self.assertEqual(codons[2].bases, "TGA")


class OrfsTest(unittest.TestCase):

    def setUp(self):
        pass

    def testGetOrfs(self):
        bases = 'ACTGTGACCTCATATTAGGGGTTT'
        codons = makeCodons(bases)
        orfEnds = getOrfEnds(codons)
        self.assertEqual(4, orfEnds[0][1] - orfEnds[0][0])
        self.assertEqual(1, orfEnds[0][0], "remember: ORF indices are 0-indexed!")

    def testGetOrfsTwoStarts(self):
        codons = makeCodons('ACTGTGACCTCATATCAGGGGTTT' + 'ACTGTGACCTCATATTAGGGGTTT' + 'ACTTTGACCTCATATTGGTGAATT')
        orfEnds = getOrfEnds(codons)
        self.assertEqual(2, len(orfEnds), 'total number of orfs')
        self.assertEqual(1, orfEnds[0][0], 'start of 1st orf')
        self.assertEqual(12, orfEnds[0][1] - orfEnds[0][0], 'length of 1st orf')
        self.assertEqual(17, orfEnds[1][0], 'start of 2nd orf')
        self.assertEqual(5, orfEnds[1][1] - orfEnds[1][0], 'length of 2nd orf')

    def testGetOrfsTwoStops(self):
        # 1, start, 7, start, 4, stop, 2, stop, 3, start, stop, 1 -> 2 ORFs
        bases = ['ACT', 'GTG', 'ACCTCATATCAGGGGTTTACT', 'GTG', 'ACCTCATATCCC', 
                 'TGA', 'CCCACT', 'TAA', 'ACCTCATAT', 'TTG', 'TGA', 'ATT']
        codons = makeCodons(''.join(bases))
        orfs = getOrfEnds(codons)
        self.assertEqual(2, len(orfs), 'total number of orfs')
        self.assertEqual(1, orfs[0][0], 'start of 1st orf')
        self.assertEqual(13, orfs[0][1] - orfs[0][0], 'length of 1st orf')
        self.assertEqual(21, orfs[1][0], 'start of 2nd orf')
        self.assertEqual(1, orfs[1][1] - orfs[1][0], 'length of 2nd orf')

    def testGetOrfsNoEndEnd(self):
        bases = ['ACC', 'GTG', 'CAC', 'TTT']
        orfs = getOrfEnds(makeCodons(''.join(bases)))
        self.assertEqual((1, None), orfs[0])
    

testClasses = [CodonsTest, OrfsTest]
