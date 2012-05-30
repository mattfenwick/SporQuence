import unittest
import model


kdIndex = { 
  'A':  1.8,
  'R': -4.5,
  'N': -3.5,
  'D': -3.5,
  'C':  2.5,
  'Q': -3.5,
  'E': -3.5,
  'G': -0.4,
  'H': -3.2,
  'I':  4.5,
  'L':  3.8,
  'K': -3.9,
  'M':  1.9,
  'F':  2.8,
  'P': -1.6,
  'S': -0.8,
  'T': -0.7,
  'W': -0.9,
  'Y': -1.3,
  'V':  4.2 
}


def getResidueScore(residue):
    '''Residue -> Float'''
    if kdIndex.has_key(residue):
        return kdIndex[residue]
    raise KeyError("missing Kyte-Doolittle index for residue <%s>" % str(residue))


def smooth(values, windowRadius):
    '''[Float] -> Int -> [Float]'''
    i = windowRadius
    smoothed = []
    width = 2 * windowRadius + 1
    while i < len(values) - windowRadius:
        start, stop = i - windowRadius, i + windowRadius + 1
        total = sum(values[start: stop])
        smoothed.append(total / width)
        i += 1
    return smoothed


def kyteDoolittle(residues, windowRadius):
    '''[Residue] -> Int -> [Float]'''
    scored = [getResidueScore(r) for r in residues]
    smoothed = smooth(scored, windowRadius)
    return smoothed


############################ translate codons to residues


residueToCodons = {
  'A': ['GCT', 'GCC', 'GCA', 'GCG'],
  'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
  'N': ['AAT', 'AAC'],
  'D': ['GAT', 'GAC'],
  'C': ['TGT', 'TGC'],
  'Q': ['CAA', 'CAG'],
  'E': ['GAA', 'GAG'],
  'G': ['GGT', 'GGC', 'GGA', 'GGG'],
  'H': ['CAT', 'CAC'],
  'I': ['ATT', 'ATC', 'ATA'],
  'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
  'K': ['AAA', 'AAG'],
  'M': ['ATG'],
  'F': ['TTT', 'TTC'],
  'P': ['CCT', 'CCC', 'CCA', 'CCG'],
  'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
  'T': ['ACT', 'ACC', 'ACA', 'ACG'],
  'W': ['TGG'],
  'Y': ['TAT', 'TAC'],
  'V': ['GTT', 'GTC', 'GTA', 'GTG']
}


def flipCode(r2cs):
    '''Map Residue [Codon] -> Map Codon Residue'''
    c2r = {}
    for (res, codons) in residueToCodons.iteritems():
        for codon in codons:
            c2r[codon] = res
    return c2r


codonToResidue = flipCode(residueToCodons)


def translateCodons(codons):
    '''[Codon] -> [Residue]'''
    return [codonToResidue[c.bases] for c in codons]



############################

class KdTest(unittest.TestCase):

    def setUp(self):
        pass

    def testIndexLength(self):
        self.assertEqual(20, len(kdIndex))

    def testGetResidueScore(self):
        self.assertEqual(1, 0)

    def testSmooth(self):
        self.assertEqual(1, 0)

    def testAlgorithm(self):
        rs = 'MATTCVGHKWERTY'
        calced = kyteDoolittle(rs, 2)
        self.assertEqual(len(calced), len(rs) - 4)


class R2CsTest(unittest.TestCase):

    def setUp(self):
        self.codons = reduce(lambda x,y: x + y, residueToCodons.values())

    def testUniqueCodons(self):
        codons = self.codons
        self.assertEqual(len(codons), len(set(codons)), 'no repeats')

    def testCodonsLength(self):
        codons = self.codons
        self.assertEqual(61, len(codons))

    def testBasesLength(self):
        self.assertEqual(20, len(residueToCodons))

    def testDNABases(self):
        for b in ''.join(self.codons):
            self.assertTrue(b in "ACGT", "base %s" % str(b))


class C2RTest(unittest.TestCase):

    def setup(self):
        pass

    def testLength(self):
        self.assertEqual(len(codonToResidue), 61)

    def testDNABases(self):
        for b in ''.join(codonToResidue.keys()):
            self.assertTrue(b in "ACGT", "base %s" % str(b))

    def testTranslateCodons(self):
        rs = translateCodons([model.Codon(seq) for seq in ['ACG', 'ATG', 'TTA']])
        self.assertEqual(len(rs), 3)
        self.assertEqual(rs[1], 'M')


def getSuite():
    suite1 = unittest.TestLoader().loadTestsFromTestCase(KdTest)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(R2CsTest)
    suite3 = unittest.TestLoader().loadTestsFromTestCase(C2RTest)
    return unittest.TestSuite([suite1, suite2, suite3])
