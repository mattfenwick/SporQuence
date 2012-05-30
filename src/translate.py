import model
import unittest


# translate codons to residues


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


def codonsToResidues(codons):
    '''[Codon] -> [Residue]'''
    return [codonToResidue[c.bases] for c in codons]


##########################################


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

    def testCodonsToResidues(self):
        rs = codonsToResidues([model.Codon(seq) for seq in ['ACG', 'ATG', 'TTA']])
        self.assertEqual(len(rs), 3)
        self.assertEqual(rs[1], 'M')


testClasses = [R2CsTest, C2RTest]
