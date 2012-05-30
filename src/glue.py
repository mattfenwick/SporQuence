import unittest
import codons as cn
import bacillussubtilis168 as bs
import kd
import translate as tr


##################### junk functions (junc-tions)

def kdFilter(kdResults):
    return len(filter(lambda x: x > 2, kdResults)) > 0

def countBig(kdResults):
    return len(filter(lambda x: x > 2, kdResults))

######################


# codons of bacillus subtilis
codons = cn.getCodons(bs.sequence)

# [ORF] :: open reading frames
orfs = cn.getOrfs(codons)

# [ORF] :: open reading frames between 50 and 80 codons long
smallOrfs = filter(lambda o: len(o.codons) >= 50 and len(o.codons) <= 80, orfs)

# [[Residue]] :: translated small ORFs
orfResidues = map(lambda orf: tr.translateCodons(orf.codons), smallOrfs)

# [[Float]] :: Kyte-Doolittles of translated small ORFs
orfKds = map(lambda o: kd.kyteDoolittle(o, 9), orfResidues)

# [[Float]] :: just the Kyte-Doolittles with a value > 2
phobicKds = filter(kdFilter, orfKds)

############################

class GlueTest(unittest.TestCase):

    def setUp(self):
        pass

    def testCodonsLength(self):
        self.assertEqual(1405202, len(codons))

    def testOrfsLength(self):
        self.assertEqual(27522, len(orfs))

    def testSmallOrfsLength(self):
        self.assertEqual(1782, len(smallOrfs))
    


testClasses = [GlueTest]
