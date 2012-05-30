import unittest
import bacillussubtilis168 as bs
import analysismodel as am


##################### junk functions (junc-tions)

def kdFilter(kdResults):
    return len(filter(lambda x: x > 2, kdResults)) > 0

def countBig(kdResults):
    return len(filter(lambda x: x > 2, kdResults))

######################


# codons of bacillus subtilis
seq = am.Sequence(bs.sequence)

# [ORF] :: open reading frames
orfs = seq.getOrfs()

# [ORF] :: open reading frames between 50 and 80 codons long
smallOrfs = filter(lambda o: len(o.getCodons()) >= 50 and len(o.getCodons()) <= 80, orfs)

# [ORF] :: just those where Kyte-Doolittles has at least one value > 2
phobicOrfs = filter(lambda o: kdFilter(o.getKyteDool(9)), smallOrfs)

############################

class GlueTest(unittest.TestCase):

    def setUp(self):
        pass

    def testCodonsLength(self):
        self.assertEqual(1405202, len(seq.getCodons()))

    def testOrfsLength(self):
        self.assertEqual(27522, len(orfs))

    def testSmallOrfsLength(self):
        self.assertEqual(1782, len(smallOrfs))
    


testClasses = [GlueTest]
