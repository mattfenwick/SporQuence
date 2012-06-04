import unittest
import bacillussubtilis168 as bs
import analysismodel as am
import yaml


######################


muncher = am.SequenceMuncher(am.Sequence(bs.sequence), kdRadius = 9, peakRadius = 8, peakHeightCutoff = 1.5)

f, r = muncher.getPromotedOrfs()['forward'], muncher.getPromotedOrfs()['reverse']

fs = [[range(len(o.getKyteDool(9))), o.getKyteDool(9), o.getTriangleKyteDool(9)] for o in f]
rs = [[range(len(o.getKyteDool(9))), o.getKyteDool(9), o.getTriangleKyteDool(9)] for o in r]



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
