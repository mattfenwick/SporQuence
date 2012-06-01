import unittest
import bacillussubtilis168 as bs
import analysismodel as am
import yaml


######################


muncher = am.SequenceMuncher(am.Sequence(bs.sequence), peakRadius = 8, peakHeightCutoff = 2.5)

f, r = muncher.getPromotedOrfs()['forward'], muncher.getPromotedOrfs()['reverse']


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
