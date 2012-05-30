import unittest
import analyze as an
import bacillussubtilis168 as bs


codons = an.getCodons(bs.sequence)

orfs = an.getOrfs(codons)

smallOrfs = filter(lambda o: len(o.codons) >= 50 and len(o.codons) <= 80, orfs)


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
    


def getSuite():
    suite1 = unittest.TestLoader().loadTestsFromTestCase(GlueTest)
    return unittest.TestSuite([suite1])
