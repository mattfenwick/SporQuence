import bacillussubtilis168 as bs
import model
import unittest



def findOrfs():
    seq = model.Sequence(bs.bases, True)
    return [seq.getOrfs(), seq.getReverseComplement().getOrfs()]


def filterByLength(orfs, low = 50, high = 80):
    return [orf for orf in orfs if low <= (len(orf.bases) / 3) <= high]


def getMediumOrfs():
    f, r = findOrfs()
    
    medFors = filterByLength(f)
    medRevs = filterByLength(r)
    
    return [medFors, medRevs]



########################################################
# unit tests
########################################################

class OrfFinderTest(unittest.TestCase):

    def setUp(self):
        pass

    def testNumOrfs(self):
        orfs = findOrfs()
        self.assertEqual(189309, sum(map(len, orfs)))

    def testFilteredOrfs(self):
        orfs = map(filterByLength, findOrfs())
        self.assertEqual(12220, sum(map(len, orfs)))

    
testClasses = [OrfFinderTest]
