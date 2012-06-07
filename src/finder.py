import bacillussubtilis168 as bs
import model
import unittest



def findOrfs():
    ''' () -> ([Orf], [Orf]) '''
    seq = model.Sequence(bs.bases, True)
    return [seq.getOrfs(), seq.getReverseComplement().getOrfs()]


def findAllOrfs():
    ''' () -> [Orf] '''
    seq = model.Sequence(bs.bases, True)
    return seq.getAllOrfs() + seq.getReverseComplement().getAllOrfs()
    


def filterByLength(orfs, low = 50, high = 80):
    return [orf for orf in orfs if low <= (len(orf.bases) / 3) <= high]


def getMediumOrfs():
    f, r = findOrfs()
    
    medFors = filterByLength(f)
    medRevs = filterByLength(r)
    
    return [medFors, medRevs]


def getAllMediumOrfs():
    '''() -> [Orf]'''
    orfs = findAllOrfs()
    return filterByLength(orfs)



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
        
    def testNumAllOrfs(self):
        orfs = findAllOrfs()
        self.assertEqual(561841, len(orfs))
        
        seq, rev = bs.bases, model.Sequence(bs.bases, True).getReverseComplement().getBases()
        starts = set(["ATG", "TTG", "CTG", "GTG"])
        
        i, j, ct = 0, 0, 0
        while i < len(seq) - 2:
            if seq[i:i + 3] in starts:
                ct += 1
            i += 1
            
        while j < len(rev) - 2:
            if rev[j:j + 3] in starts:
                ct += 1
            j += 1
            
        for x in [seq[-2:] + seq[:1], seq[-1:] + seq[:2], rev[-2:] + rev[:1], rev[-1:] + rev[:2]]:
            if x in starts:
                ct += 1
                
        self.assertEqual(ct, len(orfs))
        self.assertEqual(ct, 561841)
        

    def testFilteredAllOrfs(self):
        orfs = getAllMediumOrfs()
        self.assertEqual(42302, len(orfs))
        

    
testClasses = [OrfFinderTest]
