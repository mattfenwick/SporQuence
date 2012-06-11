import bacillussubtilis168 as bs
import model
import unittest



def find_orfs(n):
    ''' () -> ([Orf], [Orf]) '''
    seq = model.Sequence(bs.bases, True)
    return [seq.get_orfs(n), seq.get_reverse_complement().get_orfs(n)]


def find_all_orfs(n):
    ''' () -> [Orf] '''
    seq = model.Sequence(bs.bases, True)
    return seq.get_all_orfs(n) + seq.get_reverse_complement().get_all_orfs(n)
    


def filter_by_length(orfs, low = 50, high = 80):
    return [orf for orf in orfs if low <= (len(orf.bases) / 3) <= high]


def get_medium_orfs(n):
    f, r = find_orfs(n)
    
    med_Fors = filter_by_length(f)
    med_Revs = filter_by_length(r)
    
    return [med_Fors, med_Revs]


def get_all_medium_orfs(n):
    '''() -> [Orf]'''
    orfs = find_all_orfs(n)
    return filter_by_length(orfs)



########################################################
# unit tests
########################################################

class OrfFinderTest(unittest.TestCase):

    def setUp(self):
        pass

    def testNumOrfs(self):
        orfs = find_orfs(5)
        self.assertEqual(189309, sum(map(len, orfs)))

    def testFilteredOrfs(self):
        orfs = map(filter_by_length, find_orfs(5))
        self.assertEqual(12220, sum(map(len, orfs)))
        
    def testNumAllOrfs(self):
        orfs = find_all_orfs(5)
        self.assertEqual(561841, len(orfs))
        
        seq, rev = bs.bases, model.Sequence(bs.bases, True).get_reverse_complement().get_bases()
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
        orfs = get_all_medium_orfs(5)
        self.assertEqual(42302, len(orfs))
        

    
testClasses = [OrfFinderTest]
