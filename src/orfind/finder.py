import bacillussubtilis168 as bs
import unittest
import json
import fmodel






###############################

def findOrfs():
    seq = fmodel.Sequence(bs.sequence)
    return [seq.getOrfs(), seq.getReverseOrfs()]


def filterByLength(orfs, low = 50, high = 80):
    return [orf for orf in orfs if low <= (len(orf.bases) / 3) <= high]


if __name__ == "__main__":
    f, r = findOrfs()
    
    medFs = filterByLength(f)
    revFs = filterByLength(r)
    
    print json.dumps({
        'forward': [x.toJSONObject() for x in medFs],
        'reverse': [y.toJSONObject() for y in revFs]   
    })
    
    
########################################################
# unit tests
########################################################

class HuhTest(unittest.TestCase):

    def setUp(self):
        pass

    def testOne(self):
        self.assertEqual(0, 1)

    

testClasses = [HuhTest]

