import bacillussubtilis168 as bs
import sequence
import unittest



forward = bs.bases
reverse = sequence.reverseComplement(forward)



def find_once(sear, seq):
    f, r = seq.find(sear), seq.rfind(sear)
    if f == r:
        return f
    else:
        return 'No!'
    
    
    
class PetersTest(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def test_length(self):
        self.assertEqual(len(forward), len(reverse))
        
    def test_find_once(self):
        self.assertEqual('No!', find_once('A', forward))
    
    
        
testClasses = [PetersTest]