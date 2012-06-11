import finder
import filterer
import model
import unittest



# upstream, sequence, downstream, start index (0-indexed from forward strand), whether Orf is in forward strand
positives = [
    ('ATAAAG', 'ATGAAA', 'TAATTC', 3697805, False),
    ('TTTTTT', 'ATGCTT', 'TAGTTT', 1907012, True),
    ('GTCAGC', 'ATGATC', 'TGAACA', 847481,  False),
    ('TGAAGC', 'ATGAAG', 'TAAAGG', 419983,  False),
    ('GATTCT', 'ATGAAA', 'TGAATG', 1120827, False),
    ('TCGTAG', 'ATGTAT', 'TAGGCT', 738416,  False)
]

        

class ControlTest(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def test_find_positives(self):
        orfs = finder.get_all_medium_orfs(100)
        oas = model.OrfCollection(orfs)
        for (up, seq, down, ix, is_sense) in positives:
            print "trying", up, seq, down, ix
            matches = oas.filter(filterer.matchBases(seq=seq, up=up, down=down)).get_orfs()
            self.assertEqual(len(matches), 1)
            self.assertEqual(matches[0].start, ix)
            self.assertEqual(matches[0].is_sense, is_sense)
        
        
testClasses = [ControlTest]