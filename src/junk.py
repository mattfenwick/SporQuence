import finder
import model
import json
import unittest



def load_orfs(path):
    infile = open(path, 'r')
    orfs = json.loads(infile.read())
    infile.close()
    return model.OrfCollection([model.Orf.from_JSON_object(o) for o in orfs])
    
    
def find_all_medium_orfs():
    return [x.to_JSON_object() for x in finder.get_all_medium_orfs(100)]


if __name__ == "__main__":
    print json.dumps(find_all_medium_orfs())
    
    
    
class JunkTest(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def test_find(self):
        self.assertEqual(len(find_all_medium_orfs()), 42302)
        
    def test_load(self):
        path = '../allmediumorfs.txt'
        orfColl = load_orfs(path)
        self.assertEqual(len(orfColl.get_orfs()), 42302)
        
        
        
testClasses = [JunkTest]
    
    
