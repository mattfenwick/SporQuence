import finder
import filterer as fo
import model
import json
import unittest



def loadOrfAnals(path):
    infile = open(path, 'r')
    orfs = json.loads(infile.read())
    infile.close()
    return model.OACollection([model.OrfAnalysis(o) for o in orfs])
    
    
def findAllOrfs():
    return [x.toJSONObject(100) for x in finder.getAllMediumOrfs()]


if __name__ == "__main__":
    print json.dumps(findAllOrfs())
    
    
