import finder
import filterer as fo
import json
import unittest



def loadOrfs(path):
    infile = open(path, 'r')
    orfs = json.loads(infile.read())
    infile.close()
    return orfs


def getOrfAnalyses(orfs):
    return {
        'forward': fo.OACollection([fo.OrfAnalysis(o) for o in orfs['forward']]),
        'reverse': fo.OACollection([fo.OrfAnalysis(o) for o in orfs['reverse']])
    }


def findOrfs():
    [medForFs, medRevFs] = finder.getMediumOrfs()
    
    return {
        'forward': [x.toJSONObject(100) for x in medForFs],
        'reverse': [y.toJSONObject(100) for y in medRevFs]   
    }
    
    
def findAllOrfs():
    return [x.toJSONObject(100) for x in finder.getAllMediumOrfs()]


if __name__ == "__main__":
    print json.dumps(findAllOrfs())
    
    
