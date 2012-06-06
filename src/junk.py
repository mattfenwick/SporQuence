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
        'forward': [x.toJSONObject() for x in medForFs],
        'reverse': [y.toJSONObject() for y in medRevFs]   
    }


if __name__ == "__main__":
    print json.dumps(findOrfs())
    
    


