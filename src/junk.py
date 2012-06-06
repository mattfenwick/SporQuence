import finder
import filterOrfs as fo
import json
import unittest



def loadOrfs(path):
    infile = open(path, 'r')
    orfs = json.loads(infile.read())
    infile.close()
    return orfs


def getOrfAnalyses(orfs):
    return {
        'forward': fo.ORFCollection([fo.ORFAnalysis(o) for o in orfs['forward']]),
        'reverse': fo.ORFCollection([fo.ORFAnalysis(o) for o in orfs['reverse']])
    }


def findOrfs():
    [medForFs, medRevFs] = finder.getMediumOrfs()
    
    return {
        'forward': [x.toJSONObject() for x in medForFs],
        'reverse': [y.toJSONObject() for y in medRevFs]   
    }


if __name__ == "__main__":
    print json.dumps(findOrfs())
    
    


