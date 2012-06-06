import json
import filterOrfs as fo



def getOrfs(path):
    infile = open(path, 'r')
    orfs = json.loads(infile.read())
    infile.close()
    return orfs


def getOrfAnalyses(orfs):
    return {
        'forward': fo.ORFCollection([fo.ORFAnalysis(o) for o in orfs['forward']]),
        'reverse': fo.ORFCollection([fo.ORFAnalysis(o) for o in orfs['reverse']])
    }
    