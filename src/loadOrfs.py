import json



def getOrfs(path):
    infile = open(path, 'r')
    orfs = json.loads(infile.read())
    infile.close()
    return orfs