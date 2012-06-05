import orfind.finder as of
import unittest
import json




if __name__ == "__main__":
    [medForFs, medRevFs] = of.getMediumOrfs()
    
    print json.dumps({
        'forward': [x.toJSONObject() for x in medForFs],
        'reverse': [y.toJSONObject() for y in medRevFs]   
    })
    
    


