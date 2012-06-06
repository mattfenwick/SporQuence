import unittest

import fmodel
#import orfind.finder as fi
#import orfind.sequence as sq
#import orfind.bacillussubtilis168 as bs
#import orfind.model as ormodel
#import glue
#import kd
#import translate as tr
#import analysismodel as am
#import peaks



def getSuite(mod):
    subSuites = [unittest.TestLoader().loadTestsFromTestCase(c) for c in mod.testClasses]
    return unittest.TestSuite(subSuites)


def runTests():
    testModules = [fmodel]#sq, fi, ormodel, bs]#[sequence, kd, tr, am, peaks]#, bs] # glue, model]
    mySuite = unittest.TestSuite([getSuite(module) for module in testModules])
    unittest.TextTestRunner(verbosity=2).run(mySuite)
#    unittest.main()

if __name__ == "__main__":
    runTests()
