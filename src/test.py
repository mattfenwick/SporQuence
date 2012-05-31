import unittest

import model
import sequence
import bacillussubtilis168 as bs
#import glue
import kd
import translate as tr
import analysismodel as am
import peaks



def getSuite(mod):
    subSuites = [unittest.TestLoader().loadTestsFromTestCase(c) for c in mod.testClasses]
    return unittest.TestSuite(subSuites)


def runTests():
    testModules = [sequence, kd, tr, am, peaks]#, bs] # glue, model]
    mySuite = unittest.TestSuite([getSuite(module) for module in testModules])
    unittest.TextTestRunner(verbosity=2).run(mySuite)
#    unittest.main()

if __name__ == "__main__":
    runTests()
