import unittest

import model
import finder
import sequence
import bacillussubtilis168 as bs
import kd
import translate as tr
import peaks
import filterer



def getSuite(mod):
    subSuites = [unittest.TestLoader().loadTestsFromTestCase(c) for c in mod.testClasses]
    return unittest.TestSuite(subSuites)


def runTests():
    testModules = [model, sequence, bs, kd, tr, peaks, filterer, finder]
    mySuite = unittest.TestSuite([getSuite(module) for module in testModules])
    unittest.TextTestRunner(verbosity=2).run(mySuite)
#    unittest.main()

if __name__ == "__main__":
    runTests()
