import unittest

import model
import analyze as am
import bacillussubtilis168 as bs
import glue
import kd
import translate as tr


def getSuite(mod):
    subSuites = [unittest.TestLoader().loadTestsFromTestCase(c) for c in mod.testClasses]
    return unittest.TestSuite(subSuites)


def runTests():
    testModules = [am, glue, kd, tr, bs]
    mySuite = unittest.TestSuite([getSuite(module) for module in testModules])
    unittest.TextTestRunner(verbosity=2).run(mySuite)
#    unittest.main()

if __name__ == "__main__":
    runTests()
