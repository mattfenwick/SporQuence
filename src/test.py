import unittest
import model
import analyze as am
#import bacillussubtilis168 as bs
import glue
import kd


def runTests():
    testModules = [am, glue, kd]#, bs]
    mySuite = unittest.TestSuite([module.getSuite() for module in testModules])
    unittest.TextTestRunner(verbosity=2).run(mySuite)
#    unittest.main()

if __name__ == "__main__":
    runTests()
