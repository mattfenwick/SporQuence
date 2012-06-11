import unittest
import sys

import model
import finder
import sequence
import bacillussubtilis168 as bs
import kd
import translate as tr
import peaks
import filterer
import controls
import peters
import junk



_SHORTS = [model, sequence, kd, tr, peaks, filterer, peters]

_LONGS = [junk, controls, bs, finder]



def getSuite(mod):
    subSuites = [unittest.TestLoader().loadTestsFromTestCase(c) for c in mod.testClasses]
    return unittest.TestSuite(subSuites)


def runTests(test_modules):
    mySuite = unittest.TestSuite([getSuite(module) for module in test_modules])
    unittest.TextTestRunner(verbosity=2).run(mySuite)
#    unittest.main()


if __name__ == "__main__":
    try:
        switch = sys.argv[1]
        if switch == '-s':
            test_modules = _SHORTS
        elif switch == '-a':
            test_modules = _SHORTS + _LONGS
        else:
            raise ValueError('bad command-line option')
        runTests(test_modules)
    except Exception:
        print "options are -a (run all tests) and -s (run short tests)"
        raise
