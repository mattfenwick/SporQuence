import unittest
import model


kdIndex = { 
  'A':  1.8,
  'R': -4.5,
  'N': -3.5,
  'D': -3.5,
  'C':  2.5,
  'Q': -3.5,
  'E': -3.5,
  'G': -0.4,
  'H': -3.2,
  'I':  4.5,
  'L':  3.8,
  'K': -3.9,
  'M':  1.9,
  'F':  2.8,
  'P': -1.6,
  'S': -0.8,
  'T': -0.7,
  'W': -0.9,
  'Y': -1.3,
  'V':  4.2 
}


def getResidueScore(residue):
    '''Residue -> Float'''
    if kdIndex.has_key(residue):
        return kdIndex[residue]
    raise KeyError("missing Kyte-Doolittle index for residue <%s>" % str(residue))


def smooth(values, windowRadius):
    '''[Float] -> Int -> [Float]'''
    i = windowRadius
    smoothed = []
    width = 2 * windowRadius + 1
    while i < len(values) - windowRadius:
        start, stop = i - windowRadius, i + windowRadius + 1
        total = sum(values[start: stop])
        smoothed.append(total / width)
        i += 1
    return smoothed


def kyteDoolittle(residues, windowRadius):
    '''[Residue] -> Int -> [Float]'''
    scored = [getResidueScore(r) for r in residues]
    smoothed = smooth(scored, windowRadius)
    return smoothed



############################

class KdTest(unittest.TestCase):

    def setUp(self):
        pass

    def testIndexLength(self):
        self.assertEqual(20, len(kdIndex))

    def testGetResidueScore(self):
        self.assertEqual(1, 0)

    def testSmooth(self):
        self.assertEqual(1, 0)

    def testAlgorithm(self):
        rs = 'MATTCVGHKWERTY'
        calced = kyteDoolittle(rs, 2)
        self.assertEqual(len(calced), len(rs) - 4)


testClasses = [KdTest]
