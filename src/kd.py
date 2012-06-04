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
    width = float(2 * windowRadius + 1)
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


def triangleSmooth(values, windowRadius):
    '''[Float] -> Int -> [Float]'''
    i = windowRadius

    fWs, rWs = range(1, windowRadius + 1), range(windowRadius, 0, -1)
    weights = fWs + [windowRadius + 1] + rWs
    totalWeight = float(sum(weights))
    
    smoothed = []
    width = 2 * windowRadius + 1
    while i < len(values) - windowRadius:
        start, stop = i - windowRadius, i + windowRadius + 1
        vals = values[start: stop]
        
        assert len(vals) == len(weights), "need equal number of values and weights"
        
        total = sum([v * w for (v, w) in zip(vals, weights)])
        smoothed.append(total / totalWeight)
        i += 1
    return smoothed


def triangleKyteDoolittle(residues, windowRadius):
    '''[Residue] -> Int -> [Float]'''
    scored = [getResidueScore(r) for r in residues]
    smoothed = triangleSmooth(scored, windowRadius)
    return smoothed



############################

class KdTest(unittest.TestCase):

    def setUp(self):
        pass

    def testIndexLength(self):
        self.assertEqual(20, len(kdIndex))

    def testGetResidueScore(self):
        self.assertEqual((2.8, -0.7), (getResidueScore('F'), getResidueScore('T')))

    def testScoredResiduesSize(self):
        self.assertEqual(20, len(kdIndex))

    def testSmooth(self):
        self.assertEqual(1, 0)

    def testSmoothLength(self):
        smoothed = smooth(range(100), 10)
        self.assertEqual(len(smoothed), 80)

    def testAlgorithm(self):
        rs = 'MATTCV'
        c1, c2 = kyteDoolittle(rs, 1), kyteDoolittle(rs, 2)
        self.assertEqual(1.0, c1[0])
        self.assertAlmostEqual(.4 / 3, c1[1])
        self.assertAlmostEqual(0.96, c2[0])

    def testAlgorithmLength(self):
        rs = 'MATTCVGHKWERTY'
        calced = kyteDoolittle(rs, 2)
        self.assertEqual(len(calced), len(rs) - 4)

    def testTriangleKD(self):
        rs = 'MATTCV'
        c1, c2 = triangleKyteDoolittle(rs, 1), triangleKyteDoolittle(rs, 2)
        self.assertEqual((4, 2), (len(c1), len(c2)))
        self.assertEqual(1.2, c1[0])
        self.assertAlmostEqual(-0.075, c1[1])

    def testTriangleSmooth(self):
        calced = triangleSmooth([1,-1,3,4,5], 1)
        self.assertEqual(len(calced), 5 - 2)
        self.assertEqual(calced, [0.5, 2.25, 4])
        
        
        
        


testClasses = [KdTest]
