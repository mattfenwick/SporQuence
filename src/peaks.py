import unittest




def find1DPeaks(data, radius):
    '''[Float] -> Int -> [Peak]'''
    peaks = []
    i = radius
    while i < len(data) - radius:
        pt = data[i]
        neighbors = data[i - radius : i + radius + 1]
        if all([pt >= n for n in neighbors]):
            peaks.append({'index': i, 'height': pt})
            i += radius     # so that nearby points of equal height don't get counted twice ... ???
        i += 1
    return peaks


###########################################################################


class PeaksTest(unittest.TestCase):

    def setUp(self):
        pass

    def testNoPeaks(self):
        self.assertEqual(0, len(find1DPeaks(range(20), 5)))
    


testClasses = [PeaksTest]