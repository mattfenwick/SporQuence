import kd
import sequence as sq
import translate
import peaks

import unittest
import yaml



# analysis objects (i.e. ugly)

class ORF(object):

    def __init__(self, startIndex, stopIndex, sequence):
        self._startIndex = startIndex
        self._stopIndex = stopIndex
        self._sequence = sequence
        self._codons = None
        self._kyteDool = {}
        # should ORFs keep the residues too?

    ######################

    def toJSON(self):
        return {
            'startIndex': self.getStartIndex(),
            'codons': [c.toJSON() for c in self.getCodons()],
            'stopIndex': self.getStopIndex()
        }

    def __repr__(self):
        return yaml.dump(self.toJSON())

    def getStartIndex(self):
        return self._startIndex

    def getStopIndex(self):
        return self._stopIndex

    ######################

    def getResidues(self):
        return translate.codonsToResidues(self.getCodons())

    def getCodons(self):
        if self._codons is None:
            seqCodons = self._sequence.getCodons()
            start, stop = self._startIndex, self._stopIndex

            # typical case
            if stop > start:
                self._codons = seqCodons[start:stop]

            # boundary condition:  circular genome, ORF wraps around
            else:
                self._codons = seqCodons[start:] + seqCodons[:stop]

        return self._codons

    def getBases(self):
        return ''.join([c.bases for c in self.getCodons()])

    def getStopCodon(self):
        return self._sequence.getCodons()[self._stopIndex]

    def getUpstreamCodons(self, n):
        return self._sequence.getCodons()[self._startIndex - n : self._startIndex]

    def getDownstreamCodons(self, n):
        return self._sequence.getCodons()[self._stopIndex : self._stopIndex + n]

    def getKyteDool(self, windowRadius):
        if not self._kyteDool.has_key(windowRadius):
            residues = self.getResidues()
            self._kyteDool[windowRadius] = kd.kyteDoolittle(residues, windowRadius)
        return self._kyteDool[windowRadius]

    def getKDPeaks(self, kdRadius, peakRadius):
        return peaks.find1DPeaks(self.getKyteDool(kdRadius), peakRadius)



class Sequence(object):

    def __init__(self, bases):
        self._bases = bases
        self._codons = sq.makeCodons(bases)

    def getBases(self):
        return self._bases

    def getCodons(self):
        return self._codons

    def getOrfs(self):
        orfEnds = sq.getOrfEndsCircular(self._codons)
        return [ORF(start, stop, self) for (start, stop) in orfEnds]


class OrfSelector(object):

    def __init__(self, sequence):
        self._forward = sequence
        self._reverse = Sequence(sq.reverseComplement(sequence.getBases()))
        self._forwardOrfs = self._forward.getOrfs()
        self._reverseOrfs = self._reverse.getOrfs()

    def getOrfs(self):
        return {
            'forward': self._forwardOrfs,
            'reverse': self._reverseOrfs
        }

    def getSmallOrfs(self, low, high):
        def filterer(orf):
            cs = len(orf.codons)
            return low <= cs <= high
        orfs = self.getOrfs()
        return {
            'forward': filter(filterer, orfs['forward']),
            'reverse': filter(filterer, orfs['reverse'])
        }

#    def getPhobicOrfs(self, 


########################################################################################


class OrfTest(unittest.TestCase):

    def setUp(self):
        pass

    def testCodons(self):
        orf = ORF(1, 3, Sequence('ACGTATCAGCAA'))
        self.assertEqual('TATCAG', orf.getBases())

    def testCodonsLength(self):
        orf = ORF(3, 7, Sequence('ACGTGTGGCTAGCTAGCATCAGCA'))
        self.assertEqual(4, len(orf.getCodons()))
        self.assertEqual('CTA', orf.getCodons()[1].bases)

    def testCodonsWraparound(self):
        orf = ORF(4, 2, Sequence('ACGTGT' + 'GGCTAG' + 'CTAGCATCAGCA'))
        self.assertEqual('CTAGCATCAGCAACGTGT', orf.getBases())

    def testKyteDool(self):
        orf = ORF(2, 7, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual(3, len(orf.getKyteDool(1)))
        self.assertEqual(5, len(orf.getKyteDool(0)))
        self.assertEqual(orf.getKyteDool(1), orf.getKyteDool(1))

    def testGetUpstreamCodons(self):
        orf = ORF(4, 7, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual('GGC', orf.getUpstreamCodons(2)[0].bases)

    def testGetDownstreamCodons(self):
        self.assertTrue(False)

    def testGetStopCodon(self):
        orf = ORF(4, 7, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual('GCA', orf.getStopCodon().bases)

    def testGetPeaks(self):
        self.assertTrue(False)


class SequenceTest(unittest.TestCase):

    def setUp(self):
        pass

    def testCodonsLength(self):
        seq = Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA')
        self.assertEqual(9, len(seq.getCodons()))

    def testGetOrfs(self):
        seq = Sequence('ACG' + 'TTGGGCTAA' + 'ATGGCATCAGCATAA')
        orfs = seq.getOrfs()
        self.assertEqual(2, len(orfs))
        self.assertEqual(set([(1, 3), (4, 8)]), set(map(lambda o: (o.getStartIndex(), o.getStopIndex()), orfs)))
        self.assertEqual(set([2, 4]), set(map(lambda o: len(o.getCodons()), orfs)))

    def testGetOrfsWraparound(self):
        seq = Sequence('ACGTTTTAA' + 'AAAAAA' + 'TTGGGCTAT')
        orfs = seq.getOrfs()
        self.assertEqual(1, len(orfs))
        self.assertEqual('TTGGGCTATACGTTT', orfs[0].getBases())
    


testClasses = [OrfTest, SequenceTest]
