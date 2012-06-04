import kd
import sequence as sq
import translate
import peaks

import unittest
import yaml



# analysis objects (i.e. ugly)

class ORF(object):

    def __init__(self, startIndex, stopIndex, sequence):
        assert (stopIndex - startIndex) % 3 == 0, "ORF length must be multiple of 3"
        self._startIndex = startIndex
        self._stopIndex = stopIndex
        self._sequence = sequence
        self._codons = None
        self._kyteDool = {}
        # should ORFs keep the residues too?

    ######################

    def toJSON(self):
        return {
            'start'       : self.getStartIndex(),
            'sequence'    : self.getBases(),
            'stop'        : self.getStopIndex(),
            'translation' : self.getResidues(),
            'upstream'    : self.getUpstreamBases(100),
            'downstream'  : self.getDownstreamBases(100)
        }

    def __repr__(self):
        return yaml.dump(self.toJSON())

    def getStartIndex(self):
        return self._startIndex

    def getStopIndex(self):
        return self._stopIndex

    ######################

    def getBases(self):
        start, stop = self.getStartIndex(), self.getStopIndex()
        bases = self._sequence.getBases()

        # typical case
        if stop > start:
            self._bases = bases[start:stop]

        # boundary condition:  circular genome, ORF wraps around
        else:
            self._bases = bases[start:] + bases[:stop]

        return self._bases

    def getUpstreamBases(self, n):
        start = self.getStartIndex()
        return self._sequence.getBases()[start - n : start]

    def getDownstreamBases(self, n):
        '''includes stop sequence'''
        stop = self.getStopIndex()
        return self._sequence.getBases()[stop : stop + n]

    def getCodons(self):
        return sq.makeCodons(self.getBases())

    def getResidues(self):
        return ''.join(translate.codonsToResidues(self.getCodons()))

    def getKyteDool(self, windowRadius):
        if not self._kyteDool.has_key(windowRadius):
            residues = self.getResidues()
            self._kyteDool[windowRadius] = kd.kyteDoolittle(residues, windowRadius)
        return self._kyteDool[windowRadius]

    def getTriangleKyteDool(self, windowRadius):
        residues = self.getResidues()
        return kd.triangleKyteDoolittle(residues, windowRadius)

    def getKDPeaks(self, kdRadius, peakRadius):
        return peaks.find1DPeaks(self.getKyteDool(kdRadius), peakRadius)

    def getTriangleKDPeaks(self, kdRadius, peakRadius):
        return peaks.find1DPeaks(self.getTriangleKyteDool(kdRadius), peakRadius)


class Sequence(object):

    def __init__(self, bases):
        self._bases = bases
        self._codons = [None, None, None]

    def getBases(self):
        return self._bases

    def getCodons(self, n):
        assert n in [0, 1, 2], "codon alignment must be 0, 1, or 2"
        bases = self.getBases()
        if self._codons[n] is None:
            self._codons[n] = sq.makeCodons(bases[n:] + bases[:n])
        return self._codons[n]

    def getOrfs(self):
        '''finds for ORFs in all 3 alignments'''
        orfs = []
        for n in range(3): # [0, 1, 2]
            orfEnds = sq.getOrfEndsCircular(self.getCodons(n))
            orfs += [ORF(start * 3 + n, stop * 3 + n, self) for (start, stop) in orfEnds]
        return orfs



class SequenceMuncher(object):

    def __init__(self, sequence, low = 50, high = 80, kdRadius = 4, peakRadius = 4, peakHeightCutoff = 2.0, upstreamSequence = 'GG'):
        self._forward = sequence
        self._reverse = Sequence(sq.reverseComplement(sequence.getBases()))

        # configuration variables
        self.low, self.high = 50, 80
        self.kdRadius, self.peakRadius = kdRadius, peakRadius
        self.peakHeightCutoff, self.upstreamSequence = peakHeightCutoff, upstreamSequence

        self._orfs, self._smallOrfs = None, None

    def getOrfs(self):
        if self._orfs is None:
            self._orfs = {
                'forward': self._forward.getOrfs(),
                'reverse': self._reverse.getOrfs()
            }
        return self._orfs

    def getSmallOrfs(self):
        if self._smallOrfs is None:
            def filterer(orf):
                cs = len(orf.getCodons())
                return self.low <= cs <= self.high
            orfs = self.getOrfs()
            self._smallOrfs = {
                'forward': filter(filterer, orfs['forward']),
                'reverse': filter(filterer, orfs['reverse'])
            }
        return self._smallOrfs

    def getPhobicOrfs(self):
        def kdFilter(orf):
            kdResults = orf.getTriangleKyteDool(self.kdRadius)
            return len(filter(lambda x: x > self.peakHeightCutoff, kdResults)) > 0
        orfs = self.getSmallOrfs()
        return {
            'forward': filter(kdFilter, orfs['forward']),
            'reverse': filter(kdFilter, orfs['reverse'])
        }

    def getTwoPeakOrfs(self):
        def peakFilter(orf):
            kdPeaks = orf.getTriangleKDPeaks(self.kdRadius, self.peakRadius)
            highPeaks = filter(lambda p: p['height'] >= self.peakHeightCutoff, kdPeaks)
            return len(highPeaks) == 2
        orfs = self.getPhobicOrfs()
        return {
            'forward': filter(peakFilter, orfs['forward']),
            'reverse': filter(peakFilter, orfs['reverse'])
        }

    def getPromotedOrfs(self):
        def upFilter(orf):
            promoterRegion = orf.getUpstreamBases(15)[:10] # want -15 to -5 region
            return promoterRegion.find('GG') >= 0          # >= 0 means it found a match ... right?
        orfs = self.getTwoPeakOrfs()
        return {
            'forward': filter(upFilter, orfs['forward']),
            'reverse': filter(upFilter, orfs['reverse'])
        }


########################################################################################


class OrfTest(unittest.TestCase):

    def setUp(self):
        pass

    def testBases(self):
        orf = ORF(3, 9, Sequence('ACGTATCAGCAA'))
        self.assertEqual('TATCAG', orf.getBases())

    def testCodonsLength(self):
        orf = ORF(9, 21, Sequence('ACGTGTGGCTAGCTAGCATCAGCA'))
        self.assertEqual(4, len(orf.getCodons()))
        self.assertEqual('CTA', orf.getCodons()[1].bases)

    def testBasesWraparound(self):
        orf = ORF(12, 6, Sequence('ACGTGT' + 'GGCTAG' + 'CTAGCATCAGCA'))
        self.assertEqual('CTAGCATCAGCAACGTGT', orf.getBases())

    def testKyteDool(self):
        orf = ORF(6, 21, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual(3, len(orf.getKyteDool(1)))
        self.assertEqual(5, len(orf.getKyteDool(0)))
        self.assertEqual(orf.getKyteDool(1), orf.getKyteDool(1))

    def testGetUpstreamBases(self):
        orf = ORF(12, 21, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual('GGCGAT', orf.getUpstreamBases(6))

    def testGetDownstreamBases(self):
        orf = ORF(6, 15, Sequence('ACGTGT' + 'GGCGATCAA' + 'GCATCAGCAAAA'))
        self.assertEqual('GCATCAGCA', orf.getDownstreamBases(9))

    def testGetKdPeaks(self):
        self.assertTrue(False)

    def testStartStopIndex(self):
        orf = ORF(12, 21, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual(12, orf.getStartIndex())
        self.assertEqual(21, orf.getStopIndex())


class SequenceTest(unittest.TestCase):

    def setUp(self):
        pass

    def testCodonsLength(self):
        seq = Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA')
        self.assertEqual([9, 9, 9], [len(seq.getCodons(0)), len(seq.getCodons(1)), len(seq.getCodons(2))])

    def testGetOrfs(self):
        seq = Sequence('ACG' + 'TTGGGCTAA' + 'ATGGCATCAGCATAA')
        orfs = seq.getOrfs()
        self.assertEqual(2, len(orfs))
        self.assertEqual(set([(3, 9), (12, 24)]), set(map(lambda o: (o.getStartIndex(), o.getStopIndex()), orfs)))
        self.assertEqual(set([2, 4]), set(map(lambda o: len(o.getCodons()), orfs)))

    def testGetOrfsWraparound(self):
        seq = Sequence('ACGTTTTAA' + 'AAAAAA' + 'TTGGGCTAT')
        orfs = seq.getOrfs()
        self.assertEqual(1, len(orfs))
        self.assertEqual('TTGGGCTATACGTTT', orfs[0].getBases())
    


testClasses = [OrfTest, SequenceTest]
