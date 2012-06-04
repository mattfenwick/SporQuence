import unittest
import json
import sequence



class ORF(object):

    def __init__(self, startIndex, stopIndex, sequence):
        assert (stopIndex - startIndex) % 3 == 0, "ORF length must be multiple of 3"
        self.startIndex = startIndex
        self.stopIndex = stopIndex
        self.bases = self._getBases(sequence)
        self.upstream = self._getUpstream(sequence)
        self.downstream = self._getDownstream(sequence)

    ######################

    def toJSONObject(self):
        return {
            'start'       : self.startIndex,
            'sequence'    : self.bases,
            'stop'        : self.stopIndex,
            'upstream'    : self.upstream,
            'downstream'  : self.downstream
        }

    def __repr__(self):
        return json.dumps(self.toJSONObject())

    ######################

    def _getBases(self, sequence):
        start, stop = self.startIndex, self.stopIndex
        bases = sequence.getBases()

        # typical case
        if stop > start:
            myBases = bases[start:stop]

        # boundary condition:  circular genome, ORF wraps around
        else:
            myBases = bases[start:] + bases[:stop]

        return myBases

    def _getUpstream(self, sequence, n = 100):
        start = self.startIndex
        return sequence.getBases()[start - n : start]

    def _getDownstream(self, sequence, n = 100):
        '''includes stop sequence'''
        stop = self.stopIndex
        return sequence.getBases()[stop : stop + n]



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
            self._codons[n] = sequence.makeCodons(bases[n:] + bases[:n])
        return self._codons[n]

    def getOrfs(self):
        '''finds for ORFs in all 3 alignments'''
        orfs = []
        for n in range(3): # [0, 1, 2]
            orfEnds = sequence.getOrfEndsCircular(self.getCodons(n))
            orfs += [ORF(start * 3 + n, stop * 3 + n, self) for (start, stop) in orfEnds]
        return orfs
    
    def getReverseOrfs(self):
        rSeq = Sequence(sequence.reverseComplement(self._bases))
        return rSeq.getOrfs()
            
            
########################################################
# unit tests
########################################################

class OrfTest(unittest.TestCase):

    def setUp(self):
        pass

    def testOne(self):
        self.assertEqual(0, 1)


class SequenceTest(unittest.TestCase):

    def setUp(self):
        pass

    def testOne(self):
        self.assertEqual(0, 1)

    

testClasses = [OrfTest, SequenceTest]
