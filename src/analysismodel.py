import kd
import codons as cn
import translate
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

    def getCodons(self):
        if self._codons is None:
            start, stop = self._startIndex, self._stopIndex
            self._codons = self._sequence.getCodons()[start:stop]
        return self._codons

    def getStopCodon(self):
        return self._sequence.getCodons()[self._stopIndex]

    def getUpstreamCodons(self, n):
        return self._sequence.getCodons()[self._startIndex - n : self._startIndex]

    def getKyteDool(self, windowRadius):
        if not self._kyteDool.has_key(windowRadius):
            residues = translate.codonsToResidues(self.getCodons())
            self._kyteDool[windowRadius] = kd.kyteDoolittle(residues, windowRadius)
        return self._kyteDool[windowRadius]
        


class Sequence(object):

    def __init__(self, bases):
        self._codons = cn.makeCodons(bases)

    def getCodons(self):
        return self._codons

    def getOrfs(self):
        orfEnds = cn.getOrfEnds(self._codons)
        return [ORF(start, stop, self) for (start, stop) in orfEnds]


########################################################################################


class OrfTest(unittest.TestCase):

    def setUp(self):
        pass

    def testCodonsLength(self):
        orf = ORF(3, 7, Sequence('ACGTGTGGCTAGCTAGCATCAGCA'))
        self.assertEqual(4, len(orf.getCodons()))
        self.assertEqual('CTA', orf.getCodons()[1].bases)

    def testKyteDool(self):
        orf = ORF(2, 7, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual(3, len(orf.getKyteDool(1)))
        self.assertEqual(5, len(orf.getKyteDool(0)))
        self.assertEqual(orf.getKyteDool(1), orf.getKyteDool(1))

    def testGetUpstreamCodons(self):
        orf = ORF(4, 7, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual('GGC', orf.getUpstreamCodons(2)[0].bases)

    def testGetStopCodon(self):
        orf = ORF(4, 7, Sequence('ACGTGTGGCGATCAAGCATCAGCAAAA'))
        self.assertEqual('GCA', orf.getStopCodon().bases)


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
        self.assertEqual(1, orfs[0].getStartIndex())
        self.assertEqual(4, orfs[1].getStartIndex())
        self.assertEqual(8, orfs[1].getStopIndex())
        self.assertEqual(4, len(orfs[1].getCodons()))
    


testClasses = [OrfTest, SequenceTest]
