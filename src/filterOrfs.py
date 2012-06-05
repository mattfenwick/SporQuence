import orfilter.kd         as kd
import orfilter.translate  as translate
import orfilter.peaks      as peaks

import orfind.sequence     as sq

import unittest





def bothKds(residues, kdRadius):
    return [kd.kyteDoolittle(residues, kdRadius), kd.triangleKyteDoolittle(residues, kdRadius)]



# analysis objects (i.e. ugly)

class ORFAnalysis(object):

    def __init__(self, orf):
        self._orf = orf

    def getCodons(self):
        return sq.makeCodons(self._orf['sequence'])

    def getResidues(self):
        return ''.join(translate.codonsToResidues(self.getCodons()))




##########################
## filters


def hasHighKdPeak(residues, 
                  kdRadius = 9, 
                  minHeight = 2.0, 
                  algorithm = kd.triangleKyteDoolittle):
    kdResults = algorithm(residues, kdRadius)
    return len(filter(lambda x: x > minHeight, kdResults)) > 0
  
    
def hasNPhobicKdPeaks(residues, 
                      numPeaks = 2, 
                      minHeight = 1.5, 
                      peakRadius = 9, 
                      kdRadius = 9, 
                      algorithm = kd.triangleKyteDoolittle):
    hydroPhobicity = algorithm(residues, kdRadius)
    kdPeaks = peaks.find1DPeaks(hydroPhobicity, peakRadius)
    highPeaks = filter(lambda p: p['height'] >= minHeight, kdPeaks)
    return len(highPeaks) == numPeaks


def hasUpstreamSequence(orf, 
                        upstreamSequence = 'GG'):
    upstream = orf['upstream']

    promoterRegion = upstream[-15:-5]                  # want -15 to -5 region
    return promoterRegion.find(upstreamSequence) >= 0  # >= 0 means it found a match ... right?


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


    


testClasses = [OrfTest]
