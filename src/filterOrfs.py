import kd         as kd
import translate  as translate
import peaks      as peaks
import sequence   as sq

import unittest





def bothKds(residues, kdRadius):
    return [kd.kyteDoolittle(residues, kdRadius), kd.triangleKyteDoolittle(residues, kdRadius)]



##########################
## filters


def hasHighPhobicPeak(orfAnal, 
                  kdRadius = 9, 
                  minHeight = 2.0, 
                  algorithm = kd.triangleKyteDoolittle):
    residues = orfAnal.getResidues()
    kdResults = algorithm(residues, kdRadius)
    return len(filter(lambda x: x > minHeight, kdResults)) > 0
  
    
def hasNPhobicPeaks(orfAnal, 
                      numPeaks = 2, 
                      minHeight = 1.5, 
                      peakRadius = 9, 
                      kdRadius = 9, 
                      algorithm = kd.triangleKyteDoolittle):
    residues = orfAnal.getResidues()
    hydroPhobicity = algorithm(residues, kdRadius)
    kdPeaks = peaks.find1DPeaks(hydroPhobicity, peakRadius)
    highPeaks = filter(lambda p: p['height'] >= minHeight, kdPeaks)
    return len(highPeaks) == numPeaks


def hasUpstreamSequence(orfAnal, 
                        upstreamSequence = 'GG'):
    upstream = orfAnal._orf['upstream']
    promoterRegion = upstream[-15:-5]                  # want -15 to -5 region
    return promoterRegion.find(upstreamSequence) >= 0  # >= 0 means it found a match ... right?


########################
# analysis objects

algorithms = {
    'kyteDool': kd.kyteDoolittle,
    'triKyteDool': kd.triangleKyteDoolittle              
}
    
class ORFAnalysis(object):

    def __init__(self, orf):
        self.orf = orf

    def getCodons(self):
        return sq.makeCodons(self.orf['sequence'])

    def getResidues(self):
        return ''.join(translate.codonsToResidues(self.getCodons()))     

    def getPhobicity(self, algorithm, windowRadius):
        alg = algorithms[algorithm]
        residues = self.getResidues()
        return alg(residues, windowRadius)          
                              
    def getPhobicityPeaks(self, algorithm, windowRadius, peakRadius):          
        return peaks.find1DPeaks(self.getPhobicity(algorithm, windowRadius), peakRadius)


filters = {
    'nPhobicPeaks': hasNPhobicPeaks,
    'upstream': hasUpstreamSequence
}


class ORFCollection(object):
    
    def __init__(self, orfAnals):
        self._orfAnals = orfAnals
        
    def getOrfAnals(self):
        return self._orfAnals
    
    def orfFilter(self, filterName, *args, **kwargs):
        filterFunc = filters[filterName]
        filtered = [orfAnal for orfAnal in self.getOrfAnals() if filterFunc(orfAnal, *args, **kwargs)]
        return ORFCollection(filtered)
    
    def findOrf(self, 
                start = None, 
                stop = None, 
                startFilter = None, 
                stopFilter = None, 
                orfFilter = None,
                *args, **kwargs):
        orfAnals = self.getOrfAnals()
        if start is not None:
            orfAnals = filter(lambda o: o.orf['start'] == start, orfAnals)
        if stop is not None:
            orfAnals = filter(lambda o: o.orf['stop'] == stop, orfAnals)
        if startFilter is not None:
            orfAnals = filter(lambda o: startFilter(o.orf['start']), orfAnals)
        if stopFilter is not None:
            orfAnals = filter(lambda o: stopFilter(o.orf['stop']), orfAnals)
        if orfFilter is not None:
            filterFunc = filters[orfFilter]
            orfAnals = filter(lambda o: filterFunc(o, *args, **kwargs), orfAnals)
        return ORFCollection(orfAnals)
        


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
