import kd         as kd
import translate
import peaks
import sequence
import unittest





def bothKds(residues, kdRadius):
    return [kd.kyteDoolittle(residues, kdRadius), kd.triangleKyteDoolittle(residues, kdRadius)]



##########################
## filters


def hasHighPhobicPeak(orf, 
                      windowRadius  = 9, 
                      minHeight = 2.0, 
                      algorithm = kd.triangleKyteDoolittle):
    residues = translate.codonsToResidues(sequence.makeCodons(orf.bases))
    phobicity = algorithm(residues, windowRadius)
    return len(filter(lambda x: x > minHeight, phobicity)) > 0
  
    
def hasNPhobicPeaks(orf, 
                    numPeaks   = 2, 
                    minHeight  = 1.5, 
                    peakRadius = 9, 
                    kdRadius   = 9, 
                    algorithm  = kd.triangleKyteDoolittle):
    residues = translate.codonsToResidues(sequence.makeCodons(orf.bases))
    pks = peaks.get1DPeaks(algorithm(residues, windowRadius, peakRadius))
    highPeaks = filter(lambda p: p['height'] >= minHeight, pks)
    return len(highPeaks) == numPeaks


def hasUpstreamSequence(orf, 
                        upstreamSequence = 'GG'):
    promoterRegion = orf.upstream[-15:-5]              # want -15 to -5 region
    return promoterRegion.find(upstreamSequence) >= 0  # >= 0 means it found a match ... right?


def matchBases(seq, up, down):
    def isMatch(orf):
        a, b, c = orf.bases, orf.upstream, orf.downstream
        return a.startswith(seq) and b.endswith(up) and c.startswith(down)
    return isMatch




########################################################################################


class FiltersTest(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def testOne(self):
        self.assertTrue(False)


    


testClasses = [FiltersTest]
