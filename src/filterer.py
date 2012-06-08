import kd         as kd
import unittest





def bothKds(residues, kdRadius):
    return [kd.kyteDoolittle(residues, kdRadius), kd.triangleKyteDoolittle(residues, kdRadius)]



##########################
## filters


def hasHighPhobicPeak(orfAnal, 
                      kdRadius  = 9, 
                      minHeight = 2.0, 
                      algorithm = kd.triangleKyteDoolittle):
    phobicity = orfAnal.getPhobicity(algorithm)
    return len(filter(lambda x: x > minHeight, phobicity)) > 0
  
    
def hasNPhobicPeaks(orfAnal, 
                    numPeaks   = 2, 
                    minHeight  = 1.5, 
                    peakRadius = 9, 
                    kdRadius   = 9, 
                    algorithm  = kd.triangleKyteDoolittle):
    peaks = orfAnal.getPhobicityPeaks(algorithm, windowRadius, peakRadius)
    highPeaks = filter(lambda p: p['height'] >= minHeight, peaks)
    return len(highPeaks) == numPeaks


def hasUpstreamSequence(orfAnal, 
                        upstreamSequence = 'GG'):
    upstream = orfAnal.orf['upstream']
    promoterRegion = upstream[-15:-5]                  # want -15 to -5 region
    return promoterRegion.find(upstreamSequence) >= 0  # >= 0 means it found a match ... right?


def matchBases(seq, up, down):
    def isMatch(orfAnal):
        orf = orfAnal.orf
        a, b, c = orf['bases'], orf['upstream'], orf['downstream']
        return a.startswith(seq) and b.endswith(up) and c.startswith(down)
    return isMatch




########################################################################################





    


testClasses = []
