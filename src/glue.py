import unittest
import bacillussubtilis168 as bs
import analysismodel as am
import yaml


##################### junk functions (junc-tions)

def kdFilter(kdResults):
    return len(filter(lambda x: x > 2, kdResults)) > 0

def countBig(kdResults):
    return len(filter(lambda x: x > 2, kdResults))

def peakFilter(orf):
    kdPeaks = orf.getKDPeaks(4, 4)
    highPeaks = filter(lambda p: p['height'] >= 2, kdPeaks)
    return len(highPeaks) == 2

def dumpOrfAndContext(orf):
    return {
        'start': orf.getStartIndex() * 3 + 1, # 3 for (codon -> base), 1 to correct for 0-indexing
        'upstream': ''.join(map(lambda c: c.bases, orf.getUpstreamCodons(34))),
        'downstream': ''.join(map(lambda c: c.bases, orf.getDownstreamCodons(34))),
        'sequence': ''.join(map(lambda c: c.bases, orf.getCodons()))
    }

######################


# codons of bacillus subtilis
seq, rev = am.Sequence(bs.sequence), am.Sequence(bs.reverseComplement)

# [ORF] :: open reading frames
orfs, revOrfs = seq.getOrfs(), rev.getOrfs()

# [ORF] :: open reading frames between 50 and 80 codons long
smallOrfs = filter(lambda o: len(o.getCodons()) >= 50 and len(o.getCodons()) <= 80, orfs)
smallRevOrfs = filter(lambda o: len(o.getCodons()) >= 50 and len(o.getCodons()) <= 80, revOrfs)

# [ORF] :: just those where Kyte-Doolittles has at least one value > 2
phobicOrfs = filter(lambda o: kdFilter(o.getKyteDool(9)), smallOrfs)
phobicRevOrfs = filter(lambda o: kdFilter(o.getKyteDool(9)), smallRevOrfs)

orfsTwoPeaks = filter(peakFilter, phobicOrfs)
revOrfsTwoPeaks = filter(peakFilter, phobicRevOrfs)

dumped = yaml.dump(map(dumpOrfAndContext, orfsTwoPeaks))
revDumped = yaml.dump(map(dumpOrfAndContext, revOrfsTwoPeaks))

#transMembraneCounts = map(lambda ps: len(filter(lambda y: y['height'] >= 2, ps)), pks)

############################

class GlueTest(unittest.TestCase):

    def setUp(self):
        pass

    def testCodonsLength(self):
        self.assertEqual(1405202, len(seq.getCodons()))

    def testOrfsLength(self):
        self.assertEqual(27522, len(orfs))

    def testSmallOrfsLength(self):
        self.assertEqual(1782, len(smallOrfs))
    


testClasses = [GlueTest]
