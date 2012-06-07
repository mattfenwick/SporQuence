import sequence
import translate
import peaks
import unittest




########################################################################
# for finding Orfs

START_STOP_ERROR = "start and stop must be between 0 and sequence length"


class Orf(object):
    
    '''public properties:
        - start
        - stop
        - bases
        
       public methods:
        - getNormStart
        - getNormStop
        - getUpstream(n)
          - n: number of upstream bases
        - getDownstream(n)
          - n: number of downstream bases
        - toJSONObject(n)
          - n: number of up-/down-stream bases
    '''

    def __init__(self, start, stop, sequence):
        length = len(sequence.getBases())
        assert 0 <= start <= length, START_STOP_ERROR
        assert 0 <= stop <= length, START_STOP_ERROR
        
        self.start = start
        self.stop = stop
        self._sequence = sequence

        self.bases = self._getBases()
        assert len(self.bases) % 3 == 0, "Orf length must be multiple of 3 <%s>" % str(self.bases)

    ######################

    def _getBases(self):
        start, stop, seq = self.start, self.stop, self._sequence
        bases = seq.getBases()

        # typical case
        if stop > start:
            myBases = bases[start:stop]

        # boundary condition:  circular genome, Orf wraps around
        else:
            myBases = bases[start:] + bases[:stop]

        return myBases
        
    def _isSense(self):
        return self._sequence.isSense()
    
    def _getNormIndex(self, index):
        seq = self._sequence
        length = len(seq.getBases())
        assert 0 <= index < length
        if self._isSense():
            return index
        else:
            return length - index - 1
        
    #########################
        
    def getNormStart(self):
        return self._getNormIndex(self.start)
        
    def getNormStop(self):
        return self._getNormIndex(self.stop)

    def getUpstream(self, n):
        start, seq = self.start, self._sequence
        bases = seq.getBases()
        length = len(bases)

        assert n >= 0, "number of upstream bases must be positive"
        assert n <= length, "number of upstream bases requested greater than sequence length"
        
        if start - n >= 0:
            return bases[start - n : start]
        else:
            index = length + start - n
            return bases[index:] + bases[:start]

    def getDownstream(self, n):
        '''includes stop sequence'''
        stop, seq = self.stop, self._sequence
        bases = seq.getBases()
        length = len(bases)

        assert n >= 0, "number of upstream bases must be positive"
        assert n <= length, "number of upstream bases requested greater than sequence length"
        
        if stop + n <= length:
            return bases[stop : stop + n]
        else:
            index = stop + n - length
            return bases[stop:] + bases[:index]

    ######################

    def toJSONObject(self, n):
        return {
            'start'       : self.getNormStart(),
            'stop'        : self.getNormStop(),
            'bases'       : self.bases,
            'upstream'    : self.getUpstream(n),
            'downstream'  : self.getDownstream(n),
            'isSense'     : self._isSense()
        }
        



class Sequence(object):
    
    '''public properties:
        
       public methods:
        - getBases
        - isSense
        - getReverseComplement
        - getOrfs
    '''

    def __init__(self, bases, isSense):
        self._bases = bases
        self._isSense = isSense
        self._codons = [None, None, None]
        self._reverse = None
    
    def _getReverseBases(self):
        if self._reverse is None:
            self._reverse = sequence.reverseComplement(self.getBases())
        return self._reverse

    def _getCodons(self, n):
        assert n in [0, 1, 2], "codon alignment must be 0, 1, or 2"
        bases = self.getBases()
        if self._codons[n] is None:
            self._codons[n] = sequence.makeCodons(bases[n:] + bases[:n])
        return self._codons[n]
    
    def _getOrfs(self, algorithm):
        orfs = []
        for n in range(3): # [0, 1, 2]
            orfEnds = algorithm(self._getCodons(n))
            orfs += [Orf(start * 3 + n, stop * 3 + n, self) for (start, stop) in orfEnds]
        return orfs
    
    ##############################

    def getBases(self):
        return self._bases
    
    def isSense(self):
        return self._isSense
        
    def getOrfs(self):
        '''finds *only leftmost, longest* Orfs in all 3 alignments'''
        return self._getOrfs(sequence.getOrfEndsCircular)
    
    def getAllOrfs(self):
        '''finds Orfs of all sizes (including overlapping) in all 3 alignments'''
        return self._getOrfs(sequence.getAllOrfEndsCircular)
    
    def getReverseComplement(self):
        rSeq = Sequence(self._getReverseBases(), not self.isSense())
        return rSeq




########################################################################
# for analysis:  filtering and analyzing Orfs

    
class OrfAnalysis(object):

    def __init__(self, orf):
        self.orf = orf

    def getCodons(self):
        return sequence.makeCodons(self.orf['bases'])

    def getResidues(self):
        return ''.join(translate.codonsToResidues(self.getCodons()))     

    def getPhobicity(self, algorithm, windowRadius):
        residues = self.getResidues()
        return algorithm(residues, windowRadius)          
                              
    def getPhobicityPeaks(self, algorithm, windowRadius, peakRadius):          
        return peaks.find1DPeaks(self.getPhobicity(algorithm, windowRadius), peakRadius)


class OACollection(object):
    
    def __init__(self, orfAnals):
        self._orfAnals = orfAnals
        
    def getOrfAnals(self):
        return self._orfAnals
    
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
            orfAnals = filter(lambda o: orfFilter(o, *args, **kwargs), orfAnals)
        return OACollection(orfAnals)
        
     
            
########################################################
# unit tests
########################################################

class OrfTest(unittest.TestCase):

    ''' need to test: 1) forward, normal; 2) forward wrap; 3) reverse normal; 4) reverse wrap'''

    def setUp(self):
        self.forNorm = Orf(2, 14, Sequence('GA' + 'GCTAGCATCGAT'+ 'TCGAT', True))
        self.forWrap = Orf(8, 3,  Sequence('GAG' + 'CTAGC' + 'ATCGATCGA',  True))
        self.revNorm = Orf(4, 13, Sequence('GAGC' + 'TAGCATCGA' + 'TCGAT', False))
        self.revWrap = Orf(11, 6, Sequence('GAGCTA' + 'GCATC' + 'GATCGA',  False))

    def testBases(self):
        self.assertEqual('GCTAGCATCGAT', self.forNorm.bases)
        self.assertEqual('ATCGATCGAGAG', self.forWrap.bases)
        self.assertEqual('TAGCATCGA', self.revNorm.bases)
        self.assertEqual('GATCGAGAGCTA', self.revWrap.bases)
        
    def testGetUpstream(self):
        self.assertEqual('GATGA',   self.forNorm.getUpstream(5))
        self.assertEqual('CTAGC',   self.forWrap.getUpstream(5))
        self.assertEqual('GATGAGC', self.revNorm.getUpstream(7))
        self.assertEqual('TC',      self.revWrap.getUpstream(2))
        
    def testGetDownstream(self):
        self.assertEqual('TCGATGA',       self.forNorm.getDownstream(7))
        self.assertEqual('CTA',           self.forWrap.getDownstream(3))
        self.assertEqual('TCGAT',         self.revNorm.getDownstream(5))
        self.assertEqual('GCATCGATCGAG',  self.revWrap.getDownstream(12))
    
    def testGetNormStartStop(self):
        self.assertEqual((2, 14), (self.forNorm.getNormStart(), self.forNorm.getNormStop()))
        self.assertEqual((8, 3),  (self.forWrap.getNormStart(), self.forWrap.getNormStop()))
        self.assertEqual((13, 4), (self.revNorm.getNormStart(), self.revNorm.getNormStop()))
        self.assertEqual((5, 10), (self.revWrap.getNormStart(), self.revWrap.getNormStop()))
    
    def testStartStop(self):
        self.assertEqual((2, 14), (self.forNorm.start, self.forNorm.stop))
        self.assertEqual((8, 3),  (self.forWrap.start, self.forWrap.stop))
        self.assertEqual((4, 13), (self.revNorm.start, self.revNorm.stop))
        self.assertEqual((11, 6), (self.revWrap.start, self.revWrap.stop))
        
    def testToJSONObject(self):
        self.assertEqual({
            'start': 2, 
            'stop': 14, 
            'bases': 'GCTAGCATCGAT', 
            'upstream': 'TGA', 
            'downstream': 'TCG',
            'isSense': True
        }, self.forNorm.toJSONObject(3))
        self.assertEqual({
            'start': 5, 
            'stop': 10, 
            'bases': 'GATCGAGAGCTA', 
            'upstream': 'ATC', 
            'downstream': 'GCA',
            'isSense': False
        }, self.revWrap.toJSONObject(3))


class SequenceTest(unittest.TestCase):

    def setUp(self):
        self.seq = Sequence("ACGTAA" + "CCC" + "CTGAAAGGGTAG" + "ATGTTTTAC", True)
        self.orfs = self.seq.getOrfs()

    def testGetBases(self):
        self.assertEqual("ACGTAACCCCTGAAAGGGTAGATGTTTTAC", self.seq.getBases())

    def testIsSense(self):
        self.assertTrue(self.seq.isSense())
        self.assertFalse(self.seq.getReverseComplement().isSense())

    def testGetReverseComplement(self):
        self.assertEqual("GTAA" + "AACATCTACCCTTTCAGGGGTTAC" + "GT", self.seq.getReverseComplement().getBases())

    def testGetOrfs(self):
        ofs, ors = self.orfs, self.seq.getReverseComplement().getOrfs()
        self.assertEqual(2, len(ofs))
        self.assertEqual(1, len(ors))
        self.assertEqual('AACCC', ofs[0].getUpstream(5))
        self.assertEqual('TAACC', ofs[1].getDownstream(5))
        self.assertEqual('GTTAC', ors[0].getUpstream(5))
        self.assertEqual(1, ors[0].getNormStart())

    def testGetAllOrfs(self):
        seq = Sequence('AATTAAAATAGA' + 'ATGGTGTGCTGC', False)
        # 'GCAGCACACCAT' + 'TCTATTTTAATT'
        rs, fs = seq.getAllOrfs(), seq.getReverseComplement().getAllOrfs()
        self.assertEqual(4, len(rs))
        self.assertEqual(1, len(fs))
        
        f = fs[0]
        self.assertEqual('TAAT',  f.getDownstream(4))
        self.assertEqual('TTAA',  f.getUpstream(4))
        self.assertEqual('TTGCA', f.bases[:5])
        self.assertEqual((22, 19), (f.getNormStart(), f.getNormStop()))
        
        myRs = set([(11, 20), (8, 20), (6, 15), (3, 15)])
        self.assertEqual(myRs, set([(r.getNormStart(), r.getNormStop()) for r in rs]))
        


class OrfAnalysisTest(unittest.TestCase):

    def setUp(self):
        self.orf = Orf(2, 14, Sequence('TTACGCTACCTTTCGCC', True))

    def testCodonsLength(self):
        orfA = OrfAnalysis(self.orf.toJSONObject(4))
        self.assertEqual(4, len(orfA.getCodons()))
        self.assertEqual('CTA', orfA.getCodons()[1])

    def testPhobicity(self):
        self.assertTrue(False)

    def testPhobicityPeaks(self):
        self.assertTrue(False)
        


class OACollectionTest(unittest.TestCase):

    def setUp(self):
        pass
    
    def testOne(self):
        self.assertFalse(True)



testClasses = [OrfTest, SequenceTest, OrfAnalysisTest, OACollectionTest]
