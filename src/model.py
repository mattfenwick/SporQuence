import sequence
import translate
import peaks
import unittest




########################################################################
# for finding Orfs


class Orf(object):
    
    '''public (read-only) properties:
        - start (normalized by direction)
        - stop (normalized by direction)
        - bases
        - upstream
        - downstream
        - is_sense
        
       public methods:
        - toJSONObject()
        
       unsure methods:
        - get_codons
        - get_residues
        - get_phobicity
        - get_phobicity_peaks
    '''

    def __init__(self, start, stop, bases, upstream, downstream, is_sense):
        
        self.start = start
        self.stop = stop

        self.bases = bases
        assert len(self.bases) % 3 == 0, "Orf length must be multiple of 3 <%s>" % str(self.bases)
        
        self.upstream = upstream
        self.downstream = downstream
        self.is_sense = is_sense


    def to_JSON_object(self):
        return {
            'start'       : self.start,
            'stop'        : self.stop,
            'bases'       : self.bases,
            'upstream'    : self.upstream,
            'downstream'  : self.downstream,
            'is_sense'    : self.is_sense
        }
        
    @staticmethod
    def from_JSON_object(obj):
        return Orf(
            obj['start'],
            obj['stop'],
            obj['bases'],
            obj['upstream'],
            obj['downstream'],
            obj['is_sense']       
        )
        
    ##############################

    def get_codons(self):
        return sequence.makeCodons(self.bases)

    def get_residues(self):
        return ''.join(translate.codonsToResidues(self.get_codons()))     

    def get_phobicity(self, algorithm, windowRadius):
        residues = self.get_residues()
        return algorithm(residues, windowRadius)          
                              
    def get_phobicity_peaks(self, algorithm, windowRadius, peakRadius):          
        return peaks.find1DPeaks(self.get_phobicity(algorithm, windowRadius), peakRadius)
        



_START_STOP_ERROR = "start and stop must be between 0 and sequence length"

class Sequence(object):
    
    '''public (read-only) properties:
        
       public methods:
        - get_bases
        - is_sense
        - get_reverse_complement
        - get_orfs
        - get_all_orfs
        - build_orf(start, stop)
    '''

    def __init__(self, bases, is_sense):
        self._bases = bases
        self._is_sense = is_sense
        self._codons = [None, None, None]
        self._reverse = None
        
    ########################
    
    def _get_reverse_bases(self):
        if self._reverse is None:
            self._reverse = sequence.reverseComplement(self.get_bases())
        return self._reverse


    def _get_codons(self, n):
        assert n in [0, 1, 2], "codon alignment must be 0, 1, or 2"
        
        bases = self.get_bases()
        if self._codons[n] is None:
            self._codons[n] = sequence.makeCodons(bases[n:] + bases[:n])
        return self._codons[n]
    
    
    def build_orf(self, bstart, bstop, n):
        '''Build an Orf given the start and stop indices.
        
        Parameters:
         - bstart:  base-pair start index
         - bstop:  base-pair stop index
        '''
        length = len(self.get_bases())
        nstart, nstop = [self._get_norm_index(ix) for ix in [bstart, bstop]]
        bases = self.get_bases(bstart, bstop)
        up = self.get_bases((bstart - n) % length, bstart)
        down = self.get_bases(bstop, (bstop + n) % length)
        return Orf(nstart, nstop, bases, up, down, self.is_sense())
  
    
    def _get_orfs(self, algorithm, width):
        orfs = []
        for n in range(3): # [0, 1, 2]
            orfEnds = algorithm(self._get_codons(n))
            for ends in orfEnds:
                bstart, bstop = [cix * 3 + n for cix in ends]
                new_orf = self.build_orf(bstart, bstop, width)
                orfs.append(new_orf)
        return orfs
 
    
    def _get_norm_index(self, index):
        length = len(self.get_bases())
        assert 0 <= index < length

        if self.is_sense():
            return index
        else:
            return length - index - 1
    
    ##############################

    def get_bases(self, start=None, stop=None):
        if start is None or stop is None:
            return self._bases
        
        assert None not in [start, stop], "start and stop can not be None"
        length = len(self._bases)
        assert 0 <= start < length, _START_STOP_ERROR
        assert 0 <= stop < length, _START_STOP_ERROR
        
        bases = self._bases

        # typical case
        if stop > start:
            my_bases = bases[start:stop]

        # boundary condition:  circular genome, Orf wraps around
        else:
            my_bases = bases[start:] + bases[:stop]
            
        return my_bases
 
    
    def is_sense(self):
        return self._is_sense
 
        
    def get_orfs(self, n):
        '''finds *only leftmost, longest* Orfs in all 3 alignments'''
        return self._get_orfs(sequence.getOrfEndsCircular, n)

    
    def get_all_orfs(self, n):
        '''finds Orfs of all sizes (including overlapping) in all 3 alignments'''
        return self._get_orfs(sequence.getAllOrfEndsCircular, n)

    
    def get_reverse_complement(self):
        rSeq = Sequence(self._get_reverse_bases(), not self.is_sense())
        return rSeq




########################################################################
# for analysis:  filtering and analyzing Orfs


class OrfCollection(object):
    
    def __init__(self, orfs):
        self._orfs = orfs
        
    def get_orfs(self):
        return self._orfs
    
    def filter(self, f):
        orfs = filter(f, self.get_orfs())
        return OrfCollection(orfs)
        
     
            
########################################################
# unit tests
########################################################

class OrfTest(unittest.TestCase):

    def setUp(self):
        self.orf = Orf(2, 14, 'ACGCTACCTTTCGCC', 'TAATAA', 'CTCA', True)
        self.obj = {
            'start': 2, 
            'stop': 14, 
            'bases': 'ACGCTACCTTTCGCC', 
            'upstream': 'TAATAA', 
            'downstream': 'CTCA',
            'is_sense': True
        }
        
    def test_start(self):
        self.assertEqual(2, self.orf.start)
        
    def test_stop(self):
        self.assertEqual(14, self.orf.stop)
        
    def test_bases(self):
        self.assertEqual('ACGCTACCTTTCGCC', self.orf.bases)
        
    def test_upstream(self):
        self.assertEqual('TAATAA', self.orf.upstream)
        
    def test_downstream(self):
        self.assertEqual('CTCA', self.orf.downstream)
        
    def test_is_sense(self):
        self.assertEqual(True, self.orf.is_sense)
        
    def test_to_JSON_object(self):
        self.assertEqual(self.obj, self.orf.to_JSON_object())
        
    def test_from_JSON_object(self):
        self.assertEqual(self.obj['start'], self.orf.start)
        self.assertEqual(self.obj['stop'], self.orf.stop)
        self.assertEqual(self.obj['bases'], self.orf.bases)
        self.assertEqual(self.obj['upstream'], self.orf.upstream)
        self.assertEqual(self.obj['downstream'], self.orf.downstream)
        self.assertEqual(self.obj['is_sense'], self.orf.is_sense)
        
        
    ################
    # test 'questionable' methods

    def testCodonsLength(self):
        self.assertEqual(5, len(self.orf.get_codons()))
        self.assertEqual('CTA', self.orf.get_codons()[1])

    def testPhobicity(self):
        self.assertTrue(False)

    def testPhobicityPeaks(self):
        self.assertTrue(False)


class SequenceTest(unittest.TestCase):

    ''' need to test: 1) forward, normal; 2) forward wrap; 3) reverse normal; 4) reverse wrap'''

    def setUp(self):
        self.seq = Sequence("ACGTAA" + "CCC" + "CTGAAAGGGTAG" + "ATGTTTTAC", True)
        self.orfs = self.seq.get_orfs(5)
        
        self.for_norm = Sequence('GA' + 'GCTAGCATCGAT' + 'TCGAT', True).build_orf(2, 14, 12)
        self.for_wrap = Sequence('CTA' + 'GCGAG' + 'CTAGCATCGATCGAA', True).build_orf(8, 3, 12)
        self.rev_norm = Sequence('GAGC' + 'TAGCATCGA' + 'TCGAT', False).build_orf(4, 13, 12)
        self.rev_wrap = Sequence('GAGCTA' + 'GCATC' + 'GATCGA',  False).build_orf(11, 6, 12)
        
    ########

    def test_get_bases(self):
        self.assertEqual("ACGTAACCCCTGAAAGGGTAGATGTTTTAC", self.seq.get_bases())
        self.assertEqual('AACCC', self.seq.get_bases(4, 9))
        self.assertEqual('TGTTTTACACG', self.seq.get_bases(22, 3))

    def test_is_sense(self):
        self.assertTrue(self.seq.is_sense())
        self.assertFalse(self.seq.get_reverse_complement().is_sense())

    def test_get_reverse_complement(self):
        self.assertEqual("GTAA" + "AACATCTACCCTTTCAGGGGTTAC" + "GT", self.seq.get_reverse_complement().get_bases())

    def test_get_orfs(self):
        ofs, ors = self.orfs, self.seq.get_reverse_complement().get_orfs(5)
        self.assertEqual(2, len(ofs))
        self.assertEqual(1, len(ors))
        self.assertEqual('AACCC', ofs[0].upstream)
        self.assertEqual('TAACC', ofs[1].downstream)
        self.assertEqual('GTTAC', ors[0].upstream)
        self.assertEqual(1, ors[0].start)

    def test_get_all_orfs(self):
        seq = Sequence('AATTAAAATAGA' + 'ATGGTGTGCTGC', False)
        # 'GCAGCACACCAT' + 'TCTATTTTAATT'
        rs, fs = seq.get_all_orfs(5), seq.get_reverse_complement().get_all_orfs(6)
        self.assertEqual(4, len(rs))
        self.assertEqual(1, len(fs))
        
        f = fs[0]
        self.assertEqual('TAATTG',  f.downstream)
        self.assertEqual('TTTTAA',  f.upstream)
        self.assertEqual('TTGCA', f.bases[:5])
        self.assertEqual((22, 19), (f.start, f.stop))
        
        myRs = set([(11, 20), (8, 20), (6, 15), (3, 15)])
        self.assertEqual(myRs, set([(r.start, r.stop) for r in rs]))
        
    ### from Orf

    def test_bases(self):
        self.assertEqual('GCTAGCATCGAT',       self.for_norm.bases)
        self.assertEqual('CTAGCATCGATCGAACTA', self.for_wrap.bases)
        self.assertEqual('TAGCATCGA',          self.rev_norm.bases)
        self.assertEqual('GATCGAGAGCTA',       self.rev_wrap.bases)
        
    def test_orf_upstream(self):
        self.assertEqual('GATGA',   self.for_norm.upstream[-5:])
        self.assertEqual('GCGAG',   self.for_wrap.upstream[-5:])
        self.assertEqual('GATGAGC', self.rev_norm.upstream[-7:])
        self.assertEqual('TC',      self.rev_wrap.upstream[-2:])
        
    def test_orf_downstream(self):
        self.assertEqual('TCGATGA',       self.for_norm.downstream[:7])
        self.assertEqual('GCG',           self.for_wrap.downstream[:3])
        self.assertEqual('TCGAT',         self.rev_norm.downstream[:5])
        self.assertEqual('GCATCGATCGAG',  self.rev_wrap.downstream[:12])
    
    def test_start_stop(self):
        self.assertEqual((2, 14), (self.for_norm.start, self.for_norm.stop))
        self.assertEqual((8, 3),  (self.for_wrap.start, self.for_wrap.stop))
        self.assertEqual((13, 4), (self.rev_norm.start, self.rev_norm.stop))
        self.assertEqual((5, 10), (self.rev_wrap.start, self.rev_wrap.stop))
        


class OrfCollectionTest(unittest.TestCase):

    def setUp(self):
        self.oc = OrfCollection([
            Orf(14, 27, 'ACGGGGTTTCCC', 'CCC', 'ATT', True),
            Orf(21, 3, 'CGAGAATAG', 'GGG', '', False),
            Orf(18, 900, 'ACGGGGTTTCCC', '', '', True),
            Orf(45, 54, 'ACGGGGTTTCCC', '', '', False)
        ])
        
    def test_get_orfs(self):
        self.assertEqual(4, len(self.oc.get_orfs()))
    
    def test_filter(self):
        c = self.oc.filter
        cs = c(lambda o: o.start < 20), c(lambda o: o.bases == 'CGAGAATAG'), c(lambda o: o.is_sense)
        self.assertEqual([2, 1, 2], map(lambda oc: len(oc.get_orfs()), cs))
        



testClasses = [OrfTest, SequenceTest, OrfCollectionTest]
