import unittest




def makeCodons(bases):
    '''[Base] -> [Codon]'''
    if len(bases) % 3 != 0:
        raise ValueError("number of bases must be divisible by 3")

    codons = []
    x = 0
    while x < len(bases):
        codons.append(bases[x:x + 3])
        x += 3

    assert len(bases) / 3 == len(codons)

    return codons



STARTS = ["GTG", "ATG", "TTG", "CTG"]
STOPS = ["TAA", "TAG", "TGA"]

def getOrfEndsLinear(codons):
    '''[Codon] -> [(Int, Int)]
    
    Finds open reading frames in a codon sequence.
    Assumes that the sequence is linear.
    Finds longest possible ORFS -- i.e. in 
    'ATGATGCCCTAA', would find one ORF spanning
    the entire sequence.
    
    '''
    inORF, orfEnds = False, []
    ix, start = 0, None

    for c in codons:
        if not inORF and c in STARTS:
            start = ix
            inORF = True
        if inORF and c in STOPS:
            inORF = False
            orfEnds.append((start, ix))
        ix += 1

    if inORF:
        orfEnds.append((start, None))

    return orfEnds


def getOrfEndsCircular(codons):
    '''[Codon] -> [(Int, Int)]
    
    Finds open reading frames in a codon sequence.
    Assumes that the sequence is circular.
    Finds longest possible ORFS -- i.e. in 
    'TAAATGATGCCC', would find one ORF spanning
    the entire sequence.
    
    '''
    i, firstStop, firstStart, length, orfEnds = 0, None, None, len(codons), []

    # find first stop
    while i < length:
        if codons[i] in STOPS:
            firstStop = i
            break
        elif codons[i] in STARTS:
            firstStart = i
        i += 1

    # if there aren't any stops, just return
    if firstStop is None:
        # log this: assert firstStart is not None, "there are no stops so there should also be no starts"
        return []

    # begin right after the first stop
    i, start, inORF = firstStop + 1, None, False

    # iterate through the codons:   stop when reach marker again
    while i != firstStop:
        # if not in an ORF:
        #   if a start is found, then change state to 'inORF'
        #   else do nothing
        if not inORF and codons[i] in STARTS:
            start, inORF = i, True

        # if in an ORF
        #   if a stop is found, then change state to 'not inORF'
        #   else no change
        if inORF and codons[i] in STOPS:
            inORF = False
            orfEnds.append((start, i))

        #   if index overflows, start it back at 0
        i += 1
        i %= length

    if inORF:
        orfEnds.append((start, i))

    return orfEnds




def readToNextStop(codons, index):
    '''[Codons] -> Int -> (Int, Int)'''
    assert codons[index] in STARTS
        
    i, stop, length = index + 1, None, len(codons)
    while i != index:
        i %= length
        if codons[i] in STOPS:
            return (index, i)
        i += 1
        
    raise ValueError("no stop codon found")


def getAllOrfEndsCircular(codons):
    '''[Codon] -> [(Int, Int)]
    
    Finds open reading frames in codon sequence.
    Assumes that the sequence is circular.
    Finds all possible ORFS -- i.e. in 
    'TAAATGATGCCC', would find two ORFs.
    
    '''
    i, length, orfEnds = 0, len(codons), []
    while i < length:
        if codons[i] in STARTS:
            orfEnds.append(readToNextStop(codons, i))
        i += 1
    return orfEnds
    



COMPLEMENTS = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

def reverseComplement(bases):
    i, length = 0, len(bases)
    seq = [None] * length
    while i < length:
        seq[-i - 1] = COMPLEMENTS[bases[i]]
        i += 1
    return ''.join(seq)


########################################################
# unit tests
########################################################

class CodonsTest(unittest.TestCase):

    def setUp(self):
        pass

    @unittest.expectedFailure
    def testGetCodonsBadLength(self):
        bases = 'ACGTCCTG'
        codons = makeCodons(bases)

    def testGetCodons(self):
        bases = 'ACGTCCTGA'
        codons = makeCodons(bases)
        self.assertEqual(codons[0], "ACG")
        self.assertEqual(codons[2], "TGA")


class ComplementTest(unittest.TestCase):

    def setUp(self):
        pass

    def testRC(self):
        bases = 'ACGTCCTGA'
        self.assertEqual('TCAGGACGT', reverseComplement(bases))


class LinearOrfsTest(unittest.TestCase):

    def setUp(self):
        pass

    def testGetOrfs(self):
        bases = 'ACTGTGACCTCATATTAGGGGTTT'
        codons = makeCodons(bases)
        orfEnds = getOrfEndsLinear(codons)
        self.assertEqual(4, orfEnds[0][1] - orfEnds[0][0])
        self.assertEqual(1, orfEnds[0][0], "remember: ORF indices are 0-indexed!")

    def testGetOrfsTwoStarts(self):
        codons = makeCodons('ACTGTGACCTCATATCAGGGGTTT' + 'ACTGTGACCTCATATTAGGGGTTT' + 'ACTTTGACCTCATATTGGTGAATT')
        orfEnds = getOrfEndsLinear(codons)
        self.assertEqual(2, len(orfEnds), 'total number of orfs')
        self.assertEqual(1, orfEnds[0][0], 'start of 1st orf')
        self.assertEqual(12, orfEnds[0][1] - orfEnds[0][0], 'length of 1st orf')
        self.assertEqual(17, orfEnds[1][0], 'start of 2nd orf')
        self.assertEqual(5, orfEnds[1][1] - orfEnds[1][0], 'length of 2nd orf')

    def testGetOrfsTwoStops(self):
        # 1, start, 7, start, 4, stop, 2, stop, 3, start, stop, 1 -> 2 ORFs
        bases = ['ACT', 'GTG', 'ACCTCATATCAGGGGTTTACT', 'GTG', 'ACCTCATATCCC', 
                 'TGA', 'CCCACT', 'TAA', 'ACCTCATAT', 'TTG', 'TGA', 'ATT']
        codons = makeCodons(''.join(bases))
        orfs = getOrfEndsLinear(codons)
        self.assertEqual(2, len(orfs), 'total number of orfs')
        self.assertEqual(1, orfs[0][0], 'start of 1st orf')
        self.assertEqual(13, orfs[0][1] - orfs[0][0], 'length of 1st orf')
        self.assertEqual(21, orfs[1][0], 'start of 2nd orf')
        self.assertEqual(1, orfs[1][1] - orfs[1][0], 'length of 2nd orf')

    def testGetOrfsNoEndEnd(self):
        bases = ['ACC', 'GTG', 'CAC', 'TTT']
        orfs = getOrfEndsLinear(makeCodons(''.join(bases)))
        self.assertEqual((1, None), orfs[0])


class CircularOrfsTest(unittest.TestCase):

    def setUp(self):
        pass

    def testGetOrfs(self):
        bases = 'ACT' + 'GTGACCTCATATTAG' + 'GGGTTT'
        codons = makeCodons(bases)
        orfEnds = getOrfEndsCircular(codons)
        self.assertEqual(4, orfEnds[0][1] - orfEnds[0][0])
        self.assertEqual(1, orfEnds[0][0], "remember: ORF indices are 0-indexed!")

    def testGetOrfsTwoStarts(self):
        codons = makeCodons('ACT' + 'GTGACCTCATATCAGGGGTTTACT' + 'GTGACCTCATATTAG' + 'GGGTTTACT' + 'TTGACCTCATATTGGTGA' + 'ATT')
        orfEnds = getOrfEndsCircular(codons)
        self.assertEqual(set(orfEnds), set([(1, 13), (17, 22)]))

    def testGetOrfsTwoStops(self):
        # 1, start, 7, start, 4, stop, 2, stop, 3, start, stop, 1 -> 2 ORFs
        bases = ['ACT', 'GTG', 'ACCTCATATCAGGGGTTTACT', 'GTG', 'ACCTCATATCCC', 
                 'TGA', 'CCCACT', 'TAA', 'ACCTCATAT', 'TTG', 'TGA', 'ATT']
        codons = makeCodons(''.join(bases))
        orfs = getOrfEndsCircular(codons)
        self.assertEqual(set(orfs), set([(1, 14), (21, 22)]))

    def testGetOrfsWrapAround(self):
        bases = 'ACCTGA' + 'CCGCAC' + 'TTGTTT'
        orfs = getOrfEndsCircular(makeCodons(bases))
        self.assertEqual((4, 1), orfs[0])
        
        
class CircularAllOrfsTest(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def testGetOrfs(self):
        codons = makeCodons('ACTGTGACCTTGTTTTGAACT')
        orfEnds = getAllOrfEndsCircular(codons)
        self.assertEqual(set(orfEnds), set([(1, 5), (3, 5)]))
    
    def testWraparound(self):
        codons = makeCodons('TAGATTATGGTGCTG')
        orfEnds = getAllOrfEndsCircular(codons)
        self.assertEqual(set(orfEnds), set([(2, 0), (4, 0), (3, 0)]))
    
    @unittest.expectedFailure
    def testNostop(self):
        codons = makeCodons('ACTATGGTGCTG')
        orfEnds = getAllOrfEndsCircular(codons)
        self.assertEqual(len(orfEnds), 0)
    
    def testNostartNostop(self):
        codons = makeCodons('TAGATTCCCCTCTAT')
        orfEnds = getAllOrfEndsCircular(codons)
        self.assertEqual(orfEnds, [])
        
    def testAllOrfsAllAlignments(self):
        bases = 'AAATAAAATAGA' + 'ATGGTGTGCTGC'
        o1, o2, o3 = map(lambda x: getAllOrfEndsCircular(makeCodons(x)), [bases, bases[1:] + bases[:1], bases[2:] + bases[:2]])
        self.assertEqual([2, 0, 2], map(len, [o1, o2, o3]))
        self.assertEqual((4, 1), o1[0])
        self.assertEqual((5, 2), o3[0])
        self.assertEqual((6, 2), o3[1])
    

testClasses = [CodonsTest, ComplementTest, LinearOrfsTest, CircularOrfsTest, CircularAllOrfsTest]
