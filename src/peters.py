import orfind.bacillussubtilis168 as bs
import orfind.sequence as sequence


forward = bs.sequence
reverse = sequence.reverseComplement(forward)



def findOnce(sear, seq):
    f, r = seq.find(sear), seq.rfind(sear)
    if f == r:
        return f
    else:
        return 'No!'