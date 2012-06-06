import bacillussubtilis168 as bs
import sequence


forward = bs.bases
reverse = sequence.reverseComplement(forward)



def findOnce(sear, seq):
    f, r = seq.find(sear), seq.rfind(sear)
    if f == r:
        return f
    else:
        return 'No!'