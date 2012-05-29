import model


path = 'sequences/al009126/cleaned/seq.txt'


def getCodons(inpath):
    f = open(inpath, 'r')
    codons = model.CodonSequence(f.read())
    f.close()
    return codons


def makeCodons(bases):
    '''[Base] -> [Codons]'''
    assert len(bases) % 3 == 0

    codons = []
    x = 0
    while x < len(bases):
        codons.append(bases[x:x + 3])
        x += 3

    assert len(bases) / 3 == len(codons)

    return codons


def getCodons2(inpath):
    '''String -> [Codons]'''
    f = open(inpath, 'r')
    codons = makeCodons(f.read())
    f.close()
    return codons


def getOrfs(codons):
    '''[Codons] -> [ORFs]'''
    inORF = False
    orfs = []
    ix = 0
    for c in codons:
        if not inORF and c in ["GTG", "ATG", "TTG"]:
            orfs.append(model.ORF(ix, [], None))
            inORF = True
        if inORF and c in ["TAA", "TAG", "TGA"]:
            inORF = False
            orfs[-1].stop = c
        if inORF:
            orfs[-1].codons.append(c)
        ix += 1
    return orfs
