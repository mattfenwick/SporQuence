import model


def read(name):
    '''String -> String'''
    f = open(name, 'r')
    contents = f.read()
    f.close()
    return contents


def parse(fasta):
    '''String -> Sequence'''
    lines = fasta.split('\n')

    assert lines[0][0] == '>'  # fasta format first line ...

    cleaned = ''.join(lines[1:]) # ... is just metadata

    for c in cleaned:
        assert c in set("ACGT")  # wonder if 'set' makes a performance difference

    return model.Sequence(cleaned)


def writeCleaned(infile, outfile):
    seq = parse(read(infile))
    outf = open(outfile, 'w')
    outf.write(seq.bases)
    outf.close()


def getCodons(inpath):
    '''String -> [Codons]'''
    f = open(inpath, 'r')
    codons = makeCodons(f.read())
    f.close()
    return codons

###########################################
