


class Codon(object):

    def __init__(self, bases):
        assert len(bases) == 3, "codon needs 3 bases"
        self.bases = bases

    def toJSON(self):
        return {
            'bases': self.bases
        }

    def __repr__(self):
        return str(self.toJSON())


class CodonSequence(object):

    def __init__(self, bases):
        assert (len(bases) % 3) == 0, "sequence length must be divisible by 3"

        self.codons = []
        x = 0
        while x < len(bases):
            self.codons.append(Codon(bases[x:x+3]))
            x += 3

        assert len(self.codons) == len(bases) / 3

    def toJSON(self):
        return {
            'bases': self.bases
        }

    def __repr__(self):
        return str(self.toJSON())


class ORF(object):

    def __init__(self, index, codons, stop):
        self.index = index
        self.codons = codons
        self.stop = stop

    def toJSON(self):
        return {
            'index': self.index,
            'codons': self.codons,
            'stop': self.stop
        }

    def __repr__(self):
        return str(self.toJSON())
