import yaml


class Codon(object):

    def __init__(self, bases):
        if len(bases) != 3:
            raise ValueError("codon needs 3 bases")
        self.bases = bases

    def toJSON(self):
        return {
            'bases': self.bases
        }

    def __repr__(self):
        return yaml.dump(self.toJSON())


class ORF(object):

    def __init__(self, index, codons, stop):
        self.index = index
        self.codons = codons
        self.stop = stop

    def toJSON(self):
        return {
            'index': self.index,
            'codons': [c.toJSON() for c in self.codons],
            'stop': self.stop.toJSON()
        }

    def __repr__(self):
        return yaml.dump(self.toJSON())
