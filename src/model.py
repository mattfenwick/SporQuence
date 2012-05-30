import yaml


# domain objects

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
