import finder
import filterer
import model


positives = [
    ('ATAAAG', 'ATGAAA', 'TAATTC', 3697805),
    ('TTTTTT', 'ATGCTT', 'TAGTTT', 1907012),
    ('GTCAGC', 'ATGATC', 'TGAACA', 847481),
    ('TGAAGC', 'ATGAAG', 'TAAAGG', 419983),
    ('GATTCT', 'ATGAAA', 'TGAATG', 1120827),
    ('TAATAG', 'ATGAGT', 'TAGAAA', 435214) # can't find
]


def findPositives():
    orfs = finder.findAllOrfs()
    oas = model.OACollection([model.OrfAnalysis(o.toJSONObject(30)) for o in orfs])
    for up, seq, down, ix in positives:
        print "trying", up, seq, down, ix
        matches = oas.findOrf(orfFilter = filterer.matchBases(seq = seq, up = up, down = down)).getOrfAnals()
        assert len(matches) == 1, "matched <%s> Orfs" % str(len(matches))
        assert matches[0].orf['start'] == ix, "start was actually <%s>" % str(matches[0].orf['start'])
        