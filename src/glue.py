import unittest
import bacillussubtilis168 as bs




orfs = getOrfs(makeCodons(bs.sequence))

myorfs = filter(lambda o: len(o.codons) >= 50 and len(o.codons) <= 80, orfs)