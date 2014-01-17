import sys
from Bio import Motif
from Bio.Alphabet import IUPAC

motif = sys.argv[1]

m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
from Bio.Seq import Seq

with open("matches.txt", "r") as inFile:
	for line in inFile:
		m.add_instance(Seq(line[:-1], m.alphabet))

p = m.pwm()
with open("outputs/logos/"+motif+"_pwm.txt", "w") as pwmFile:
	print >> pwmFile, "  G     A     T     C"
	for pos in range(len(m)):
		print >> pwmFile, "%.3f %.3f %.3f %.3f" % (p[pos]["G"], p[pos]["A"], p[pos]["T"], p[pos]["C"])

with open("outputs/logos/"+motif+"_transfac.txt", "w") as transFile:
	print >> transFile, m.format("transfac")

m.weblogo("outputs/logos/"+motif+"logo.png")

