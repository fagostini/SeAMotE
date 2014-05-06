import sys
from Bio import Motif
from Bio.Alphabet import IUPAC

motif = sys.argv[1]

m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
ml = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
from Bio.Seq import Seq

myLimit = 0
with open("matches.txt", "r") as inFile:
	for line in inFile:
		m.add_instance(Seq(line[:-1], m.alphabet))
		if myLimit < 35000:	
			ml.add_instance(Seq(line[:-1], ml.alphabet))
		myLimit += 1

p = m.pwm()
with open("logos/"+motif+"_pwm.txt", "w") as pwmFile:
	print >> pwmFile, "  A     C     G     T"
	for pos in range(len(m)):
		print >> pwmFile, "%.6f %.6f %.6f %.6f" % (p[pos]["A"], p[pos]["C"], p[pos]["G"], p[pos]["T"])

pl = ml.pwm()
ml.weblogo("logos/"+motif+"_logo.png")

