CC = gcc
CCFLAGS = -ansi -O3 -pthread -Wall -pedantic
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CCFLAGS += -D LINUX -lm
endif
ifeq ($(UNAME_S),Darwin)
	CCFLAGS += -D OSX
endif
UNAME_P := $(shell uname -m)
ifeq ($(UNAME_P),x86_64)
	CCFLAGS += -m64
endif

OBJ = main.o
LIBS = my_library.lib RNA_lib.lib DNA_lib.lib
BIN = main
LOG = test.log
QUICKFILE = quick_posi.oneline quick_nega.oneline
RNAFILE = posi10_mRNA.oneline nega10_mRNA.oneline
REALRNA = positive_mRNA.oneline negative_mRNA.oneline
DNAFILE = posi10_DNA.oneline nega10_DNA.oneline
REALDNA = positive_DNA.oneline negative_DNA.oneline

# PROTFILE = posiPROT10.oneline negaPROT10.oneline
# REALPROT = positive_protein.oneline negative_protein.oneline

all: $(BIN)

$(BIN):	$(OBJ)
	$(CC) $(OBJ) $(CCFLAGS) -o $(BIN)

main.o:	main.c main.h $(LIBS)
	$(CC) main.c -c

clean:
	@rm -f $(OBJ) $(BIN) $(LOG)

runRNA:	$(BIN) $(REALRNA)
	@./$(BIN) -p positive_mRNA.oneline -n negative_mRNA.oneline -t 0.9 > $(LOG)
	@echo "The algorithm RNA PASSED" || echo "The algorithm RNA FAILED"

runDNA:	$(BIN) $(REALDNA)
	@./$(BIN) -p positive_DNA.oneline -n negative_DNA.oneline -t 0.9 > $(LOG)
	@echo "The algorithm DNA PASSED" || echo "The algorithm DNA FAILED"

#runPROT:	$(BIN) $(REALPROT)
#	@./$(BIN) $(REALPROT) > $(LOG)
#	@echo "The algorithm PASSED" || echo "The algorithm FAILED"
	
tests: testRNA testDNA

quick:	$(BIN) $(QUICKFILE)
	@./$(BIN) -p quick_posi.oneline -n quick_nega.oneline -t 0.8 -a rna > $(LOG)
	@echo "Test Quick PASSED" || echo "Test Quick FAILED"
	
testRNA:	$(BIN) $(RNAFILE)
	@./$(BIN) -p posi10_mRNA.oneline -n nega10_mRNA.oneline -t 0.95 -a rna > $(LOG)
	@echo "Test RNA PASSED" || echo "Test RNA FAILED"

testDNA:	$(BIN) $(DNAFILE)
	@./$(BIN) -p posi10_DNA.oneline -n nega10_DNA.oneline -t 0.95 -a rna > $(LOG)
	@echo "Test DNA PASSED" || echo "Test DNA FAILED"

#testPROT:	$(BIN) $(PROTFILE)
#	@./$(BIN) $(PROTFILE) > $(LOG) 
#	@echo "Test PROTEIN PASSED" || echo "Test PROTEIN FAILED"
