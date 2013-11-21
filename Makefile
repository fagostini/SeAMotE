# Compiler
CC = gcc
# Compiler flags
CCFLAGS = -ansi
# Object flags
OBFLAGS = -O3 -Wall -pedantic
# System-specific flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	LNFLAGS += -D LINUX -lm -lpthread
endif
ifeq ($(UNAME_S),Darwin)
	LNFLAGS += -D OSX -pthread
endif
# Arch-specific flags
UNAME_P := $(shell uname -m)
ifeq ($(UNAME_P),x86_64)
	CCFLAGS += -m64
endif

# Log file
LOG = test.log
# Sources
OLD_SRC = main.c
SRC = motifSearch.c
# Objects
OLD_OBJ = main.o
OBJ = motifSearch.o
# Headers
OLD_HDR = main.h
HDR = motifSearch.h
# Libraries
LIBS = my_library.lib RNA_lib.lib DNA_lib.lib
# Binaries
OLD_BIN = main
BIN = motifSearch
# Files
QUICKFILE = quick_posi.oneline quick_nega.oneline
RNAFILE = posi10_mRNA.oneline nega10_mRNA.oneline
REALRNA = positive_mRNA.oneline negative_mRNA.oneline
DNAFILE = posi10_DNA.oneline nega10_DNA.oneline
REALDNA = positive_DNA.oneline negative_DNA.oneline

# Make options
all: $(BIN) $(OLD_BIN)

$(BIN):	$(OBJ)
	$(CC) $(OBJ) $(OBFLAGS) -o $(BIN) $(LNFLAGS)

$(OBJ):	$(SRC) $(HDR) $(LIBS)
	$(CC) $(SRC) $(CCFLAGS) -c

$(OLD_BIN):	$(OLD_OBJ)
	$(CC) $(OLD_OBJ) $(OBFLAGS) -o $(OLD_BIN) $(LNFLAGS)

$(OLD_OBJ):	$(OLD_SRC) $(OLD_HDR) $(LIBS)
	$(CC) $(OLD_SRC) $(CCFLAGS) -c

clean:
	@rm -f $(OBJ) $(BIN) $(OLD_OBJ) $(OLD_BIN) $(LOG)

# Testing options
tests: testRNA testDNA

old:	$(OLD_BIN) $(QUICKFILE)
	@./$(OLD_BIN) -p quick_posi.oneline -n quick_nega.oneline -t 0.8 -a rna > $(LOG)
	@echo "Test Quick PASSED" || echo "Test Quick FAILED"

quick:	$(BIN) $(QUICKFILE)
	@./$(BIN) -p quick_posi.oneline -n quick_nega.oneline -t 0.8 -a rna > $(LOG)
	@echo "Test Quick PASSED" || echo "Test Quick FAILED"
	
testRNA:	$(BIN) $(RNAFILE)
	@./$(BIN) -p posi10_mRNA.oneline -n nega10_mRNA.oneline -t 0.95 -a rna > $(LOG)
	@echo "Test RNA PASSED" || echo "Test RNA FAILED"

testDNA:	$(BIN) $(DNAFILE)
	@./$(BIN) -p posi10_DNA.oneline -n nega10_DNA.oneline -t 0.95 -a rna > $(LOG)
	@echo "Test DNA PASSED" || echo "Test DNA FAILED"

runRNA:	$(BIN) $(REALRNA)
	@./$(BIN) -p positive_mRNA.oneline -n negative_mRNA.oneline -t 0.9 > $(LOG)
	@echo "The algorithm RNA PASSED" || echo "The algorithm RNA FAILED"

runDNA:	$(BIN) $(REALDNA)
	@./$(BIN) -p positive_DNA.oneline -n negative_DNA.oneline -t 0.9 > $(LOG)
	@echo "The algorithm DNA PASSED" || echo "The algorithm DNA FAILED"
