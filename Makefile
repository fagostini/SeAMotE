# Compiler
CC = gcc
# Compiler flags
CCFLAGS = -ansi -pthread
# Linker flags
LNFLAGS = -O3 -Wall -pedantic
# System-specific flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CCFLAGS += -D LINUX -lm
endif
ifeq ($(UNAME_S),Darwin)
	CCFLAGS += -D OSX
endif
# Arch-specific flags
UNAME_P := $(shell uname -m)
ifeq ($(UNAME_P),x86_64)
	CCFLAGS += -m64
endif

# Log file
LOG = test.log
# Sources
SRC = main.c
# Objects
OBJ = main.o
# Headers
HDR = main.h
# Libraries
LIBS = my_library.lib RNA_lib.lib DNA_lib.lib
# Binaries
BIN = main
# Files
QUICKFILE = quick_posi.oneline quick_nega.oneline
RNAFILE = posi10_mRNA.oneline nega10_mRNA.oneline
REALRNA = positive_mRNA.oneline negative_mRNA.oneline
DNAFILE = posi10_DNA.oneline nega10_DNA.oneline
REALDNA = positive_DNA.oneline negative_DNA.oneline

# Make options
all: $(BIN)

$(BIN):	$(OBJ)
	$(CC) $(OBJ) $(LNFLAGS) -o $(BIN)

$(OBJ):	$(SRC) $(HDR) $(LIBS)
	$(CC) $(CCFLAGS) $(SRC) -c

clean:
	@rm -f $(OBJ) $(BIN) $(LOG)

# Testing options
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

runRNA:	$(BIN) $(REALRNA)
	@./$(BIN) -p positive_mRNA.oneline -n negative_mRNA.oneline -t 0.9 > $(LOG)
	@echo "The algorithm RNA PASSED" || echo "The algorithm RNA FAILED"

runDNA:	$(BIN) $(REALDNA)
	@./$(BIN) -p positive_DNA.oneline -n negative_DNA.oneline -t 0.9 > $(LOG)
	@echo "The algorithm DNA PASSED" || echo "The algorithm DNA FAILED"
