CC=g++
CFLAGS=-O3
LDFLAGS=-static

PROGRAM=obo2csv backpropagate curate_GAF fasta2tsv tsv2fasta fasta2miss subset_GOA subset_tsv cd-hit calculate_ic

all: ${PROGRAM}

obo2csv: obo2csv.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

backpropagate: backpropagate.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

curate_GAF: curate_GAF.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

curate_sifts_go: curate_sifts_go.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

curate_ppi_go: curate_ppi_go.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2tsv: fasta2tsv.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

tsv2fasta: tsv2fasta.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2miss: fasta2miss.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cif2fasta: cif2fasta.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cif2chain: cif2chain.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

subset_GOA: subset_GOA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

filter_GOA: filter_GOA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

subset_tsv: subset_tsv.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

subset_ppi: subset_ppi.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

calculate_ic: calculate_ic.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

blast2homo: blast2homo.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cd-hit: cdhit_src/cdhit.c++ cdhit_src/cdhit-common.c++ cdhit_src/cdhit-utility.c++ cdhit_src/cdhit-common.h cdhit_src/cdhit-utility.h
	cd cdhit_src; make cd-hit; cd .. ; cp cdhit_src/cd-hit cd-hit

clean:
	rm ${PROGRAM}
