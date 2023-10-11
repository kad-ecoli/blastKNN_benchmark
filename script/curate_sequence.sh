#!/bin/bash
# download sequence
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

cd $dbdir
$bindir/makeblastdb -in   $dbdir/curated/sequence/uniprot_sprot_exp.fasta -dbtype prot -title uniprot_sprot_exp -parse_seqids -logfile $dbdir/tmp/makeblastdb.log
$bindir/diamond prepdb -d $dbdir/curated/sequence/uniprot_sprot_exp.fasta
$bindir/mmseqs/bin/mmseqs createdb $dbdir/curated/sequence/uniprot_sprot_exp.fasta $dbdir/curated/sequence/uniprot_sprot_exp
$bindir/hhsuite2/scripts/kClust2db.py $dbdir/curated/sequence/uniprot_sprot_exp.fasta $dbdir/curated/sequence/uniprot_sprot_exp.hh2
