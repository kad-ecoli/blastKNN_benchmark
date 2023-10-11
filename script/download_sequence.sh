#!/bin/bash
# download sequence
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

cd $dbdir
mkdir -p $dbdir/raw/sequence/
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.212.gz -O $dbdir/raw/sequence/goa_uniprot_all.gaf.train.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz -O $dbdir/raw/sequence/uniprot_sprot.fasta.gz

mkdir -p $dbdir/curated/sequence/
mkdir -p $dbdir/tmp
zcat $dbdir/raw/sequence/goa_uniprot_all.gaf.train.gz|grep "^UniProtKB" | grep -P "(\tEXP\t)|(\tIDA\t)|(\tIPI\t)|(\tIMP\t)|(\tIGI\t)|(\tIEP\t)|(\tTAS\t)|(\tIC\t)|(\tHTP\t)|(\tHDA\t)|(\tHMP\t)|(\tHGI\t)|(\tHEP\t)" > $dbdir/tmp/goa_uniprot_exp.gaf
$bindir/curate_GAF $dbdir/tmp/goa_uniprot_exp.gaf $dbdir/curated/sequence/uniprot_sprot_exp $dbdir/curated/ontology/alt_id.csv $dbdir/curated/ontology/is_a.csv $dbdir/evidence_exp.list

zcat $dbdir/raw/sequence/uniprot_sprot.fasta.gz $dbdir/raw/cafa-5-protein-function-prediction/train_sequences.fasta.gz|cut -f1,2 -d'|'|sed 's/>sp|/>/g' > $dbdir/tmp/uniprot_sprot_exp.fasta
if [ -s $dbdir/curated/sequence/uniprot_sprot_exp.fasta ];then
   cat $dbdir/curated/sequence/uniprot_sprot_exp.fasta >> $dbdir/tmp/uniprot_sprot_exp.fasta
fi
$bindir/fasta2tsv $dbdir/tmp/uniprot_sprot_exp.fasta - |sort -k1,1 -u > $dbdir/tmp/uniprot_sprot_exp.tsv
$bindir/tsv2fasta $dbdir/tmp/uniprot_sprot_exp.tsv $dbdir/tmp/uniprot_sprot_exp.fasta
$bindir/fasta2miss $dbdir/curated/sequence/uniprot_sprot_exp.list $dbdir/tmp/uniprot_sprot_exp.fasta $dbdir/tmp/miss.list $dbdir/curated/sequence/uniprot_sprot_exp.fasta
for target in `cat $dbdir/tmp/miss.list`;do
   curl "https://rest.uniprot.org/unisave/$target?format=fasta&versions=1" |cut -f1,2 -d'|' |sed 's/>sp|/>/g' |sed 's/>tr|/>/g' 
done >> $dbdir/curated/sequence/uniprot_sprot_exp.fasta

$bindir/makeblastdb -in   $dbdir/curated/sequence/uniprot_sprot_exp.fasta -dbtype prot -title uniprot_sprot_exp -parse_seqids -logfile $dbdir/tmp/makeblastdb.log
$bindir/diamond prepdb -d $dbdir/curated/sequence/uniprot_sprot_exp.fasta

$bindir/calculate_ic $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a $dbdir/curated/ontology/naive.F $dbdir/curated/ontology/is_a.csv $dbdir/curated/ontology/name.csv
$bindir/calculate_ic $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a $dbdir/curated/ontology/naive.P $dbdir/curated/ontology/is_a.csv $dbdir/curated/ontology/name.csv
$bindir/calculate_ic $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a $dbdir/curated/ontology/naive.C $dbdir/curated/ontology/is_a.csv $dbdir/curated/ontology/name.csv

cd $dbdir
$bindir/makeblastdb -in   $dbdir/curated/sequence/uniprot_sprot_exp.fasta -dbtype prot -title uniprot_sprot_exp -parse_seqids -logfile $dbdir/tmp/makeblastdb.log
$bindir/diamond prepdb -d $dbdir/curated/sequence/uniprot_sprot_exp.fasta
$bindir/mmseqs/bin/mmseqs createdb $dbdir/curated/sequence/uniprot_sprot_exp.fasta $dbdir/curated/sequence/uniprot_sprot_exp
$bindir/hhsuite2/scripts/kClust2db.py $dbdir/curated/sequence/uniprot_sprot_exp.fasta $dbdir/curated/sequence/uniprot_sprot_exp.hh2
