#!/bin/bash
# download sequence
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

cd $dbdir
mkdir -p $dbdir/raw/sequence/
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.216.gz -O $dbdir/raw/sequence/goa_uniprot_all.gaf.validate.gz

mkdir -p $dbdir/curated/sequence/
mkdir -p $dbdir/tmp
zcat $dbdir/raw/sequence/goa_uniprot_all.gaf.validate.gz | grep "^UniProtKB" | grep -P "(\tEXP\t)|(\tIDA\t)|(\tIPI\t)|(\tIMP\t)|(\tIGI\t)|(\tIEP\t)|(\tTAS\t)|(\tIC\t)|(\tHTP\t)|(\tHDA\t)|(\tHMP\t)|(\tHGI\t)|(\tHEP\t)" > $dbdir/tmp/validate.gaf
$bindir/curate_GAF $dbdir/tmp/validate.gaf $dbdir/tmp/validate $dbdir/curated/ontology/alt_id.csv $dbdir/curated/ontology/is_a.csv $dbdir/evidence_exp.list
$bindir/make_target.py $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a $dbdir/tmp/validate.F.is_a $dbdir/curated/sequence/validate.F.is_a
$bindir/make_target.py $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a $dbdir/tmp/validate.P.is_a $dbdir/curated/sequence/validate.P.is_a
$bindir/make_target.py $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a $dbdir/tmp/validate.C.is_a $dbdir/curated/sequence/validate.C.is_a
cat $dbdir/curated/sequence/validate.F.is_a $dbdir/curated/sequence/validate.P.is_a $dbdir/curated/sequence/validate.C.is_a|cut -f1|sort|uniq > $dbdir/curated/sequence/validate.list
$bindir/combineFPC.py $dbdir/curated/sequence/validate.list $dbdir/curated/sequence/validate.F.is_a $dbdir/curated/sequence/validate.P.is_a $dbdir/curated/sequence/validate.C.is_a $dbdir/curated/sequence/validate.3.is_a

zcat $dbdir/raw/sequence/uniprot_sprot.fasta.gz $dbdir/raw/cafa-5-protein-function-prediction/train_sequences.fasta.gz|cut -f1,2 -d'|'|sed 's/>sp|/>/g' > $dbdir/tmp/uniprot_sprot_exp.fasta
cat $dbdir/curated/sequence/uniprot_sprot_exp.fasta $dbdir/curated/sequence/validate.fasta >> $dbdir/tmp/uniprot_sprot_exp.fasta
cat $dbdir/tmp/uniprot_sprot_exp.fasta | $bindir/fasta2tsv - - |sort -u -k1,1 | $bindir/tsv2fasta - - | $bindir/fasta2miss $dbdir/curated/sequence/validate.list - $dbdir/tmp/miss.list $dbdir/tmp/validate.fasta
mv $dbdir/tmp/validate.fasta $dbdir/curated/sequence/validate.fasta
for target in `cat $dbdir/tmp/miss.list`;do
   curl "https://rest.uniprot.org/unisave/$target?format=fasta&versions=1" |cut -f1,2 -d'|' |sed 's/>sp|/>/g' |sed 's/>tr|/>/g' 
done >> $dbdir/curated/sequence/validate.fasta
