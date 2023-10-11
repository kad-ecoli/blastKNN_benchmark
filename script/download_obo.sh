#!/bin/bash
# download go-basic.obo
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

cd $dbdir

mkdir -p $dbdir/raw/ontology/
#wget http://purl.obolibrary.org/obo/go/go-basic.obo -O $dbdir/raw/ontology/go-basic.obo

mkdir -p $dbdir/curated/ontology/
$bindir/obo2csv $dbdir/raw/ontology/go-basic.obo $dbdir/curated/ontology/is_a.csv $dbdir/curated/ontology/name.csv $dbdir/curated/ontology/alt_id.csv

#mkdir -p $dbdir/tmp
#$bindir/IA.py F $dbdir/curated/ontology/is_a.csv $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a $dbdir/tmp/IA.F
#$bindir/IA.py P $dbdir/curated/ontology/is_a.csv $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a $dbdir/tmp/IA.P
#$bindir/IA.py C $dbdir/curated/ontology/is_a.csv $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a $dbdir/tmp/IA.C
#cat $dbdir/tmp/IA.F $dbdir/tmp/IA.P $dbdir/tmp/IA.C |sort > $dbdir/curated/ontology/IA.txt
