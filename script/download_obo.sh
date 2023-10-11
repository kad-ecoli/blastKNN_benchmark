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
