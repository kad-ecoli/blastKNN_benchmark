#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

rm -rf   $dbdir/result/tmp
mkdir -p $dbdir/result/tmp
cd       $dbdir/result
for dataset in `echo validate`;do
    mkdir -p $dbdir/result/tmp/mmseqs
    $bindir/mmseqs/bin/mmseqs easy-search $dbdir/curated/sequence/$dataset.fasta $dbdir/curated/sequence/uniprot_sprot_exp mmseqs_$dataset.m6 $dbdir/result/tmp/mmseqs --format-output query,qlen,target,tlen,evalue,bits,alnlen,fident
    rm -rf $dbdir/result/tmp/mmseqs

    for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
        $bindir/parse_blast.py mmseqs_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a mmseqs_$dataset.F.$m.tsv $m
        $bindir/parse_blast.py mmseqs_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a mmseqs_$dataset.P.$m.tsv $m
        $bindir/parse_blast.py mmseqs_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a mmseqs_$dataset.C.$m.tsv $m
    
        cat mmseqs_$dataset.F.$m.tsv mmseqs_$dataset.P.$m.tsv mmseqs_$dataset.C.$m.tsv > mmseqs_$dataset.$m.tsv
        rm  mmseqs_$dataset.F.$m.tsv mmseqs_$dataset.P.$m.tsv mmseqs_$dataset.C.$m.tsv 

        $bindir/assess_result.py $dbdir/curated/ontology/IA.txt mmseqs_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a mmseqs_$dataset.$m.txt
    done
done
