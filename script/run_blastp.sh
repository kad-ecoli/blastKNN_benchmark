#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
cd       $dbdir/result
for dataset in `echo validate`;do
    $bindir/blastp -db $dbdir/curated/sequence/uniprot_sprot_exp.fasta -outfmt '6 qacc qlen sacc slen evalue bitscore length nident' -query $dbdir/curated/sequence/$dataset.fasta -out blastp_$dataset.m6

    for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
    	$bindir/parse_blast.py blastp_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a blastp_$dataset.F.$m.tsv $m
    	$bindir/parse_blast.py blastp_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a blastp_$dataset.P.$m.tsv $m
        $bindir/parse_blast.py blastp_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a blastp_$dataset.C.$m.tsv $m
    
        cat blastp_$dataset.F.$m.tsv blastp_$dataset.P.$m.tsv blastp_$dataset.C.$m.tsv > blastp_$dataset.$m.tsv
        rm  blastp_$dataset.F.$m.tsv blastp_$dataset.P.$m.tsv blastp_$dataset.C.$m.tsv 

    	$bindir/assess_result.py $dbdir/curated/ontology/IA.txt blastp_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a blastp_$dataset.$m.txt
    done
done
