#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
cd       $dbdir/result
for dataset in `echo validate`;do
    $bindir/ghostx aln -d $dbdir/curated/sequence/uniprot_sprot_exp.ghostx -i $dbdir/curated/sequence/$dataset.fasta -o ghostx_$dataset.ghostx
    $bindir/ghostx2tsv.py $dbdir/curated/sequence/$dataset.fasta $dbdir/curated/sequence/uniprot_sprot_exp.fasta ghostx_$dataset.ghostx ghostx_$dataset.m6
    rm ghostx_$dataset.ghostx

    for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
        $bindir/parse_blast.py ghostx_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a ghostx_$dataset.F.$m.tsv $m
        $bindir/parse_blast.py ghostx_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a ghostx_$dataset.P.$m.tsv $m
        $bindir/parse_blast.py ghostx_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a ghostx_$dataset.C.$m.tsv $m
    
        cat ghostx_$dataset.F.$m.tsv ghostx_$dataset.P.$m.tsv ghostx_$dataset.C.$m.tsv > ghostx_$dataset.$m.tsv
        rm  ghostx_$dataset.F.$m.tsv ghostx_$dataset.P.$m.tsv ghostx_$dataset.C.$m.tsv 

        $bindir/assess_result.py $dbdir/curated/ontology/IA.txt ghostx_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a ghostx_$dataset.$m.txt
    done
done
