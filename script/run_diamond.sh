#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
cd       $dbdir/result
for dataset in `echo validate`;do
    $bindir/diamond blastp --db $dbdir/curated/sequence/uniprot_sprot_exp.fasta --outfmt 6 qseqid qlen sseqid slen evalue bitscore length nident --query $dbdir/curated/sequence/$dataset.fasta --out diamond_$dataset.m6

    for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
        $bindir/parse_blast.py diamond_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a diamond_$dataset.F.$m.tsv $m
        $bindir/parse_blast.py diamond_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a diamond_$dataset.P.$m.tsv $m
        $bindir/parse_blast.py diamond_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a diamond_$dataset.C.$m.tsv $m
    
        cat diamond_$dataset.F.$m.tsv diamond_$dataset.P.$m.tsv diamond_$dataset.C.$m.tsv > diamond_$dataset.$m.tsv
        rm  diamond_$dataset.F.$m.tsv diamond_$dataset.P.$m.tsv diamond_$dataset.C.$m.tsv 

        $bindir/assess_result.py $dbdir/curated/ontology/IA.txt diamond_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a diamond_$dataset.$m.txt
    done
done

for sensitive in `echo faster fast mid-sensitive sensitive more-sensitive very-sensitive ultra-sensitive`;do
    for dataset in `echo validate`;do
        $bindir/diamond blastp --$sensitive --db $dbdir/curated/sequence/uniprot_sprot_exp.fasta --outfmt 6 qseqid qlen sseqid slen evalue bitscore length nident --query $dbdir/curated/sequence/$dataset.fasta --out ${sensitive}_$dataset.m6

        for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
            $bindir/parse_blast.py ${sensitive}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a ${sensitive}_$dataset.F.$m.tsv $m
            $bindir/parse_blast.py ${sensitive}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a ${sensitive}_$dataset.P.$m.tsv $m
            $bindir/parse_blast.py ${sensitive}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a ${sensitive}_$dataset.C.$m.tsv $m
    
            cat ${sensitive}_$dataset.F.$m.tsv ${sensitive}_$dataset.P.$m.tsv ${sensitive}_$dataset.C.$m.tsv > ${sensitive}_$dataset.$m.tsv
            rm  ${sensitive}_$dataset.F.$m.tsv ${sensitive}_$dataset.P.$m.tsv ${sensitive}_$dataset.C.$m.tsv 

            $bindir/assess_result.py $dbdir/curated/ontology/IA.txt ${sensitive}_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a ${sensitive}_$dataset.$m.txt
        done
    done
done
