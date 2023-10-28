#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
cd       $dbdir/result

for dataset in `echo validate`;do
    for evalue in `echo 0.0001 0.001 0.01 0.1 1 10 100`;do
        for max_target_seqs in `echo 10 25 50 100 200 300 500 1000 2000`;do
            if [ -s "diamond_${evalue}_${max_target_seqs}_$dataset.txt" ];then
                continue
            fi
            
            echo diamond_${evalue}_${max_target_seqs}_$dataset
            $bindir/diamond blastp --ultra-sensitive --db $dbdir/curated/sequence/uniprot_sprot_exp.fasta --outfmt 6 qseqid qlen sseqid slen evalue bitscore length nident --query $dbdir/curated/sequence/$dataset.fasta --out diamond_${evalue}_${max_target_seqs}_$dataset.m6 --max-target-seqs $max_target_seqs --evalue $evalue

            $bindir/parse_blast.py diamond_${evalue}_${max_target_seqs}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a diamond_${evalue}_${max_target_seqs}_$dataset.F.tsv 2
            $bindir/parse_blast.py diamond_${evalue}_${max_target_seqs}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a diamond_${evalue}_${max_target_seqs}_$dataset.P.tsv 2
            $bindir/parse_blast.py diamond_${evalue}_${max_target_seqs}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a diamond_${evalue}_${max_target_seqs}_$dataset.C.tsv 2

            cat diamond_${evalue}_${max_target_seqs}_$dataset.F.tsv diamond_${evalue}_${max_target_seqs}_$dataset.P.tsv diamond_${evalue}_${max_target_seqs}_$dataset.C.tsv > diamond_${evalue}_${max_target_seqs}_$dataset.tsv
            rm diamond_${evalue}_${max_target_seqs}_$dataset.F.tsv diamond_${evalue}_${max_target_seqs}_$dataset.P.tsv diamond_${evalue}_${max_target_seqs}_$dataset.C.tsv

            $bindir/assess_result.py $dbdir/curated/ontology/IA.txt diamond_${evalue}_${max_target_seqs}_$dataset.tsv $dbdir/curated/sequence/$dataset.3.is_a diamond_${evalue}_${max_target_seqs}_$dataset.txt
        done
    done
done
