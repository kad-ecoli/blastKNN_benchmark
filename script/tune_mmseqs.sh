#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

rm -rf   $dbdir/result/tmp
mkdir -p $dbdir/result/tmp
cd       $dbdir/result
for dataset in `echo validate`;do
    for s in `echo 1.0 4.0 5.7 7.5 10`;do
        for max_seqs in `echo 10 25 50 100 200 300 500 1000`;do
            if [ -s "mmseqs_${s}_${max_seqs}_$dataset.txt" ];then
                continue
            fi
            mkdir -p $dbdir/result/tmp/$s
            $bindir/mmseqs/bin/mmseqs easy-search $dbdir/curated/sequence/$dataset.fasta $dbdir/curated/sequence/uniprot_sprot_exp mmseqs_${s}_${max_seqs}_$dataset.m6 $dbdir/result/tmp/$s --format-output query,qlen,target,tlen,evalue,bits,alnlen,fident -s $s --max-seqs $max_seqs
            rm -rf $dbdir/result/tmp/$s

            $bindir/parse_blast.py mmseqs_${s}_${max_seqs}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a mmseqs_${s}_${max_seqs}_$dataset.F.tsv 2
            $bindir/parse_blast.py mmseqs_${s}_${max_seqs}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a mmseqs_${s}_${max_seqs}_$dataset.P.tsv 2
            $bindir/parse_blast.py mmseqs_${s}_${max_seqs}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a mmseqs_${s}_${max_seqs}_$dataset.C.tsv 2

            cat mmseqs_${s}_${max_seqs}_$dataset.F.tsv mmseqs_${s}_${max_seqs}_$dataset.P.tsv mmseqs_${s}_${max_seqs}_$dataset.C.tsv > mmseqs_${s}_${max_seqs}_$dataset.tsv
            rm  mmseqs_${s}_${max_seqs}_$dataset.F.tsv mmseqs_${s}_${max_seqs}_$dataset.P.tsv mmseqs_${s}_${max_seqs}_$dataset.C.tsv

            $bindir/assess_result.py $dbdir/curated/ontology/IA.txt mmseqs_${s}_${max_seqs}_$dataset.tsv $dbdir/curated/sequence/$dataset.3.is_a mmseqs_${s}_${max_seqs}_$dataset.txt
        done
    done
done


for dataset in `echo validate`;do
    for s in `echo 7.5`;do
        for max_seqs in `echo 25 50 100`;do
            ln -s mmseqs_${s}_${max_seqs}_$dataset.txt mmseqs_${s}_${max_seqs}_1_$dataset.txt
            for num_iterations in `echo 2 3 4`;do
                #if [ -s "mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.txt" ];then
                if [ -s "mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.tsv" ];then
                    continue
                fi
                mkdir -p $dbdir/result/tmp/$s
                $bindir/mmseqs/bin/mmseqs easy-search $dbdir/curated/sequence/$dataset.fasta $dbdir/curated/sequence/uniprot_sprot_exp mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.m6 $dbdir/result/tmp/$s --format-output query,qlen,target,tlen,evalue,bits,alnlen,fident -s $s --max-seqs $max_seqs --num-iterations $num_iterations
                rm -rf $dbdir/result/tmp/$s

                $bindir/parse_blast.py mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.F.tsv 2
                $bindir/parse_blast.py mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.P.tsv 2
                $bindir/parse_blast.py mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.C.tsv 2

                cat mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.F.tsv mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.P.tsv mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.C.tsv > mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.tsv
                rm  mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.F.tsv mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.P.tsv mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.C.tsv

                $bindir/assess_result.py $dbdir/curated/ontology/IA.txt mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.tsv $dbdir/curated/sequence/$dataset.3.is_a mmseqs_${s}_${max_seqs}_${num_iterations}_$dataset.txt
            done
        done
    done
done
