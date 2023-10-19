#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
mkdir -p $dbdir/tmp/hh2
mkdir -p $dbdir/slurm
cd       $dbdir/result
for dataset in `echo validate`;do
    for target in `$bindir/fasta2tsv $dbdir/curated/sequence/$dataset.fasta|cut -f1`;do
        if [ -s "$dbdir/tmp/hh2/$target.m6" ];then
            continue
        fi
        echo $target
        sequence=`$bindir/fasta2tsv $dbdir/curated/sequence/$dataset.fasta |grep $target|cut -f2`
        echo ">$target" > $dbdir/tmp/hh2/$target.fasta
        echo $sequence >> $dbdir/tmp/hh2/$target.fasta
        echo "#!/bin/bash" > $dbdir/slurm/$target.sh
        echo "#SBATCH --mem=5gb" >> $dbdir/slurm/$target.sh
        echo "#SBATCH -t 24:00:00" >> $dbdir/slurm/$target.sh
        echo "#SBATCH -o $dbdir/slurm/$target.out" >> $dbdir/slurm/$target.sh
        echo "#SBATCH -p sigbio" >> $dbdir/slurm/$target.sh
        echo "$bindir/hhsuite2/bin/hhblits -i $dbdir/tmp/hh2/$target.fasta -d $dbdir/curated/sequence/uniprot_sprot_exp.hh2 -o $dbdir/tmp/hh2/$target.hhr" >> $dbdir/slurm/$target.sh
        echo "$bindir/hhblits2tsv.py $dbdir/tmp/hh2/$target.fasta $dbdir/curated/sequence/uniprot_sprot_exp.fasta $dbdir/tmp/hh2/$target.hhr $dbdir/tmp/hh2/$target.m6" >> $dbdir/slurm/$target.sh
        chmod +x $dbdir/slurm/$target.sh
        #$dbdir/slurm/$target.sh
        sbatch   $dbdir/slurm/$target.sh
    done
    
    for target in `$bindir/fasta2tsv $dbdir/curated/sequence/$dataset.fasta|cut -f1`;do
        cat $dbdir/tmp/hh2/$target.m6
    done > $dbdir/result/hhblits_$dataset.m6

    for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
        $bindir/parse_blast.py hhblits_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a hhblits_$dataset.F.$m.tsv $m
        $bindir/parse_blast.py hhblits_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a hhblits_$dataset.P.$m.tsv $m
        $bindir/parse_blast.py hhblits_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a hhblits_$dataset.C.$m.tsv $m
    
        cat hhblits_$dataset.F.$m.tsv hhblits_$dataset.P.$m.tsv hhblits_$dataset.C.$m.tsv > hhblits_$dataset.$m.tsv
        rm  hhblits_$dataset.F.$m.tsv hhblits_$dataset.P.$m.tsv hhblits_$dataset.C.$m.tsv

        $bindir/assess_result.py $dbdir/curated/ontology/IA.txt hhblits_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a hhblits_$dataset.$m.txt
    done
done
