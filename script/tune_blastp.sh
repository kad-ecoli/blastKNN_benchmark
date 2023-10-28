#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
cd       $dbdir/result
mkdir -p $dbdir/slurm

for dataset in `echo validate`;do
    for evalue in `echo 0.0001 0.001 0.01 0.1 1 10 100`;do
        for max_target_seqs in `echo 10 25 50 100 200 300 500 1000`;do
            prefix=blastp_${evalue}_${max_target_seqs}_$dataset
            if [ -s "$prefix.txt" ];then
                continue
            fi
            
            echo $prefix

            echo "#!/bin/bash" > $dbdir/slurm/$prefix.sh
            echo "#SBATCH --mem=5gb" >> $dbdir/slurm/$prefix.sh
            echo "#SBATCH -t 48:00:00" >> $dbdir/slurm/$prefix.sh
            echo "#SBATCH -o $dbdir/slurm/$prefix.out" >> $dbdir/slurm/$prefix.sh
            echo "#SBATCH -p sigbio" >> $dbdir/slurm/$prefix.sh
            echo "cd $dbdir/result" >> $dbdir/slurm/$prefix.sh
            echo "$bindir/blastp -db $dbdir/curated/sequence/uniprot_sprot_exp.fasta -outfmt '6 qacc qlen sacc slen evalue bitscore length nident' -query $dbdir/curated/sequence/$dataset.fasta -evalue $evalue -max_target_seqs $max_target_seqs -out $dbdir/result/$prefix.m6" >> $dbdir/slurm/$prefix.sh
            echo "$bindir/parse_blast.py $dbdir/result/$prefix.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a $dbdir/result/$prefix.F.tsv 2"  >> $dbdir/slurm/$prefix.sh
            echo "$bindir/parse_blast.py $dbdir/result/$prefix.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a $dbdir/result/$prefix.P.tsv 2"  >> $dbdir/slurm/$prefix.sh
            echo "$bindir/parse_blast.py $dbdir/result/$prefix.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a $dbdir/result/$prefix.C.tsv 2"  >> $dbdir/slurm/$prefix.sh
            echo "cat $dbdir/result/$prefix.F.tsv $dbdir/result/$prefix.P.tsv $dbdir/result/$prefix.C.tsv > $prefix.tsv" >> $dbdir/slurm/$prefix.sh
            echo "rm  $dbdir/result/$prefix.F.tsv $dbdir/result/$prefix.P.tsv $dbdir/result/$prefix.C.tsv" >> $dbdir/slurm/$prefix.sh
            echo "$bindir/assess_result.py $dbdir/curated/ontology/IA.txt $prefix.tsv $dbdir/curated/sequence/$dataset.3.is_a $prefix.txt"  >> $dbdir/slurm/$prefix.sh

            chmod +x $dbdir/slurm/$prefix.sh
            #$dbdir/slurm/$prefix.sh
            sbatch   $dbdir/slurm/$prefix.sh
        done
    done
done
