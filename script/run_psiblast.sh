#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
mkdir -p $dbdir/tmp
mkdir -p $dbdir/slurm
cd       $dbdir/result
for dataset in `echo validate`;do
    $bindir/fasta2tsv $dbdir/curated/sequence/$dataset.fasta > $dbdir/tmp/$dataset.tsv
    split -l 100  $dbdir/tmp/$dataset.tsv $dbdir/tmp/psiblast_${dataset}
    
    for filename in `ls $dbdir/tmp/|grep psiblast_|grep -vF '.'`;do
        if [ -s "$dbdir/tmp/$filename.m6" ];then
            continue
        fi
        echo "#!/bin/bash" > $dbdir/slurm/$filename.sh
        echo "#SBATCH --mem=5gb" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -t 168:00:00" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -o $dbdir/slurm/$filename.out" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -p sigbio" >> $dbdir/slurm/$filename.sh
        echo "$bindir/tsv2fasta $dbdir/tmp/$filename $dbdir/tmp/$filename.fasta" >> $dbdir/slurm/$filename.sh
        echo "$bindir/psiblast -db $dbdir/curated/sequence/uniprot_sprot_exp.fasta -outfmt '6 qacc qlen sacc slen evalue bitscore length nident' -query $dbdir/tmp/$filename.fasta -out $dbdir/tmp/$filename.m6 -num_iterations 3" >> $dbdir/slurm/$filename.sh
        chmod +x $dbdir/slurm/$filename.sh
        #$dbdir/slurm/$filename.sh
        sbatch   $dbdir/slurm/$filename.sh
    done
    
    for filename in `ls $dbdir/tmp/|grep psiblast_|grep -vF '.'`;do
        cat $dbdir/tmp/$filename.m6
    done | grep -P '\t' > psiblast_$dataset.m6

    for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
        $bindir/parse_blast.py psiblast_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a psiblast_$dataset.F.$m.tsv $m
        $bindir/parse_blast.py psiblast_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a psiblast_$dataset.P.$m.tsv $m
        $bindir/parse_blast.py psiblast_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a psiblast_$dataset.C.$m.tsv $m
    
        cat psiblast_$dataset.F.$m.tsv psiblast_$dataset.P.$m.tsv psiblast_$dataset.C.$m.tsv > psiblast_$dataset.$m.tsv
        rm  psiblast_$dataset.F.$m.tsv psiblast_$dataset.P.$m.tsv psiblast_$dataset.C.$m.tsv 

        $bindir/assess_result.py $dbdir/curated/ontology/IA.txt psiblast_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a psiblast_$dataset.$m.txt
    done
done
