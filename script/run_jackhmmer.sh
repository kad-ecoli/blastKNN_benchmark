#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
mkdir -p $dbdir/tmp
mkdir -p $dbdir/slurm
cd       $dbdir/result
for dataset in `echo validate`;do
    #rm $dbdir/tmp/jackhmmer_$dataset*
    $bindir/fasta2tsv $dbdir/curated/sequence/$dataset.fasta > $dbdir/tmp/$dataset.tsv
    split -l 100  $dbdir/tmp/$dataset.tsv $dbdir/tmp/jackhmmer_${dataset}
    
    for filename in `ls $dbdir/tmp/|grep jackhmmer_|grep -vF '.'`;do
        if [ -s "$dbdir/tmp/$filename.m6" ];then
            continue
        fi
        echo "#!/bin/bash" > $dbdir/slurm/$filename.sh
        echo "#SBATCH --mem=10gb" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -t 168:00:00" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -o $dbdir/slurm/$filename.out" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -p sigbio" >> $dbdir/slurm/$filename.sh
        echo "$bindir/tsv2fasta $dbdir/tmp/$filename $dbdir/tmp/$filename.fasta" >> $dbdir/slurm/$filename.sh
        echo "$bindir/jackhmmer -A $dbdir/tmp/$filename.sto --tblout $dbdir/tmp/$filename.tblout --noali $dbdir/tmp/$filename.fasta $dbdir/curated/sequence/uniprot_sprot_exp.fasta > /dev/null"  >> $dbdir/slurm/$filename.sh
        echo "$bindir/esl-reformat pfam $dbdir/tmp/$filename.sto  | grep -vP '^#=G' > $dbdir/tmp/$filename.pfam"  >> $dbdir/slurm/$filename.sh
        echo "$bindir/jackhmmer2tsv.py $dbdir/tmp/$filename.fasta $dbdir/curated/sequence/uniprot_sprot_exp.fasta $dbdir/tmp/$filename.pfam $dbdir/tmp/$filename.tblout $dbdir/tmp/$filename.m6"  >> $dbdir/slurm/$filename.sh
        chmod +x $dbdir/slurm/$filename.sh
        #$dbdir/slurm/$filename.sh
        sbatch   $dbdir/slurm/$filename.sh
    done
    
    for filename in `ls $dbdir/tmp/|grep jackhmmer_|grep -vF '.'`;do
        cat $dbdir/tmp/$filename.m6
    done > jackhmmer_$dataset.m6

    for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
        $bindir/parse_blast.py jackhmmer_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a jackhmmer_$dataset.F.$m.tsv $m
        $bindir/parse_blast.py jackhmmer_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a jackhmmer_$dataset.P.$m.tsv $m
        $bindir/parse_blast.py jackhmmer_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a jackhmmer_$dataset.C.$m.tsv $m
    
        cat jackhmmer_$dataset.F.$m.tsv jackhmmer_$dataset.P.$m.tsv jackhmmer_$dataset.C.$m.tsv > jackhmmer_$dataset.$m.tsv
        rm  jackhmmer_$dataset.F.$m.tsv jackhmmer_$dataset.P.$m.tsv jackhmmer_$dataset.C.$m.tsv 

        $bindir/assess_result.py $dbdir/curated/ontology/IA.txt jackhmmer_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a jackhmmer_$dataset.$m.txt
    done
done
