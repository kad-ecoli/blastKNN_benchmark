#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
dbdir=`dirname $bindir`

mkdir -p $dbdir/result
mkdir -p $dbdir/tmp
mkdir -p $dbdir/slurm
cd       $dbdir/result
for dataset in `echo validate`;do
    #rm $dbdir/tmp/phmmer_$dataset*
    $bindir/fasta2tsv $dbdir/curated/sequence/$dataset.fasta > $dbdir/tmp/$dataset.tsv
    split -l 100  $dbdir/tmp/$dataset.tsv $dbdir/tmp/phmmer_${dataset}
    
    for filename in `ls $dbdir/tmp/|grep phmmer_|grep -vF '.'`;do
        if [ -s "$dbdir/tmp/$filename.m6" ];then
            continue
        fi
        echo "#!/bin/bash" > $dbdir/slurm/$filename.sh
        echo "#SBATCH --mem=5gb" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -t 168:00:00" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -o $dbdir/slurm/$filename.out" >> $dbdir/slurm/$filename.sh
        echo "#SBATCH -p sigbio" >> $dbdir/slurm/$filename.sh
        echo "$bindir/fasta2tsv $dbdir/curated/sequence/uniprot_sprot_exp.fasta > $dbdir/tmp/$filename.db" >> $dbdir/slurm/$filename.sh
        echo "cat $dbdir/tmp/$filename >> $dbdir/tmp/$filename.db" >> $dbdir/slurm/$filename.sh
        echo "cat $dbdir/tmp/$filename.db | sort | uniq | grep -v '^\\s*$' | $bindir/tsv2fasta - $dbdir/tmp/$filename.db.fasta" >> $dbdir/slurm/$filename.sh
        echo "$bindir/tsv2fasta $dbdir/tmp/$filename $dbdir/tmp/$filename.fasta" >> $dbdir/slurm/$filename.sh
        echo "$bindir/phmmer --incE 10 -A $dbdir/tmp/$filename.sto --tblout $dbdir/tmp/$filename.tblout --noali $dbdir/tmp/$filename.fasta $dbdir/tmp/$filename.db.fasta > /dev/null"  >> $dbdir/slurm/$filename.sh
        echo "$bindir/esl-reformat pfam $dbdir/tmp/$filename.sto  | grep -vP '^#=G' > $dbdir/tmp/$filename.pfam"  >> $dbdir/slurm/$filename.sh
        echo "$bindir/jackhmmer2tsv.py $dbdir/tmp/$filename.fasta $dbdir/tmp/$filename.db.fasta $dbdir/tmp/$filename.pfam $dbdir/tmp/$filename.tblout $dbdir/tmp/$filename.m6"  >> $dbdir/slurm/$filename.sh
        chmod +x $dbdir/slurm/$filename.sh
        #$dbdir/slurm/$filename.sh
        sbatch   $dbdir/slurm/$filename.sh
    done
    
    for filename in `ls $dbdir/tmp/|grep phmmer_|grep -vF '.'`;do
        cat $dbdir/tmp/$filename.m6
    done > phmmer_$dataset.m6

    for m in `echo 1 2 3 4 5 6 7 13 14 15 16`;do
        $bindir/parse_blast.py phmmer_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.F.is_a phmmer_$dataset.F.$m.tsv $m
        $bindir/parse_blast.py phmmer_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.P.is_a phmmer_$dataset.P.$m.tsv $m
        $bindir/parse_blast.py phmmer_$dataset.m6 $dbdir/curated/sequence/uniprot_sprot_exp.C.is_a phmmer_$dataset.C.$m.tsv $m
    
        cat phmmer_$dataset.F.$m.tsv phmmer_$dataset.P.$m.tsv phmmer_$dataset.C.$m.tsv > phmmer_$dataset.$m.tsv
        rm  phmmer_$dataset.F.$m.tsv phmmer_$dataset.P.$m.tsv phmmer_$dataset.C.$m.tsv 

        $bindir/assess_result.py $dbdir/curated/ontology/IA.txt phmmer_$dataset.$m.tsv $dbdir/curated/sequence/$dataset.3.is_a phmmer_$dataset.$m.txt
    done
done
