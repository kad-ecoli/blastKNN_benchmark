# blastKNN_benchmark #

A large-scale assessment of sequence database search tools for homology-based protein function prediction

#### download data ####
```bash
script/download_obo.sh
script/download_sequence.sh
script/download_validate.sh
```

The training data set are available at
`` curated/sequence/uniprot_sprot_exp* ``.
The validation data set are available at
`` curated/sequence/validate* ``.
Among these files, the one ending with ".is_a" are the ground truth labels.

#### run function prediction ####
```bash
script/run_blastp.sh   
script/run_diamond.sh
script/run_mmseqs.sh
script/run_hhblits.sh
script/run_jackhmmer.sh
script/run_psiblast.sh
script/run_phmmer.sh

script/tune_blastp.sh
script/tune_diamond.sh
script/tune_mmseqs.sh
```

#### plot function prediction result ####
```bash
script/plot_default.py
script/plot_diamond.py

script/tune_blastp.py
script/tune_diamond.py
script/tune_mmseqs.py
script/tune_mmseqs_iteration.py
script/tune.py

script/plot_best.py
```
