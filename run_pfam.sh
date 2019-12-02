#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out run_pfam.log

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi

module load hmmer
module load db-pfam

IN=M_tuberculosis.pep.fasta
DB=$PFAM_DB/Pfam-A.hmm
BASE=$(basename $IN .pep.fasta)
DOMTBL=$BASE.domtbl.tab
OUT=$BASE.hmmscan
if [ ! -f $DOMTBL ]; then
    hmmscan --cpu $CPUS --domtbl $DOMTBL --tblout $BASE.tbl.tab --cut_ga -o $OUT $DB $IN 
fi

awk '{print $1,$2,$4}' $DOMTBL | grep -v '^#' > $BASE.domtbl_simple.tab
