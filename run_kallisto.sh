#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out kallisto.log

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi

module load kallisto
SAMPLES=samples.csv
CDSFILE=M_tuberculosis.cds.fasta
FASTQDIR=fastq
INDEX=Mtub.idx
OUTDIR=kallisto_results
mkdir -p $OUTDIR
if [ ! -f $INDEX ]; then
    kallisto index -i $INDEX $CDSFILE
fi

IFS=,
tail -n +2 $SAMPLES | while read ACC CONDITION PH REPLICATE
do
    if [ ! -s $OUTDIR/${CONDITION}_${PH}_r${REPLICATE}/abundance.tsv ]; then
	kallisto quant -i $INDEX -t $CPUS -o $OUTDIR/${CONDITION}_${PH}_r${REPLICATE} -l 300 --sd 20 --single $FASTQDIR/${ACC}.fastq.gz
    fi
done
