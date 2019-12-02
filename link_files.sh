#!/usr/bin/bash
mkdir -p fastq
ln -s /bigdata/gen220/shared/data/M_tuberculosis/M_tuberculosis.cds.fasta .
ln -s /bigdata/gen220/shared/data/M_tuberculosis/M_tuberculosis.pep.fasta .
ln -s /bigdata/gen220/shared/data/M_tuberculosis/*.fastq.gz fastq/

curl -O http://current.geneontology.org/ontology/external2go/pfam2go
