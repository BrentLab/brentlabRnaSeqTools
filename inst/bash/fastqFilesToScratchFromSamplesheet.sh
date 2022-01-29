#!/bin/bash

#./fastqFilesToScratchFromSamplesheet.sh \
#    /path/to/nf_sample_sheet.tsv \ # (this is $1)
#    /lts/mblab/Crypto/rnaseq_data/lts_sequence # (this is $2)

sed 's/,/\t/g' $1 > tmp_sample_sheet1234.tsv

nrow=$(wc -l tmp_sample_sheet1234.tsv | grep -oP "^[[:digit:]]+")


for i in $(seq 1 $nrow); do

    # NOTE!! ONE OF THE FOLLOWING
    # for the novoalign samplesheet
    read -r fastq1 strandedness < <(sed -n ${i}p tmp_sample_sheet1234.tsv)
    # for the nf-co pipeline
    #read -r sample fastq1 fastq2 strandedness < <(sed -n ${i}p $1)

    dest_dir=$(dirname $fastq1)
    fastq1_basename=$(basename $fastq1)
    fastq1_dirname=$(basename $dest_dir)

    source=${2}/${fastq1_dirname}/${fastq1_basename}

    mkdir -p $dest_dir

    echo "copying $i of $nrow"

    rsync -aHv $source $dest_dir

done

rm tmp_sample_sheet1234.tsv
