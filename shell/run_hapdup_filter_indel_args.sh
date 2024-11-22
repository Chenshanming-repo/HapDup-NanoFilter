#!/bin/bash
CTG=$1
OUT_DIR="/public/home/hpc224712204/bench/test_hapdup/$CTG/${CTG}_filter_indel"
/public/home/hpc224712204/software/hapdup/hapdup.py --use-filter-indel --assembly ~/data/hg002_tag/hg002_ont/$CTG/flye_hap/assembly.fasta --bam ~/data/hg002_tag/hg002_ont/$CTG/all_${CTG}_assembly.bam --out-dir $OUT_DIR -t 40 --rtype ont

/public/home/hpc224712204/bench/test_hapdup/run_merqury.sh $OUT_DIR
