#!/bin/bash
CTG=$1
cd /public/home/hpc224712204/software/hapdup
/public/home/hpc224712204/software/hapdup/hapdup.py --assembly ~/data/hg002_tag/hg002_ont/$CTG/flye_hap/assembly.fasta --bam ~/data/hg002_tag/hg002_ont/$CTG/all_${CTG}_assembly.bam --out-dir /public/home/hpc224712204/bench/test_hapdup/$CTG/${CTG}_raw -t 40 --rtype ont

/public/home/hpc224712204/bench/test_hapdup_chr22/run_merqury.sh /public/home/hpc224712204/bench/test_hapdup/$CTG/${CTG}_raw
