#!/bin/bash
CTG="chr1"
/public/home/hpc224712204/software/hapdup/hapdup.py --use-filter --assembly ~/data/hg002_tag/hg002_ont/$CTG/flye_hap/assembly.fasta --bam ~/data/hg002_tag/hg002_ont/$CTG/all_${CTG}_assembly.bam --out-dir /public/home/hpc224712204/bench/test_hapdup/$CTG/${CTG}_filter -t 40 --rtype ont

/public/home/hpc224712204/bench/test_hapdup_chr22/run_merqury.sh /public/home/hpc224712204/bench/test_hapdup/$CTG/${CTG}_filter