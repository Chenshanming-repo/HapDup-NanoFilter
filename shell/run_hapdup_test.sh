#!/bin/bash

/public/home/hpc224712204/software/hapdup/hapdup.py --assembly ~/data/hg002_tag/hg002_ont/chr22/flye_hap/assembly.fasta --bam ~/data/hg002_tag/hg002_ont/chr22/all_chr22_assembly.bam --out-dir ~/bench/test -t 40 --rtype ont

/public/home/hpc224712204/bench/test_hapdup_chr22/run_merqury.sh /public/home/hpc224712204/bench/test
