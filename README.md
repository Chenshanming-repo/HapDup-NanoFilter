[![DOI](https://zenodo.org/badge/892602566.svg)](https://doi.org/10.5281/zenodo.16777890)
# HapDup-NanoFilter

## Overview

This software is a modification version of HapDup which integrates with NanoFilter.

## Usage

```bash
# only filter indel                 
python hapdup.py --use-multiprocess \
                 --use-filter-indel \
                 --assembly test_data/ref_test.fasta \
                 --bam test_data/all_assembly_test.bam \
                 --out-dir <output_directory> \
                 -t <num_threads> \
                 --rtype <hifi/ont/ont_r10>
# filter snv and indel
python hapdup.py --use-multiprocess \
                 --use-filter-all \
                 --assembly test_data/ref_test.fasta \
                 --bam test_data/all_assembly_test.bam \
                 --out-dir <output_directory> \
                 -t <num_threads> \
                 --rtype <hifi/ont/ont_r10>
```

## Output

The output directory will contain:
* `hapdup_dual_{1,2}.fasta` - dual assembly
* `phased_blocks_hp{1,2}.bed` - phased blocks coordinates (in dual assmeblies)
* `hapdup_phased_{1,2}.fasta` - haplotype-resolved assmebly

## Installation 

```bash
git clone https://github.com/Chenshanming-repo/HapDup-NanoFilter.git
cd HapDup-NanoFilter

conda env create -f environment.yml
conda activate hapdup

git submodule update --init --recursive

#build and install Flye
pushd submodules/Flye/ && python setup.py install && popd

#build and install Margin
pushd submodules/margin/ && mkdir build && cd build && cmake .. && make && cp ./margin $CONDA_PREFIX/bin/ && popd

#build and install PEPPER and its dependencies
pushd submodules/pepper/ && python -m pip install . && popd

# install HapCUT2
git clone https://github.com/vibansal/HapCUT2.git
cd HapCUT2
make
# add HAPCUT2 and extractHAIRS to enviroment path
export PATH=$PATH:{PATH_TO_HAPCUT2}/build
```






## Contact

For issues or questions, please open an issue on the project repository or contact the developer (224712204@csu.edu.cn).
