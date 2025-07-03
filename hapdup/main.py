#!/usr/bin/env python3

import os
import sys
import argparse
from distutils import spawn
import subprocess
import threading
import shutil
import logging
from hapdup.hapcut2_multiprocess import run_hapcut2_multiprocess, run_hapcut2_filter_multiprocess
from hapdup.find_breakpoints import find_breakpoints
from hapdup.bed_liftover import liftover_parallel
from hapdup.apply_inversions import apply_inversions
from hapdup.filter_misplaced_alignments import filter_alignments_parallel
from hapdup.cut_phased_blocks import cut_phased_blocks
from hapdup.MultiStepConsistencyAnalysisHapDup import filter_main
from hapdup.__version__ import __version__

FLYE = "flye"
SAMTOOLS = "flye-samtools"
MINIMAP = "flye-minimap2"
PEPPER_VARIANT = "pepper_variant"
EXTRACT_HAIRS="extractHAIRS"
HAPCUT2="HAPCUT2"
MARGIN = "margin"

PEPPER_MODEL_DIR = os.environ["PEPPER_MODEL_DIR"]
PEPPER_MODEL = {"hifi" : os.path.join(PEPPER_MODEL_DIR, "PEPPER_VARIANT_HIFI_V7.pkl"),
                 "ont" : os.path.join(PEPPER_MODEL_DIR, "PEPPER_VARIANT_ONT_R941_GUPPY5_SUP_V7.pkl"),
                 "ont_r10" : os.path.join(PEPPER_MODEL_DIR, "PEPPER_VARIANT_ONT_R10_Q20_V7.pkl")}
MARGIN_CONFIG_DIR = os.environ["MARGIN_CONFIG_DIR"]
#MARGIN_CONFIG = {"hifi" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.pb-hifi.hapDup.json"),
#
#MARGIN_CONFIG_DIR = "/public/home/hpc224712204/software/hapdup/submodules/margin/params/phase"
MARGIN_CONFIG = {"hifi" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.pb-hifi.hapDup.json"),
                 "ont" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.ont-r94g507.hapDup.filter.json"),
                 "ont_r10" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.ont-r94g507.hapDup.filter.json")}
MARGIN_ONLY_SNV_CONFIG = {"hifi" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.pb-hifi.hapDup.json"),
                 "ont" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.ont-r94g507.hapDup.json"),
                 "ont_r10" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.ont-r94g507.hapDup.json")}
MARGIN_FILTER_CONFIG = {"hifi" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.pb-hifi.hapDup.json"),
                 "ont" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.ont-r94g507.hapDup.filter.json"),
                 "ont_r10" : os.path.join(MARGIN_CONFIG_DIR, "allParams.haplotag.ont-r94g507.hapDup.filter.json")}

logger = logging.getLogger()


def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def _version():
    return str(__version__)

def bool_to_int(bool_var):
    return 1 if bool_var else 0

def gzip_vcf(vcf_fn):
    raw_vcf = os.path.abspath(vcf_fn)
    gz_vcf = os.path.abspath(vcf_fn + ".gz")
    subprocess.check_call(f"bgzip -c {raw_vcf} > {gz_vcf}", shell=True)
    subprocess.check_call(f"tabix -p vcf -f {gz_vcf}", shell=True)

def main():
    parser = argparse.ArgumentParser(description="Reassemble haplotypes from collapsed haploid assmebly")

    parser.add_argument("--assembly", dest="assembly",
                        metavar="path", required=True,
                        help="path to haploid assembly (contigs in fasta format)")
    parser.add_argument("--bam", dest="bam", required=True, metavar="path",
                        default=None, help="path to the alignment of reads on the assembly in bam format")
    parser.add_argument("--bam-index", dest="bam_index", required=False, metavar="path",
                        default=None, help="path to bam index (if non-standard)")
    parser.add_argument("--pepper-model", dest="pepper_model", required=False, metavar="path",
                        default=None, help="path to custom pepper model")
    parser.add_argument("--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("--overwrite", dest="overwrite",
                        default=False, action="store_true",
                        help="Do not attempt to restart from complete phases, overwrites existing results")
    parser.add_argument("--use-unphased", dest="use_unphased",
                        default=False, action="store_true",
                        help="Use unphased reads for polishing")
    parser.add_argument("--rtype", dest="rtype", default="ont", required=True, metavar="(ont|hifi|ont_r10)",
                        help="Long reads type, ONT or HiFi: (ont|hifi|ont_r10)")

    parser.add_argument("--min-aligned-length", dest="min_aligned_length", type=int,
                        default=10000, metavar="int", help="minimum aligned length for each read [10000]")
    parser.add_argument("--max-read-error", dest="max_read_error", type=float,
                        default=0.1, metavar="float", help="maximum read mapping error rate [0.1]")
    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=10, metavar="int", help="number of parallel threads [10]")
    parser.add_argument("-v", "--version", action="version", version=_version())

    # additional
    parser.add_argument("--use-filter-all", dest="use_filter_all",
                        default=False, action="store_true",
                        help="Use consistency to filter SNP and INDEL")
    parser.add_argument("--use-filter-indel", dest="use_filter_indel",
                        default=False, action="store_true",
                        help="Use consistency to filter INDEL, and remain SNP")
    parser.add_argument("--phasing", dest="phasing",
                        default=False, action="store_true",
                        help="Only phasing, not polish or assembly")
    parser.add_argument("--only-snvs", dest="only_snvs",
                        default=False, action="store_true",
                        help="Only use snv to phasing")
    parser.add_argument("--use-multiprocess", dest="use_multiprocess",
                        default=False, action="store_true",
                        help="Use multiprocess to run HapCUT2")
    # bench
    parser.add_argument("--bench-snvs", dest="bench_snvs",
                        default=False, action="store_true",
                        help="Use consistency to filter SNP and INDEL, and only SNP to phasing")
    parser.add_argument("--bench-qual", dest="bench_qual",
                        default=False, action="store_true",
                        help="Filter vcf by qual 10.")


    args = parser.parse_args()

    for e in [SAMTOOLS, FLYE, MINIMAP, PEPPER_VARIANT]:
        if not spawn.find_executable(e):
            print("Not installed: " + e, file=sys.stderr)
            return 1

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    # if args.use_filter_indel:
    #    args.use_filter = True

    if args.rtype not in ["ont", "hifi", "ont_r10"]:
        print("Reads type could be not in ['ont' 'hifi' 'ont_r10']", file=sys.stderr)
        return 1

    def file_check(path):
        if not os.path.isfile(path):
            logger.error("Missing output: %s", path)
            raise Exception("Missing output")

    hapdup_log = os.path.join(args.out_dir, "hapdup.log")
    _enable_logging(hapdup_log, debug=False, overwrite=False)

    file_check(args.bam)
    file_check(args.assembly)

    file_check(PEPPER_MODEL[args.rtype])

    tmp_dir = os.path.join(args.out_dir, "tmp")
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    filtered_bam = os.path.join(args.out_dir, "filtered.bam")
    pepper_dir = os.path.join(args.out_dir, "pepper")
    #pepper_vcf = os.path.abspath(os.path.join(pepper_dir, "PEPPER_VARIANT_SNP_OUTPUT.vcf"))
    #pepper_vcf = os.path.abspath(os.path.join(pepper_dir, "PEPPER_VARIANT_OUTPUT_PHASING.vcf"))
    pepper_vcf = os.path.abspath(os.path.join(pepper_dir, "PEPPER_VARIANT_FULL.vcf"))

    hapcut2_dir = os.path.join(args.out_dir, "hapcut2")
    margin_dir = os.path.join(args.out_dir, "margin")
    haplotagged_bam = os.path.join(margin_dir, "MARGIN_PHASED") + ".haplotagged.bam"
    phased_blocks_bed = os.path.join(margin_dir, "MARGIN_PHASED") + ".phaseset.bed"

    polished_flye_hap = {1 : os.path.join(args.out_dir, "flye_hap_1", "polished_1.fasta"),
                         2 : os.path.join(args.out_dir, "flye_hap_2", "polished_1.fasta")}

    structural_dir = os.path.join(args.out_dir, "structural")
    inversions_bed = os.path.join(structural_dir, "inversions.bed")

    final_dual_asm = {1 : os.path.join(args.out_dir, "hapdup_dual_1.fasta"),
                      2 : os.path.join(args.out_dir, "hapdup_dual_2.fasta")}

    final_phased_asm = {1 : os.path.join(args.out_dir, "hapdup_phased_1.fasta"),
                        2 : os.path.join(args.out_dir, "hapdup_phased_2.fasta")}

    #STAGE 1: filter suspicious alignments, index the resulting bam
    overwrite = args.overwrite
    if (os.path.isfile(filtered_bam) or os.path.isfile(haplotagged_bam)) and not overwrite:
        logger.info("Skipped filtering phase")
    else:
        logger.info("Filtering alignments")
        filter_alignments_parallel(args.bam, filtered_bam, min(args.threads, 30),
                                   args.min_aligned_length, args.max_read_error,
                                   args.bam_index)
        file_check(filtered_bam)

        # index_cmd = [SAMTOOLS, "index", "-@4", filtered_bam]
        index_cmd = [SAMTOOLS, "index", "-@", f"{args.threads}",filtered_bam]
        logger.info("Running: %s", " ".join(index_cmd))
        subprocess.check_call(" ".join(index_cmd), shell=True)
        overwrite = True

    #STAGE 2: Run PEPPER

    if os.path.isfile(pepper_vcf) and not overwrite:
        logger.info("Skipped pepper phase")
    else:
        logger.info(f"PEPPER Model Path: {PEPPER_MODEL_DIR}, using {PEPPER_MODEL[args.rtype]}")
        pepper_log = os.path.join(pepper_dir, "pepper.log")
        if not os.path.isdir(pepper_dir):
            os.mkdir(pepper_dir)

        if args.pepper_model is None:
            args.pepper_model = PEPPER_MODEL[args.rtype]
        model_copy = os.path.join(pepper_dir, "pepper_model.bin")
        shutil.copyfile(args.pepper_model, model_copy)

        reads_arg = None
        if args.rtype == "ont":
            reads_arg = "--ont_r9_guppy5_sup"
        elif args.rtype == "hifi":
            reads_arg = "--hifi"
        elif args.rtype == "ont_r10":
            reads_arg = "--ont_r10_q20"

        pepper_cmd = [PEPPER_VARIANT, "call_variant", "-b", os.path.abspath(filtered_bam), "-f", os.path.abspath(args.assembly),
                      "-o", os.path.abspath(pepper_dir), "-m", model_copy, "-t", str(args.threads), "-s", "Sample", reads_arg,
                      "--include-supplementary", "--no_quantized", "2>&1", "|tee", pepper_log]

        logger.info("Running: %s", " ".join(pepper_cmd))
        subprocess.check_call(" ".join(pepper_cmd), shell=True)
        subprocess.call("rm -r " + os.path.join(pepper_dir, "images*"), shell=True)
        subprocess.call("rm -r " + os.path.join(pepper_dir, "predictions*"), shell=True)
        file_check(pepper_vcf)
        overwrite = True
    # STAGE 2.1: delete indel lines in pepper_vcf if only_snvs is True
    if args.only_snvs:
        pepper_no_indel_vcf = f"{pepper_vcf}.no.indel"

        from hapdup.filter_indel_in_vcf import output_vcf_file_by_positions_only_snp
        output_vcf_file_by_positions_only_snp(pepper_vcf, pepper_no_indel_vcf)

        logger.info(f"Changing {pepper_vcf} to {pepper_no_indel_vcf}")
        pepper_vcf = pepper_no_indel_vcf
        # gzip_vcf(pepper_vcf)

    logger.info(f"gzip {pepper_vcf}")
    gzip_vcf(pepper_vcf)

    #STAGE 3: Phase with HapCUT2
    if os.path.isfile(haplotagged_bam) and not overwrite:
        logger.info("Skipped margin phase")
    else:
        hapcut2_log = os.path.join(margin_dir, "hapcut2.log")
        margin_log = os.path.join(margin_dir, "margin.log")
        if not os.path.isdir(margin_dir):
            os.mkdir(margin_dir)
        if not os.path.isdir(hapcut2_dir):
            os.mkdir(hapcut2_dir)

        if not args.use_filter_all and not args.use_filter_indel:
            margin_input_vcf = pepper_vcf
        else:
            logger.info("Use consistency filter")
            margin_input_vcf = f"{pepper_vcf}.filtered"

            # subprocess.check_call(f"rm -rf {tmp_dir}/*", shell=True)

            # Filter variant by consistency
            logger.info("Filtering (mp mode)")
            run_hapcut2_filter_multiprocess(os.path.abspath(args.assembly), os.path.abspath(filtered_bam), pepper_vcf, margin_input_vcf, tmp_dir,
                                        args.only_snvs, args.rtype, int(args.threads), FILTER_SNP=args.use_filter_all)

            if args.bench_snvs:
                logger.info("Only phasing filtered snp (bench)")
                from hapdup.filter_indel_in_vcf import output_vcf_file_by_positions_only_snp
                output_vcf_file_by_positions_only_snp(margin_input_vcf, f"{pepper_vcf}.filtered.bench.snv")

                margin_input_vcf = f"{pepper_vcf}.filtered.bench.snv"


        # gzip_vcf(margin_input_vcf)

        logger.info(f"HapCUT2 input vcf: {margin_input_vcf}")

        #  MARGIN phasing
        if args.use_filter_all and args.only_snvs:
            logger.error("args.use_filter_all and args.only_snvs all True")
        PARAM_FILE=MARGIN_CONFIG[args.rtype]
        if args.use_filter_all:
            PARAM_FILE=MARGIN_FILTER_CONFIG[args.rtype]
        if args.only_snvs:
            PARAM_FILE=MARGIN_ONLY_SNV_CONFIG[args.rtype]

        margin_cmd = [MARGIN, "phase", os.path.abspath(filtered_bam), os.path.abspath(args.assembly), margin_input_vcf,
                      PARAM_FILE,
                      "-t", str(args.threads), "-o", os.path.abspath(os.path.join(margin_dir, "MARGIN_PHASED")),
                      "2>&1", "|tee", margin_log]
        logger.info("Running: %s", " ".join(margin_cmd))
        subprocess.check_call(" ".join(margin_cmd), shell=True)
        file_check(haplotagged_bam)
        #subprocess.call("rm " + os.path.abspath(filtered_bam), shell=True)

        index_cmd = [SAMTOOLS, "index", "-@", str(args.threads), haplotagged_bam]
        logger.info("Running: %s", " ".join(index_cmd))
        subprocess.check_call(" ".join(index_cmd), shell=True)
        overwrite = True
    if args.phasing:
        sys.exit(0)
    #STAGE 4: polish haplotypes with Flye
    if all(map(os.path.isfile, polished_flye_hap.values())) and not overwrite:
        logger.info("Skipped Flye phase")
    else:
        def run_flye_hp(_hp):
            _threads = max(1, int(args.threads) // 2)
            flye_out = os.path.join(args.out_dir, "flye_hap_{0}".format(_hp))

            _reads_arg = None
            if args.rtype == "ont":
                _reads_arg = "--nano-raw"
            elif args.rtype == "hifi":
                _reads_arg = "--pacbio-hifi"
            elif args.rtype == "ont_r10":
                _reads_arg = "--nano-hq"

            hp_string = "0,{}".format(_hp) if args.use_unphased else str(_hp)

            flye_cmd = [FLYE, "--polish-target", os.path.abspath(args.assembly), _reads_arg,  haplotagged_bam, "-t", str(_threads),
                        "-o", flye_out, "--polish-haplotypes", hp_string, f"2>{args.out_dir}/flye.err.log"]
            logger.info("Running: %s", " ".join(flye_cmd))
            subprocess.check_call(" ".join(flye_cmd), shell=True)

        threads = []
        for hp in [1, 2]:
            threads.append(threading.Thread(target=run_flye_hp, args=(hp,)))
            threads[-1].start()
        for t in threads:
            t.join()
        overwrite = True

    #STAGE 5: structural polishing

    if not os.path.isdir(structural_dir):
        os.mkdir(structural_dir)
    logger.info("Finding breakpoints")
    find_breakpoints(haplotagged_bam, structural_dir, args.threads)

    for hp in [1, 2]:
        minimap_out = os.path.join(structural_dir, "liftover_hp{0}.bam".format(hp))
        minimap_cmd = [MINIMAP, "-ax", "asm5", "-t", str(args.threads), "-K", "5G", "-I", "64G",
                       args.assembly, polished_flye_hap[hp], "2>/dev/null", "|",
                       SAMTOOLS, "sort", "-m", "4G", "-@4", ">", minimap_out]
        logger.info("Running: %s", " ".join(minimap_cmd))
        subprocess.check_call(" ".join(minimap_cmd), shell=True)
        subprocess.check_call(SAMTOOLS + " index -@ 4 {0}".format(minimap_out), shell=True)

        inversions_hp = os.path.join(structural_dir, "inversions_hp{0}.bed".format(hp))
        #bed_liftover(inversions_bed, minimap_out, open(inversions_hp, "w"))
        liftover_parallel(inversions_bed, minimap_out, open(inversions_hp, "w"), False, args.threads)

        phased_blocks_hp = os.path.join(args.out_dir, "phased_blocks_hp{0}.bed".format(hp))
        #bed_liftover(phased_blocks_bed, minimap_out, open(phased_blocks_hp, "w"))
        liftover_parallel(phased_blocks_bed, minimap_out, open(phased_blocks_hp, "w"), False, args.threads)

        apply_inversions(inversions_hp, polished_flye_hap[hp], final_dual_asm[hp], hp)
        cut_phased_blocks(phased_blocks_hp, final_dual_asm[hp], final_phased_asm[hp])


if __name__ == "__main__":
    main()
