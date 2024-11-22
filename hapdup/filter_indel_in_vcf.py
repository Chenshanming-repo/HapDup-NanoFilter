from hapdup.MultiStepConsistencyAnalysisHapDup import read_pepper_vcf
from hapdup.VariantRecord import VariantType

def output_vcf_file_filter_by_qual(vcf_fn, new_vcf_fn, qual=100000):
    clair3_vcf_records_list = read_pepper_vcf(vcf_fn)
    clair3_vcf_records_dict = { vcf_record.key : vcf_record for vcf_record in clair3_vcf_records_list }


    with open(new_vcf_fn, "w") as f_clair3:
        with open(vcf_fn, "r") as clair3:
            for line in clair3.readlines():
                if line.startswith("#"):
                    f_clair3.write(f"{line.strip()}\n")
                    continue
                items = line.split("\t")
                key = f"{items[0]}:{int(items[1])}"

                if int(clair3_vcf_records_dict[key].qual) < qual:
                    continue

                f_clair3.write(f"{line.strip()}\n")

    read_pepper_vcf(new_vcf_fn)

def output_vcf_file_by_positions_only_snp(vcf_fn, new_vcf_fn, qual=100000):
    clair3_vcf_records_list = read_pepper_vcf(vcf_fn)
    clair3_vcf_records_dict = { vcf_record.key : vcf_record for vcf_record in clair3_vcf_records_list }


    with open(new_vcf_fn, "w") as f_clair3:
        with open(vcf_fn, "r") as clair3:
            for line in clair3.readlines():
                if line.startswith("#"):
                    f_clair3.write(f"{line.strip()}\n")
                    continue
                items = line.split("\t")
                key = f"{items[0]}:{int(items[1])}"

                if clair3_vcf_records_dict[key].is_INDEL() and int(clair3_vcf_records_dict[key].qual) < qual:
                    continue

                f_clair3.write(f"{line.strip()}\n")

    read_pepper_vcf(new_vcf_fn)

def output_vcf_file_heter(vcf_fn, new_vcf_fn):
    clair3_vcf_records_list = read_pepper_vcf(vcf_fn)
    clair3_vcf_records_dict = { vcf_record.key : vcf_record for vcf_record in clair3_vcf_records_list }


    with open(new_vcf_fn, "w") as f_clair3:
        with open(vcf_fn, "r") as clair3:
            for line in clair3.readlines():
                if line.startswith("#"):
                    f_clair3.write(f"{line.strip()}\n")
                    continue
                items = line.split("\t")
                key = f"{items[0]}:{int(items[1])}"
                if clair3_vcf_records_dict[key].variant_type.name == VariantType.VARIANT_UNKNOWN.name:
                    continue

                if clair3_vcf_records_dict[key].is_INDEL() and not clair3_vcf_records_dict[key].is_heter():
                    continue



                f_clair3.write(f"{line.strip()}\n")

    read_pepper_vcf(new_vcf_fn)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="vcf",  dest="vcf_fn", type=str)
    parser.add_argument("-o", help="output vcf",  dest="output_vcf_fn", type=str)
    args = parser.parse_args()
    output_vcf_file_heter(args.vcf_fn, args.output_vcf_fn)
    # output_vcf_file_by_positions_only_snp(args.vcf_fn, args.output_vcf_fn, qual=15)
    # output_vcf_file_by_positions_only_snp(args.vcf_fn, args.output_vcf_fn)
    # output_vcf_file_filter_by_qual(args.vcf_fn, args.output_vcf_fn, qual=15)
