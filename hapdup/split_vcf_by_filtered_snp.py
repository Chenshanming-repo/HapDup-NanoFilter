from hapdup.MultiStepConsistencyAnalysisHapDup import read_pepper_vcf
from hapdup.VariantRecord import VariantType

def main(vcf_raw, vcf_filtered, out_dir):
    vcf_records_list1 = read_pepper_vcf(vcf_raw)
    vcf_records_dict1 = { vcf_record.key : vcf_record for vcf_record in vcf_records_list1 }
    vcf_records_list2 = read_pepper_vcf(vcf_filtered)
    vcf_records_dict2 = { vcf_record.key : vcf_record for vcf_record in vcf_records_list2 }

    with open(f"{out_dir}/filter_snp.vcf", "w") as f_clair3:
        with open(vcf_raw, "r") as clair3:
            for line in clair3.readlines():
                if line.startswith("#"):
                    f_clair3.write(f"{line.strip()}\n")
                    continue
                items = line.split("\t")
                key = f"{items[0]}:{int(items[1])}"
                if vcf_records_dict1[key].is_INDEL():
                    continue
                if not key in vcf_records_dict2:
                    f_clair3.write(f"{line.strip()}\n")
    with open(f"{out_dir}/remain_snp.vcf", "w") as f_clair3:
        with open(vcf_raw, "r") as clair3:
            for line in clair3.readlines():
                if line.startswith("#"):
                    f_clair3.write(f"{line.strip()}\n")
                    continue
                items = line.split("\t")
                key = f"{items[0]}:{int(items[1])}"
                if vcf_records_dict1[key].is_INDEL():
                    continue
                if key in vcf_records_dict2:
                    f_clair3.write(f"{line.strip()}\n")


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
    parser.add_argument("--raw", help="raw vcf",  dest="raw_vcf_fn", type=str)
    parser.add_argument("--filter", help="filtered vcf",  dest="filtered_vcf_fn", type=str)
    parser.add_argument("--out_dir", help="output dir",  dest="output_dir", type=str)
    args = parser.parse_args()
    # output_vcf_file_heter(args.vcf_fn, args.output_vcf_fn)
    main(args.raw_vcf_fn, args.filtered_vcf_fn, args.output_dir)
    # output_vcf_file_by_positions_only_snp(args.vcf_fn, args.output_vcf_fn)
    # output_vcf_file_filter_by_qual(args.vcf_fn, args.output_vcf_fn, qual=10)
