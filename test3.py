import vcf
from vcf import utils
import Bio.bgzf as bgzf
import os

def add_parent_variants(parent_file, child_file, out_file):
    vcf_reader_parent = vcf.Reader(filename=parent_file)
    vcf_reader_child = vcf.Reader(filename=child_file)

    vcf_writer = vcf.Writer(bgzf.BgzfWriter(filename=out_file), vcf_reader_child)

    iterate_simulatnously = utils.walk_together(vcf_reader_parent,vcf_reader_child)
    for readers in iterate_simulatnously:
        if readers[1]:
            vcf_writer.write_record(readers[1])
        elif (readers[0]):
            vcf_writer.write_record(readers[0])


# def add_parent_variants(parent_file, child_file, out_file):
#     vcf_reader_parent = vcf.Reader(filename=parent_file)
#     vcf_reader_child = vcf.Reader(filename=child_file)
#
#     output = out_file[:out_file.rfind(".gz")] if out_file.endswith(".gz") else out_file
#     vcf_writer = vcf.Writer(open(output, 'w'), vcf_reader_child)
#
#     iterate_simulatnously = utils.walk_together(vcf_reader_parent,vcf_reader_child)
#     for readers in iterate_simulatnously:
#         if readers[1]:
#             vcf_writer.write_record(readers[1])
#         elif (readers[0]):
#             vcf_writer.write_record(readers[0])
#     if out_file.endswith(".gz"):
#         bgzf.

if __name__ == "__main__":

    print("hello")

    add_parent_variants("test_arabidopsis\output_all_chroms_Internal_2_golden.vcf.gz","test_arabidopsis\output_all_chroms_B_golden.vcf.gz",'test_arabidopsis\output_all_chroms_B_golden_final.vcf.gz')

    print("world")
    print(os.path.realpath(os.path.dirname(__file__)))