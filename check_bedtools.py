import sys

import pybedtools
from pybedtools import BedTool

def load_mutation_regions(annotations_file):
    if annotations_file is not None:
        try:
            annotations = BedTool(annotations_file)
            annotation_types = ['exon', 'CDS', 'five_prime_UTR', 'mRNA', 'three_prime_UTR', 'gene', 'ncRNA_gene',
                                'lnc_RNA', 'tRNA', 'ncRNA', 'miRNA', 'snoRNA', 'snRNA', 'rRNA']
            for chrom in annotations.filter(lambda x: x.fields[2] == 'chromosome'):
                print('----------------------------------------------------------------------')
                print(f'chromosome {chrom.name} length is {chrom.length}')
                for current_type in annotation_types:
                    total_current_type_length = 0
                    for element in annotations.filter(lambda x: x.fields[2] == current_type and x.chrom == chrom.chrom):
                        total_current_type_length += element.length
                    print(f'{current_type}\t{total_current_type_length}')
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")
            sys.exit(1)

def cds_exon_intersection(annotations_file):
    if annotations_file is not None:
        try:
            annotations = pybedtools.example_bedtool(annotations_file)
            for chrom in annotations.filter(lambda x: x.fields[2] == 'chromosome'):
                print(f'chromosome {chrom.name} length is {chrom.length}')
                print('----------------------------- First try -----------------------------')
                exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                # total_length = 0
                # for element in exon_elements.intersect(cds_elements):
                #     # print(element)
                #     total_length += element.length
                # print(f'exon+CDS\t{total_length}')
                exon_intersect_cds = exon_elements.intersect(cds_elements)
                exon_intersect_cds_df = exon_intersect_cds.to_dataframe()
                exon_intersect_cds_df.to_csv('exon_intersect_cds.csv')

                # -v
                exon_non_cds = exon_elements.intersect(cds_elements, v=True)
                exon_non_cds_df = exon_non_cds.to_dataframe()
                exon_non_cds_df.to_csv('exon_non_cds_df')

                # -u
                exon_intersect_cds_unique = exon_elements.intersect(cds_elements, u=True)
                exon_intersect_cds_unique_df = exon_intersect_cds_unique.to_dataframe()
                exon_intersect_cds_unique_df.to_csv('exon_intersect_cds_unique_df')


                print('----------------------------- Second try -----------------------------')
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                # total_length = 0
                # for element in cds_elements.intersect(exon_elements):
                #     # print(element)
                #     total_length += element.length
                # print(f'exon+CDS\t{total_length}')
                cds_intersect_exon = cds_elements.intersect(exon_elements)
                cds_intersect_exon_df = cds_intersect_exon.to_dataframe()
                cds_intersect_exon_df.to_csv('cds_intersect_exon.csv')

                # -v
                cds_non_exon = cds_elements.intersect(exon_elements, v=True)
                cds_non_exon_df = cds_non_exon.to_dataframe()
                cds_non_exon_df.to_csv('cds_non_exon_df')

                # -u
                cds_intersect_exon_unique = cds_elements.intersect(exon_elements, u=True)
                cds_intersect_exon_unique_df = cds_intersect_exon_unique.to_dataframe()
                cds_intersect_exon_unique_df.to_csv('cds_intersect_exon_unique_df')

                break # TODO remove - currently we're checking only the first chromosome
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")
            sys.exit(1)


if __name__ == "__main__":
    cds_exon_intersection(sys.argv[1])