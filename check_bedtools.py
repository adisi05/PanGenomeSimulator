import sys

import pybedtools
from pybedtools import BedTool
import pandas as pd

def load_mutation_regions(annotations_file, squash_annotations=False):
    if annotations_file is not None:
        try:
            annotations = pybedtools.example_bedtool(annotations_file)
            annotations = annotations.sort() #TODO needed indeed?
            annotation_types = ['exon', 'CDS', 'five_prime_UTR', 'mRNA', 'three_prime_UTR', 'gene', 'ncRNA_gene',
                                'lnc_RNA', 'tRNA', 'ncRNA', 'miRNA', 'snoRNA', 'snRNA', 'rRNA']
            for chrom in annotations.filter(lambda x: x.fields[2] == 'chromosome'):
                print('----------------------------------------------------------------------')
                print(f'chromosome {chrom.name} length is {chrom.length}')
                for strand in ['+', '-']:
                    print(f'strand {strand} :')
                    for current_type in annotation_types:
                        total_current_type_length = 0
                        annotations_of_type = annotations.filter(lambda x: x.fields[2] == current_type and x.chrom == chrom.chrom and x.strand == strand)
                        annotations_of_type = annotations_of_type.saveas(f'annotations_of_type_{current_type}_strand_{strand}.bed') #, trackline='track name="a and b"')
                        annotations_of_type = annotations_of_type.sort() #TODO needed indeed?
                        if squash_annotations:
                            annotations_of_type = annotations_of_type.merge()
                        for element in annotations_of_type:
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
                exon_intersect_cds_df = exon_intersect_cds_df.drop_duplicates()
                exon_intersect_cds_df.to_csv('exon_intersect_cds.csv')

                # -v
                exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                # exon_non_cds = exon_elements.intersect(cds_elements, v=True)
                exon_non_cds = exon_elements.subtract(cds_elements)
                exon_non_cds_df = exon_non_cds.to_dataframe()
                exon_non_cds_df.to_csv('exon_non_cds.csv')

                # -u
                exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                exon_intersect_cds_unique = exon_elements.intersect(cds_elements, u=True)
                exon_intersect_cds_unique_df = exon_intersect_cds_unique.to_dataframe()
                exon_intersect_cds_unique_df.to_csv('exon_intersect_cds_unique.csv')


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
                cds_intersect_exon_df = cds_intersect_exon_df.drop_duplicates()
                cds_intersect_exon_df.to_csv('cds_intersect_exon.csv')

                # -v
                exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                # cds_non_exon = cds_elements.intersect(exon_elements, v=True)
                cds_non_exon = cds_elements.subtract(exon_elements)
                cds_non_exon_df = cds_non_exon.to_dataframe()
                cds_non_exon_df.to_csv('cds_non_exon.csv')


                # -u
                exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                cds_intersect_exon_unique = cds_elements.intersect(exon_elements, u=True)
                cds_intersect_exon_unique_df = cds_intersect_exon_unique.to_dataframe()
                cds_intersect_exon_unique_df.to_csv('cds_intersect_exon_unique.csv')

                break # TODO remove - currently we're checking only the first chromosome
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")
            sys.exit(1)


def seperate_exons_genes_intergenics(annotations_file):
    if annotations_file is not None:
        try:
            annotations = pybedtools.example_bedtool(annotations_file)
            for chrom in annotations.filter(lambda x: x.fields[2] == 'chromosome'):
                exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                exon_elements = exon_elements.sort()
                exon_elements = exon_elements.merge()
                exon_elements = exon_elements.saveas(f'exon_elements_chrom_{chrom.chrom}.bed')
                exon_elements_df = exon_elements.to_dataframe()
                exon_elements_df['feature'] = 'exon'
                exon_elements_df = exon_elements_df.loc[:,['feature','start','end']]
                exon_elements_df = exon_elements_df.drop_duplicates()
                exon_elements_df.to_csv(f'exon_elements_chrom_{chrom.chrom}.csv')

                # exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene' and x.chrom == chrom.chrom)
                gene_elements = gene_elements.sort()
                gene_elements = gene_elements.saveas(f'gene_elements_chrom_{chrom.chrom}.bed')
                intron_elements = gene_elements.subtract(exon_elements)
                intron_elements = intron_elements.sort()
                intron_elements = intron_elements.merge()
                intron_elements_df = intron_elements.to_dataframe()
                intron_elements_df['feature'] = 'intron'
                intron_elements_df = intron_elements_df.loc[:,['feature','start','end']]
                intron_elements_df = intron_elements_df.drop_duplicates()
                intron_elements_df.to_csv(f'intron_elements_chrom_{chrom.chrom}.csv')


                # gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene' and x.chrom == chrom.chrom)
                all_elements = annotations.filter(lambda x: x.chrom == chrom.chrom)
                all_elements = all_elements.sort()
                all_elements = all_elements.saveas(f'all_elements_chrom_{chrom.chrom}.bed')
                intergenic_elements = all_elements.subtract(gene_elements).subtract(exon_elements)
                intergenic_elements = intergenic_elements.sort()
                intergenic_elements = intergenic_elements.merge()
                intergenic_elements_df = intergenic_elements.to_dataframe()
                intergenic_elements_df['feature'] = 'intergenic'
                intergenic_elements_df = intergenic_elements_df.loc[:,['feature','start','end']]
                intergenic_elements_df = intergenic_elements_df.drop_duplicates()
                intergenic_elements_df.to_csv(f'intergenic_elements_chrom_{chrom.chrom}.csv')

                exons_genes_intergenics = pd.concat([exon_elements_df, intron_elements_df, intergenic_elements_df], ignore_index=True)
                exons_genes_intergenics.sort_values(by=['start', 'end'], inplace=True)
                exons_genes_intergenics = exons_genes_intergenics.reset_index(drop=True)
                exons_genes_intergenics.to_csv(f'exons_genes_intergenics_chrom_{chrom.chrom}.csv')
                for i, row in exons_genes_intergenics.iterrows():
                    if i == 0:
                        continue
                    prev = exons_genes_intergenics.iloc[i-1,:]
                    curr = row
                    if prev['end'] != curr['start']:
                        print(f'chromosome {chrom.chrom}: index {i}')


        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")
            sys.exit(1)


if __name__ == "__main__":
    # args: /groups/itay_mayrose/adisivan/arabidopsis/ensemblgenomes/gff3/Arabidopsis_thaliana.TAIR10.53.gff3
    # load_mutation_regions(sys.argv[1], squash_annotations=True)
    # cds_exon_intersection(sys.argv[1])
    seperate_exons_genes_intergenics(sys.argv[1])