import os
import sys
from os.path import exists

import pybedtools
import pandas as pd

#TODO??? -> annotations_df['track_len'] = annotations_df.end - annotations_df.start + 1
def seperate_exons_genes_intergenics(annotations_file, working_dir):
    if working_dir:
        os.chdir(working_dir)
    #TODO just for debug. Remove later
    if exists('all_chroms_annotaions.csv'):
        all_chroms_annotaions = pd.read_csv('all_chroms_annotaions.csv')
        return all_chroms_annotaions

    if annotations_file is not None:
        try:
            annotations = pybedtools.example_bedtool(annotations_file)
            chroms_annotaions = []
            for chrom in annotations.filter(lambda x: x.fields[2] == 'chromosome'):
                exon_elements = annotations.filter(lambda x: x.fields[2] == 'exon' and x.chrom == chrom.chrom)
                exon_elements = exon_elements.sort()
                exon_elements = exon_elements.merge()
                exon_elements_df = exon_elements.to_dataframe()
                exon_elements_df['feature'] = 'exon'
                exon_elements_df = exon_elements_df.loc[:,['feature','start','end']]
                exon_elements_df = exon_elements_df.drop_duplicates()

                gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene' and x.chrom == chrom.chrom)
                gene_elements = gene_elements.sort()
                intron_elements = gene_elements.subtract(exon_elements)
                intron_elements = intron_elements.sort()
                intron_elements = intron_elements.merge()
                intron_elements_df = intron_elements.to_dataframe()
                intron_elements_df['feature'] = 'intron'
                intron_elements_df = intron_elements_df.loc[:,['feature','start','end']]
                intron_elements_df = intron_elements_df.drop_duplicates()

                all_elements = annotations.filter(lambda x: x.chrom == chrom.chrom)
                all_elements = all_elements.sort()
                intergenic_elements = all_elements.subtract(gene_elements).subtract(exon_elements)
                intergenic_elements = intergenic_elements.sort()
                intergenic_elements = intergenic_elements.merge()
                intergenic_elements_df = intergenic_elements.to_dataframe()
                intergenic_elements_df['feature'] = 'intergenic'
                intergenic_elements_df = intergenic_elements_df.loc[:,['feature','start','end']]
                intergenic_elements_df = intergenic_elements_df.drop_duplicates()

                exons_genes_intergenics = pd.concat([exon_elements_df, intron_elements_df, intergenic_elements_df], ignore_index=True)
                exons_genes_intergenics.sort_values(by=['start', 'end'], inplace=True)
                exons_genes_intergenics = exons_genes_intergenics.reset_index(drop=True)
                exons_genes_intergenics.insert(0,'chrom', str(chrom.chrom))
                exons_genes_intergenics[['chrom', 'feature']] = exons_genes_intergenics[['chrom', 'feature']].astype(str)

                # exons_genes_intergenics[['chrom', 'feature']] = exons_genes_intergenics[['chrom', 'feature']].astype(str)

                chroms_annotaions.append(exons_genes_intergenics)
            #TODO remove index column
            all_chroms_annotaions = pd.concat(chroms_annotaions)
            all_chroms_annotaions.to_csv('all_chroms_annotaions.csv')
            return all_chroms_annotaions
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")

if __name__ == "__main__":
    dir = '/groups/itay_mayrose/adisivan/PanGenomeSimulator/test_arabidopsis'
    seperate_exons_genes_intergenics(sys.argv[1], dir)