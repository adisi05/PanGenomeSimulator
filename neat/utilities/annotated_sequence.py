import os
from os.path import exists

import numpy as np
import pybedtools
import pandas as pd

from neat.utilities.common_data_structues import Region

#TODO check it's working, especially now with CDS..
def to_annotations_df(args, working_dir):
    # Process bed/gff file
    annotations_df = None
    if args.b:
        print('Processing bed/gff file...')
        try:
            annotations_df = seperate_cds_genes_intergenics(args.b, working_dir)
            annotations_df[['chrom', 'feature']] = annotations_df[['chrom', 'feature']].astype(str)
            # TODO validate chromosome length are same as in reference
        except ValueError:
            print('Problem parsing bed file. Ensure bed file is tab separated, standard bed format')
    return annotations_df

#TODO??? -> annotations_df['track_len'] = annotations_df.end - annotations_df.start + 1
def seperate_cds_genes_intergenics(annotations_file, working_dir):
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
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                cds_elements = cds_elements.sort()
                cds_elements = cds_elements.merge()
                cds_elements_df = cds_elements.to_dataframe()
                cds_elements_df['feature'] = 'CDS'
                cds_elements_df = cds_elements_df.loc[:,['feature','start','end']]
                cds_elements_df = cds_elements_df.drop_duplicates()

                gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene' and x.chrom == chrom.chrom)
                gene_elements = gene_elements.sort()
                intron_elements = gene_elements.subtract(cds_elements)
                intron_elements = intron_elements.sort()
                intron_elements = intron_elements.merge()
                intron_elements_df = intron_elements.to_dataframe()
                intron_elements_df['feature'] = 'intron'
                intron_elements_df = intron_elements_df.loc[:,['feature','start','end']]
                intron_elements_df = intron_elements_df.drop_duplicates()

                all_elements = annotations.filter(lambda x: x.chrom == chrom.chrom)
                all_elements = all_elements.sort()
                intergenic_elements = all_elements.subtract(gene_elements).subtract(cds_elements)
                intergenic_elements = intergenic_elements.sort()
                intergenic_elements = intergenic_elements.merge()
                intergenic_elements_df = intergenic_elements.to_dataframe()
                intergenic_elements_df['feature'] = 'intergenic'
                intergenic_elements_df = intergenic_elements_df.loc[:,['feature','start','end']]
                intergenic_elements_df = intergenic_elements_df.drop_duplicates()

                cds_genes_intergenics = pd.concat([cds_elements_df, intron_elements_df, intergenic_elements_df], ignore_index=True)
                cds_genes_intergenics.sort_values(by=['start', 'end'], inplace=True)
                cds_genes_intergenics = cds_genes_intergenics.reset_index(drop=True)
                cds_genes_intergenics.insert(0,'chrom', str(chrom.chrom))
                cds_genes_intergenics[['chrom', 'feature']] = cds_genes_intergenics[['chrom', 'feature']].astype(str)

                # cds_genes_intergenics[['chrom', 'feature']] = cds_genes_intergenics[['chrom', 'feature']].astype(str)

                chroms_annotaions.append(cds_genes_intergenics)
            #TODO remove index column
            all_chroms_annotaions = pd.concat(chroms_annotaions)
            all_chroms_annotaions.to_csv('all_chroms_annotaions.csv')
            return all_chroms_annotaions
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")


class AnnotatedSeqence:
    _sequence_per_chrom = {}
    _code_to_annotation = {2:Region.CDS.value, 1:Region.INTRON.value, 0:Region.INTERGENIC.value}
    _annotation_to_code = {Region.CDS.value:2, Region.INTRON.value:1, Region.INTERGENIC.value:0}
    _relevant_regions = []

    def __init__(self, annotations_df: pd.DataFrame):
        if annotations_df is None or annotations_df.empty:
            self._sequence_per_chrom = None
            self._relevant_regions.append(Region.ALL)
            return

        for i, annotation in annotations_df.iterrows():
            current_chrom = annotation['chrom']
            if not current_chrom in self._sequence_per_chrom:
                self._sequence_per_chrom[current_chrom] = np.array([])
            region = Region(annotation['feature'])
            if region not in self._relevant_regions:
                self._relevant_regions(region)
            annotation_length = annotation['end'] - annotation['start']
            current_sequence = [self._annotation_to_code[region.value]] * annotation_length
            self._sequence_per_chrom[current_chrom] = \
                np.concatenate(self._sequence_per_chrom[current_chrom], current_sequence)
        if len(self._relevant_regions) == 0:
            self._relevant_regions.append(Region.ALL)

    def get_regions(self):
        return self._relevant_regions

    def get_annotation(self, chrom, index):
        if not self._sequence_per_chrom:
            return Region.ALL
        return self._code_to_annotation[self._sequence_per_chrom[chrom][index]]

    def get_annotation_start_end(self, chrom, index):
        # TODO implement
        pass

    def get_nucleotides_counts_per_region(self, chrom, start=-1, end=-1):
        start = start if start != -1 else 0
        end = end if end != -1 else len(self._sequence_per_chrom[chrom])
        window = self._sequence_per_chrom[chrom][start:end]
        counts_per_region = {}
        for region in self.get_regions():
            counts_per_region[region] = window.count(str(self._annotation_to_code(region)))
        return counts_per_region

    def get_mask_in_window_of_region(self, region, start, end):
        # check for cache
        if start != self._recent_start or end != self._recent_end:
            self._recent_start = start
            self._recent_end = end
            self._mask_in_window_per_region = {}
        if region in self._mask_in_window_per_region:
            return self._mask_in_window_per_region[region]

        # compute if no cache
        window_sequence = self._sequence_per_chrom[self.chromosome_name][start:end]
        relevant_code = self._annotation_to_code[region.value]
        self._mask_in_window_per_region[region] = np.where(window_sequence == relevant_code, 1, 0)
        return self._mask_in_window_per_region[region]