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
                cds_elements_df['feature'] = Region.CDS.value
                cds_elements_df = cds_elements_df.loc[:,['feature','start','end']]
                cds_elements_df = cds_elements_df.drop_duplicates()

                gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene' and x.chrom == chrom.chrom)
                gene_elements = gene_elements.sort()
                intron_elements = gene_elements.subtract(cds_elements)
                intron_elements = intron_elements.sort()
                intron_elements = intron_elements.merge()
                intron_elements_df = intron_elements.to_dataframe()
                intron_elements_df['feature'] = Region.INTRON.value
                intron_elements_df = intron_elements_df.loc[:,['feature','start','end']]
                intron_elements_df = intron_elements_df.drop_duplicates()

                all_elements = annotations.filter(lambda x: x.chrom == chrom.chrom)
                all_elements = all_elements.sort()
                intergenic_elements = all_elements.subtract(gene_elements).subtract(cds_elements)
                intergenic_elements = intergenic_elements.sort()
                intergenic_elements = intergenic_elements.merge()
                intergenic_elements_df = intergenic_elements.to_dataframe()
                intergenic_elements_df['feature'] = Region.INTERGENIC.value
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


class AnnotatedSequence:
    _sequence_per_chrom = {}
    _code_to_annotation = {2:Region.CDS.value, 1:Region.INTRON.value, 0:Region.INTERGENIC.value}
    _annotation_to_code = {Region.CDS.value:2, Region.INTRON.value:1, Region.INTERGENIC.value:0}
    _relevant_regions = []
    _annotations_df = None

    def __init__(self, annotations_df: pd.DataFrame):
        if annotations_df is None or annotations_df.empty:
            self._sequence_per_chrom = None
            self._relevant_regions.append(Region.ALL)
            return

        self._annotations_df = annotations_df
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

    def _get_annotation_by_position(self, chrom, pos) -> (pd.DataFrame, int):
        condition = (self._annotations_df['chrom'] == chrom) & \
                    (self._annotations_df['start'] <= pos) & (self._annotations_df['end'] >= pos)
        annotation_indices = self._annotations_df.index[condition].tolist()
        if len(annotation_indices) != 1:
            raise Exception(f"Can't determine annotation for chromosome {chrom} index {pos}")
        return self._annotations_df.iloc[annotation_indices[0]], annotation_indices[0]

    def get_regions(self):
        return self._relevant_regions

    def get_region_by_index(self, chrom, index) -> Region:
        # if not self._sequence_per_chrom:
        #     return Region.ALL
        # return self._code_to_annotation[self._sequence_per_chrom[chrom][index]]
        if not self._annotations_df:
            return Region.ALL
        annotation, _ = self._get_annotation_by_position(chrom, index)
        return annotation.iloc[0, annotation.columns.get_loc('feature')].item()

    def get_annotation_start_end(self, chrom, index) -> (int, int):
        annotation, _ = self._get_annotation_by_position(chrom, index)

        return annotation.iloc[0, annotation.columns.get_loc('start')].item(),\
               annotation.iloc[0, annotation.columns.get_loc('end')].item()

    def get_encapsulating_trinuc_positions(self, chrom, index):
        cds, _ = self._get_annotation_by_position(chrom, index)
        if cds.iloc[0, cds.columns.get_loc('feature')].item() != Region.CDS.value:
            raise Exception("Trinuc is only relevant for CDS region")
        cds_start = cds.iloc[0, cds.columns.get_loc('start')].item()
        cds_end = cds.iloc[0, cds.columns.get_loc('end')].item()
        reading_frame_offset = cds.iloc[0, cds.columns.get_loc('reading_frame_offset')].item() + (index - cds_start)
        reading_frame_offset = reading_frame_offset % 3
        first = index - reading_frame_offset
        second = index + 1 - reading_frame_offset
        third = index + 2 - reading_frame_offset

        if first < cds_start or second < cds_start:
            prev_cds = self._annotations_df[(self._annotations_df['chrom'] == chrom) &
                                            (self._annotations_df['feature'] == Region.CDS.value) &
                                            (self._annotations_df['end'] < cds_start)]
            if len(prev_cds.index) != 1:
                raise Exception(f"Can't determine previous CDS for chromosome {chrom} index {index}")
            prev_cds_end = prev_cds.iloc[0, cds.columns.get_loc('end')].item()
            if second < cds_start:
                second = prev_cds_end
                first =  prev_cds_end - 1
            elif first < cds_start:
                first = prev_cds_end

        if cds_end < second or cds_end < third:
            next_cds = self._annotations_df[(self._annotations_df['chrom'] == chrom) &
                                            (self._annotations_df['feature'] == Region.CDS.value) &
                                            (self._annotations_df['start'] > index)]
            if len(next_cds.index) != 1:
                raise Exception(f"Can't determine next CDS for chromosome {chrom} index {index}")
            next_cds_start = prev_cds.iloc[0, cds.columns.get_loc('start')].item()
            if cds_end < second:
                second = next_cds_start
                third = next_cds_start + 1
            elif cds_end < third:
                third = next_cds_start

        return first, second, third


    def get_involved_annotations(self, chrom, start, end) -> pd.DataFrame:
        annotations = self._annotations_df[(self._annotations_df['chrom'] == chrom) &
                                          (start <= self._annotations_df['end'] ) &
                                          (self._annotations_df['start'] <= end)]
        return annotations

    def mute_encapsulating_gene(self, chrom, index):
        annotation, annotation_index = self._get_annotation_by_position(chrom, index)
        gene = annotation.iloc[0, annotation.columns.get_loc('gene')].item()
        gene_annotations_indices = self._annotations_df.index[(self._annotations_df['chrom'] == chrom) &
                                                              (self._annotations_df['gene'] == gene)].tolist()
        first_index = gene_annotations_indices[0]
        last_index = gene_annotations_indices[-1]

        # melt with previous intergenic region if exists
        if first_index != 0 and \
            self._annotations_df.iloc[first_index - 1]['feature'].item() == Region.INTERGENIC.value:
            first_index -= 1

        # melt with next intergenic region if exists
        if last_index + 1 != len(self._annotations_df) and \
            self._annotations_df.iloc[last_index + 1]['feature'].item() == Region.INTERGENIC.value:
            last_index += 1

        # create new intergenic region
        intergenic_start = self._annotations_df.iloc[first_index]['start'].item()
        intergenic_end = self._annotations_df.iloc[last_index]['end'].item()
        new_intergenic = pd.DataFrame({'chrom': chrom,
                                       'start': intergenic_start,
                                       'end': intergenic_end,
                                       'feature': Region.INTERGENIC.value}) #TODO add index?

        # insert new intergenic region instead of muted gene
        dfs_to_concat = []
        if(first_index != 0):
            dfs_to_concat.append(self._annotations_df.iloc[:first_index - 1])
        dfs_to_concat.append(new_intergenic)
        if(last_index + 1 != len(self._annotations_df)):
            dfs_to_concat.append(self._annotations_df.iloc[last_index + 1:])
        self._annotations_df = pd.concat(dfs_to_concat).reset_index(drop=True)


    def handle_insertion(self, chrom, index, insertion_len):
        # TODO implement - elongate sequence and shift downstream genes, including current (if not intergenic region)
        pass


    def handle_deletion(self, chrom, start, end):
        # TODO implement - shorten sequence and shift downstream genes, including current (if not intergenic region)
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