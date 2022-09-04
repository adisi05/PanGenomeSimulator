import os
from os.path import exists
from typing import Optional

import numpy as np
import pybedtools
import pandas as pd

from neat.utilities.common_data_structues import Region, Strand


def to_annotations_df(args, working_dir):
    # Process bed/gff file
    annotations_df = None
    if args.b:
        print('Processing bed/gff file...')
        try:
            annotations_df = separate_cds_genes_intergenics(args.b, working_dir)
            annotations_df[['chrom', 'region']] = annotations_df[['chrom', 'region']].astype(str)
            # TODO validate chromosome length are same as in reference
        except ValueError:
            print('Problem parsing bed file. Ensure bed file is tab separated, standard bed format')
    return annotations_df

def separate_cds_genes_intergenics(annotations_file, working_dir):
    # GENE_ID = 'gene_id'
    # def extract_gene_id(attributes: str):
    #     return next(filter(lambda x: x.startswith(GENE_ID), attributes.split(';')), None).split('=')[1]

    def assign_gene(start: int, end: int, gene_elements_df: pd.DataFrame) -> (int, int):
        gene_ids = gene_elements_df.index[
            (gene_elements_df['start'] - 1 <= start) & (end <= gene_elements_df['end'])].tolist()
        # gene_elements_df['start'] - 1 correction is because of some bug in pybedtools

        if len(gene_ids) < 1:
            return 0, Strand.UNKNOWN.value
        if len(gene_ids) > 1:
            # We assume genes from different strands don't overlap, but if it happens nonetheless - ignore
            return 0, Strand.UNKNOWN.value
        strand = Strand(gene_elements_df.iloc[gene_ids[0]]['strand'])
        return gene_ids[0] + 1, strand.value # indecies in pandas are 0-based

    if working_dir:
        os.chdir(working_dir)
    #TODO just for debug. Remove later
    if exists('all_chroms_annotaions.csv'):
        all_chroms_annotaions = pd.read_csv('all_chroms_annotaions.csv')
        return all_chroms_annotaions

    if annotations_file is not None:
        try:
            annotations = pybedtools.example_bedtool(annotations_file)
            all_chroms_annotaions = []
            for chrom in annotations.filter(lambda x: x.fields[2] == 'chromosome'):
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                cds_elements = cds_elements.sort()
                cds_elements = cds_elements.merge()
                cds_elements_df = cds_elements.to_dataframe()
                cds_elements_df['region'] = Region.CDS.value
                del cds_elements_df['chrom']
                cds_elements_df = cds_elements_df.drop_duplicates()

                gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene' and x.chrom == chrom.chrom)
                gene_elements = gene_elements.sort()
                non_coding_gene_elements = gene_elements.subtract(cds_elements)
                non_coding_gene_elements = non_coding_gene_elements.sort()
                non_coding_gene_elements = non_coding_gene_elements.merge()
                non_coding_gene_elements_df = non_coding_gene_elements.to_dataframe()
                non_coding_gene_elements_df['region'] = Region.NON_CODING_GENE.value
                del non_coding_gene_elements_df['chrom']
                non_coding_gene_elements_df = non_coding_gene_elements_df.drop_duplicates()

                all_elements = annotations.filter(lambda x: x.chrom == chrom.chrom)
                all_elements = all_elements.sort()
                intergenic_elements = all_elements.subtract(gene_elements).subtract(cds_elements)
                intergenic_elements = intergenic_elements.sort()
                intergenic_elements = intergenic_elements.merge()
                intergenic_elements_df = intergenic_elements.to_dataframe()
                intergenic_elements_df['region'] = Region.INTERGENIC.value
                del intergenic_elements_df['chrom']
                intergenic_elements_df = intergenic_elements_df.drop_duplicates()

                cds_genes_intergenics = pd.concat([cds_elements_df, non_coding_gene_elements_df, intergenic_elements_df], ignore_index=True)
                cds_genes_intergenics.sort_values(by=['start', 'end'], inplace=True)
                cds_genes_intergenics.insert(0,'chrom', str(chrom.chrom))
                cds_genes_intergenics[['chrom', 'region']] = cds_genes_intergenics[['chrom', 'region']].astype(str)

                # cds_genes_intergenics[['chrom', 'region']] = cds_genes_intergenics[['chrom', 'region']].astype(str)

                gene_elements_df = gene_elements.to_dataframe()
                cds_genes_intergenics[['gene','strand']] = cds_genes_intergenics.apply(
                    lambda x: assign_gene(x['start'], x['end'], gene_elements_df), axis=1, result_type='expand')

                # add_reading_frames_test(cds_genes_intergenics, str(chrom.chrom))

                all_chroms_annotaions.append(cds_genes_intergenics)
            all_chroms_annotaions = pd.concat(all_chroms_annotaions).reset_index(drop=True)
            all_chroms_annotaions.to_csv('all_chroms_annotaions.csv')
            return all_chroms_annotaions
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")


# TODO for test. remove later
# def add_reading_frames_test(annotations_df, chromosome):
#     if annotations_df is None or annotations_df.empty:
#         return
#
#     for gene in annotations_df.gene.unique():
#         if gene == 0:
#             continue # intergenic region
#
#         cds_annotations = annotations_df[(annotations_df['gene'] == gene) &
#                                             (annotations_df['region'] == Region.CDS.value)]
#         cds_total_len = cds_annotations['end'].sum()-cds_annotations['start'].sum()
#         if cds_total_len % 3 != 0:
#             #TODO raise Exception?
#             print(f'Total length of CDS elements in gene {gene} in chromosome {chromosome} '
#                             f'cannot be divided by 3')
#         reading_offset = 0
#         for _, annotation in cds_annotations.iterrows():
#             annotation['reading_offset'] = reading_offset
#             reading_offset += annotation['end'] - annotation['start']
#             reading_offset = reading_offset % 3




class AnnotatedSequence:
    _relevant_regions = []
    _chromosome = None
    """
    annotations are 0-based, the length of each is end-start
    (meaning end is not included in the annotation positions, it's the position right after)
    """
    _annotations_df = None

    """ Cache parameters for recent window """
    _cached_start = None
    _cached_end = None
    _cached_mask_in_window_per_region = None

    def __init__(self, annotations_df: pd.DataFrame, chromosome : str, is_sorted : bool =False):
        self._chromosome = chromosome
        if annotations_df is None or annotations_df.empty:
            self._relevant_regions.append(Region.ALL)
            return

        if 'chrom' in annotations_df.columns:
            self._annotations_df = annotations_df[annotations_df['chrom'] == chromosome].copy()
            del self._annotations_df['chrom']
            self._annotations_df.reset_index(drop=True, inplace=True)
        else:
            self._annotations_df = annotations_df.copy()

        if not is_sorted:
            self._annotations_df.sort_values('start', inplace=True)
            self._annotations_df.reset_index(drop=True, inplace=True)


        relevant_region_names = self._annotations_df['region'].unique()
        self._relevant_regions = [Region(name) for name in relevant_region_names]

        if len(self._relevant_regions) == 0:
            self._relevant_regions.append(Region.ALL)

        self._add_reading_frames()

    def len(self):
        if self._annotations_df is None or self._annotations_df.empty:
            return None
        return self._annotations_df.iloc[-1]['end'].item()

    def _get_annotation_by_position(self, pos) -> (pd.DataFrame, int):
        condition = (self._annotations_df['start'] <= pos) & (pos < self._annotations_df['end'])
        annotation_indices = self._annotations_df.index[condition].tolist()
        if len(annotation_indices) != 1:
            raise Exception(f"Can't determine annotation for chromosome {self._chromosome} position {pos}")
        return self._annotations_df.iloc[annotation_indices[0]], annotation_indices[0]

    def get_regions(self) -> list[Region]:
        return self._relevant_regions

    def get_region_by_position(self, pos) -> (Region, Strand):
        if not self._annotations_df:
            return Region.ALL
        annotation, _ = self._get_annotation_by_position(pos)
        return annotation.iloc[0]['region'].item(), Strand(annotation.iloc[0]['strand'].item())

    def get_annotation_start_end(self, pos) -> (int, int):
        annotation, _ = self._get_annotation_by_position(pos)

        return annotation.iloc[0]['start'].item(),\
               annotation.iloc[0]['end'].item()

    def get_last_coding_position_of_encapsulating_gene(self, pos) -> (Optional[int], Strand):
        annotation, annotation_index = self._get_annotation_by_position(pos)

        gene = annotation.iloc[0]['gene'].item()
        if not gene or gene == 0:
            raise Exception("Coding positions are only relevant for genes")

        strand = Strand(annotation.iloc[0]['strand'].item())
        if strand == Strand.UNKNOWN:
            return None, strand

        gene_cds_annotations = self._annotations_df[(self._annotations_df['gene'] == gene) &
                                                    (self._annotations_df['region'] == Region.CDS.value)]

        if strand == Strand.FORWARD:
            last_cds = gene_cds_annotations[-1]
            return last_cds['end'].item() - 1, strand
        else:
            last_cds = gene_cds_annotations[0]
            return last_cds['start'].item(), strand

    def get_encapsulating_codon_positions(self, pos) -> (int, int, int):
        cds, _ = self._get_annotation_by_position(pos)
        if cds.iloc[0]['region'].item() != Region.CDS.value:
            raise Exception("Codon is only relevant for CDS region")
        cds_start = cds.iloc[0]['start'].item()
        cds_end = cds.iloc[0]['end'].item()
        reading_offset = cds.iloc[0]['reading_offset'].item() + (pos - cds_start)
        reading_offset = reading_offset % 3
        first = pos - reading_offset
        second = pos + 1 - reading_offset
        third = pos + 2 - reading_offset

        if first < cds_start or second < cds_start:
            prev_cds = self._annotations_df[(self._annotations_df['region'] == Region.CDS.value) &
                                            (self._annotations_df['end'] <= cds_start)]
            if len(prev_cds.index) != 1:
                raise Exception(f"Can't determine previous CDS for chromosome {self._chromosome} index {pos}")
            prev_cds_end = prev_cds.iloc[0]['end'].item()
            if second < cds_start:
                second = prev_cds_end - 1
                first =  prev_cds_end - 2
            elif first < cds_start:
                first = prev_cds_end - 1

        if cds_end <= second or cds_end <= third:
            next_cds = self._annotations_df[(self._annotations_df['region'] == Region.CDS.value) &
                                            (cds_end <= self._annotations_df['start'])]
            if len(next_cds.index) != 1:
                raise Exception(f"Can't determine next CDS for chromosome {self._chromosome} position {pos}")
            next_cds_start = next_cds.iloc[0]['start'].item()
            if cds_end <= second:
                second = next_cds_start
                third = next_cds_start + 1
            elif cds_end < third:
                third = next_cds_start

        return first, second, third


    def get_annotations_in_range(self, start, end) -> pd.DataFrame:
        annotations = self._annotations_df[(start < self._annotations_df['end'] ) &
                                           (self._annotations_df['start'] < end)]
        return annotations

    def mute_encapsulating_gene(self, pos) -> None:
        annotation, annotation_index = self._get_annotation_by_position(pos)
        gene = annotation.iloc[0]['gene'].item()
        gene_annotations_indices = self._annotations_df.index[(self._annotations_df['gene'] == gene)].tolist()
        first_index = gene_annotations_indices[0]
        last_index = gene_annotations_indices[-1]

        # melt with previous intergenic region if exists
        if first_index != 0 and \
            self._annotations_df.iloc[first_index - 1]['region'].item() == Region.INTERGENIC.value:
            first_index -= 1

        # melt with next intergenic region if exists
        if last_index + 1 != len(self._annotations_df) and \
            self._annotations_df.iloc[last_index + 1]['region'].item() == Region.INTERGENIC.value:
            last_index += 1

        # create new intergenic region
        intergenic_start = self._annotations_df.iloc[first_index]['start'].item()
        intergenic_end = self._annotations_df.iloc[last_index]['end'].item()
        new_intergenic = pd.DataFrame({'start': intergenic_start,
                                       'end': intergenic_end,
                                       'region': Region.INTERGENIC.value},
                                      index=[first_index])

        # insert new intergenic region instead of muted gene
        dfs_to_concat = []
        if(first_index != 0):
            dfs_to_concat.append(self._annotations_df.iloc[:first_index])
        dfs_to_concat.append(new_intergenic)
        if(last_index + 1 != len(self._annotations_df)):
            dfs_to_concat.append(self._annotations_df.iloc[last_index + 1:])
        self._annotations_df = pd.concat(dfs_to_concat).reset_index(drop=True)


    def handle_insertion(self, pos, insertion_len) -> None:
        _, index = self._get_annotation_by_position(pos)
        self._annotations_df.iloc[index]['end'] += insertion_len
        if index + 1 != len(self._annotations_df):
            self._annotations_df.iloc[index+1:]['start'] += insertion_len
            self._annotations_df.iloc[index+1:]['end'] += insertion_len


    def handle_deletion(self, pos, deletion_len) -> None:
        """
        Delete right after pos, sequence at the length deletion_len
        :param start:
        :param end:
        :return:
        """
        annotation, index = self._get_annotation_by_position(pos)
        annotation_residue = self._annotations_df.iloc[index]['end'].item() - pos - 1
        deleted_already = min(annotation_residue, deletion_len)
        annotation['end'] = annotation['end'] - deleted_already
        index += 1

        annotations_to_delete = []
        # deletion involves more than one annotation
        while deleted_already < deletion_len and index < len(self._annotations_df):
            annotation_len = self._annotations_df.iloc[index]['end'].item() - \
                             self._annotations_df.iloc[index]['start'].item()
            if annotation_len <= deletion_len - deleted_already:
                annotations_to_delete.append(index)
                deleted_already += annotation_len
            else:
                self._annotations_df.iloc[index]['start'] = pos + 1
                self._annotations_df.iloc[index]['end'] = (pos + annotation_len) - (deletion_len - deleted_already)
                deleted_already = deletion_len
            index += 1

        if index != len(self._annotations_df):
            self._annotations_df.iloc[index:]['start'] -= deletion_len
            self._annotations_df.iloc[index:]['end'] -= deletion_len

        if len(annotations_to_delete):
            self._annotations_df = self._annotations_df.drop(annotations_to_delete).reset_index(drop=True)

    def get_nucleotides_counts_per_region(self, start=-1, end=-1) -> dict:
        start = start if start != -1 else 0
        end = end if end != -1 else self.len()
        relevant_annotations = self.get_annotations_in_range(start, end)
        counts_per_region = {}
        for _, annotation in relevant_annotations.iterrows():
            region = annotation['region'].item()
            if region not in counts_per_region:
                counts_per_region[region] = 0
            region_start = max(start, annotation['start'].item())
            region_end = min(end, annotation['end'].item())
            counts_per_region[region] += region_end - region_start
        return counts_per_region

    def get_mask_in_window_of_region(self, relevant_region, start, end) -> dict:
        # if cache is valid and contains the region - return it
        if start == self._cached_start and end == self._cached_end and relevant_region in self._cached_mask_in_window_per_region:
            return self._cached_mask_in_window_per_region[relevant_region]

        # else - initialize cache
        self._cached_start = start
        self._cached_end = end
        self._cached_mask_in_window_per_region = {region: np.array([]) for region in self.get_regions()}

        # compute if no cache
        relevant_annotations = self.get_annotations_in_range(start, end)
        for _, annotation in relevant_annotations.iterrows():
            annotation_start = max(start, annotation['start'].item())
            annotation_end = min(end, annotation['end'].item())
            annotation_length = annotation_end - annotation_start
            for region in self._cached_mask_in_window_per_region.keys():
                mask = 1 if region == relevant_region else 0
                self._cached_mask_in_window_per_region[region] = np.concatenate(
                    self._cached_mask_in_window_per_region[region], mask * annotation_length)

        return self._cached_mask_in_window_per_region[relevant_region]

    def _add_reading_frames(self): #TODO move to creation of dataframe?
        if self._annotations_df is None or self._annotations_df.empty:
            return

        for gene in self._annotations_df.gene.unique():
            if gene == 0:
                continue  # intergenic region

            cds_annotations = self._annotations_df[(self._annotations_df['gene'] == gene) &
                                             (self._annotations_df['region'] == Region.CDS.value)]
            cds_total_len = cds_annotations['end'].sum() - cds_annotations['start'].sum()
            if cds_total_len % 3 != 0:
                # TODO raise Exception?
                print(f'Total length of CDS elements in gene {gene} in chromosome {self._chromosome} '
                      f'cannot be divided by 3')
            reading_offset = 0
            for _, annotation in cds_annotations.iterrows():
                annotation['reading_offset'] = reading_offset
                reading_offset += annotation['end'] - annotation['start']
                reading_offset = reading_offset % 3


if __name__ == "__main__":
    file = '/groups/itay_mayrose/adisivan/arabidopsis/ensemblgenomes/gff3/Arabidopsis_thaliana.TAIR10.53_1000.gff3'
    dir = '/groups/itay_mayrose/adisivan/PanGenomeSimulator/test_arabidopsis'
    separate_cds_genes_intergenics(file, dir)

    #ID=gene:AT1G01100;Name=RPP1A;biotype=protein_coding;description=60S acidic ribosomal protein P1-1 [Source:UniProtKB/Swiss-Prot%3BAcc:Q8LCW9];gene_id=AT1G01100;logic_name=araport11
