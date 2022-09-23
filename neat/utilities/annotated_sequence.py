import os
from os.path import exists
from typing import Optional, List

import numpy as np
import pybedtools
import pandas as pd

from common_data_structues import Region, Strand


def to_annotations_df(file_path, output_dir=None):
    # Process bed/gff file
    annotations_df = None
    if file_path:
        print('Processing bed/gff file...')
        try:
            annotations_df = separate_cds_genes_intergenics(file_path, output_dir)
            annotations_df[['chrom', 'region']] = annotations_df[['chrom', 'region']].astype(str)
            # TODO validate chromosome length are same as in reference
        except ValueError:
            print('Problem parsing bed file. Ensure bed file is tab separated, standard bed format')
    return annotations_df


def separate_cds_genes_intergenics(annotations_file, output_dir=None):
    # GENE_ID = 'gene_id'
    # def extract_gene_id(attributes: str):
    #     return next(filter(lambda x: x.startswith(GENE_ID), attributes.split(';')), None).split('=')[1]

    # check if CSV exists
    _, extension = os.path.splitext(annotations_file)
    if extension == '.csv':
        csv_file = annotations_file
    else:
        # TODO remove?
        csv_file = 'all_chroms_annotations.csv'
        if output_dir:
            csv_file = os.path.join(output_dir, csv_file)
    if exists(csv_file):
        all_chroms_annotations = pd.read_csv(csv_file)
        return all_chroms_annotations

    if annotations_file is not None:
        try:
            annotations = pybedtools.example_bedtool(annotations_file)
            all_chroms_annotations = []
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

                cds_genes_intergenics = pd.concat([cds_elements_df, non_coding_gene_elements_df,
                                                   intergenic_elements_df], ignore_index=True)
                cds_genes_intergenics.sort_values(by=['start', 'end'], inplace=True)
                cds_genes_intergenics.insert(0, 'chrom', str(chrom.chrom))
                cds_genes_intergenics[['chrom', 'region']] = cds_genes_intergenics[['chrom', 'region']].astype(str)

                # cds_genes_intergenics[['chrom', 'region']] = cds_genes_intergenics[['chrom', 'region']].astype(str)

                gene_elements_df = gene_elements.to_dataframe()
                gene_elements_df = gene_elements_df[['strand', 'start', 'end']]
                cds_genes_intergenics[['gene', 'strand']] = cds_genes_intergenics.apply(
                    lambda x: assign_gene(x['start'], x['end'], gene_elements_df), axis=1, result_type='expand')

                # add_reading_frames_test(cds_genes_intergenics, str(chrom.chrom))

                all_chroms_annotations.append(cds_genes_intergenics)
            all_chroms_annotations = pd.concat(all_chroms_annotations).reset_index(drop=True)
            all_chroms_annotations.to_csv(csv_file)
            return all_chroms_annotations
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")


def assign_gene(start: int, end: int, gene_elements_df: pd.DataFrame) -> (int, int):
    gene_ids = gene_elements_df.index[
        (gene_elements_df['start'] - 1 <= start) & (end <= gene_elements_df['end'])].tolist()
    # gene_elements_df['start'] - 1 correction is because of some bug in pybedtools

    if len(gene_ids) < 1:
        return 0, Strand.UNKNOWN.value
    if len(gene_ids) > 1:
        # We assume genes from different strands don't overlap, but if it happens nonetheless - ignore
        return 0, Strand.UNKNOWN.value
    strand = Strand(gene_elements_df.loc[gene_ids[0], 'strand'])
    return gene_ids[0] + 1, strand.value  # indices in pandas are 0-based


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

    def __init__(self, annotations_df: pd.DataFrame, chromosome: str, is_sorted: bool = False):
        self._chromosome = chromosome
        if annotations_df is None or annotations_df.empty:
            self._relevant_regions.append(Region.ALL)
            return

        self._annotations_df = annotations_df.copy()
        if 'chrom' in annotations_df.columns:
            self._annotations_df = self._annotations_df[annotations_df['chrom'] == chromosome]
        if not is_sorted:
            self._annotations_df.sort_values('start', inplace=True)
        self._annotations_df = self._annotations_df[['start', 'end', 'region', 'gene', 'strand']]
        self._annotations_df.reset_index(drop=True, inplace=True)

        relevant_region_names = self._annotations_df['region'].unique()
        self._relevant_regions = [Region(name) for name in relevant_region_names]

        if len(self._relevant_regions) == 0:
            self._relevant_regions.append(Region.ALL)

        self._update_reading_frames()
        self._discard_genes_without_strand()

    def len(self):
        if self._annotations_df is None or self._annotations_df.empty:
            return None
        return self._annotations_df.iloc[-1]['end'].item()

    def _get_annotation_by_position(self, pos) -> (pd.Series, int):
        condition = (self._annotations_df['start'] <= pos) & (pos < self._annotations_df['end'])
        annotation_indices = self._annotations_df.index[condition].tolist()
        if len(annotation_indices) != 1:
            raise Exception(f"Can't determine annotation for chromosome {self._chromosome} position {pos}")
        return self._annotations_df.loc[annotation_indices[0]], annotation_indices[0]

    def get_regions(self) -> List[Region]:
        return self._relevant_regions

    def get_region_by_position(self, pos) -> (Region, Strand):
        if self._annotations_df is None or self._annotations_df.empty:
            return Region.ALL
        annotation, _ = self._get_annotation_by_position(pos)
        return Region(annotation['region']), Strand(annotation['strand'])

    def get_annotation_start_end(self, pos) -> (int, int):
        annotation, _ = self._get_annotation_by_position(pos)

        return annotation['start'], annotation['end']

    def get_last_coding_position_of_encapsulating_gene(self, pos: int) -> (Optional[int], Strand):
        annotation, annotation_index = self._get_annotation_by_position(pos)

        gene = annotation['gene']
        if not gene or gene == 0:
            raise Exception("Coding positions are only relevant for genes")

        strand = Strand(annotation.iloc[0]['strand'])
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

    def get_encapsulating_codon_positions(self, pos: int) -> (int, int, int):
        cds, _ = self._get_annotation_by_position(pos)
        if cds['region'] != Region.CDS.value:
            raise Exception("Codon is only relevant for CDS region")
        cds_start = cds['start']
        cds_end = cds['end']
        reading_offset = cds['reading_offset'] + (pos - cds_start)
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
                first = prev_cds_end - 2
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

    def get_annotations_in_range(self, start: int, end: int) -> pd.DataFrame:
        annotations = self._annotations_df[(start < self._annotations_df['end']) &
                                           (self._annotations_df['start'] < end)]
        return annotations

    def mute_gene(self, position_on_gene: int = -1, gene_id: int = 0) -> None:
        if 0 <= position_on_gene < self.len():
            annotation, annotation_index = self._get_annotation_by_position(position_on_gene)
            gene_id = annotation['gene']
        if gene_id == 0:
            return
        gene_annotations_indices = self._annotations_df.index[(self._annotations_df['gene'] == gene_id)].tolist()
        first_index = gene_annotations_indices[0]
        last_index = gene_annotations_indices[-1]

        # melt with previous intergenic region if exists
        if first_index != 0 and \
                self._annotations_df.loc[first_index - 1, 'region'] == Region.INTERGENIC.value:
            first_index -= 1

        # melt with next intergenic region if exists
        if last_index + 1 != len(self._annotations_df) and \
                self._annotations_df.loc[last_index + 1, 'region'] == Region.INTERGENIC.value:
            last_index += 1

        # create new intergenic region
        intergenic_start = self._annotations_df.loc[first_index, 'start'].item()
        intergenic_end = self._annotations_df.loc[last_index, 'end'].item()
        new_intergenic = pd.DataFrame({'start': intergenic_start,
                                       'end': intergenic_end,
                                       'region': Region.INTERGENIC.value,
                                       'gene': 0,
                                       'strand': Strand.UNKNOWN.value},
                                      index=[first_index])

        # insert new intergenic region instead of muted gene
        dfs_to_concat = []
        if first_index != 0:
            dfs_to_concat.append(self._annotations_df.loc[:first_index - 1])
        dfs_to_concat.append(new_intergenic)
        if last_index + 1 != len(self._annotations_df):
            dfs_to_concat.append(self._annotations_df.loc[last_index + 1:])
        self._annotations_df = pd.concat(dfs_to_concat).reset_index(drop=True)

    def handle_insertion(self, pos: int, insertion_len: int) -> None:
        _, index = self._get_annotation_by_position(pos)
        self._annotations_df.loc[index, 'end'] += insertion_len
        if index + 1 != len(self._annotations_df):
            self._annotations_df.loc[index+1:, 'start'] += insertion_len
            self._annotations_df.loc[index+1:, 'end'] += insertion_len

    def handle_deletion(self, pos: int, deletion_len: int) -> None:
        """
        Delete right after pos, sequence at the length deletion_len
        """
        annotation, index = self._get_annotation_by_position(pos)
        annotation_residue = self._annotations_df.loc[index, 'end'].item() - pos - 1
        deleted_already = min(annotation_residue, deletion_len)
        annotation['end'] -= deleted_already
        index += 1

        annotations_to_delete = []
        # deletion involves more than one annotation
        while deleted_already < deletion_len and index < len(self._annotations_df):
            annotation_len = self._annotations_df.loc[index, 'end'].item() - \
                             self._annotations_df.loc[index, 'start'].item()
            if annotation_len <= deletion_len - deleted_already:
                annotations_to_delete.append(index)
                deleted_already += annotation_len
            else:
                self._annotations_df.loc[index, 'start'] = pos + 1
                self._annotations_df.loc[index, 'end'] = (pos + annotation_len) - (deletion_len - deleted_already)
                deleted_already = deletion_len
            index += 1

        if index != len(self._annotations_df):
            self._annotations_df.loc[index:, 'start'] -= deletion_len
            self._annotations_df.loc[index:, 'end'] -= deletion_len

        if len(annotations_to_delete):
            self._annotations_df = self._annotations_df.drop(annotations_to_delete).reset_index(drop=True)

    def get_nucleotides_counts_per_region(self, start=-1, end=-1) -> dict:
        start = start if start != -1 else 0
        end = end if end != -1 else self.len()
        relevant_annotations = self.get_annotations_in_range(start, end)
        counts_per_region = {}
        for _, annotation in relevant_annotations.iterrows():
            region = Region(annotation['region'])
            if region not in counts_per_region:
                counts_per_region[region] = 0
            region_start = max(start, annotation['start'])
            region_end = min(end, annotation['end'])
            counts_per_region[region] += region_end - region_start
        return counts_per_region

    def get_mask_in_window_of_region(self, relevant_region: Region, start: int, end: int) -> dict:
        # if cache is valid and contains the region - return it
        if start == self._cached_start and end == self._cached_end and \
                relevant_region.value in self._cached_mask_in_window_per_region:
            return self._cached_mask_in_window_per_region[relevant_region.value]

        # else - initialize cache
        self._cached_start = start
        self._cached_end = end
        self._cached_mask_in_window_per_region = {region: [] for region in self.get_regions()}

        # compute if no cache
        relevant_annotations = self.get_annotations_in_range(start, end)
        for _, annotation in relevant_annotations.iterrows():
            annotation_start = max(start, annotation['start'])
            annotation_end = min(end, annotation['end'])
            annotation_length = annotation_end - annotation_start
            for region in self._cached_mask_in_window_per_region.keys():
                mask = 1 if region == relevant_region.value else 0
                if len(self._cached_mask_in_window_per_region[region]) > 0:
                    self._cached_mask_in_window_per_region[region] = np.concatenate(
                        [self._cached_mask_in_window_per_region[region], np.full(annotation_length, mask)])
                else:
                    self._cached_mask_in_window_per_region[region] = np.full(annotation_length, mask)

        return self._cached_mask_in_window_per_region[relevant_region.value]

    def _update_reading_frames(self):
        if self._annotations_df is None or self._annotations_df.empty:
            return

        self._annotations_df['reading_offset'] = np.nan
        gene_ids = self._annotations_df.gene.unique()
        for gene_id in gene_ids:
            if gene_id == 0:
                continue  # intergenic region

            cds_annotations = self._annotations_df[(self._annotations_df['gene'] == gene_id) &
                                                   (self._annotations_df['region'] == Region.CDS.value)]
            cds_total_len = cds_annotations['end'].sum() - cds_annotations['start'].sum()

            if cds_total_len % 3 == 0:
                reading_offset = 0
                for idx, annotation in cds_annotations.iterrows():
                    self._annotations_df.loc[idx, 'reading_offset'] = reading_offset
                    reading_offset += annotation['end'] - annotation['start']
                    reading_offset = reading_offset % 3

            else:
                # TODO raise Exception?
                print(f'Total length of CDS elements in gene {gene_id} in chromosome {self._chromosome} '
                      f'cannot be divided by 3')
                self.mute_gene(gene_id=gene_id)

    def _discard_genes_without_strand(self):
        gene_ids = self._annotations_df.gene.unique()
        for gene_id in gene_ids:
            if gene_id == 0:
                continue  # intergenic region

            gene_annotations = self._annotations_df[(self._annotations_df['gene'] == gene_id)].reset_index(drop=True)
            if gene_annotations is None or gene_annotations.empty:
                continue
            strand = Strand(gene_annotations.iloc[0]['strand'])
            if strand == Strand.UNKNOWN:
                self.mute_gene(gene_id=gene_id)


if __name__ == "__main__":
    filename = '/groups/itay_mayrose/adisivan/arabidopsis/ensemblgenomes/gff3/Arabidopsis_thaliana.TAIR10.53_1000.gff3'
    dirname = '/groups/itay_mayrose/adisivan/PanGenomeSimulator/test_arabidopsis'
    separate_cds_genes_intergenics(filename, dirname)

    # ID=gene:AT1G01100;Name=RPP1A;biotype=protein_coding;description=60S acidic ribosomal protein P1-1
    # [Source:UniProtKB/Swiss-Prot%3BAcc:Q8LCW9];gene_id=AT1G01100;logic_name=araport11
