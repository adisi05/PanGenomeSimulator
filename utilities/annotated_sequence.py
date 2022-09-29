import time
from typing import Optional, List

import numpy as np
import pandas as pd

from common_data_structues import Region, Strand


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

    def __init__(self, annotations_df: pd.DataFrame, chromosome: str, is_sorted: bool = False, debug: bool = False):
        self.debug = debug
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

        strand = Strand(annotation['strand'])
        if strand.value == Strand.UNKNOWN.value:
            return None, strand

        gene_cds_annotations = self._annotations_df[(self._annotations_df['gene'] == gene) &
                                                    (self._annotations_df['region'] == Region.CDS.value)]

        if strand.value == Strand.FORWARD.value:
            last_cds = gene_cds_annotations.iloc[-1]
            return last_cds['end'].item() - 1, strand
        else:
            last_cds = gene_cds_annotations.iloc[0]
            return last_cds['start'].item(), strand

    def get_encapsulating_codon_positions(self, pos: int) -> (int, int, int):
        cds, _ = self._get_annotation_by_position(pos)
        if cds['region'] != Region.CDS.value:
            raise Exception("Codon is only relevant for CDS region")
        cds_start = cds['start']
        cds_end = cds['end']
        gene = cds['gene']
        reading_offset = int(cds['reading_offset']) + (pos - cds_start)
        reading_offset = reading_offset % 3
        first = pos - reading_offset
        second = pos + 1 - reading_offset
        third = pos + 2 - reading_offset

        if first < cds_start or second < cds_start:
            prev_cds = self._annotations_df[(self._annotations_df['gene'] == gene) &
                                            (self._annotations_df['region'] == Region.CDS.value) &
                                            (self._annotations_df['end'] <= cds_start)]
            prev_cds_end = prev_cds.iloc[-1]['end'].item()
            if second < cds_start:
                second = prev_cds_end - 1
                first = prev_cds_end - 2
            elif first < cds_start:
                first = prev_cds_end - 1

        if cds_end <= second or cds_end <= third:
            next_cds = self._annotations_df[(self._annotations_df['gene'] == gene) &
                                            (self._annotations_df['region'] == Region.CDS.value) &
                                            (cds_end <= self._annotations_df['start'])]
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

        start = time.time()
        if self.debug:
            print(f"Trying to mute gene {gene_id}")

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

        if self.debug:
            end = time.time()
            print(f"Gene muting took {int(end - start)} seconds.")

    def handle_insertion(self, pos: int, insertion_len: int) -> None:
        _, index = self._get_annotation_by_position(pos)
        self._annotations_df.loc[index, 'end'] += insertion_len
        if index + 1 != len(self._annotations_df):
            self._annotations_df.loc[index + 1:, 'start'] += insertion_len
            self._annotations_df.loc[index + 1:, 'end'] += insertion_len

    def handle_deletion(self, pos: int, deletion_len: int) -> None:
        """
        Delete right after pos, sequence at the length deletion_len
        """
        _, index = self._get_annotation_by_position(pos)
        annotation_residue = self._annotations_df.loc[index, 'end'].item() - pos - 1
        deleted_already = min(annotation_residue, deletion_len)
        self._annotations_df.loc[index, 'end'] -= deleted_already
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

    def get_nucleotides_counts_per_region(self, start: int = -1, end: int = -1) -> dict:
        start = start if start != -1 else 0
        end = end if end != -1 else self.len()
        relevant_annotations = self.get_annotations_in_range(start, end)
        counts_per_region = {region.value: 0 for region in self.get_regions()}
        for _, annotation in relevant_annotations.iterrows():
            region = Region(annotation['region'])
            region_start = max(start, annotation['start'])
            region_end = min(end, annotation['end'])
            counts_per_region[region.value] += region_end - region_start
        return counts_per_region

    def get_mask_in_window_of_region(self, relevant_region: Region, start: int, end: int) -> dict:
        # if cache is valid and contains the region - return it
        if start == self._cached_start and end == self._cached_end and \
                relevant_region.value in self._cached_mask_in_window_per_region:
            return self._cached_mask_in_window_per_region[relevant_region.value]

        # else - initialize cache
        self._cached_start = start
        self._cached_end = end
        self._cached_mask_in_window_per_region = {}

        # compute if no cache
        self.compute_mask_in_window(start, end)

        return self._cached_mask_in_window_per_region[relevant_region.value]

    def compute_mask_in_window(self, start: int, end: int):
        relevant_annotations = self.get_annotations_in_range(start, end)
        for _, annotation in relevant_annotations.iterrows():
            annotation_start = max(start, annotation['start'])
            annotation_end = min(end, annotation['end'])
            annotation_length = annotation_end - annotation_start
            annotation_region = Region(annotation['region'])
            for region in self.get_regions():
                mask = 1 if region.value == annotation_region.value else 0
                if region.value in self._cached_mask_in_window_per_region:
                    self._cached_mask_in_window_per_region[region.value] = np.concatenate(
                        [self._cached_mask_in_window_per_region[region.value],
                         np.full(annotation_length, mask)])
                else:
                    self._cached_mask_in_window_per_region[region.value] = np.full(annotation_length, mask)

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
            if strand.value == Strand.UNKNOWN.value:
                self.mute_gene(gene_id=gene_id)