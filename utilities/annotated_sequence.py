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
        self._genes_presence_absence_dict = {}

        if annotations_df is None or annotations_df.empty:
            self._relevant_regions.append(Region.ALL)
            return

        self._annotations_df = annotations_df.copy()
        if 'chrom' in annotations_df.columns:

            # if chromosome name is not the exactly the name in the dataframe, try to find it
            if self._chromosome not in self._annotations_df.chrom.unique():
                split = self._chromosome.split(' ')
                for item in split:
                    if item in self._annotations_df.chrom.unique():
                        self._chromosome = item
                        print(f"Chromosome {chromosome} was not found in the annotation dataframe. "
                              f"Instead, found {self._chromosome}, so using it.")
                        break

            self._annotations_df = self._annotations_df[annotations_df['chrom'] == self._chromosome]
            if self._annotations_df.empty:
                print(f"Chromosome {self._chromosome} was not found in the annotation dataframe."
                      f"Therefore, ignoring the annotations for this chromosome.")
            del self._annotations_df['chrom']
        if not is_sorted:
            self._annotations_df.sort_values('start', inplace=True)
        self._annotations_df.reset_index(drop=True, inplace=True)

        relevant_region_names = self._annotations_df['region'].unique()
        self._relevant_regions = [Region(name) for name in relevant_region_names]
        if len(self._relevant_regions) == 0:
            self._relevant_regions.append(Region.ALL)

        relevant_genes = self._annotations_df.gene_id.unique()
        self._genes_presence_absence_dict = {gene_id: True for gene_id in relevant_genes}

    def len(self):
        if self._annotations_df is None or self._annotations_df.empty:
            return None
        return self._annotations_df.iloc[-1]['end']

    def get_annotations_df(self) -> pd.DataFrame:
        annotations_df = self._annotations_df.copy()
        if 'chrom' not in annotations_df.columns:
            annotations_df['chrom'] = self._chromosome
        return annotations_df

    def get_genes_presence_absence_dict(self) -> dict:
        return self._genes_presence_absence_dict

    def get_regions(self) -> List[Region]:
        return self._relevant_regions

    def get_region_by_position(self, pos) -> (Region, Strand):
        if self._annotations_df is None or self._annotations_df.empty:
            return Region.ALL
        annotation, _ = self._get_annotation_by_position(pos)
        return Region(annotation['region']), Strand(annotation['strand'])

    def get_annotation_start_end(self, pos) -> (int, int):
        if self._annotations_df is None or self._annotations_df.empty:
            return None, None

        annotation, _ = self._get_annotation_by_position(pos)

        return annotation['start'], annotation['end']

    def get_last_coding_position_of_encapsulating_gene(self, pos: int) -> (Optional[int], Strand):
        annotation, annotation_index = self._get_annotation_by_position(pos)

        gene = annotation['gene_id']
        if not gene or gene == 0:
            raise Exception("Coding positions are only relevant for genes")

        strand = Strand(annotation['strand'])
        if strand.value == Strand.UNKNOWN.value:
            return None, strand

        gene_cds_annotations = self._annotations_df[(self._annotations_df['gene_id'] == gene) &
                                                    (self._annotations_df['region'] == Region.CDS.value)]

        if strand.value == Strand.FORWARD.value:
            last_cds = gene_cds_annotations.iloc[-1]
            return last_cds['end'] - 1, strand
        else:
            last_cds = gene_cds_annotations.iloc[0]
            return last_cds['start'], strand

    def get_encapsulating_codon_positions(self, pos: int) -> (int, int, int):
        if self.debug:
            print(f"Trying to find codon positions for nucleotide at position {pos}.")
        cds, _ = self._get_annotation_by_position(pos)
        if cds['region'] != Region.CDS.value:
            raise Exception("Codon is only relevant for CDS region")
        cds_start = cds['start']
        cds_end = cds['end']
        cds_gene = cds['gene_id']
        cds_strand = cds['strand']
        if self.debug:
            print(f"CDS details: start={cds_start}, CDS end={cds_end}, gene={cds_gene}, strand={cds_strand}.")

        # find positions of 3 nucleotides
        if cds_strand == Strand.FORWARD.value:
            pos_frame = (int(cds['frame']) + (pos - cds_start)) % 3
            first = pos - pos_frame
            second = pos + 1 - pos_frame
            third = pos + 2 - pos_frame
        elif cds_strand == Strand.REVERSE.value:
            # initially, find codons in reversed read
            reversed_cds_start = cds_end - 1
            pos_frame = (int(cds['frame']) + (reversed_cds_start - pos)) % 3
            reversed_first = pos + pos_frame
            reversed_second = pos - 1 + pos_frame
            reversed_third = pos - 2 + pos_frame
            # then convert to forward read
            first = reversed_third
            second = reversed_second
            third = reversed_first
        else:
            raise Exception("Unknown strand for CDS element")
        if self.debug:
            print(f"Found this positions: {first}-{third}.")

        # get previous CDS (which is actually next if reverse strand...) if needed
        if first < cds_start or second < cds_start:
            if self.debug:
                print(f"Getting previous CDS because current CDS start is {cds_start}.")
            prev_cds = self._annotations_df[(self._annotations_df['gene_id'] == cds_gene) &
                                            (self._annotations_df['region'] == Region.CDS.value) &
                                            (self._annotations_df['end'] <= cds_start)]
            if len(prev_cds) == 0 and self.debug:
                print(f"No previous CDS element is available.")
            prev_cds_end = prev_cds.iloc[-1]['end']
            if second < cds_start:
                second = prev_cds_end - 1
                first = prev_cds_end - 2
            elif first < cds_start:
                first = prev_cds_end - 1

        # get next CDS (which is actually previous if reverse strand...) if needed
        if cds_end <= second or cds_end <= third:
            if self.debug:
                print(f"Getting next CDS because current CDS end is {cds_end}.")
            next_cds = self._annotations_df[(self._annotations_df['gene_id'] == cds_gene) &
                                            (self._annotations_df['region'] == Region.CDS.value) &
                                            (cds_end <= self._annotations_df['start'])]
            if len(next_cds) == 0 and self.debug:
                print(f"No next CDS element is available.")
            next_cds_start = next_cds.iloc[0]['start']
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
            gene_id = annotation['gene_id']
        if gene_id == 0:
            return

        start = time.time()
        if self.debug:
            print(f"Trying to mute gene {gene_id}")

        gene_annotations_indices = self._annotations_df.index[(self._annotations_df['gene_id'] == gene_id)].tolist()
        first_index = gene_annotations_indices[0]
        last_index = gene_annotations_indices[-1]

        # melt with previous intergenic region if exists
        if first_index != 0 and \
                self._annotations_df.iloc[first_index - 1]['region'] == Region.INTERGENIC.value:
            first_index -= 1

        # melt with next intergenic region if exists
        if last_index + 1 != len(self._annotations_df) and \
                self._annotations_df.iloc[last_index + 1]['region'] == Region.INTERGENIC.value:
            last_index += 1

        # create new intergenic region
        intergenic_start = self._annotations_df.iloc[first_index]['start']
        intergenic_end = self._annotations_df.iloc[last_index]['end']
        new_intergenic = pd.DataFrame({'start': intergenic_start,
                                       'end': intergenic_end,
                                       'region': Region.INTERGENIC.value,
                                       'gene_id': 0,
                                       'strand': Strand.UNKNOWN.value,
                                       'frame': np.nan},
                                      index=[first_index])

        # insert new intergenic region instead of muted gene
        dfs_to_concat = []
        if first_index != 0:
            dfs_to_concat.append(self._annotations_df.iloc[:first_index, :])
        dfs_to_concat.append(new_intergenic)
        if last_index + 1 != len(self._annotations_df):
            dfs_to_concat.append(self._annotations_df.iloc[last_index + 1:, :])
        self._annotations_df = pd.concat(dfs_to_concat).reset_index(drop=True)

        self._genes_presence_absence_dict[gene_id] = False

        if self.debug:
            end = time.time()
            print(f"Gene muting took {0:.3f} seconds.".format(end-start))
            print(f"New intergenic region is now from start={intergenic_start} to end={intergenic_end}.")
            print(f"Merged the (current) annotation indexes: first_index={first_index} to last_index={last_index}.")
            print(f"The result: lost {last_index-first_index} annotations.")

    def handle_insertion(self, pos: int, insertion_len: int) -> None:
        start = time.time()

        _, index = self._get_annotation_by_position(pos)
        self._annotations_df.iloc[index, self._annotations_df.columns.get_loc('end')] += insertion_len
        if index + 1 != len(self._annotations_df):
            self._annotations_df.iloc[index + 1:, self._annotations_df.columns.get_loc('start')] += insertion_len
            self._annotations_df.iloc[index + 1:, self._annotations_df.columns.get_loc('end')] += insertion_len
            if self.debug:
                print(f"Shifted forward annotations from index {index + 1} by {insertion_len}.")

        if self.debug:
            end = time.time()
            print(f"handle_insertion took {0:.3f} seconds.".format(end-start))

    def handle_deletion(self, pos: int, deletion_len: int) -> None:
        """
        Delete right after pos, sequence at the length deletion_len
        """
        start = time.time()

        _, index = self._get_annotation_by_position(pos)
        annotation_residue = self._annotations_df.iloc[index]['end'] - pos - 1
        deleted_already = min(annotation_residue, deletion_len)
        self._annotations_df.iloc[index, self._annotations_df.columns.get_loc('end')] -= deleted_already
        index += 1

        annotations_to_delete = []
        # deletion involves more than one annotation
        while deleted_already < deletion_len and index < len(self._annotations_df):
            annotation_len = self._annotations_df.iloc[index]['end'] - \
                             self._annotations_df.iloc[index]['start']

            if annotation_len <= (deletion_len - deleted_already):
                annotations_to_delete.append(index)
                deleted_already += annotation_len
            else:
                self._annotations_df.iloc[index, self._annotations_df.columns.get_loc('start')] = pos + 1
                self._annotations_df.iloc[index, self._annotations_df.columns.get_loc('end')] = \
                    (pos + 1 + annotation_len) - (deletion_len - deleted_already)
                deleted_already = deletion_len

            index += 1

        if index != len(self._annotations_df):
            self._annotations_df.iloc[index:, self._annotations_df.columns.get_loc('start')] -= deletion_len
            self._annotations_df.iloc[index:, self._annotations_df.columns.get_loc('end')] -= deletion_len
            if self.debug:
                print(f"Shifted backwards annotations from index {index + 1} by {deletion_len}.")

        if len(annotations_to_delete):
            self._annotations_df = self._annotations_df.drop(annotations_to_delete).reset_index(drop=True)

        if self.debug:
            end = time.time()
            print(f"handle_deletion took {0:.3f} seconds.".format(end-start))

    def get_nucleotides_counts_per_region(self, start: int = -1, end: int = -1) -> dict:
        if self._annotations_df is None or self._annotations_df.empty:
            return {Region.ALL.value: end - start}

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
        self._compute_mask_in_window(start, end)

        return self._cached_mask_in_window_per_region[relevant_region.value]

    def _compute_mask_in_window(self, start: int, end: int):
        if self._annotations_df is None or self._annotations_df.empty:
            self._cached_mask_in_window_per_region = {Region.ALL.value: np.full(end - start, 1)}
            return

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

    def _get_annotation_by_position(self, pos) -> (pd.Series, int):
        condition = (self._annotations_df['start'] <= pos) & (pos < self._annotations_df['end'])
        annotation_indices = self._annotations_df.index[condition].tolist()
        if len(annotation_indices) != 1:
            raise Exception(f"Can't determine annotation for chromosome {self._chromosome} position {pos}.\n"
                            f"Number of proposed annotations is {len(annotation_indices)}")
        return self._annotations_df.iloc[annotation_indices[0]], annotation_indices[0]
