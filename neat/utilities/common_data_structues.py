import numpy as np
import pandas as pd
from enum import Enum


class MutType(Enum):
    SNP = 'snp'
    INDEL = 'indel'
    SV = 'SV'

class Region(Enum):
    CDS = 'CDS'
    INTRON = 'intron'
    INTERGENIC = 'intergenic'
    ALL = 'all'
    #TODO see how we pass through evey Region which is not ALL, or only through ALL, when inserting mutations.
    # need some kind of "strategy" solution


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
