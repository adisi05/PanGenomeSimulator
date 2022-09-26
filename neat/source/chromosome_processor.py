import copy
import pickle
import random
import sys
from dataclasses import dataclass
from typing import Optional, List, Union, Tuple

import numpy as np
from Bio import Seq
from Bio.Seq import MutableSeq

from neat.utilities.annotated_sequence import AnnotatedSequence
from neat.utilities.common_data_structues import Region, MutType, Strand
from neat.utilities.gen_mut_model import Stats
from probability import DiscreteDistribution, poisson_list

import pandas as pd


"""
Constants needed for analysis
"""
MAX_ATTEMPTS = 100  # max attempts to insert a mutation into a valid position
MAX_MUTFRAC = 0.3  # TODO rethink it!!! the maximum percentage of a window/sequence that can contain mutations

NUCL = ['A', 'C', 'G', 'T']
TRI_IND = {'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3, 'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
           'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11, 'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15}
NUC_IND = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
ALL_TRI = [NUCL[i] + NUCL[j] + NUCL[k] for i in range(len(NUCL)) for j in range(len(NUCL)) for k in range(len(NUCL))]
ALL_IND = {ALL_TRI[i]: i for i in range(len(ALL_TRI))}

# DEBUG
IGNORE_TRINUC = False

# percentile resolution used for fraglen quantizing
COV_FRAGLEN_PERCENTILE = 10.
LARGE_NUMBER = 9999999999

"""
DEFAULT MUTATION MODELS
"""

DEFAULT_1_OVERALL_MUT_RATE = 0.001
DEFAULT_1_HOMOZYGOUS_FREQ = 0.010
DEFAULT_1_INDEL_FRACTION = 0.05
DEFAULT_1_INS_VS_DEL = 0.6
DEFAULT_1_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_1_INS_LENGTH_WEIGHTS = [0.4, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.034, 0.033, 0.033]
DEFAULT_1_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5]
DEFAULT_1_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1]
example_matrix_1 = [[0.0, 0.15, 0.7, 0.15],
                    [0.15, 0.0, 0.15, 0.7],
                    [0.7, 0.15, 0.0, 0.15],
                    [0.15, 0.7, 0.15, 0.0]]
DEFAULT_1_TRI_FREQS = [copy.deepcopy(example_matrix_1) for _ in range(16)]
DEFAULT_1_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
DEFAULT_MODEL_1 = [DEFAULT_1_OVERALL_MUT_RATE,
                   DEFAULT_1_HOMOZYGOUS_FREQ,
                   DEFAULT_1_INDEL_FRACTION,
                   DEFAULT_1_INS_VS_DEL,
                   DEFAULT_1_INS_LENGTH_VALUES,
                   DEFAULT_1_INS_LENGTH_WEIGHTS,
                   DEFAULT_1_DEL_LENGTH_VALUES,
                   DEFAULT_1_DEL_LENGTH_WEIGHTS,
                   DEFAULT_1_TRI_FREQS,
                   DEFAULT_1_TRINUC_BIAS]

DEFAULT_2_OVERALL_MUT_RATE = 0.002
DEFAULT_2_HOMOZYGOUS_FREQ = 0.200
DEFAULT_2_INDEL_FRACTION = 0.1
DEFAULT_2_INS_VS_DEL = 0.3
DEFAULT_2_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_2_INS_LENGTH_WEIGHTS = [0.1, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
# noinspection DuplicatedCode
DEFAULT_2_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5]
DEFAULT_2_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1]
example_matrix_2 = [[0.0, 0.15, 0.7, 0.15],
                    [0.15, 0.0, 0.15, 0.7],
                    [0.7, 0.15, 0.0, 0.15],
                    [0.15, 0.7, 0.15, 0.0]]
DEFAULT_2_TRI_FREQS = [copy.deepcopy(example_matrix_2) for _ in range(16)]
DEFAULT_2_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
DEFAULT_MODEL_2 = [DEFAULT_2_OVERALL_MUT_RATE,
                   DEFAULT_2_HOMOZYGOUS_FREQ,
                   DEFAULT_2_INDEL_FRACTION,
                   DEFAULT_2_INS_VS_DEL,
                   DEFAULT_2_INS_LENGTH_VALUES,
                   DEFAULT_2_INS_LENGTH_WEIGHTS,
                   DEFAULT_2_DEL_LENGTH_VALUES,
                   DEFAULT_2_DEL_LENGTH_WEIGHTS,
                   DEFAULT_2_TRI_FREQS,
                   DEFAULT_2_TRINUC_BIAS]


# see https://stackoverflow.com/questions/68840767/what-is-the-purpose-of-sort-index-of-a-dataclass-in-python
@dataclass(order=True)
@dataclass
class Mutation:
    position: int  # = field(init=False, repr=False)
    ref_nucl: str  # = field(init=False, repr=False)
    new_nucl: str
    mut_type: MutType

    def get_offset_change(self):
        return len(self.new_nucl) - len(self.ref_nucl)


class WindowUnit:
    start_offset: int = 0
    start: int = 0
    end: int = 0
    original_end: int = 0
    # Blocklist explanation:
    # blocklist[pos] = 0		safe to insert variant here
    # blocklist[pos] = 1		indel inserted here
    # blocklist[pos] = 2		snp inserted here
    # blocklist[pos] = 3		invalid position for various processing reasons
    blocklist: dict  # TODO re-consider !!!

    def initialize(self, start=-1, end=-1, final_start_end=False):
        # TODO sanity check about start, end. Maybe consider N regions?
        self.original_end = end
        self.blocklist = {}  # np.zeros(end-start, dtype='<i4')
        return self.update(start=start, end=end, final_start_end=final_start_end)

    def update(self, start=-1, end=-1, final_start_end=False, end_shift=0):
        # TODO sanity check about start, end. Maybe consider N regions?
        changed = False
        if start != -1 and end != -1:
            if not final_start_end:
                start += self.start_offset
                end += self.start_offset
            if self.start != start or self.end != end:
                self.start = start
                self.end = end
                changed = True
        elif end_shift != 0:
            self.end += end_shift
            changed = True
        return changed

    def finalize(self):
        self.start_offset += (self.end - self.original_end)
        # TODO add the feature to calculate the "bite" from the next window to come?

    def update_blocklist(self, mutation):
        for k in range(mutation.position, mutation.position + len(mutation.ref_nucl)):  # TODO continue here
            self.blocklist[k] = 1 if mutation.mut_type.value == MutType.INDEL.value else 2

    def check_blocklist(self, mutation):
        for k in range(mutation.position, mutation.position + len(mutation.ref_nucl)):
            if self.blocklist.get(k, 0):
                return False
        return True


class RandomMutationPool:
    def __init__(self, indels_per_region: dict, snps_per_region: dict, max_mutations_in_window: int,
                 sv_list: list = [], debug: bool = False):
        # TODO add SVs
        self.indels_per_region = indels_per_region
        self.snps_per_region = snps_per_region
        self.options = {}
        for region_name, count in indels_per_region.items():
            if count != 0:
                self.options[(MutType.INDEL.value, region_name)] = count
        for region_name, count in snps_per_region.items():
            if count != 0:
                self.options[(MutType.SNP.value, region_name)] = count
        self.overall_count = min(max_mutations_in_window, sum(self.options.values()))
        self.debug = debug

        if self.debug:
            print(f"Created random mutations pool")
            print(f"Overall planned random mutations within window: {self.overall_count}")
            print(f"Mutations distribution: {list(self.options.items())}")


    def has_next(self) -> bool:
        return self.overall_count > 0

    def get_next(self) -> (MutType, Region):
        if self.overall_count <= 0:
            return None
        option_with_count_list = self.options.items()
        # https://pynative.com/python-weighted-random-choices-with-probability/
        choice = random.choices([opt[0] for opt in option_with_count_list],
                                weights=[opt[1] for opt in option_with_count_list],
                                k=1)[0]
        self.overall_count -= 1
        self.options[choice] -= 1
        if self.options[choice] == 0:
            del self.options[choice]
        mut_type_name, region_name = choice

        if self.debug:
            print(f"Remained mutations count within window: {self.overall_count}")

        return MutType(mut_type_name), Region(region_name)


class ChromosomeProcessor:
    """
    Container for reference sequences, applies mutations
    """

    def __init__(self, chromosome_name: str, chromosome_sequence: Seq, annotations_df: pd.DataFrame,
                 annotations_sorted: bool = False, mut_models=None, mut_rate=None, dist=None, debug: bool = False):
        self.chromosome_name = chromosome_name
        self.chromosome_sequence = MutableSeq(str(chromosome_sequence))
        # TODO consider using Seq class to benefit from the class-supported methods
        self.seq_len = len(chromosome_sequence)
        self.annotated_seq = AnnotatedSequence(annotations_df, chromosome_name, is_sorted=annotations_sorted)
        self._update_mut_models(mut_models, mut_rate, dist)
        self.window_unit = WindowUnit()
        self.debug = debug

    def _update_mut_models(self, model_data: dict, mut_rate, dist):
        # TODO figure out: called one time? or a few times, for each window?
        if not model_data:
            self.model_data = {region.value: copy.deepcopy(DEFAULT_MODEL_1) for region in
                               self.annotated_seq.get_regions()}
        else:
            self.model_data = copy.deepcopy(model_data)

        # do we need to rescale mutation frequencies?
        self.mut_rescale = mut_rate
        if self.mut_rescale is None:
            self.mut_scalar_per_region = {region.value: 1.0 for region in self.annotated_seq.get_regions()}
        else:
            self.mut_scalar_per_region = \
                {region.value: float(self.mut_rescale) // (self.model_data[region.value][0] /
                                                           float(len(self.model_data[region.value])))
                 for region in self.annotated_seq.get_regions()}
        if dist:
            self.mut_scalar_per_region = {region.value: self.mut_scalar_per_region[region.value] * dist for region in
                                          self.annotated_seq.get_regions()}

        # init mutation models
        #
        # self.models_per_region[region.value][0] = average mutation rate
        # self.models_per_region[region.value][1] = p(mut is homozygous | mutation occurs)
        # self.models_per_region[region.value][2] = p(mut is indel | mut occurs)
        # self.models_per_region[region.value][3] = p(insertion | indel occurs)
        # self.models_per_region[region.value][4] = distribution of insertion lengths
        # self.models_per_region[region.value][5] = distribution of deletion lengths
        # self.models_per_region[region.value][6] = distribution of trinucleotide SNP transitions
        # self.models_per_region[region.value][7] = p(trinuc mutates)
        self.model_per_region = {}
        for region in self.annotated_seq.get_regions():
            data = self.model_data[region.value]
            self.model_per_region[region.value] = \
                [self.mut_scalar_per_region[region.value] * data[0], data[1], data[2], data[3],
                 DiscreteDistribution(data[5], data[4]), DiscreteDistribution(data[7], data[6]), []]
            for m in data[8]:
                # noinspection PyTypeChecker
                self.model_per_region[region.value][6].append(
                    [DiscreteDistribution(m[0], NUCL), DiscreteDistribution(m[1], NUCL),
                     DiscreteDistribution(m[2], NUCL), DiscreteDistribution(m[3], NUCL)])
            self.model_per_region[region.value].append([m for m in data[9]])

    def get_window_mutations(self) -> RandomMutationPool:
        # NOTE: window can be also a whole non-N region or the entire chromosome
        indels_to_add_window_per_region, snps_to_add_window_per_region = \
            self.get_planned_snps_and_indels_in_window_per_region()
        max_mutations_in_window = round(
            MAX_MUTFRAC * (self.window_unit.end - self.window_unit.start))  # TODO rethink it
        return RandomMutationPool(indels_to_add_window_per_region, snps_to_add_window_per_region,
                                  max_mutations_in_window, debug=self.debug)

    def get_planned_snps_and_indels_in_window_per_region(self) -> (dict, dict):
        indel_poisson_per_region = self.init_poisson(type_is_indel=True)
        snp_poisson_per_region = self.init_poisson(type_is_indel=False)
        indels_to_add_window_per_region = {region.value: indel_poisson_per_region[region.value].sample()
                                           for region in self.annotated_seq.get_regions()}
        snps_to_add_window_per_region = {region.value: snp_poisson_per_region[region.value].sample()
                                         for region in self.annotated_seq.get_regions()}
        return indels_to_add_window_per_region, snps_to_add_window_per_region

    def _update_trinuc_bias_of_window(self):
        # initialize/update trinuc snp bias
        # compute mutation positional bias given trinucleotide strings of the sequence (ONLY AFFECTS SNPs)
        #
        # note: since indels are added before snps, it's possible these positional biases aren't correctly utilized
        #       at positions affected by indels. At the moment I'm going to consider this negligible.
        window_seq_len = self.window_unit.end - self.window_unit.start
        trinuc_snp_bias_of_window_per_region = {region.value: [0. for _ in range(window_seq_len)] for region in
                                                self.annotated_seq.get_regions()}
        self.trinuc_bias_per_region = {region.value: None for region in self.annotated_seq.get_regions()}
        for region in self.annotated_seq.get_regions():  # TODO switch order between loops?
            region_mask = self.annotated_seq.get_mask_in_window_of_region(region, self.window_unit.start,
                                                                          self.window_unit.end)
            for i in range(0 + 1, window_seq_len - 1):
                codon = str(self.chromosome_sequence[self.window_unit.start + i - 1:self.window_unit.start + i + 2])
                trinuc_snp_bias_of_window_per_region[region.value][i] = \
                    region_mask[i] * self.model_per_region[region.value][7][ALL_IND[codon]]
            self.trinuc_bias_per_region[region.value] = \
                DiscreteDistribution(trinuc_snp_bias_of_window_per_region[region.value][0 + 1:window_seq_len - 1],
                                     range(0 + 1, window_seq_len - 1))

    def init_poisson(self, type_is_indel: bool = True):
        list_per_region = {}
        poisson_per_region = {}
        nucleotides_counts_per_region = self.annotated_seq.get_nucleotides_counts_per_region(self.window_unit.start,
                                                                                             self.window_unit.end)
        for region in self.annotated_seq.get_regions():
            param = self.model_per_region[region.value][2] if type_is_indel else (1. - self.model_per_region[region.value][2])
            list_per_region[region.value] = \
                nucleotides_counts_per_region[region.value] * param * self.model_per_region[region.value][2]
            k_range = range(int(nucleotides_counts_per_region[region.value] * MAX_MUTFRAC))
            poisson_per_region[region.value] = poisson_list(k_range, list_per_region[region.value])
            # TODO validate this. How does this distribution work? should we really multiply by MAX_MUTFRAC?
        return poisson_per_region

    def insert_given_mutations(self, vars_in_current_window, start: int = -1, end: int = -1,
                               final_start_end: bool = False) -> List[Tuple]:  # , use_sequence_offset=True, offset=0):
        self.window_unit.initialize(start=start, end=end, final_start_end=final_start_end)

        mutations_to_insert = self.validate_given_mutations_list(vars_in_current_window)

        # TODO actually insert mutations!!!
        #  consider current_offset and self.window_unit.update(end_shift=mutation.get_offset_change())
        vcf_mutations = self.prepare_mutations_to_vcf(mutations_to_insert, mutations_already_inserted=False)
        self.window_unit.finalize()
        return vcf_mutations

    def validate_given_mutations_list(self, input_df: pd.DataFrame):
        inserted_mutations = []
        for _, row in input_df.iterrows():
            ref_nucl = row['allele']
            new_nucl = row['alternatives'][0]  # take the first alternative
            mut_type = MutType.SNP if len(ref_nucl) == 1 and len(new_nucl) == 1 else MutType.INDEL
            mutation = Mutation(row['pos'] + self.window_unit.start_offset, ref_nucl, new_nucl, mut_type)
            if self.validate_given_mutation(mutation):
                self.window_unit.update_blocklist(mutation)
                inserted_mutations.append(mutation)

        return inserted_mutations

    def validate_given_mutation(self, mutation: Mutation):
        if mutation.position < self.window_unit.start or mutation.position >= self.window_unit.end:
            print('\nError: Attempting to insert variant out of window bounds.')
            print(f'mutation.position={mutation.position}')
            print(f'self.window_unit.start={self.window_unit.start}')
            print(f'self.window_unit.end={self.window_unit.end}')

            sys.exit(1)
        # TODO - use this instead? : ref_len = max([len(input_variable[1]), len(my_alt)])
        if mutation.position + len(mutation.ref_nucl) >= self.window_unit.end:
            # TODO mind that by using self.window_unit.end and not self.seq_len + self.sequence_offset + current_offset
            #  we don't allow deletions to take a "bite" form the next window. Should we aloow that?
            return False
        return self.window_unit.check_blocklist(mutation)

    def generate_random_mutations(self, start: int, end: int) -> List[Tuple]:
        self.window_unit.initialize(start=start, end=end)
        random_mutations_pool = self.get_window_mutations()
        if not IGNORE_TRINUC:
            self._update_trinuc_bias_of_window()
        # TODO consider ? random_snps_minus_inserted = max(self.snps_to_add[i] - len(self.snp_list[i]), 0)
        # TODO consider ? random_indels_minus_inserted = max(self.indels_to_add[i] - len(self.indel_list[i]), 0)

        inserted_mutations = []
        while random_mutations_pool.has_next():
            mut_type, region = random_mutations_pool.get_next()
            # TODO add a check to see if the mutation was really inserted? something like status code?
            inserted_mutation, window_shift = self.insert_random_mutation(mut_type, region)
            if not inserted_mutation:
                continue
            inserted_mutations.append(inserted_mutation)
            # check window
            if window_shift != 0:
                self.window_unit.update(end_shift=window_shift)
            # check annotations
            annotation_changed = self.handle_annotations_after_mutated_sequence(inserted_mutation)
            # if at least one has changed - sample mutations again
            if annotation_changed or window_shift != 0:
                if not IGNORE_TRINUC:
                    self._update_trinuc_bias_of_window()

        vcf_mutations = self.prepare_mutations_to_vcf(inserted_mutations, mutations_already_inserted=True)
        self.window_unit.finalize()
        return vcf_mutations

    def insert_random_mutation(self, mut_type: MutType, region: Region) -> (Optional[Mutation], int):
        window_shift = 0

        position = self.find_position_for_mutation(mut_type, region)
        if position == -1:
            return None, window_shift

        inserted_mutation = None
        if mut_type.value == MutType.SNP.value:
            inserted_mutation = self.insert_snp(position, region)

        elif mut_type.value == MutType.INDEL.value:
            inserted_mutation, window_shift = self.insert_indel(position, region)

        # elif mut_type == MutType.SV:
        #   ...

        return inserted_mutation, window_shift

    def find_position_for_mutation(self, mut_type: MutType, region: Region) -> int:
        # TODO use blocklist?
        region_mask = self.annotated_seq.get_mask_in_window_of_region(region, self.window_unit.start,
                                                                      self.window_unit.end)
        if 1 not in region_mask:
            return -1  # current annotation doesn't exist in window
        for attempt in range(MAX_ATTEMPTS):
            if mut_type.value == MutType.INDEL.value or IGNORE_TRINUC:
                # -2 because we don't allow SNP in the window start/end
                k = self.window_unit.end - self.window_unit.start - 2
                if k < 1:
                    return -1
                event_pos = random.choices(
                    range(self.window_unit.start + 1, self.window_unit.end - 1),
                    weights=region_mask[self.window_unit.start + 1:self.window_unit.end - 1],
                    k=1)[0]
                # https://pynative.com/python-weighted-random-choices-with-probability/
                # TODO if event_pos is ok return it, otherwise keep trying
                return event_pos
            else:
                event_pos = self.trinuc_bias_per_region[region.value].sample()
                # TODO if event_pos is ok return it, otherwise keep trying
                return event_pos
        return -1

    def insert_snp(self, position: int, region: Region):
        snp = self.get_specific_snp(position, region)
        self.mutate_sequence(snp)
        return snp

    def get_specific_snp(self, position: int, region: Region):
        ref_nucl = self.chromosome_sequence[position]
        context = str(self.chromosome_sequence[position - 1]) + str(self.chromosome_sequence[position + 1])
        # sample from tri-nucleotide substitution matrices to get SNP alt allele
        new_nucl = self.model_per_region[region.value][6][TRI_IND[context]][NUC_IND[ref_nucl]].sample()
        snp = Mutation(position, ref_nucl, new_nucl, MutType.SNP)
        # self.blocklist[snp.position] = 2  # TODO use blocklist?
        return snp

    def insert_indel(self, position: int, region: Region):
        indel = self.get_specific_indel(position, region)
        window_shift = self.mutate_sequence(indel)
        return indel, window_shift

    def get_specific_indel(self, position: int, region: Region):
        # insertion
        if random.random() <= self.model_per_region[region.value][3]:
            indel_len = self.model_per_region[region.value][4].sample()
            # sequence content of random insertions is uniformly random (change this later, maybe)
            indel_seq = ''.join([random.choice(NUCL) for _ in range(indel_len)])
            ref_nucl = self.chromosome_sequence[position]
            indel = Mutation(position, ref_nucl, ref_nucl + indel_seq, MutType.INDEL)

        # deletion
        else:
            indel_len = self.model_per_region[region.value][5].sample()

            # skip if deletion too close to boundary
            if position + indel_len + 1 >= self.window_unit.end:
                # TODO mind that by using self.window_unit.end,
                #  and not self.seq_len + self.sequence_offset + current_offset,
                #  we don't allow deletions to take a "bite" form the next window.
                #  - should we allow that?
                indel_len = self.window_unit.end - 2 - position

            # forbid deletion to crossing the boundary of annotation!
            _, annotation_end = self.annotated_seq.get_annotation_start_end(position)
            if position + indel_len > annotation_end:
                indel_len = annotation_end - position
            if indel_len < 1:
                # shoudn't occure
                raise Exception("deletion length is less than 0, program error")

            if indel_len == 1:
                indel_seq = self.chromosome_sequence[position + 1]
            else:
                indel_seq = str(self.chromosome_sequence[position + 1:position + indel_len + 1])
            ref_nucl = self.chromosome_sequence[position]
            indel = Mutation(position, ref_nucl + indel_seq, ref_nucl, MutType.INDEL)

        # TODO use blocklist?
        # for k in range(position, position + indel_len + 1):
        #     self.blocklist[k] = 1

        return indel

    def mutate_sequence(self, mutation: Mutation):
        ref_start = mutation.position
        ref_end = ref_start + len(mutation.ref_nucl)
        window_shift = mutation.get_offset_change()

        if mutation.ref_nucl != str(self.chromosome_sequence[ref_start:ref_end]):
            print('\nError: Something went wrong!\n', mutation, [ref_start, ref_end],
                  str(self.chromosome_sequence[ref_start:ref_end]), '\n')
            sys.exit(1)
        else:
            # alter reference sequence
            self.chromosome_sequence = self.chromosome_sequence[:ref_start] + MutableSeq(mutation.new_nucl) + \
                                       self.chromosome_sequence[ref_end:]

        print(f"Inserted mutation of type {mutation.mut_type.value} at position {mutation.position}. "
              f"Window shift is {window_shift}")
        return window_shift

    def prepare_mutations_to_vcf(self, inserted_mutations: List[Mutation], mutations_already_inserted: bool) -> List[Tuple]:
        inserted_mutations.sort()
        vcf_mutations = []
        mutations_affected_offset = 0
        for mutation in inserted_mutations:
            vcf_position = mutation.position - self.window_unit.start_offset
            if mutations_already_inserted:
                vcf_position -= mutations_affected_offset
                mutations_affected_offset += mutation.get_offset_change()
            vcf_mutations.append(tuple([vcf_position, mutation.ref_nucl, mutation.new_nucl, 1, 'WP=1']))
        # TODO: combine multiple variants that happened to occur at same position into single vcf entry?
        #       reconsider blocklist for that!
        return vcf_mutations

    def handle_annotations_after_mutated_sequence(self, inserted_mutation: Mutation) -> bool:
        if self.annotated_seq.get_regions() == [Region.ALL]:
            return False  # sequence is not annotated

        # SNP
        if inserted_mutation.mut_type.value == MutType.SNP.value:
            return self.handle_annotations_after_snp(inserted_mutation)
        # INDEL - insertion
        elif inserted_mutation.mut_type.value == MutType.INDEL.value \
                and len(inserted_mutation.ref_nucl) < len(inserted_mutation.new_nucl):
            self.handle_annotations_after_small_insertion(inserted_mutation)
            return True
        # INDEL - deletion
        elif inserted_mutation.mut_type.value == MutType.INDEL.value \
                and len(inserted_mutation.ref_nucl) >= len(inserted_mutation.new_nucl):
            self.handle_annotations_after_small_deletion(inserted_mutation)
            return True
        # SVs and other - currently not supported
        else:
            raise Exception("currently supporting only SNPs and INDELs")

    def handle_annotations_after_snp(self, inserted_mutation: Mutation) -> bool:
        """
        Handle annotations after SNP
        :param inserted_mutation:
        :return: True if annotation has changed somehow, False otherwise
        """
        region, strand = self.annotated_seq.get_region_by_position(inserted_mutation.position)

        if region.value == Region.INTERGENIC.value or region.value == Region.NON_CODING_GENE.value:
            return False  # current behaviour is not to check for start/stop codon in these regions.

        elif region.value == Region.CDS.value:
            start, _, end = self.annotated_seq.get_encapsulating_codon_positions(inserted_mutation.position)
            # assuming sequence has been mutated already!
            codon = self.chromosome_sequence[start:end + 1]
            is_stop = is_stop_codon(codon, strand)
            is_last_codon = self.is_last_coding_position(start, end)
            if (is_stop and not is_last_codon) or (not is_stop and is_last_codon):
                self.annotated_seq.mute_gene(position_on_gene=inserted_mutation.position)
                return True
            return False
        else:
            raise Exception("unknown annotation")

    def handle_annotations_after_small_insertion(self, inserted_mutation: Mutation):
        if self.should_mute_gene_after_small_insertion(inserted_mutation):
            self.annotated_seq.mute_gene(position_on_gene=inserted_mutation.position)
        self.annotated_seq.handle_insertion(inserted_mutation.position, len(inserted_mutation.new_nucl) - 1)

    def should_mute_gene_after_small_insertion(self, inserted_mutation: Mutation) -> bool:
        region, strand = self.annotated_seq.get_region_by_position(inserted_mutation.position)

        if region.value == Region.INTERGENIC.value or region.value == Region.NON_CODING_GENE.value:
            return False

        elif region.value == Region.CDS.value:
            frameshift = (len(inserted_mutation.new_nucl) - 1) % 3 != 0
            if frameshift:
                return True
            else:
                first_codon_start, _, first_codon_end = \
                    self.annotated_seq.get_encapsulating_codon_positions(inserted_mutation.position)
                added_codons = int((len(inserted_mutation.new_nucl) - 1) / 3)
                for i in range(added_codons):
                    # assuming sequence has been mutated already!
                    codon = self.chromosome_sequence[first_codon_start + (3 * i): first_codon_end + 1 + (3 * i)]
                    if is_stop_codon(codon, strand):
                        return True
            return False

        else:
            raise Exception("unknown annotation")

    def handle_annotations_after_small_deletion(self, inserted_mutation: Mutation):
        mut_start = inserted_mutation.position
        mut_end = inserted_mutation.position + len(inserted_mutation.ref_nucl) - 1
        involved_annotations = self.annotated_seq.get_annotations_in_range(mut_start, mut_end)
        if len(involved_annotations) > 1:
            raise Exception("currently not supporting large deletions (SVs)")
            # how many genes involved? get names by order
            # if 0 genes - must be intergenic region (should be handled above - 1 notification)
            #    delete - shorten intergenic region and shift downstream genes
            # elif 1 gene
            #     get all cds involved. was there a frameshift? have we gotten stop codon?
            #     mute if needed
            #     shorten cds + non_coding_gene regions as needed
            # elif 2 genes
            #     ...
            # else - 3 and more genes
            #     ...

        if self.should_mute_gene_after_small_deletion(inserted_mutation):
            self.annotated_seq.mute_gene(position_on_gene=inserted_mutation.position)
        self.annotated_seq.handle_deletion(inserted_mutation.position, len(inserted_mutation.new_nucl) - 1)

    def should_mute_gene_after_small_deletion(self, inserted_mutation: Mutation) -> bool:
        region, strand = self.annotated_seq.get_region_by_position(inserted_mutation.position)

        if region.value == Region.INTERGENIC.value or region.value == Region.NON_CODING_GENE.value:
            return False

        elif region.value == Region.CDS.value:
            frameshift = (len(inserted_mutation.ref_nucl) - 1) % 3 != 0
            if frameshift:
                return True
            else:
                start, _, end = self.annotated_seq.get_encapsulating_codon_positions(inserted_mutation.position)
                # assuming sequence has been mutated already!
                codon = self.chromosome_sequence[start:end + 1]
                is_stop = is_stop_codon(codon, strand)
                is_last_codon = self.is_last_coding_position(start, end)
                if (is_stop and not is_last_codon) or (not is_stop and is_last_codon):
                    return True
            return False

        else:
            raise Exception("unknown annotation")

    def is_last_coding_position(self, start: int, end: int):
        last_coding_position, strand = self.annotated_seq.get_last_coding_position_of_encapsulating_gene(end)
        is_last_codon = (last_coding_position == end and strand.value == Strand.FORWARD.value) or \
                        (last_coding_position == start and strand.value == Strand.REVERSE.value)
        return is_last_codon


def is_stop_codon(codon: str, strand: Strand) -> bool:
    if strand.value == Strand.FORWARD.value:
        return codon in ['TAG', 'TAA', 'TGA']

    if strand.value == Strand.REVERSE.value:
        return codon in ['CTA', 'TTA', 'TCA']

    return False


# TODO use self.annotated_seq.get_regions()?
# parse mutation model pickle file
def parse_input_mutation_model(model_file: str = None, which_default: int = 1):
    if which_default == 1:
        mut_model = {region.value: [copy.deepcopy(n) for n in DEFAULT_MODEL_1] for region in Region}
    elif which_default == 2:
        mut_model = {region.value: [copy.deepcopy(n) for n in DEFAULT_MODEL_2] for region in Region}
    else:
        print('\nError: Unknown default mutation model specified\n')
        sys.exit(1)

    if model_file is not None:
        pickle_dict = pickle.load(open(model_file, "rb"))
        region_list = list(Region)
        region_list.append('')
        region_names_found = []
        for region in region_list:
            if pickel_key(region, Stats.AVG_MUT_RATE) not in pickle_dict:
                continue

            region_name = region.value if region != '' else Region.ALL.value
            region_names_found.append(region_name)

            mut_model[region_name][0] = pickle_dict[pickel_key(region, Stats.AVG_MUT_RATE)]
            mut_model[region_name][2] = 1. - pickle_dict[pickel_key(region, Stats.SNP_FREQ)]

            ins_list = pickle_dict[pickel_key(region, Stats.INDEL_FREQ)]
            if len(ins_list):
                ins_count = sum([ins_list[k] for k in ins_list.keys() if k >= 1])
                del_count = sum([ins_list[k] for k in ins_list.keys() if k <= -1])
                ins_vals = [k for k in sorted(ins_list.keys()) if k >= 1]
                ins_weight = [ins_list[k] / float(ins_count) for k in ins_vals]
                del_vals = [k for k in sorted([abs(k) for k in ins_list.keys() if k <= -1])]
                del_weight = [ins_list[-k] / float(del_count) for k in del_vals]
            else:  # degenerate case where no indel stats are provided
                ins_count = 1
                del_count = 1
                ins_vals = [1]
                ins_weight = [1.0]
                del_vals = [1]
                del_weight = [1.0]
            mut_model[region_name][3] = ins_count / float(ins_count + del_count)
            mut_model[region_name][4] = ins_vals
            mut_model[region_name][5] = ins_weight
            mut_model[region_name][6] = del_vals
            mut_model[region_name][7] = del_weight

            trinuc_trans_prob = pickle_dict[pickel_key(region, Stats.TRINUC_TRANS_PROBS)]
            for k in sorted(trinuc_trans_prob.keys()):
                my_ind = TRI_IND[k[0][0] + k[0][2]]
                (k1, k2) = (NUC_IND[k[0][1]], NUC_IND[k[1][1]])
                mut_model[region_name][8][my_ind][k1][k2] = trinuc_trans_prob[k]
            for i in range(len(mut_model[region_name][8])):
                for j in range(len(mut_model[region_name][8][i])):
                    for k in range(len(mut_model[region_name][8][i][j])):
                        # if trinuc not present in input mutation model, assign it uniform probability
                        if float(sum(mut_model[region_name][8][i][j])) < 1e-12:
                            mut_model[region_name][8][i][j] = [0.25, 0.25, 0.25, 0.25]
                        else:
                            mut_model[region_name][8][i][j][k] /= float(sum(mut_model[region_name][8][i][j]))

            trinuc_mut_prob = pickle_dict[pickel_key(region, Stats.TRINUC_MUT_PROB)]
            which_have_we_seen = {n: False for n in ALL_TRI}
            trinuc_mean = np.mean(list(trinuc_mut_prob.values()))
            for trinuc in trinuc_mut_prob.keys():
                mut_model[region_name][9][ALL_IND[trinuc]] = trinuc_mut_prob[trinuc]
                which_have_we_seen[trinuc] = True
            for trinuc in which_have_we_seen.keys():
                if not which_have_we_seen[trinuc]:
                    mut_model[region_name][9][ALL_IND[trinuc]] = trinuc_mean
        print(f'found the next regions in the model: {region_names_found}')
    return mut_model


def pickel_key(region: Union[Region, str], stats: Stats) -> str:
    return f'{region.value}.{stats.value}' if region != '' else stats.value
