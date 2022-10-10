import random
import sys
import time
from dataclasses import dataclass
from typing import Optional, List, Tuple

import copy
from Bio import Seq
from Bio.Seq import MutableSeq

from mutation_model import MODEL_AVG_MUT_RATE, MODEL_P_INDEL, MODEL_P_INSERTION, MODEL_INS_LEN_DSTRBTN, \
    MODEL_DEL_LEN_DSTRBTN, MODEL_TRINUC_TRANS_DSTRBTN, MODEL_P_TRINUC_MUT
from utilities.annotated_sequence import AnnotatedSequence
from utilities.common_data_structues import Region, MutType, Strand, STOP_CODONS_FORWARD_STRAND, \
    STOP_CODONS_REVERSE_STRAND, VALID_NUCL, TRI_IND, NUC_IND, ALL_IND
from utilities.probability import DiscreteDistribution, poisson_list

import pandas as pd

"""
Constants needed for analysis
"""
MAX_ATTEMPTS = 100  # max attempts to insert a mutation into a valid position
MAX_MUTFRAC = 0.3
# TODO rethink: 1. Do we need it? 2. Do we want to change the scalar? 3. Do we want to make one for each region?

# DEBUG
IGNORE_TRINUC = False

# percentile resolution used for fraglen quantizing
COV_FRAGLEN_PERCENTILE = 10.
LARGE_NUMBER = 9999999999

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
    def __init__(self, debug: bool = False):

        self.windows_start_offset: int = 0
        self.start: int = 0
        self.end: int = 0
        self.original_end: int = 0
        self.debug = debug
        # Blocklist explanation:
        # blocklist[pos] = 0		safe to insert variant here
        # blocklist[pos] = 1		indel inserted here
        # blocklist[pos] = 2		snp inserted here
        # blocklist[pos] = 3		invalid position for various processing reasons
        self.blocklist = {}  # TODO re-consider !!!

    def next_window(self, new_start: int = -1, new_end: int = -1):
        # TODO sanity check about start, end. Maybe consider N regions?
        last_window_offset = self.end - self.original_end
        new_start += (self.windows_start_offset + last_window_offset)
        new_end += (self.windows_start_offset + last_window_offset)
        if new_start < 0 or new_end < 0 or not (self.end <= new_start < new_end):
            raise Exception(f'Illegal window positions. Start and end should be positive integers. '
                            f'The new window should come after the previous one and not overlap with it. '
                            f'Got the next values: start={new_start}, end={new_end}, '
                            f'while the previous values are: start={self.start}, end={self.end}')
        self.windows_start_offset += last_window_offset
        self.start = new_start
        self.end = new_end
        self.original_end = self.end
        self.blocklist = {}  # np.zeros(end-start, dtype='<i4')
        if self.debug:
            print(f"Updated window. Overall windows offset is {self.windows_start_offset}"
                  f" and therefore the adjusted parameters are window: start={self.start}, end={self.end}")

    def adjust_window(self, end_shift: int = 0):
        # TODO sanity check about start, end. Maybe consider N regions?
        if self.debug:
            print(f"Adjusting window, end shift is {end_shift}")
        self.end += end_shift

    def update_blocklist(self, mutation: Mutation):
        for k in range(mutation.position, mutation.position + len(mutation.ref_nucl)):  # TODO continue here
            self.blocklist[k] = 1 if mutation.mut_type.value == MutType.INDEL.value else 2

    def check_blocklist(self, mutation: Mutation):
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
                 annotations_sorted: bool = False, mut_model=None, mut_rate=None, dist=None, debug: bool = False):
        self.debug = debug
        self.chrom_name = chromosome_name
        self.chrom_sequence = MutableSeq(str(chromosome_sequence))
        # TODO consider using Seq class to benefit from the class-supported methods
        self.seq_len = len(chromosome_sequence)
        self.annotated_seq = AnnotatedSequence(annotations_df, chromosome_name, is_sorted=annotations_sorted,
                                               debug=self.debug)
        if self.annotated_seq.len() and len(chromosome_sequence) != self.annotated_seq.len():
            print(f'Chromosome {chromosome_name} in the reference is {len(chromosome_sequence)},'
                  f' while in the annotations file it is {self.annotated_seq.len()}.\n'
                  f'Currently continue normally and hope for the best...')
        self._load_mutation_model(mut_model, mut_rate, dist)
        self.window_unit = WindowUnit(self.debug)

    def _load_mutation_model(self, mut_model: dict, mut_rate: float, dist: float):
        self.model_per_region = copy.deepcopy(mut_model)

        model_regions = [region_name for region_name in self.model_per_region.keys()]
        for region in self.annotated_seq.get_regions():
            if region.value not in model_regions:
                self.model_per_region[region.value] = self.model_per_region[Region.ALL.value]
            if dist:
                self.model_per_region[region.value][MODEL_AVG_MUT_RATE] *= dist
            if mut_rate:  # TODO remove this feature?
                self.model_per_region[region.value][MODEL_AVG_MUT_RATE] = mut_rate

        relevant_regions = [region.value for region in self.annotated_seq.get_regions()]
        redundant_regions = []
        for region_name in self.model_per_region.keys():
            if region_name not in relevant_regions:
                redundant_regions.append(region_name)
        for region_name in redundant_regions:
            del self.model_per_region[region_name]

    def get_window_mutations(self) -> RandomMutationPool:
        # NOTE: window can be also a whole non-N region or the entire chromosome
        indels_to_add_window_per_region, snps_to_add_window_per_region = \
            self.get_planned_snps_and_indels_in_window_per_region()
        max_mutations_in_window = round(MAX_MUTFRAC * (self.window_unit.end - self.window_unit.start))
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

        start = time.time()
        window_seq_len = self.window_unit.end - self.window_unit.start
        trinuc_snp_bias_of_window_per_region = {region.value: [0. for _ in range(window_seq_len)] for region in
                                                self.annotated_seq.get_regions()}
        self.trinuc_bias_per_region = {region.value: None for region in self.annotated_seq.get_regions()}
        for region in self.annotated_seq.get_regions():  # TODO switch order between loops?
            region_mask = self.annotated_seq.get_mask_in_window_of_region(region, self.window_unit.start,
                                                                          self.window_unit.end)
            for i in range(0 + 1, window_seq_len - 1):
                codon = str(self.chrom_sequence[self.window_unit.start + i - 1:self.window_unit.start + i + 2])
                trinuc_snp_bias_of_window_per_region[region.value][i] = \
                    region_mask[i] * self.model_per_region[region.value][MODEL_P_TRINUC_MUT][ALL_IND[codon]]
            self.trinuc_bias_per_region[region.value] = \
                DiscreteDistribution(trinuc_snp_bias_of_window_per_region[region.value], range(window_seq_len))
            # from initialization, the probability of the first and the last element is 0
        end = time.time()
        if self.debug:
            print(f"Updating window trinuc bias took {0:.3f} seconds.".format(end-start))

    def init_poisson(self, type_is_indel: bool = True):
        poisson_per_region = {}
        window_counts_per_region = self.annotated_seq.get_nucleotides_counts_per_region(self.window_unit.start,
                                                                                        self.window_unit.end)
        for region in self.annotated_seq.get_regions():
            param = self.model_per_region[region.value][MODEL_P_INDEL] if type_is_indel else (
                    1. - self.model_per_region[region.value][MODEL_P_INDEL])
            poisson_lambda = \
                window_counts_per_region[region.value] * param * self.model_per_region[region.value][MODEL_AVG_MUT_RATE]
            k_range = range(int(window_counts_per_region[region.value] * MAX_MUTFRAC))
            poisson_per_region[region.value] = poisson_list(k_range, poisson_lambda)
        return poisson_per_region

    def insert_given_mutations(self, vars_in_current_window) -> List[Tuple]:
        mutations_to_insert = self.validate_given_mutations_list(vars_in_current_window)

        # TODO actually insert mutations!!!
        #  consider current_offset and self.window_unit.update(end_shift=mutation.get_offset_change())
        vcf_mutations = self.prepare_mutations_to_vcf(mutations_to_insert, mutations_already_inserted=False)
        return vcf_mutations

    def validate_given_mutations_list(self, input_df: pd.DataFrame):
        inserted_mutations = []
        for _, row in input_df.iterrows():
            ref_nucl = row['allele']
            new_nucl = row['alternatives'][0]  # take the first alternative
            mut_type = MutType.SNP if len(ref_nucl) == 1 and len(new_nucl) == 1 else MutType.INDEL
            mutation = Mutation(row['pos'] + self.window_unit.windows_start_offset, ref_nucl, new_nucl, mut_type)
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

    def next_window(self, start: int, end: int):
        self.window_unit.next_window(start, end)

    def generate_random_mutations(self) -> List[Tuple]:
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
                self.window_unit.adjust_window(end_shift=window_shift)
            # check annotations
            annotation_changed = self.handle_annotations_after_mutated_sequence(inserted_mutation)
            # if at least one has changed - sample mutations again
            if annotation_changed or window_shift != 0:
                if not IGNORE_TRINUC:
                    self._update_trinuc_bias_of_window()

        vcf_mutations = self.prepare_mutations_to_vcf(inserted_mutations, mutations_already_inserted=True)
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
                    range(self.window_unit.start + 1, self.window_unit.end - 1), weights=region_mask[1:-1], k=1)[0]
                # https://pynative.com/python-weighted-random-choices-with-probability/
                # TODO if event_pos is ok return it, otherwise keep trying
                return event_pos
            else:
                event_pos = self.window_unit.start + self.trinuc_bias_per_region[region.value].sample()
                if event_pos <= self.window_unit.start or self.window_unit.end - 1 <= event_pos:
                    continue
                # TODO if event_pos is ok return it, otherwise keep trying
                return event_pos
        return -1

    def insert_snp(self, position: int, region: Region):
        snp = self.get_specific_snp(position, region)
        self.mutate_sequence(snp)
        return snp

    def get_specific_snp(self, position: int, region: Region):
        ref_nucl = self.chrom_sequence[position]
        context = str(self.chrom_sequence[position - 1]) + str(self.chrom_sequence[position + 1])
        # sample from tri-nucleotide substitution matrices to get SNP alt allele
        new_nucl = self.model_per_region[region.value][MODEL_TRINUC_TRANS_DSTRBTN][TRI_IND[context]][NUC_IND[ref_nucl]].sample()
        snp = Mutation(position, ref_nucl, new_nucl, MutType.SNP)
        # self.blocklist[snp.position] = 2  # TODO use blocklist?
        return snp

    def insert_indel(self, position: int, region: Region):
        indel = self.get_specific_indel(position, region)
        window_shift = self.mutate_sequence(indel)
        return indel, window_shift

    def get_specific_indel(self, position: int, region: Region):
        # insertion
        if random.random() <= self.model_per_region[region.value][MODEL_P_INSERTION]:
            indel_len = self.model_per_region[region.value][MODEL_INS_LEN_DSTRBTN].sample()
            # sequence content of random insertions is uniformly random (change this later, maybe)
            indel_seq = ''.join([random.choice(VALID_NUCL) for _ in range(indel_len)])
            ref_nucl = self.chrom_sequence[position]
            indel = Mutation(position, ref_nucl, ref_nucl + indel_seq, MutType.INDEL)

        # deletion
        else:
            indel_len = self.model_per_region[region.value][MODEL_DEL_LEN_DSTRBTN].sample()

            # skip if deletion too close to boundary
            if position + indel_len >= self.window_unit.end:
                # TODO mind that by using self.window_unit.end,
                #  and not self.seq_len + self.sequence_offset + current_offset,
                #  we don't allow deletions to take a "bite" form the next window.
                #  - should we allow that?
                indel_len = self.window_unit.end - 1 - position

            # forbid deletion to crossing the boundary of annotation!
            _, annotation_end = self.annotated_seq.get_annotation_start_end(position)
            if annotation_end and position + indel_len >= annotation_end:
                indel_len = annotation_end - 1 - position

            if indel_len < 0:
                # shouldn't occure
                indel_len = 0
                if self.debug:
                    print(f"Deletion length is less than 0, changing to 0. "
                          f"This is a program error, and need further investigation")

            if indel_len == 0:
                indel_seq = ''
            else:
                indel_seq = str(self.chrom_sequence[position + 1:position + indel_len + 1])
            ref_nucl = self.chrom_sequence[position]
            indel = Mutation(position, ref_nucl + indel_seq, ref_nucl, MutType.INDEL)

        # TODO use blocklist?
        # for k in range(position, position + indel_len + 1):
        #     self.blocklist[k] = 1

        return indel

    def mutate_sequence(self, mutation: Mutation):
        print(f"Trying to insert mutation of type {mutation.mut_type.value} at position {mutation.position}.")

        ref_start = mutation.position
        ref_end = ref_start + len(mutation.ref_nucl) #TODO
        window_shift = mutation.get_offset_change()

        if mutation.ref_nucl != str(self.chrom_sequence[ref_start:ref_end]):
            print('\nError: Something went wrong!\n', mutation, [ref_start, ref_end],
                  str(self.chrom_sequence[ref_start:ref_end]), '\n')
            sys.exit(1)
        else:
            # alter reference sequence
            self.chrom_sequence = self.chrom_sequence[:ref_start] + MutableSeq(mutation.new_nucl) + \
                                  self.chrom_sequence[ref_end:]

        print(f"Inserted mutation successfully. Window shift is {window_shift}")
        return window_shift

    def prepare_mutations_to_vcf(self, inserted_mutations: List[Mutation], mutations_already_inserted: bool) \
            -> List[Tuple]:
        inserted_mutations.sort()
        vcf_mutations = []
        mutations_affected_offset = 0
        for mutation in inserted_mutations:
            vcf_position = mutation.position - self.window_unit.windows_start_offset
            if mutations_already_inserted:
                vcf_position -= mutations_affected_offset
                mutations_affected_offset += mutation.get_offset_change()
            vcf_mutations.append(tuple([vcf_position, mutation.ref_nucl, mutation.new_nucl, 1, 'WP=1']))
        # TODO: combine multiple variants that happened to occur at same position into single vcf entry?
        #       reconsider blocklist for that!
        return vcf_mutations

    def handle_annotations_after_mutated_sequence(self, inserted_mutation: Mutation) -> bool:
        if [region.value for region in self.annotated_seq.get_regions()] == [Region.ALL.value]:
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
            codon = self.chrom_sequence[start:end + 1]
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
        self.annotated_seq.handle_insertion(inserted_mutation.position,
                                            len(inserted_mutation.new_nucl) - len(inserted_mutation.ref_nucl))

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
                    codon = self.chrom_sequence[first_codon_start + (3 * i): first_codon_end + 1 + (3 * i)]
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
        self.annotated_seq.handle_deletion(inserted_mutation.position,
                                           len(inserted_mutation.ref_nucl) - len(inserted_mutation.new_nucl))

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
                codon = self.chrom_sequence[start:end + 1]
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


def is_stop_codon(codon: str, strand: Strand, debug: bool = False) -> bool:
    stop = False

    if strand.value == Strand.FORWARD.value:
        stop = codon in STOP_CODONS_FORWARD_STRAND

    if strand.value == Strand.REVERSE.value:
        stop = codon in STOP_CODONS_REVERSE_STRAND

    if debug:
        print("Encountered stop codon")

    return stop
