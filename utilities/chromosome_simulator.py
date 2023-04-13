import random
import sys
import time
from typing import Optional, List, Tuple
import pandas as pd
import copy
from Bio import Seq
from Bio.Seq import MutableSeq

from utilities.annotated_sequence import AnnotatedSequence
from utilities.common_data_structues import Region, MutType, Strand, STOP_CODONS_FORWARD_STRAND, \
    STOP_CODONS_REVERSE_STRAND, VALID_NUCL, TRI_IND, NUC_IND, ALL_IND, Mutation, ModelKeys
from utilities.io.logger import Logger
from utilities.probability import DiscreteDistribution, poisson_list
from utilities.random_mutations_pool import RandomMutationsPool
from utilities.simulation_window import SimulationWindow

"""
Constants needed for analysis
"""
MAX_ATTEMPTS = 100  # max attempts to insert a mutation into a valid position
MAX_MUTFRAC = 0.3

# DEBUG
IGNORE_TRINUC = False


class ChromosomeSimulator:
    """
    Container for reference sequences, applies mutations
    """

    def __init__(self, chromosome_name: str, chromosome_sequence: Seq, annotations_df: pd.DataFrame,
                 annotations_sorted: bool = False, mut_model=None, mut_scalar=None, dist=None, logger: Logger = None):
        self.logger = logger if logger else Logger()
        self.chrom_name = chromosome_name
        # TODO consider using Seq class to benefit from the class-supported methods:
        self.chrom_sequence = MutableSeq(str(chromosome_sequence))
        self.seq_len = len(chromosome_sequence)
        self.annotated_seq = AnnotatedSequence(annotations_df, chromosome_name, is_sorted=annotations_sorted,
                                               logger=self.logger)
        if self.annotated_seq.len() and len(chromosome_sequence) != self.annotated_seq.len():
            self.logger.message(f'Chromosome {chromosome_name} in the reference is {len(chromosome_sequence)},'
                                f' while in the annotations file it is {self.annotated_seq.len()}.\n'
                                f'Currently continue normally and hope for the best...')
        self._load_mutation_model(mut_model, mut_scalar, dist)
        self.window_unit = SimulationWindow(self.logger)

    def get_annotations_df(self) -> pd.DataFrame:
        return self.annotated_seq.get_annotations_df()

    def next_window(self, start: int, end: int):
        self.window_unit.next_window(start, end)

    def generate_random_mutations(self) -> (List[Tuple], int, int):
        random_mutations_pool = self._get_window_mutations()
        if not IGNORE_TRINUC:
            self._init_trinuc_bias_of_window()

        inserted_mutations = []
        snps_count = 0
        indels_count = 0
        while random_mutations_pool.has_next():
            mut_type, region = random_mutations_pool.get_next()
            inserted_mutation, window_shift = self._insert_random_mutation(mut_type, region)
            if not inserted_mutation:
                continue
            inserted_mutations.append(inserted_mutation)
            if inserted_mutation.mut_type.value == MutType.SNP.value:
                snps_count += 1
            else:
                indels_count += 1
            # handle window
            if window_shift != 0:
                self.window_unit.adjust_window(end_shift=window_shift)
            # handle annotations
            self._handle_annotations_after_mutated_sequence(inserted_mutation)
            # handle trinuc bias
            if not IGNORE_TRINUC:
                self._update_trinuc_bias_of_window(inserted_mutation)

        vcf_mutations = self._prepare_mutations_to_vcf(inserted_mutations, mutations_already_inserted=True)
        return vcf_mutations, snps_count, indels_count

    def _load_mutation_model(self, mut_model: dict, mut_scalar: float, dist: float):
        self.model_per_region = copy.deepcopy(mut_model)

        model_regions = [region_name for region_name in self.model_per_region.keys()]
        for region in self.annotated_seq.get_regions():
            if region.value not in model_regions:
                self.model_per_region[region.value] = copy.deepcopy(self.model_per_region[Region.ALL.value])
            if dist:
                self.model_per_region[region.value][ModelKeys.AVG_MUT_RATE] *= dist
            if mut_scalar:
                self.model_per_region[region.value][ModelKeys.AVG_MUT_RATE] *= mut_scalar

        relevant_regions = [region.value for region in self.annotated_seq.get_regions()]
        redundant_regions = []
        for region_name in self.model_per_region.keys():
            if region_name not in relevant_regions:
                redundant_regions.append(region_name)
        for region_name in redundant_regions:
            del self.model_per_region[region_name]

    def _get_window_mutations(self) -> RandomMutationsPool:
        # NOTE: window can be also a whole non-N region or the entire chromosome
        indels_to_add_window_per_region, snps_to_add_window_per_region = \
            self._get_planned_snps_and_indels_in_window_per_region()
        max_mutations_in_window = round(MAX_MUTFRAC * (self.window_unit.end - self.window_unit.start))
        return RandomMutationsPool(indels_to_add_window_per_region, snps_to_add_window_per_region,
                                   max_mutations_in_window, logger=self.logger)

    def _get_planned_snps_and_indels_in_window_per_region(self) -> (dict, dict):
        indel_poisson_per_region = self._init_poisson(type_is_indel=True)
        snp_poisson_per_region = self._init_poisson(type_is_indel=False)
        indels_to_add_window_per_region = {region.value: indel_poisson_per_region[region.value].sample()
                                           for region in self.annotated_seq.get_regions()}
        snps_to_add_window_per_region = {region.value: snp_poisson_per_region[region.value].sample()
                                         for region in self.annotated_seq.get_regions()}
        return indels_to_add_window_per_region, snps_to_add_window_per_region

    def _init_poisson(self, type_is_indel: bool = True):
        poisson_per_region = {}
        window_counts_per_region = self.annotated_seq.get_nucleotides_counts_per_region(self.window_unit.start,
                                                                                        self.window_unit.end)
        for region in self.annotated_seq.get_regions():
            param = self.model_per_region[region.value][ModelKeys.P_INDEL] if type_is_indel else (
                    1. - self.model_per_region[region.value][ModelKeys.P_INDEL])
            poisson_lambda = \
                window_counts_per_region[region.value] * param * self.model_per_region[region.value][ModelKeys.AVG_MUT_RATE]
            k_range = range(int(window_counts_per_region[region.value] * MAX_MUTFRAC))
            poisson_per_region[region.value] = poisson_list(k_range, poisson_lambda)
        return poisson_per_region

    def _init_trinuc_bias_of_window(self):
        """
        initialize trinuc snp bias
        compute mutation positional bias given trinucleotide strings of the sequence (ONLY AFFECTS SNPs)
        :return:
        """
        t_start = time.time()
        window_seq_len = self.window_unit.end - self.window_unit.start
        trinuc_snp_bias_of_window_per_region = {region.value: [0. for _ in range(window_seq_len)] for region in
                                                self.annotated_seq.get_regions()}
        self.trinuc_bias_per_region = {region.value: None for region in self.annotated_seq.get_regions()}
        for region in self.annotated_seq.get_regions():
            region_mask = self.annotated_seq.get_mask_in_window_of_region(region, self.window_unit.start,
                                                                          self.window_unit.end)
            for i in range(0 + 1, window_seq_len - 1):
                codon = str(self.chrom_sequence[self.window_unit.start + i - 1:self.window_unit.start + i + 2])
                trinuc_snp_bias_of_window_per_region[region.value][i] = \
                    region_mask[i] * self.model_per_region[region.value][ModelKeys.P_TRINUC_MUT][ALL_IND[codon]]
            self.trinuc_bias_per_region[region.value] = \
                DiscreteDistribution(trinuc_snp_bias_of_window_per_region[region.value], range(window_seq_len))
            # from initialization, the probability of the first and the last element is 0
        t_end = time.time()
        self.logger.debug_message(f"Initializing window trinuc bias took {0:.3f} seconds.".format(t_end-t_start))

    def _update_trinuc_bias_of_window(self, mut: Mutation):
        """
        update trinuc snp bias after mutation has occurred
        :return:
        """
        t_start = time.time()

        window_seq_len = self.window_unit.end - self.window_unit.start
        forbidden_zone_start = max(mut.position - 1, self.window_unit.start) - self.window_unit.start
        forbidden_zone_end = min(mut.position + len(mut.new_nucl), self.window_unit.end - 1) - self.window_unit.start
        forbidden_positions = range(forbidden_zone_start, forbidden_zone_end + 1)

        for region in self.annotated_seq.get_regions():
            if self.trinuc_bias_per_region[region.value].degenerate is None:
                new_weights = self.trinuc_bias_per_region[region.value].weights
                # INSERTION
                if mut.get_offset_change() > 0:
                    offset = mut.get_offset_change()
                    for i in range(offset):
                        new_weights.insert(mut.position + 1 + i - self.window_unit.start, 0.)
                # DELETION
                elif mut.get_offset_change() < 0:
                    offset = -mut.get_offset_change()
                    del new_weights[mut.position + 1 - self.window_unit.start:
                                    mut.position + offset + 1 - self.window_unit.start]
            else:
                new_weights = [self.trinuc_bias_per_region[region.value].degenerate for _ in range(window_seq_len)]

            for pos in forbidden_positions:
                new_weights[pos] = 0.

            self.trinuc_bias_per_region[region.value] = DiscreteDistribution(new_weights, range(window_seq_len))

        t_end = time.time()
        self.logger.debug_message(f"Updating window trinuc bias took {0:.3f} seconds.".format(t_end-t_start))

    def _insert_random_mutation(self, mut_type: MutType, region: Region) -> (Optional[Mutation], int):
        window_shift = 0

        position = self._find_position_for_mutation(mut_type, region)
        if position == -1:
            return None, window_shift

        inserted_mutation = None
        if mut_type.value == MutType.SNP.value:
            inserted_mutation = self._insert_snp(position, region)

        elif mut_type.value == MutType.INDEL.value:
            inserted_mutation, window_shift = self._insert_indel(position, region)

        # elif mut_type == MutType.SV:
        #   ...

        return inserted_mutation, window_shift

    def _find_position_for_mutation(self, mut_type: MutType, region: Region) -> int:
        region_mask = self.annotated_seq.get_mask_in_window_of_region(region, self.window_unit.start,
                                                                      self.window_unit.end)
        if 1 not in region_mask:
            return -1  # current annotation doesn't exist in window
        for _ in range(MAX_ATTEMPTS):
            if mut_type.value == MutType.INDEL.value or IGNORE_TRINUC:
                # -2 because we don't allow SNP in the window start/end
                k = self.window_unit.end - self.window_unit.start - 2
                if k < 1:
                    return -1
                event_pos = random.choices(
                    range(self.window_unit.start + 1, self.window_unit.end - 1), weights=region_mask[1:-1], k=1)[0]
                # https://pynative.com/python-weighted-random-choices-with-probability/
                return event_pos
            else:
                event_pos = self.window_unit.start + self.trinuc_bias_per_region[region.value].sample()
                if event_pos <= self.window_unit.start or self.window_unit.end - 1 <= event_pos:
                    continue
                return event_pos
        return -1

    def _insert_snp(self, position: int, region: Region):
        snp = self._get_specific_snp(position, region)
        self._mutate_sequence(snp)
        return snp

    def _get_specific_snp(self, position: int, region: Region):
        ref_nucl = self.chrom_sequence[position]
        context = str(self.chrom_sequence[position - 1]) + str(self.chrom_sequence[position + 1])
        # sample from tri-nucleotide substitution matrices to get SNP alt allele
        new_nucl = self.model_per_region[region.value][ModelKeys.TRINUC_TRANS_DSTRBTN][TRI_IND[context]][NUC_IND[ref_nucl]]\
            .sample()
        snp = Mutation(position, ref_nucl, new_nucl, MutType.SNP)
        return snp

    def _insert_indel(self, position: int, region: Region):
        indel = self._get_specific_indel(position, region)
        window_shift = self._mutate_sequence(indel)
        return indel, window_shift

    def _get_specific_indel(self, position: int, region: Region):
        # insertion
        if random.random() <= self.model_per_region[region.value][ModelKeys.P_INSERTION]:
            indel_len = self.model_per_region[region.value][ModelKeys.INS_LEN_DSTRBTN].sample()
            # sequence content of random insertions is uniformly random (change this later, maybe)
            indel_seq = ''.join([random.choice(VALID_NUCL) for _ in range(indel_len)])
            ref_nucl = self.chrom_sequence[position]
            indel = Mutation(position, ref_nucl, ref_nucl + indel_seq, MutType.INDEL)

        # deletion
        else:
            indel_len = self.model_per_region[region.value][ModelKeys.DEL_LEN_DSTRBTN].sample()

            # skip if deletion too close to boundary
            if position + indel_len >= self.window_unit.end:
                indel_len = self.window_unit.end - 1 - position
                # - mind that by using self.window_unit.end,
                #  and not self.seq_len + self.sequence_offset + current_offset,
                #  we don't allow deletions to take a "bite" form the next window.
                #  Should we allow that?

            # forbid deletion to crossing the boundary of annotation!
            _, annotation_end = self.annotated_seq.get_annotation_start_end(position)
            if annotation_end and position + indel_len >= annotation_end:
                indel_len = annotation_end - 1 - position

            if indel_len < 0:
                # shouldn't occure
                indel_len = 0
                self.logger.debug_message(f"Deletion length is less than 0, changing to 0. "
                                          f"This is a program error, and need further investigation")

            if indel_len == 0:
                indel_seq = ''
            else:
                indel_seq = str(self.chrom_sequence[position + 1:position + indel_len + 1])
            ref_nucl = self.chrom_sequence[position]
            indel = Mutation(position, ref_nucl + indel_seq, ref_nucl, MutType.INDEL)

        return indel

    def _mutate_sequence(self, mutation: Mutation):

        ref_start = mutation.position
        ref_end = ref_start + len(mutation.ref_nucl)
        window_shift = mutation.get_offset_change()

        if mutation.ref_nucl != str(self.chrom_sequence[ref_start:ref_end]):
            self.logger.message(f'Error: Failed to insert the next mutation: {mutation}\n'
                                f'[ref_start, ref_end]='
                                f'{[ref_start, ref_end]}, {str(self.chrom_sequence[ref_start:ref_end])}')
            sys.exit(1)
        else:
            # alter reference sequence
            self.chrom_sequence = self.chrom_sequence[:ref_start] + MutableSeq(mutation.new_nucl) + \
                                  self.chrom_sequence[ref_end:]

        self.logger.message(f"Inserted mutation: "
                            f"{mutation.mut_type.value} at position {mutation.position}. "
                            f"Window shift: {window_shift}")
        return window_shift

    def _handle_annotations_after_mutated_sequence(self, inserted_mutation: Mutation):
        if [region.value for region in self.annotated_seq.get_regions()] == [Region.ALL.value]:
            return

        # SNP
        if inserted_mutation.mut_type.value == MutType.SNP.value:
            self._handle_annotations_after_snp(inserted_mutation)
        # INDEL - insertion
        elif inserted_mutation.mut_type.value == MutType.INDEL.value \
                and len(inserted_mutation.ref_nucl) < len(inserted_mutation.new_nucl):
            self._handle_annotations_after_small_insertion(inserted_mutation)
        # INDEL - deletion
        elif inserted_mutation.mut_type.value == MutType.INDEL.value \
                and len(inserted_mutation.ref_nucl) >= len(inserted_mutation.new_nucl):
            self._handle_annotations_after_small_deletion(inserted_mutation)
        # SVs and other - currently not supported
        else:
            raise Exception("currently supporting only SNPs and INDELs")

    def _handle_annotations_after_snp(self, inserted_mutation: Mutation):
        """
        Handle annotations after SNP
        :param inserted_mutation:
        :return: True if annotation has changed somehow, False otherwise
        """
        region, strand = self.annotated_seq.get_region_by_position(inserted_mutation.position)

        if region.value == Region.INTERGENIC.value or region.value == Region.NON_CODING_GENE.value:
            return  # current behaviour is not to check for start/stop codon in these regions.

        elif region.value == Region.CDS.value:
            start, _, end = self.annotated_seq.get_encapsulating_codon_positions(inserted_mutation.position)
            # assuming sequence has been mutated already!
            codon = self.chrom_sequence[start:end + 1]
            is_stop = is_stop_codon(codon, strand)
            is_last_codon = self._is_last_coding_position(start, end)
            if (is_stop and not is_last_codon) or (not is_stop and is_last_codon):
                if is_stop and not is_last_codon:
                    self.logger.message('SNP caused stop codon gain')
                elif not is_stop and is_last_codon:
                    self.logger.message('SNP caused stop codon loss')
                self.annotated_seq.mute_gene(position_on_gene=inserted_mutation.position)
        else:
            raise Exception("unknown annotation")

    def _handle_annotations_after_small_insertion(self, inserted_mutation: Mutation):
        if self._should_mute_gene_after_small_insertion(inserted_mutation):
            self.annotated_seq.mute_gene(position_on_gene=inserted_mutation.position)
        self.annotated_seq.handle_insertion(inserted_mutation.position,
                                            len(inserted_mutation.new_nucl) - len(inserted_mutation.ref_nucl))

    def _should_mute_gene_after_small_insertion(self, inserted_mutation: Mutation) -> bool:
        region, strand = self.annotated_seq.get_region_by_position(inserted_mutation.position)

        if region.value == Region.INTERGENIC.value or region.value == Region.NON_CODING_GENE.value:
            return False

        elif region.value == Region.CDS.value:
            frameshift = (len(inserted_mutation.new_nucl) - 1) % 3 != 0
            if frameshift:
                self.logger.message('Insertion caused frame shift')
                return True
            else:
                first_codon_start, _, first_codon_end = \
                    self.annotated_seq.get_encapsulating_codon_positions(inserted_mutation.position)
                added_codons = int((len(inserted_mutation.new_nucl) - 1) / 3)
                for i in range(added_codons):
                    # assuming sequence has been mutated already!
                    codon = self.chrom_sequence[first_codon_start + (3 * i): first_codon_end + 1 + (3 * i)]
                    if is_stop_codon(codon, strand):
                        self.logger.message('Insertion caused stop codon gain')
                        return True
            return False

        else:
            raise Exception("unknown annotation")

    def _handle_annotations_after_small_deletion(self, inserted_mutation: Mutation):
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

        if self._should_mute_gene_after_small_deletion(inserted_mutation):
            self.annotated_seq.mute_gene(position_on_gene=inserted_mutation.position)
        self.annotated_seq.handle_deletion(inserted_mutation.position,
                                           len(inserted_mutation.ref_nucl) - len(inserted_mutation.new_nucl))

    def _should_mute_gene_after_small_deletion(self, inserted_mutation: Mutation) -> bool:
        region, strand = self.annotated_seq.get_region_by_position(inserted_mutation.position)

        if region.value == Region.INTERGENIC.value or region.value == Region.NON_CODING_GENE.value:
            return False

        elif region.value == Region.CDS.value:
            frameshift = (len(inserted_mutation.ref_nucl) - 1) % 3 != 0
            if frameshift:
                self.logger.message('Deletion caused frame shift')
                return True
            else:
                start, _, end = self.annotated_seq.get_encapsulating_codon_positions(inserted_mutation.position)
                # assuming sequence has been mutated already!
                codon = self.chrom_sequence[start:end + 1]
                is_stop = is_stop_codon(codon, strand)
                is_last_codon = self._is_last_coding_position(start, end)
                if is_stop and not is_last_codon:
                    self.logger.message('Deletion caused stop codon gain')
                    return True
                elif not is_stop and is_last_codon:
                    self.logger.message('Deletion caused stop codon loss')
                    return True
            return False

        else:
            raise Exception("unknown annotation")

    def _is_last_coding_position(self, start: int, end: int):
        last_coding_position, strand = self.annotated_seq.get_last_coding_position_of_encapsulating_gene(end)
        is_last_codon = (last_coding_position == end and strand.value == Strand.FORWARD.value) or \
                        (last_coding_position == start and strand.value == Strand.REVERSE.value)
        return is_last_codon

    def _prepare_mutations_to_vcf(self, inserted_mutations: List[Mutation], mutations_already_inserted: bool) \
            -> List[Tuple]:
        inserted_mutations.sort()  # sort by positions, which derives from Mutation class definition
        vcf_mutations = []
        mutations_affected_offset = 0
        for mutation in inserted_mutations:
            vcf_position = mutation.position - self.window_unit.windows_start_offset
            if mutations_already_inserted:
                vcf_position -= mutations_affected_offset
                mutations_affected_offset += mutation.get_offset_change()
            vcf_mutations.append(tuple([vcf_position, mutation.ref_nucl, mutation.new_nucl]))
        return vcf_mutations

    def get_genes(self) -> set:
        return self.annotated_seq.get_genes()


def is_stop_codon(codon: str, strand: Strand) -> bool:
    stop = False

    if strand.value == Strand.FORWARD.value:
        stop = codon in STOP_CODONS_FORWARD_STRAND

    if strand.value == Strand.REVERSE.value:
        stop = codon in STOP_CODONS_REVERSE_STRAND
    return stop
