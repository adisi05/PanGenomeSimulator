import random
import copy
import pathlib
import bisect
import pickle
import sys
from enum import Enum

import numpy as np
import pandas as pd
from Bio.Seq import Seq, MutableSeq

from probability import DiscreteDistribution, poisson_list

"""
Constants needed for analysis
"""
MAX_ATTEMPTS = 100  # max attempts to insert a mutation into a valid position
MAX_MUTFRAC = 0.3  # the maximum percentage of a window that can contain mutations

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


class SequenceContainer:
    """
    Container for reference sequences, applies mutations
    """

    def __init__(self, x_offset, sequence, ploidy, window_overlap, mut_models=None, mut_rate=None, dist=None):

        self.update_sequence(x_offset, sequence, ploidy, window_overlap, mut_models, mut_rate, dist)

    def update_sequence(self, x_offset, sequence, ploidy, window_overlap, mut_models=None, mut_rate=None, dist=None):
        # initialize mutation models or reinitialize if changed
        if ploidy != self.ploidy or mut_rate != self.mut_rescale or mut_models is not None:
            self.ploidy = ploidy
            self.update_mut_models(mut_models, mut_rate, dist)

        # initialize poisson attributes or redo them if sequence length is different than previous window
        if len(sequence) != self.seq_len:
            self.seq_len = len(sequence)
            self.indel_poisson_per_region, self.snp_poisson_per_region = self.init_poisson()

        # basic vars
        self.update_basic_vars(x_offset, sequence, window_overlap)

        # sample the number of variants that will be inserted into each ploid
        self.indels_to_add_per_region =  {region: [n.sample() for n in self.indel_poisson_per_region[region]] for region in Regions}
        self.snps_to_add_per_region = {region: [n.sample() for n in self.snp_poisson_per_region[region]] for region in Regions}

        # initialize trinuc snp bias
        if not IGNORE_TRINUC:
            self.update_trinuc_bias()


    def update_mut_models(self, mut_models, mut_rate, dist):
        if not mut_models:
            single_ploid_model = {region: copy.deepcopy(DEFAULT_MODEL_1) for region in Regions}
            default_model = [single_ploid_model for _ in range(self.ploidy)]
            self.model_data = default_model[:self.ploidy]
        else:
            if len(mut_models) != self.ploidy:
                print('\nError: Number of mutation models received is not equal to specified ploidy\n')
                sys.exit(1)
            self.model_data = copy.deepcopy(mut_models)

        # do we need to rescale mutation frequencies?
        mut_rate_sum_per_region = {region: sum([n[region][0] for n in self.model_data]) for region in Regions}
        self.mut_rescale = mut_rate
        if self.mut_rescale is None:
            self.mut_scalar_per_region = {region: 1.0 for region in Regions}
        else:
            self.mut_scalar_per_region = {region: float(self.mut_rescale) //
                         (mut_rate_sum_per_region[region] / float(len(self.model_data))) for region in Regions}
        if dist:
            self.mut_scalar_per_region = {region: self.mut_scalar_per_region[region] * dist for region in Regions}

        # how are mutations spread to each ploid, based on their specified mut rates?
        self.ploid_mut_frac_per_region = \
            {region: [float(n[0]) / mut_rate_sum_per_region[region] for n in self.model_data] for region in Regions}
        self.ploid_mut_prior_per_region = \
            {region: DiscreteDistribution(self.ploid_mut_frac_per_region[region], range(self.ploidy))
             for region in Regions}

        # init mutation models
        #
        # self.models_per_region[region][ploid][0] = average mutation rate
        # self.models_per_region[region][ploid][1] = p(mut is homozygous | mutation occurs)
        # self.models_per_region[region][ploid][2] = p(mut is indel | mut occurs)
        # self.models_per_region[region][ploid][3] = p(insertion | indel occurs)
        # self.models_per_region[region][ploid][4] = distribution of insertion lengths
        # self.models_per_region[region][ploid][5] = distribution of deletion lengths
        # self.models_per_region[region][ploid][6] = distribution of trinucleotide SNP transitions
        # self.models_per_region[region][ploid][7] = p(trinuc mutates)
        self.models_per_region = {region: [] for region in Regions}
        for region in Regions:
            for n in self.model_data:
                self.models_per_region[region].append(
                    [self.mut_scalar_per_region[region] * n[0], n[1], n[2], n[3], DiscreteDistribution(n[5], n[4]),
                     DiscreteDistribution(n[7], n[6]), []])
                for m in n[8]:
                    # noinspection PyTypeChecker
                    self.models_per_region[region][-1][6].append(
                        [DiscreteDistribution(m[0], NUCL), DiscreteDistribution(m[1], NUCL),
                         DiscreteDistribution(m[2], NUCL), DiscreteDistribution(m[3], NUCL)])
                self.models_per_region[region][-1].append([m for m in n[9]])


    def update_basic_vars(self, x_offset, sequence, window_overlap):
        self.x = x_offset
        self.sequences = [Seq(str(sequence)) for _ in range(self.ploidy)]
        self.indel_list = [[] for _ in range(self.ploidy)]
        self.snp_list = [[] for _ in range(self.ploidy)]

        # Blacklist explanation:
        # black_list[ploid][pos] = 0		safe to insert variant here
        # black_list[ploid][pos] = 1		indel inserted here
        # black_list[ploid][pos] = 2		snp inserted here
        # black_list[ploid][pos] = 3		invalid position for various processing reasons
        self.black_list = [np.zeros(self.seq_len, dtype='<i4') for _ in range(self.ploidy)]

        # disallow mutations to occur on window overlap points
        self.win_buffer = window_overlap
        for p in range(self.ploidy):
            self.black_list[p][-self.win_buffer] = 3
            self.black_list[p][-self.win_buffer - 1] = 3


    def update_trinuc_bias(self):
        # initialize/update trinuc snp bias
        # compute mutation positional bias given trinucleotide strings of the sequence (ONLY AFFECTS SNPs)
        #
        # note: since indels are added before snps, it's possible these positional biases aren't correctly utilized
        #       at positions affected by indels. At the moment I'm going to consider this negligible.
        trinuc_snp_bias = [[0. for _ in range(self.seq_len)] for _ in range(self.ploidy)]
        self.trinuc_bias_per_region = {region: [None for _ in range(self.ploidy)] for region in Regions}
        for p in range(self.ploidy):
            for region in Regions:
                for i in range(self.win_buffer + 1, self.seq_len - 1):
                    trinuc_snp_bias[p][i] = self.models_per_region[region][p][7][ALL_IND[str(self.sequences[p][i - 1:i + 2])]]
                self.trinuc_bias_per_region[region][p] = DiscreteDistribution(trinuc_snp_bias[p][self.win_buffer + 1:self.seq_len - 1],
                                                           range(self.win_buffer + 1, self.seq_len - 1))

    def init_poisson(self):
        ind_l_list_per_region = {}
        snp_l_list_per_region = {}
        for region in Regions:
            ind_l_list_per_region[region] = [self.seq_len * self.models_per_region[region][i][0] * self.models_per_region[region][i][2] * self.ploid_mut_frac_per_region[region][i] for i in
                          range(len(self.models_per_region[region]))]
            snp_l_list_per_region[region] = [self.seq_len * self.models_per_region[region][i][0] * (1. - self.models_per_region[region][i][2]) * self.ploid_mut_frac_per_region[region][i] for i in
                          range(len(self.models_per_region[region]))]
        k_range = range(int(self.seq_len * MAX_MUTFRAC))
        # return (indel_poisson, snp_poisson)
        # TODO These next two lines are really slow. Maybe there's a better way
        return {region:
                    [poisson_list(k_range, ind_l_list_per_region[region][n])
                     for n in range(len(self.models_per_region[region]))]
                for region in Regions}, \
               {region:
                    [poisson_list(k_range, snp_l_list_per_region[region][n])
                     for n in range(len(self.models_per_region[region]))]
                for region in Regions}



    def insert_given_mutations(self, input_list):
        for input_variable in input_list:
            which_alts, which_ploids = self.determine_given_mutation_ploids(input_variable)

            for i in range(len(which_ploids)):
                p = which_ploids[i]
                my_alt = input_variable[2][which_alts[i]]
                my_var = (input_variable[0] - self.x, input_variable[1], my_alt)
                # This is a potential fix implemented by Zach in a previous commit. He left the next line in.
                # in_len = max([len(input_variable[1]), len(my_alt)])
                in_len = len(input_variable[1])

                if my_var[0] < 0 or my_var[0] >= len(self.black_list[p]):
                    print('\nError: Attempting to insert variant out of window bounds:')
                    print(my_var, '--> blackList[0:' + str(len(self.black_list[p])) + ']\n')
                    sys.exit(1)
                if len(input_variable[1]) == 1 and len(my_alt) == 1:
                    if self.black_list[p][my_var[0]]:
                        continue
                    self.snp_list[p].append(my_var)
                    self.black_list[p][my_var[0]] = 2
                else:
                    indel_failed = False
                    for k in range(my_var[0], my_var[0] + in_len):
                        if k >= len(self.black_list[p]):
                            indel_failed = True
                            continue
                        if self.black_list[p][k]:
                            indel_failed = True
                            continue
                    if indel_failed:
                        continue
                    for k in range(my_var[0], my_var[0] + in_len):
                        self.black_list[p][k] = 1
                    self.indel_list[p].append(my_var)

    def determine_given_mutation_ploids(self, input_variable):
        which_ploids = []
        wps = input_variable[4][0]
        # if no genotype given, assume heterozygous and choose a single ploid based on their mut rates
        if wps is None:
            region = self.get_region()
            which_ploids.append(self.ploid_mut_prior_per_region[region].sample())
            which_alts = [0]
        else:
            if '/' in wps or '|' in wps:
                if '/' in wps:
                    splt = wps.split('/')
                else:
                    splt = wps.split('|')
                which_ploids = []
                for i in range(len(splt)):
                    if splt[i] == '1':
                        which_ploids.append(i)
                # assume we're just using first alt for inserted variants?
                which_alts = [0 for _ in which_ploids]
            # otherwise assume monoploidy
            else:
                which_ploids = [0]
                which_alts = [0]
        # ignore invalid ploids
        for i in range(len(which_ploids) - 1, -1, -1):
            if which_ploids[i] >= self.ploidy:
                del which_ploids[i]
        return which_alts, which_ploids

    def insert_random_mutations(self, start, end):

        all_indels = self.pick_random_indels()

        all_snps = self.pick_random_snps()

        self.insert_snps(all_snps)

        self.insert_indels(all_indels)

        output_variants = self.sum_up_variants(all_indels, all_snps)
        return output_variants

    def sum_up_variants(self, all_indels, all_snps):
        # tally up all the variants we handled...
        count_dict = {}
        all_variants = [sorted(all_snps[i] + all_indels[i]) for i in range(self.ploidy)]
        for i in range(len(all_variants)):
            for j in range(len(all_variants[i])):
                all_variants[i][j] = tuple([all_variants[i][j][0] + self.x]) + all_variants[i][j][1:]
                t = tuple(all_variants[i][j])
                if t not in count_dict:
                    count_dict[t] = []
                count_dict[t].append(i)
        # TODO: combine multiple variants that happened to occur at same position into single vcf entry?
        output_variants = []
        for k in sorted(count_dict.keys()):
            output_variants.append(k + tuple([len(count_dict[k]) / float(self.ploidy)]))
            ploid_string = ['0' for _ in range(self.ploidy)]
            for k2 in [n for n in count_dict[k]]:
                ploid_string[k2] = '1'
            output_variants[-1] += tuple(['WP=' + '/'.join(ploid_string)])
        return output_variants

    def insert_snps(self, all_snps):# TODO here should update annotation
        # combine random snps with inserted snps, remove any snps that overlap indels
        for p in range(len(all_snps)):
            all_snps[p].extend(self.snp_list[p])
            all_snps[p] = [n for n in all_snps[p] if self.black_list[p][n[0]] != 1]
        # MODIFY REFERENCE STRING: SNPS
        for i in range(len(all_snps)):
            temp = MutableSeq(self.sequences[i])
            for j in range(len(all_snps[i])):
                v_pos = all_snps[i][j][0]

                if all_snps[i][j][1] != temp[v_pos]:
                    print('\nError: Something went wrong!\n', all_snps[i][j], temp[v_pos], '\n')
                    print(all_snps[i][j])
                    print(self.sequences[i][v_pos])
                    sys.exit(1)
                else:
                    temp[v_pos] = all_snps[i][j][2]
            self.sequences[i] = Seq(temp)

    def insert_indels(self, all_indels): # TODO here should update annotation?
        # organize the indels we want to insert
        for i in range(len(all_indels)):
            all_indels[i].extend(self.indel_list[i])
        all_indels_ins = [sorted([list(m) for m in n]) for n in all_indels]
        # MODIFY REFERENCE STRING: INDELS
        for i in range(len(all_indels_ins)):
            rolling_adj = 0

            for j in range(len(all_indels_ins[i])):
                v_pos = all_indels_ins[i][j][0] + rolling_adj
                v_pos2 = v_pos + len(all_indels_ins[i][j][1])
                indel_length = len(all_indels_ins[i][j][2]) - len(all_indels_ins[i][j][1])
                rolling_adj += indel_length

                if all_indels_ins[i][j][1] != str(self.sequences[i][v_pos:v_pos2]):
                    print('\nError: Something went wrong!\n', all_indels_ins[i][j], [v_pos, v_pos2],
                          str(self.sequences[i][v_pos:v_pos2]), '\n')
                    sys.exit(1)
                else:
                    # alter reference sequence
                    self.sequences[i] = self.sequences[i][:v_pos] + Seq(all_indels_ins[i][j][2]) + \
                                        self.sequences[i][v_pos2:]


    def pick_random_snps(self):
        # add random snps
        all_snps = [[] for _ in self.sequences]
        for region in Regions:
            for i in range(self.ploidy):
                random_snps_minus_inserted = max(self.snps_to_add[i] - len(self.snp_list[i]), 0)
                for j in range(random_snps_minus_inserted):
                    which_ploids = self.determine_random_mutation_ploids(region, i)

                    event_pos = self.find_position_snp(region, which_ploids)
                    if event_pos == -1:
                        continue

                    ref_nucl = self.sequences[i][event_pos]
                    context = str(self.sequences[i][event_pos - 1]) + str(self.sequences[i][event_pos + 1])
                    # sample from tri-nucleotide substitution matrices to get SNP alt allele
                    new_nucl = self.models_per_region[region][i][6][TRI_IND[context]][NUC_IND[ref_nucl]].sample()
                    my_snp = (event_pos, ref_nucl, new_nucl)

                    for p in which_ploids:
                        all_snps[p].append(my_snp)
                        self.black_list[p][my_snp[0]] = 2
        return all_snps

    # TODO consider regions
    def find_position_snp(self, region, which_ploid):
        # try to find suitable places to insert snps
        event_pos = ???
        if IGNORE_TRINUC:
            event_pos = random.randint(self.win_buffer + 1,
                                       self.seq_len - 2)
        else:
            ploid_to_use = which_ploid[random.randint(0, len(which_ploid) - 1)]
            event_pos = self.trinuc_bias_per_region[region][
                ploid_to_use].sample()
        for p in which_ploid:
            if self.black_list[p][event_pos]:
                event_pos = -1
        for attempt in range(MAX_ATTEMPTS):
            # based on the mutation model for the specified ploid, choose a SNP location based on trinuc bias
            # (if there are multiple ploids, choose one at random)
            if event_pos != -1:
                break
        return event_pos

    def pick_random_indels(self):
        # add random indels
        all_indels = [[] for _ in self.sequences]
        for region in Regions:
            for i in range(self.ploidy):
                random_indels_minus_inserted = max(self.indels_to_add[i] - len(self.indel_list[i]), 0)
                for j in range(random_indels_minus_inserted):
                    which_ploid = self.determine_random_mutation_ploids(region, i)

                    event_pos = self.find_position_indel(which_ploid)
                    if event_pos == -1:
                        continue

                    # insertion
                    if random.random() <= self.models_per_region[region][i][3]:
                        in_len = self.models_per_region[region][i][4].sample()
                        # sequence content of random insertions is uniformly random (change this later, maybe)
                        in_seq = ''.join([random.choice(NUCL) for _ in range(in_len)])
                        ref_nucl = self.sequences[i][event_pos]
                        my_indel = (event_pos, ref_nucl, ref_nucl + in_seq)
                    # deletion
                    else:
                        in_len = self.models_per_region[region][i][5].sample()
                        # skip if deletion too close to boundary
                        if event_pos + in_len + 1 >= len(self.sequences[i]):
                            continue
                        if in_len == 1:
                            in_seq = self.sequences[i][event_pos + 1]
                        else:
                            in_seq = str(self.sequences[i][event_pos + 1:event_pos + in_len + 1])
                        ref_nucl = self.sequences[i][event_pos]
                        my_indel = (event_pos, ref_nucl + in_seq, ref_nucl)

                    # if event too close to boundary, skip. if event conflicts with other indel, skip.
                    skip_event = False
                    if event_pos + len(my_indel[1]) >= self.seq_len - self.win_buffer - 1:
                        skip_event = True
                    if skip_event:
                        continue
                    for p in which_ploid:
                        for k in range(event_pos, event_pos + in_len + 1):
                            if self.black_list[p][k]:
                                skip_event = True
                    if skip_event:
                        continue

                    for p in which_ploid:
                        for k in range(event_pos, event_pos + in_len + 1):
                            self.black_list[p][k] = 1
                        all_indels[p].append(my_indel)

        return all_indels

    def determine_random_mutation_ploids(self, region, i):
        # insert homozygous indel
        if random.random() <= self.models_per_region[region][i][1]:
            which_ploid = range(self.ploidy)
        # insert heterozygous indel
        else:
            which_ploid = [self.ploid_mut_prior_per_region[region].sample()]
        return which_ploid

    # TODO consider regions
    def find_position_indel(self, which_ploid):
        # try to find suitable places to insert indels
        event_pos = -1 ???
        for attempt in range(MAX_ATTEMPTS):
            event_pos = random.randint(self.win_buffer, self.seq_len - 1)
            for p in which_ploid:
                if self.black_list[p][event_pos]:
                    event_pos = -1
            if event_pos != -1:
                break
        return event_pos

# parse mutation model pickle file
def parse_input_mutation_model(model=None, which_default=1):
    if which_default == 1:
        out_model = {region: [copy.deepcopy(n) for n in DEFAULT_MODEL_1] for region in Regions}
    elif which_default == 2:
        out_model = {region: [copy.deepcopy(n) for n in DEFAULT_MODEL_2] for region in Regions}
    else:
        print('\nError: Unknown default mutation model specified\n')
        sys.exit(1)

    if model is not None:
        pickle_dict = pickle.load(open(model, "rb"))
        for region in Regions:
            out_model[region][0] = pickle_dict[f'{region.value}.AVG_MUT_RATE']
            out_model[region][2] = 1. - pickle_dict[f'{region.value}.SNP_FREQ']

            ins_list = pickle_dict[f'{region.value}.INDEL_FREQ']
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
            out_model[region][3] = ins_count / float(ins_count + del_count)
            out_model[region][4] = ins_vals
            out_model[region][5] = ins_weight
            out_model[region][6] = del_vals
            out_model[region][7] = del_weight

            trinuc_trans_prob = pickle_dict[f'{region.value}.TRINUC_TRANS_PROBS']
            for k in sorted(trinuc_trans_prob.keys()):
                my_ind = TRI_IND[k[0][0] + k[0][2]]
                (k1, k2) = (NUC_IND[k[0][1]], NUC_IND[k[1][1]])
                out_model[region][8][my_ind][k1][k2] = trinuc_trans_prob[k]
            for i in range(len(out_model[region][8])):
                for j in range(len(out_model[region][8][i])):
                    for l in range(len(out_model[region][8][i][j])):
                        # if trinuc not present in input mutation model, assign it uniform probability
                        if float(sum(out_model[region][8][i][j])) < 1e-12:
                            out_model[region][8][i][j] = [0.25, 0.25, 0.25, 0.25]
                        else:
                            out_model[region][8][i][j][l] /= float(sum(out_model[region][8][i][j]))

            trinuc_mut_prob = pickle_dict[f'{region.value}.TRINUC_MUT_PROB']
            which_have_we_seen = {n: False for n in ALL_TRI}
            trinuc_mean = np.mean(list(trinuc_mut_prob.values()))
            for trinuc in trinuc_mut_prob.keys():
                out_model[region][9][ALL_IND[trinuc]] = trinuc_mut_prob[trinuc]
                which_have_we_seen[trinuc] = True
            for trinuc in which_have_we_seen.keys():
                if not which_have_we_seen[trinuc]:
                    out_model[region][9][ALL_IND[trinuc]] = trinuc_mean

    return out_model

#####################################
#    Supporting Data Structures     #
#####################################

class Regions(Enum):
    EXON = 'exon'
    INTRON = 'intron'
    INTERGENIC = 'intergenic'
    ALL = 'all'


class Stats(Enum):
    # how many times do we observe each trinucleotide in the reference (and input bed region, if present)?
    TRINUC_REF_COUNT = 'TRINUC_REF_COUNT'
    # [(trinuc_a, trinuc_b)] = # of times we observed a mutation from trinuc_a into trinuc_b
    TRINUC_TRANSITION_COUNT = 'TRINUC_TRANSITION_COUNT'
    # total count of SNPs
    SNP_COUNT = 'SNP_COUNT'
    # overall SNP transition probabilities
    SNP_TRANSITION_COUNT = 'SNP_TRANSITION_COUNT'
    # total count of indels, indexed by length
    INDEL_COUNT = 'INDEL_COUNT'
    # tabulate how much non-N reference sequence we've eaten through
    TOTAL_REFLEN = 'TOTAL_REFLEN'
    # detect variants that occur in a significant percentage of the input samples (pos,ref,alt,pop_fraction)
    COMMON_VARIANTS = 'COMMON_VARIANTS'
    # identify regions that have significantly higher local mutation rates than the average
    HIGH_MUT_REGIONS = 'HIGH_MUT_REGIONS'
    # list to be used for counting variants that occur multiple times in file (i.e. in multiple samples)
    VDAT_COMMON = 'VDAT_COMMON'
    SNP_FREQ = 'SNP_FREQ'
    AVG_INDEL_FREQ = 'AVG_INDEL_FREQ'
    INDEL_FREQ = 'INDEL_FREQ'
    AVG_MUT_RATE = 'AVG_MUT_RATE'


class AnnotatedSeqence:
    _chrom_seqeunces = {}
    _code_to_annotation = {0:Regions.EXON.value, 1:Regions.INTRON.value, 2:Regions.INTERGENIC.value}
    _annotation_to_code = {Regions.EXON.value:0, Regions.INTRON.value:1, Regions.INTERGENIC.value:2}

    def __init__(self, annotations_df: pd.DataFrame):
        if annotations_df is None or annotations_df.empty:
            self._chrom_seqeunces = None
            return

        for i, annotation in annotations_df.iterrows():
            if not annotation['chrom'] in self.chrom_seqeunces:
                self._chrom_seqeunces[annotation['chrom']] = []
            annotation_length = annotation['end'] - annotation['start']
            current_sequence = [self._annotation_to_code[annotation['feature']]] * annotation_length
            self._chrom_seqeunces[annotation['chrom']] = self._chrom_seqeunces[annotation['chrom']] + current_sequence

    def get_annotation(self, chrom, index):
        if not self._chrom_seqeunces:
            return Regions.ALL
        return self._code_to_annotation[self.chrom_seqeunces[chrom][index]]


class RegionStats:
    def __init__(self, annotations_df = None):
        self._regions_stats = {Regions.ALL: self.create_stats_dict()}
        _annotated_sequence = AnnotatedSeqence(None)
        if annotations_df:
            self._regions_stats[Regions.EXON] = self.create_stats_dict()
            self._regions_stats[Regions.INTRON] = self.create_stats_dict()
            self._regions_stats[Regions.INTERGENIC] = self.create_stats_dict()
            _annotated_sequence = AnnotatedSeqence(annotations_df)

    def get_region(self, chrom, index):
        return self._annotated_sequence.get_annotation(chrom, index)

    def get_stat_by_region(self, region_name, stat_name):
        return self._regions_stats[region_name][stat_name]

    def get_stat_by_location(self, chrom, index, stat_name):
        region_name = self.get_region(chrom, index)
        return self.get_stat_by_region(region_name, stat_name)

    def get_all_stats(self):
        return self._regions_stats

    @staticmethod
    def create_stats_dict():
        return {
            Stats.TRINUC_REF_COUNT: {},
            Stats.TRINUC_TRANSITION_COUNT: {},
            Stats.SNP_COUNT: [0],
            Stats.SNP_TRANSITION_COUNT: {},
            Stats.INDEL_COUNT: {},
            Stats.TOTAL_REFLEN: [0],
            Stats.COMMON_VARIANTS: [],
            Stats.HIGH_MUT_REGIONS: []
        }