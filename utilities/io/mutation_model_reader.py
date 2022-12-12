import pickle
import sys
from typing import Union

import numpy as np

from utilities import DiscreteDistribution
from utilities.common_data_structues import Region, ALL_TRI, ALL_IND, TRI_IND, NUC_IND, VALID_NUCL, ModelStats

"""
MUTATION MODEL KEYS
"""
# average mutation rate
MODEL_AVG_MUT_RATE = 'AVG_MUT_RATE'
# p(mut is indel | mut occurs)
MODEL_P_INDEL = 'P_INDEL'
# p(insertion | indel occurs)
MODEL_P_INSERTION = 'P_INSERTION'
# distribution of insertion lengths
MODEL_INS_LEN_DSTRBTN = 'INS_LEN_DSTRBTN'
# distribution of deletion lengths
MODEL_DEL_LEN_DSTRBTN = 'DEL_LEN_DSTRBTN'
# distribution of trinucleotide SNP transitions
MODEL_TRINUC_TRANS_DSTRBTN = 'TRINUC_TRANS_DSTRBTN'
# p(trinuc mutates)
MODEL_P_TRINUC_MUT = 'P_TRINUC_MUT'

"""
More data keys relevant for model building
"""
INS_LENGTH_VALUES = 'INS_LENGTH_VALUES'
INS_LENGTH_WEIGHTS = 'INS_LENGTH_WEIGHTS'
DEL_LENGTH_VALUES = 'DEL_LENGTH_VALUES'
DEL_LENGTH_WEIGHTS = 'DEL_LENGTH_WEIGHTS'
TRINUC_FREQS = 'TRINUC_FREQS'
TRINUC_BIAS = 'TRINUC_BIAS'



def load_mutation_model_from_file(model_file: str = None) -> dict:
    """
    Parse mutation model from a pickle file, initialize and return it
    :param model_file: pickle file path
    :return: mutation model, ready to use
    """
    data_from_file = parse_mutation_model(model_file)
    mutation_model = init_model(data_from_file)
    return mutation_model


def parse_mutation_model(model_file: str = None):
    """
    Parse mutation model pickle file
    """

    mut_model = {region.value: {} for region in Region}
    if model_file is not None:
        pickle_dict = pickle.load(open(model_file, "rb"))
        region_names_found = []

        region_list = list(Region)
        region_list.append('')
        for region in region_list:
            if pickel_key(region, ModelStats.AVG_MUT_RATE) not in pickle_dict:  # checking that the current region exists
                continue

            region_name = region.value if region != '' else Region.ALL.value
            region_names_found.append(region_name)

            mut_model[region_name][MODEL_AVG_MUT_RATE] = pickle_dict[pickel_key(region, ModelStats.AVG_MUT_RATE)]
            mut_model[region_name][MODEL_P_INDEL] = 1. - pickle_dict[pickel_key(region, ModelStats.SNP_FREQ)]

            ins_list = pickle_dict[pickel_key(region, ModelStats.INDEL_FREQ)]
            parse_mutation_model_indel_stats(ins_list, mut_model, region_name)

            trinuc_trans_prob = pickle_dict[pickel_key(region, ModelStats.TRINUC_TRANS_PROBS)]
            parse_mutation_model_trinuc_trans_probs(mut_model, region_name, trinuc_trans_prob)

            trinuc_mut_prob = pickle_dict[pickel_key(region, ModelStats.TRINUC_MUT_PROB)]
            parse_mutation_model_trinuc_mut_probs(mut_model, region_name, trinuc_mut_prob)
        print(f'found the next regions in the model: {region_names_found}')
        for region in Region:
            if region.value not in region_names_found:
                del mut_model[region.value]

    else:
        print('\nError: No mutation model specified\n')
        sys.exit(1)
    return mut_model


def parse_mutation_model_trinuc_mut_probs(mut_model, region_name, trinuc_mut_prob):
    which_have_we_seen = {n: False for n in ALL_TRI}
    trinuc_mean = np.mean(list(trinuc_mut_prob.values()))
    mut_model[region_name][TRINUC_BIAS] = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
    for trinuc in trinuc_mut_prob.keys():
        mut_model[region_name][TRINUC_BIAS][ALL_IND[trinuc]] = trinuc_mut_prob[trinuc]
        which_have_we_seen[trinuc] = True
    for trinuc in which_have_we_seen.keys():
        if not which_have_we_seen[trinuc]:
            mut_model[region_name][TRINUC_BIAS][ALL_IND[trinuc]] = trinuc_mean


def parse_mutation_model_trinuc_trans_probs(mut_model, region_name, trinuc_trans_prob):
    mut_model[region_name][TRINUC_FREQS] = [np.zeros((4, 4)) for _ in range(16)]
    for k in sorted(trinuc_trans_prob.keys()):
        my_ind = TRI_IND[k[0][0] + k[0][2]]
        (k1, k2) = (NUC_IND[k[0][1]], NUC_IND[k[1][1]])
        mut_model[region_name][TRINUC_FREQS][my_ind][k1][k2] = trinuc_trans_prob[k]
    for i in range(len(mut_model[region_name][TRINUC_FREQS])):
        for j in range(len(mut_model[region_name][TRINUC_FREQS][i])):
            for k in range(len(mut_model[region_name][TRINUC_FREQS][i][j])):
                # if trinuc not present in input mutation model, assign it uniform probability
                if float(sum(mut_model[region_name][TRINUC_FREQS][i][j])) < 1e-12:
                    mut_model[region_name][TRINUC_FREQS][i][j] = [0.25, 0.25, 0.25, 0.25]
                else:
                    mut_model[region_name][TRINUC_FREQS][i][j][k] /= float(
                        sum(mut_model[region_name][TRINUC_FREQS][i][j]))


def parse_mutation_model_indel_stats(ins_list, mut_model, region_name):
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
    mut_model[region_name][MODEL_P_INSERTION] = ins_count / float(ins_count + del_count)
    mut_model[region_name][INS_LENGTH_VALUES] = ins_vals
    mut_model[region_name][INS_LENGTH_WEIGHTS] = ins_weight
    mut_model[region_name][DEL_LENGTH_VALUES] = del_vals
    mut_model[region_name][DEL_LENGTH_WEIGHTS] = del_weight


def pickel_key(region: Union[Region, str], stats: ModelStats) -> str:
    return f'{region.value}.{stats.value}' if region != '' else stats.value


def init_model(model_from_file: dict) -> dict:
    """
    Initialize mutation model, ready to use
    """

    model_per_region = {}
    for region_name in model_from_file.keys():
        data = model_from_file[region_name]
        model_per_region[region_name] = {}
        model_per_region[region_name][MODEL_AVG_MUT_RATE] = data[MODEL_AVG_MUT_RATE]
        model_per_region[region_name][MODEL_P_INDEL] = data[MODEL_P_INDEL]
        model_per_region[region_name][MODEL_P_INSERTION] = data[MODEL_P_INSERTION]
        model_per_region[region_name][MODEL_INS_LEN_DSTRBTN] = \
            DiscreteDistribution(data[INS_LENGTH_WEIGHTS], data[INS_LENGTH_VALUES])
        model_per_region[region_name][MODEL_DEL_LEN_DSTRBTN] = \
            DiscreteDistribution(data[DEL_LENGTH_WEIGHTS], data[DEL_LENGTH_VALUES])

        model_per_region[region_name][MODEL_TRINUC_TRANS_DSTRBTN] = []
        for m in data[TRINUC_FREQS]:
            # noinspection PyTypeChecker
            model_per_region[region_name][MODEL_TRINUC_TRANS_DSTRBTN].append(
                [DiscreteDistribution(m[0], VALID_NUCL), DiscreteDistribution(m[1], VALID_NUCL),
                 DiscreteDistribution(m[2], VALID_NUCL), DiscreteDistribution(m[3], VALID_NUCL)])
        model_per_region[region_name][MODEL_P_TRINUC_MUT] = [m for m in data[TRINUC_BIAS]]

    return model_per_region