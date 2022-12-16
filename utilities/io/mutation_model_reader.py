import pickle
import sys
from typing import Union, List

import numpy as np

from utilities import DiscreteDistribution
from utilities.common_data_structues import Region, ALL_TRI, ALL_IND, TRI_IND, NUC_IND, VALID_NUCL, \
    ModelStats, ModelKeys


def load_model_from_file(model_file: str = None) -> dict:
    """
    Parse mutation model from a pickle file, initialize and return it
    :param model_file: pickle file path
    :return: mutation model, ready to use
    """
    mut_model = {region.value: {} for region in Region}
    if model_file is not None:
        pickle_dict = pickle.load(open(model_file, "rb"))
        region_names_found = []

        region_list = list(Region)
        region_list.append('')
        for region in region_list:
            # checking that the current region exists:
            if pickel_key(region, ModelStats.AVG_MUT_RATE) not in pickle_dict:
                continue

            region_name = region.value if region != '' else Region.ALL.value
            region_names_found.append(region_name)

            mut_model[region_name][ModelKeys.AVG_MUT_RATE] = pickle_dict[pickel_key(region, ModelStats.AVG_MUT_RATE)]
            mut_model[region_name][ModelKeys.P_INDEL] = 1. - pickle_dict[pickel_key(region, ModelStats.SNP_FREQ)]

            ins_list = pickle_dict[pickel_key(region, ModelStats.INDEL_FREQ)]
            parse_indel_stats(mut_model, region_name, ins_list)

            trinuc_trans_prob = pickle_dict[pickel_key(region, ModelStats.TRINUC_TRANS_PROBS)]
            trinuc_mut_prob = pickle_dict[pickel_key(region, ModelStats.TRINUC_MUT_PROB)]
            parse_trinuc_stats(mut_model, region_name, trinuc_mut_prob, trinuc_trans_prob)

        print(f'found the next regions in the model: {region_names_found}')
        for region in Region:
            if region.value not in region_names_found:
                del mut_model[region.value]

    else:
        print('\nError: No mutation model specified\n')
        sys.exit(1)
    return mut_model


def parse_indel_stats(mut_model, region_name, ins_list):
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
    mut_model[region_name][ModelKeys.P_INSERTION] = ins_count / float(ins_count + del_count)
    mut_model[region_name][ModelKeys.INS_LEN_DSTRBTN] = \
        DiscreteDistribution(ins_weight, ins_vals)
    mut_model[region_name][ModelKeys.DEL_LEN_DSTRBTN] = \
        DiscreteDistribution(del_weight, del_vals)


def parse_trinuc_stats(mut_model, region_name, region_trinuc_mut_prob: dict, region_trinuc_trans_prob: dict):
    region_trinuc_freqs = parse_trinuc_trans_probs(region_trinuc_trans_prob)
    region_trinuc_bias = parse_trinuc_mut_probs(region_trinuc_mut_prob)
    mut_model[region_name][ModelKeys.TRINUC_TRANS_DSTRBTN] = []
    for m in region_trinuc_freqs:
        mut_model[region_name][ModelKeys.TRINUC_TRANS_DSTRBTN].append(
            [DiscreteDistribution(m[0], VALID_NUCL), DiscreteDistribution(m[1], VALID_NUCL),
             DiscreteDistribution(m[2], VALID_NUCL), DiscreteDistribution(m[3], VALID_NUCL)])
    mut_model[region_name][ModelKeys.P_TRINUC_MUT] = [m for m in region_trinuc_bias]


def parse_trinuc_mut_probs(region_trinuc_mut_prob: dict) -> List:
    which_have_we_seen = {n: False for n in ALL_TRI}
    trinuc_mean = np.mean(list(region_trinuc_mut_prob.values()))
    region_trinuc_bias = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
    for trinuc in region_trinuc_mut_prob.keys():
        region_trinuc_bias[ALL_IND[trinuc]] = region_trinuc_mut_prob[trinuc]
        which_have_we_seen[trinuc] = True
    for trinuc in which_have_we_seen.keys():
        if not which_have_we_seen[trinuc]:
            region_trinuc_bias[ALL_IND[trinuc]] = trinuc_mean
    return region_trinuc_bias


def parse_trinuc_trans_probs(region_trinuc_trans_prob: dict) -> List:
    region_trinuc_freqs = [np.zeros((4, 4)) for _ in range(16)]
    for k in sorted(region_trinuc_trans_prob.keys()):
        my_ind = TRI_IND[k[0][0] + k[0][2]]
        (k1, k2) = (NUC_IND[k[0][1]], NUC_IND[k[1][1]])
        region_trinuc_freqs[my_ind][k1][k2] = region_trinuc_trans_prob[k]
    for i in range(len(region_trinuc_freqs)):
        for j in range(len(region_trinuc_freqs[i])):
            for k in range(len(region_trinuc_freqs[i][j])):
                # if trinuc not present in input mutation model, assign it uniform probability
                if float(sum(region_trinuc_freqs[i][j])) < 1e-12:
                    region_trinuc_freqs[i][j] = [0.25, 0.25, 0.25, 0.25]
                else:
                    region_trinuc_freqs[i][j][k] /= float(
                        sum(region_trinuc_freqs[i][j]))
    return region_trinuc_freqs


def pickel_key(region: Union[Region, str], stats: ModelStats) -> str:
    return f'{region.value}.{stats.value}' if region != '' else stats.value
