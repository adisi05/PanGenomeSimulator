#!/usr/bin/env source

# Python 3 ready

import os
import re
import pickle
import argparse
import numpy as np
from Bio import SeqIO
import pandas as pd
from enum import Enum

from utilities.common_data_structues import Region, VALID_NUCL, VALID_TRINUC
from utilities.annotated_sequence import AnnotatedSequence
from utilities.genome_annotations import read_annotations_csv

# if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
VCF_DEFAULT_POP_FREQ = 0.00001
VCF_CHROM_COL = 0


# TODO move to a different file?
class Stats(Enum):
    # how many times do we observe each trinucleotide in the reference (and input bed region, if present)?
    TRINUC_REF_COUNT = 'TRINUC_REF_COUNT'
    # [(trinuc_ref, trinuc_alt)] = # of times we observed a mutation from trinuc_ref into trinuc_alt
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
    # TODO write explanation
    SNP_FREQ = 'SNP_FREQ'
    # TODO write explanation
    AVG_INDEL_FREQ = 'AVG_INDEL_FREQ'
    # TODO write explanation
    INDEL_FREQ = 'INDEL_FREQ'
    # TODO write explanation
    AVG_MUT_RATE = 'AVG_MUT_RATE'
    # frequency of snp transitions, given a snp occurs.
    SNP_TRANS_FREQ = 'SNP_TRANS_FREQ'
    # frequency that each trinuc mutated into anything else
    TRINUC_MUT_PROB = 'TRINUC_MUT_PROB'
    # frequency that a trinuc mutates into another trinuc, given that it mutated
    TRINUC_TRANS_PROBS = 'TRINUC_TRANS_PROBS'


class RegionStats:
    def __init__(self, annotations_df=None):
        self._regions_stats = {Region.ALL.value: self.create_stats_dict()}
        self._annotated_sequence_per_chrom = {}
        if annotations_df is not None and not annotations_df.empty:
            self._regions_stats[Region.CDS.value] = self.create_stats_dict()
            self._regions_stats[Region.NON_CODING_GENE.value] = self.create_stats_dict()
            self._regions_stats[Region.INTERGENIC.value] = self.create_stats_dict()
            for chrom in annotations_df.chrom.unique():
                self._annotated_sequence_per_chrom[chrom] = AnnotatedSequence(annotations_df, chrom)

    def get_region_by_chrom_and_pos(self, chrom: str, pos: int) -> Region:
        region, _ = self._annotated_sequence_per_chrom[chrom].get_region_by_position(pos)
        return region

    def get_stat_by_region(self, region: Region, stat: Stats):
        return self._regions_stats[region.value][stat.value]

    def get_stat_by_location(self, chrom: str, pos: int, stat: Stats):
        region = self.get_region_by_chrom_and_pos(chrom, pos)
        return self.get_stat_by_region(region, stat)

    def get_all_stats(self):
        return self._regions_stats

    @staticmethod
    def create_stats_dict():
        return {
            Stats.TRINUC_REF_COUNT.value: {},
            Stats.TRINUC_TRANSITION_COUNT.value: {},
            Stats.SNP_COUNT.value: [0],
            Stats.SNP_TRANSITION_COUNT.value: {},
            Stats.INDEL_COUNT.value: {},
            Stats.TOTAL_REFLEN.value: [0],
            Stats.COMMON_VARIANTS.value: [],
            Stats.HIGH_MUT_REGIONS.value: [],
            Stats.VDAT_COMMON.value: {},
            Stats.TRINUC_MUT_PROB.value: {},
            Stats.TRINUC_TRANS_PROBS.value: {},
            Stats.SNP_TRANS_FREQ.value: {}
        }

#####################################
#				main()				#
#####################################


def main():
    args = parse_arguments()
    (ref, vcf, out_pickle, save_trinuc, skip_common, annotations_file, workdir) = (
        args.r, args.m, args.o, args.save_trinuc, args.skip_common, args.b, args.w)

    chrom_sequences_dict, chrom_names = process_reference(ref)

    indices_to_indels, matching_chromosomes, matching_variants = process_vcf(chrom_names, vcf)

    if not workdir:
        workdir = os.path.dirname(out_pickle)
    annotations_df = read_annotations_csv(annotations_file, workdir)

    regions_stats = RegionStats(annotations_df)

    process_trinuc_counts(chrom_sequences_dict, matching_chromosomes, annotations_df, regions_stats)

    # Load and process variants in each reference sequence individually, for memory reasons...
    print('Creating mutational model...')
    for chrom_name in matching_chromosomes:
        collect_basic_stats_for_chrom(chrom_sequences_dict, chrom_name, indices_to_indels, matching_variants,
                                      annotations_df, regions_stats)

    compute_probabilities(regions_stats)

    save_stats_to_file(out_pickle, skip_common, regions_stats)


#########################################################
#				VARIOUS HELPER FUNCTIONS				#
#########################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description='gen_mut_model.source',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-r', type=str, required=True, metavar='/path/to/reference.fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-m', type=str, required=True, metavar='/path/to/mutations.vcf',
                        help="Mutation file for organism in VCF format")
    parser.add_argument('-o', type=str, required=True, metavar='/path/to/output/and/prefix',
                        help="Name of output file (final model will append \'.p\')")
    parser.add_argument('-b', type=str, required=False, metavar='Bed file of regions to include '
                                                                '(use bedtools complement if you have a '
                                                                'bed of exclusion areas)', default=None,
                        help="only_use_these_regions.bed")
    parser.add_argument('-w', type=str, required=False, metavar='/path/to/working/dir/', default=None,
                        help="Name of working directory to process annotations. If not given, the directory of the "
                             "output file wll be taken instead")
    parser.add_argument('--save-trinuc', required=False, action='store_true', default=False,
                        help='save trinucleotide counts for reference')
    parser.add_argument('--skip-common', required=False, action='store_true', default=False,
                        help='Do not save common snps + high mut regions')
    args = parser.parse_args()
    return args


def process_reference(ref):
    # Process reference file
    print('Processing reference...')
    reference = None
    try:
        reference = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
    except ValueError:
        print("Problems parsing reference file. Ensure reference is in proper fasta format")
        exit(1)
    chrom_sequences_dict = {}
    for chrom in reference.keys():
        chrom_split = chrom.split("|")
        chrom_sequences_dict[chrom_split[0]] = reference[chrom]
    chrom_names = list(map(str, chrom_sequences_dict.keys()))
    return chrom_sequences_dict, chrom_names


def process_vcf(chrom_list, vcf):
    # Process VCF file. First check if it's been entered as a TSV
    if vcf[-3:] == 'tsv':
        print("Warning! TSV file must follow VCF specifications.")
    # Pre-parsing to find all the matching chromosomes between ref and vcf
    print('Processing VCF file...')
    variants = []
    try:
        variants = pd.read_csv(vcf, sep='\t', comment='#', index_col=None, header=None)
        variants[VCF_CHROM_COL] = variants[VCF_CHROM_COL].map(str)
    except ValueError:
        print("VCF must be in standard VCF format with tab-separated columns")
    # Narrow chromosomes to those matching the reference
    # This is in part to make sure the names match
    variant_chroms = variants[VCF_CHROM_COL].unique() if not variants.empty else []
    matching_chromosomes = []
    for chrom_name in chrom_list:
        if chrom_name not in variant_chroms:
            continue
        else:
            matching_chromosomes.append(chrom_name)
    # Check to make sure there are some matches
    if not matching_chromosomes:
        print("Found no chromosomes in common between VCF and Fasta. Please fix the chromosome names and try again")
        exit(1)
    # Double check that there are matches
    matching_variants = []
    try:
        matching_variants = variants[variants[VCF_CHROM_COL].isin(matching_chromosomes)]
    except ValueError:
        print("Problem matching variants with reference.")
    if len(matching_variants) == 0:
        print("There is no overlap between reference and variant file. This could be a chromosome naming problem")
        exit(1)
    # Rename header in dataframe for processing
    matching_variants = matching_variants.rename(columns={0: 'CHROM', 1: 'chr_start', 2: 'ID', 3: 'REF', 4: 'ALT',
                                                          5: 'QUAL', 6: 'FILTER', 7: 'INFO'})
    # Change the indexing by -1 to match source format indexing (0-based)
    matching_variants['chr_start'] = matching_variants['chr_start'] - 1
    matching_variants['chr_end'] = matching_variants['chr_start']
    # Process the variant table
    indices_to_indels = \
        matching_variants.loc[matching_variants.ALT.apply(len) != matching_variants.REF.apply(len)].index
    # indels in vcf don't include the preserved first nucleotide, so lets trim the vcf alleles
    ref_values_to_change = matching_variants.loc[indices_to_indels, 'REF'].copy().str[1:]
    alt_values_to_change = matching_variants.loc[indices_to_indels, 'ALT'].copy().str[1:]
    matching_variants.loc[indices_to_indels, "REF"] = ref_values_to_change
    matching_variants.loc[indices_to_indels, "ALT"] = alt_values_to_change
    matching_variants.replace('', '-', inplace=True)
    # If multi-alternate alleles are present, lets just ignore this variant. I may come back and improve this later
    indices_to_ignore = matching_variants[matching_variants['ALT'].str.contains(',')].index
    matching_variants = matching_variants.drop(indices_to_ignore)
    # if we encounter a multi-np (i.e. 3 nucl --> 3 different nucl), let's skip it for now...
    # Alt and Ref contain no dashes
    no_dashes = matching_variants[
        ~matching_variants['REF'].str.contains('-') & ~matching_variants['ALT'].str.contains('-')].index
    # Alt and Ref lengths are greater than 1
    long_variants = matching_variants[
        (matching_variants['REF'].apply(len) > 1) & (matching_variants['ALT'].apply(len) > 1)].index
    complex_variants = list(set(no_dashes) & set(long_variants))
    matching_variants = matching_variants.drop(complex_variants)
    # This is solely to make regex easier later, since we can't predict where in the line a string will be
    new_info = ';' + matching_variants['INFO'].copy() + ';'
    matching_variants['INFO'] = new_info
    return indices_to_indels, matching_chromosomes, matching_variants  # TODO list(indices_to_indels) ?


def cluster_list(variants_position_list_sorted: list, delta: float) -> list:
    """
    Clusters a sorted list
    :param variants_position_list_sorted: a sorted list
    :param delta: the value to compare list items to
    :return: a clustered list of values
    """
    out_list = [[variants_position_list_sorted[0]]]
    previous_value = variants_position_list_sorted[0]
    current_index = 0
    for item in variants_position_list_sorted[1:]:
        if item - previous_value <= delta:
            out_list[current_index].append(item)
        else:
            current_index += 1
            out_list.append([])
            out_list[current_index].append(item)
        previous_value = item
    return out_list


def process_trinuc_counts(chrom_sequences_dict, matching_chromosomes, annotations_df, regions_stats):
    for chrom_name in matching_chromosomes:

        # using genome annotations
        if annotations_df is not None and not annotations_df.empty:
            chrom_annotations = annotations_df[annotations_df['chrom'] == chrom_name]
            for i, annotation in chrom_annotations.iterrows():
                sub_seq = chrom_sequences_dict[chrom_name][annotation['start']: annotation['end']].seq
                region = regions_stats.get_region_by_chrom_and_pos(chrom_name, annotation['start'])
                update_trinuc_ref_count(sub_seq, regions_stats, region)

        # without genome annotations
        else:
            sub_seq = chrom_sequences_dict[chrom_name].seq
            update_trinuc_ref_count(sub_seq, regions_stats)


def collect_basic_stats_for_chrom(chrom_sequences_dict, chrom_name, indices_to_indels, matching_variants,
                                  annotations_df,
                                  regions_stats):
    update_total_reflen(chrom_sequences_dict, chrom_name, regions_stats, annotations_df)

    # Create a view that narrows variants list to current chromosome
    variants_to_process = matching_variants[matching_variants['CHROM'] == chrom_name].copy()

    # process SNPs
    snp_df = variants_to_process[~variants_to_process.index.isin(indices_to_indels)]
    # no '-' characters allowed, and chrStart must be same as chrEnd:
    snp_df = snp_df.loc[snp_df['chr_start'] == snp_df['chr_end']]
    process_snps(snp_df, chrom_name, str(chrom_sequences_dict[chrom_name].seq), regions_stats)

    # process indels
    indel_df = variants_to_process[variants_to_process.index.isin(indices_to_indels)]
    process_indels(indel_df, chrom_name, regions_stats)

    process_common_variants(chrom_sequences_dict, chrom_name, regions_stats)


def process_snps(snp_df, chrom_name, chrom_sequence, regions_stats):
    if not snp_df.empty:
        # only consider positions where ref allele in vcf matches the nucleotide in our reference
        for index, snp in snp_df.iterrows():
            trinuc_to_analyze = str(chrom_sequence[snp.chr_start - 1: snp.chr_start + 2])
            if trinuc_to_analyze not in VALID_TRINUC:
                continue
            if snp.REF == trinuc_to_analyze[1]:
                trinuc_ref = trinuc_to_analyze
                trinuc_alt = trinuc_to_analyze[0] + snp_df.loc[index, 'ALT'] + trinuc_to_analyze[2]
                if trinuc_alt not in VALID_TRINUC:
                    continue
                region = regions_stats.get_region_by_chrom_and_pos(chrom_name, snp.chr_start)
                update_trinuc_transition_count(trinuc_ref, trinuc_alt, regions_stats, region)
                update_snp_transition_count(str(snp.REF), str(snp.ALT), regions_stats, region)
                update_snp_count(regions_stats, region)

                pop_freq = VCF_DEFAULT_POP_FREQ
                if ';CAF=' in snp_df.loc[index, 'INFO']:
                    caf_str = re.findall(r";CAF=.*?(?=;)", snp.INFO)[0]
                    if ',' in caf_str:
                        pop_freq = float(caf_str[5:].split(',')[1])
                update_vdat_common(chrom_name, snp, pop_freq, regions_stats, region)
            else:
                print('\nError: ref allele in variant call does not match reference.\n')
                exit(1)


def process_indels(indel_df, chrom_name, regions_stats):
    if not indel_df.empty:
        for index, indel in indel_df.iterrows():
            if "-" in indel.REF:
                len_ref = 0
            else:
                len_ref = len(indel.REF)
            if "-" in indel.ALT:
                len_alt = 0
            else:
                len_alt = len(indel.ALT)
            if len_ref != len_alt:
                indel_len = len_alt - len_ref
                region = regions_stats.get_region_by_chrom_and_pos(chrom_name, indel.chr_start)
                update_indel_count(indel_len, regions_stats, region)

                pop_freq = VCF_DEFAULT_POP_FREQ
                if ';CAF=' in indel.INFO:
                    caf_str = re.findall(r";CAF=.*?(?=;)", indel.INFO)[0]
                    if ',' in caf_str:
                        pop_freq = float(caf_str[5:].split(',')[1])
                update_vdat_common(chrom_name, indel, pop_freq, regions_stats, region)


def process_common_variants(chrom_sequences_dict, chrom_name, regions_stats):
    chr_start_idx = 0
    ref_idx1 = 1
    ref_idx2 = 2
    alt_idx = 3
    pop_freq_idx = 4

    for region_name, region_stats in regions_stats.get_all_stats().items():
        vdat_common_per_chrom = region_stats[Stats.VDAT_COMMON.value]
        if chrom_name not in vdat_common_per_chrom:
            vdat_common_per_chrom[chrom_name] = []
        vdat_common = vdat_common_per_chrom[chrom_name]
        if not len(vdat_common):
            print(f'Found no variants for chromosome {chrom_name} in region {region_name}.')
            continue
        # identify common mutations
        percentile_var = 95
        min_value = np.percentile([variant[pop_freq_idx] for variant in vdat_common], percentile_var)
        for variant in sorted(vdat_common):
            if variant[4] >= min_value:
                region_stats[Stats.COMMON_VARIANTS.value].append(
                    (chrom_name, variant[chr_start_idx], variant[ref_idx1], variant[alt_idx], variant[pop_freq_idx]))
        vdat_common = {(variant[chr_start_idx], variant[ref_idx1], variant[ref_idx2], variant[alt_idx]):
                           variant[pop_freq_idx]
                       for variant in vdat_common}
        # identify areas that have contained significantly higher random mutation rates
        dist_thresh = 2000
        percentile_clust = 97
        scaler = 1000
        # identify regions with disproportionately more variants in them
        variants_position_list_sorted = sorted([variant[chr_start_idx] for variant in vdat_common.keys()])
        variants_clustered_by_position = cluster_list(variants_position_list_sorted, dist_thresh)
        variants_clustered_by_len_and_pos = \
            [(len(variants_clustered_by_position[i]), min(variants_clustered_by_position[i]),
              max(variants_clustered_by_position[i]), i)
             for i in range(len(variants_clustered_by_position))]

        candidate_regions = []
        for cnddt in variants_clustered_by_len_and_pos:
            count_idx = 0
            min_indx = 1
            max_idx = 2
            bi = int((cnddt[min_indx] - dist_thresh) / float(scaler)) * scaler
            bf = int((cnddt[max_idx] + dist_thresh) / float(scaler)) * scaler
            candidate_regions.append(
                (cnddt[count_idx] / float(bf - bi), max([0, bi]), min([len(chrom_sequences_dict[chrom_name]), bf])))
        minimum_value = np.percentile([cand[count_idx] for cand in candidate_regions], percentile_clust)
        for cnddt in candidate_regions:
            if cnddt[0] >= minimum_value:
                region_stats[Stats.HIGH_MUT_REGIONS.value].append((chrom_name, cnddt[1], cnddt[2], cnddt[0]))
        # collapse overlapping regions
        for i in range(len(region_stats[Stats.HIGH_MUT_REGIONS.value]) - 1, 0, -1):
            if region_stats[Stats.HIGH_MUT_REGIONS.value][i - 1][2] >= \
                    region_stats[Stats.HIGH_MUT_REGIONS.value][i][1] and \
                    region_stats[Stats.HIGH_MUT_REGIONS.value][i - 1][0] == \
                    region_stats[Stats.HIGH_MUT_REGIONS.value][i][0]:
                # Might need to research a more accurate way to get the mutation rate for this region
                avg_mut_rate = 0.5 * region_stats[Stats.HIGH_MUT_REGIONS.value][i - 1][3] + \
                               0.5 * region_stats[Stats.HIGH_MUT_REGIONS.value][i][3]
                region_stats[Stats.HIGH_MUT_REGIONS.value][i - 1] = (
                    region_stats[Stats.HIGH_MUT_REGIONS.value][i - 1][0],
                    region_stats[Stats.HIGH_MUT_REGIONS.value][i - 1][1],
                    region_stats[Stats.HIGH_MUT_REGIONS.value][i][2],
                    avg_mut_rate)
                del region_stats[Stats.HIGH_MUT_REGIONS.value][i]


def compute_probabilities(regions_stats):
    # if for some reason we didn't find any valid input variants AT ALL, exit gracefully...
    total_var = regions_stats.get_stat_by_region(Region.ALL, Stats.SNP_COUNT)[0] + \
                sum(regions_stats.get_stat_by_region(Region.ALL, Stats.INDEL_COUNT).values())
    if total_var == 0:
        print(
            '\nError: No valid variants were found, model could not be created. '
            '(Are you using the correct reference?)\n')
        exit(1)

    for region_name, region_stats in regions_stats.get_all_stats().items():
        for trinuc in sorted(region_stats[Stats.TRINUC_REF_COUNT.value].keys()):
            my_count = 0
            for k in sorted(region_stats[Stats.TRINUC_TRANSITION_COUNT.value].keys()):
                if k[0] == trinuc:
                    my_count += region_stats[Stats.TRINUC_TRANSITION_COUNT.value][k]
            region_stats[Stats.TRINUC_MUT_PROB.value][trinuc] = \
                my_count / float(region_stats[Stats.TRINUC_REF_COUNT.value][trinuc])
            for k in sorted(region_stats[Stats.TRINUC_TRANSITION_COUNT.value].keys()):
                if k[0] == trinuc:
                    region_stats[Stats.TRINUC_TRANS_PROBS.value][k] = \
                        region_stats[Stats.TRINUC_TRANSITION_COUNT.value][k] / float(my_count)

        for n1 in VALID_NUCL:
            rolling_tot = sum([region_stats[Stats.SNP_TRANSITION_COUNT.value][(n1, n2)] for n2 in VALID_NUCL if (n1, n2)
                               in region_stats[Stats.SNP_TRANSITION_COUNT.value]])
            for n2 in VALID_NUCL:
                key2 = (n1, n2)
                if key2 in region_stats[Stats.SNP_TRANSITION_COUNT.value]:
                    region_stats[Stats.SNP_TRANS_FREQ.value][key2] = \
                        region_stats[Stats.SNP_TRANSITION_COUNT.value][key2] / float(rolling_tot)

        # compute average snp and indel frequencies
        total_var = region_stats[Stats.SNP_COUNT.value][0] + sum(region_stats[Stats.INDEL_COUNT.value].values())
        if total_var != 0:
            region_stats[Stats.SNP_FREQ.value] = region_stats[Stats.SNP_COUNT.value][0] / float(total_var)
            region_stats[Stats.AVG_INDEL_FREQ.value] = 1. - region_stats[Stats.SNP_FREQ.value]
            region_stats[Stats.INDEL_FREQ.value] = \
                {k: (region_stats[Stats.INDEL_COUNT.value][k] / float(total_var)) /
                    region_stats[Stats.AVG_INDEL_FREQ.value]
                 for k in region_stats[Stats.INDEL_COUNT.value].keys()}
        else:
            region_stats[Stats.SNP_FREQ.value] = 0.
            region_stats[Stats.AVG_INDEL_FREQ.value] = 1.
            region_stats[Stats.INDEL_FREQ.value] = {}
        region_stats[Stats.AVG_MUT_RATE.value] = total_var / region_stats[Stats.TOTAL_REFLEN.value][0]

        # if values weren't found in data, appropriately append null entries
        print_trinuc_warning = False
        for trinuc in VALID_TRINUC:
            trinuc_mut = [trinuc[0] + n + trinuc[2] for n in VALID_NUCL if n != trinuc[1]]
            if trinuc not in region_stats[Stats.TRINUC_MUT_PROB.value]:
                region_stats[Stats.TRINUC_MUT_PROB.value][trinuc] = 0.
                print_trinuc_warning = True
            for trinuc2 in trinuc_mut:
                if (trinuc, trinuc2) not in region_stats[Stats.TRINUC_TRANS_PROBS.value]:
                    region_stats[Stats.TRINUC_TRANS_PROBS.value][(trinuc, trinuc2)] = 0.
                    print_trinuc_warning = True
        if print_trinuc_warning:
            print(
                'Warning: Some trinucleotides transitions were not encountered in the input dataset, '
                'probabilities of 0.0 have been assigned to these events.')

        #
        # print some stuff
        #
        print(f'Probabilities for region {region_name}:')

        for k in sorted(region_stats[Stats.TRINUC_MUT_PROB.value].keys()):
            print('p(' + k + ' mutates) =', region_stats[Stats.TRINUC_MUT_PROB.value][k])

        for k in sorted(region_stats[Stats.TRINUC_TRANS_PROBS.value].keys()):
            print('p(' + k[0] + ' --> ' + k[1] + ' | ' + k[0] + ' mutates) =',
                  region_stats[Stats.TRINUC_TRANS_PROBS.value][k])

        for k in sorted(region_stats[Stats.INDEL_FREQ.value].keys()):
            if k > 0:
                print('p(ins length = ' + str(abs(k)) + ' | indel occurs) =',
                      region_stats[Stats.INDEL_FREQ.value][k])
            else:
                print('p(del length = ' + str(abs(k)) + ' | indel occurs) =',
                      region_stats[Stats.INDEL_FREQ.value][k])

        for k in sorted(region_stats[Stats.SNP_TRANS_FREQ.value].keys()):
            print('p(' + k[0] + ' --> ' + k[1] + ' | SNP occurs) =', region_stats[Stats.SNP_TRANS_FREQ.value][k])

        print('p(snp)   =', region_stats[Stats.SNP_FREQ.value])
        print('p(indel) =', region_stats[Stats.AVG_INDEL_FREQ.value])
        print('overall average mut rate:', region_stats[Stats.AVG_MUT_RATE.value])
        print('total variants processed:', total_var)


def update_trinuc_ref_count(sub_seq, regions_stats, region: Region = Region.ALL, also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    for trinuc in VALID_TRINUC:
        for current_region in regions_to_update:
            if trinuc not in regions_stats.get_stat_by_region(current_region, Stats.TRINUC_REF_COUNT):
                regions_stats.get_stat_by_region(current_region, Stats.TRINUC_REF_COUNT)[trinuc] = 0
            regions_stats.get_stat_by_region(current_region, Stats.TRINUC_REF_COUNT)[trinuc] += sub_seq.count_overlap(
                trinuc)


def update_total_reflen(chrom_sequences_dict, chrom_name, regions_stats, annotations_df):
    # Count the number of non-N nucleotides for the reference
    total_reflen_all_chrom = len(chrom_sequences_dict[chrom_name].seq) - chrom_sequences_dict[chrom_name].seq.count('N')
    update_total_reflen_for_region(total_reflen_all_chrom, regions_stats, Region.ALL)
    chrom_annotations = annotations_df[annotations_df['chrom'] == chrom_name]
    if annotations_df is not None and not annotations_df.empty:
        for i, annotation in chrom_annotations.iterrows():
            sub_seq = chrom_sequences_dict[chrom_name][annotation['start']: annotation['end']].seq
            reflen_annotation = len(sub_seq) - sub_seq.count('N')
            update_total_reflen_for_region(reflen_annotation, regions_stats, Region(annotation['region']))


def update_total_reflen_for_region(total_reflen_region_chrom, regions_stats, region: Region = Region.ALL):
    total_reflen_region = regions_stats.get_stat_by_region(region, Stats.TOTAL_REFLEN)
    total_reflen_region[0] += total_reflen_region_chrom


def update_trinuc_transition_count(trinuc_ref, trinuc_alt, regions_stats, region: Region = Region.ALL,
                                   also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    key = (trinuc_ref, trinuc_alt)
    for current_region in regions_to_update:
        if key not in regions_stats.get_stat_by_region(current_region, Stats.TRINUC_TRANSITION_COUNT):
            regions_stats.get_stat_by_region(current_region, Stats.TRINUC_TRANSITION_COUNT)[key] = 0
        regions_stats.get_stat_by_region(current_region, Stats.TRINUC_TRANSITION_COUNT)[key] += 1


def update_snp_transition_count(ref, alt, regions_stats, region: Region = Region.ALL, also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    key2 = (ref, alt)
    for current_region in regions_to_update:
        if key2 not in regions_stats.get_stat_by_region(current_region, Stats.SNP_TRANSITION_COUNT):
            regions_stats.get_stat_by_region(current_region, Stats.SNP_TRANSITION_COUNT)[key2] = 0
        regions_stats.get_stat_by_region(current_region, Stats.SNP_TRANSITION_COUNT)[key2] += 1


def update_snp_count(regions_stats, region: Region = Region.ALL, also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    for current_region in regions_to_update:
        regions_stats.get_stat_by_region(current_region, Stats.SNP_COUNT)[0] += 1


def update_indel_count(indel_len, regions_stats, region: Region = Region.ALL, also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    for current_region in regions_to_update:
        if indel_len not in regions_stats.get_stat_by_region(current_region, Stats.INDEL_COUNT):
            regions_stats.get_stat_by_region(current_region, Stats.INDEL_COUNT)[indel_len] = 0
        regions_stats.get_stat_by_region(current_region, Stats.INDEL_COUNT)[indel_len] += 1


def update_vdat_common(chrom_name, variant, pop_freq, regions_stats, region: Region = Region.ALL,
                       also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    for current_region in regions_to_update:
        vdat_common_per_chrom = regions_stats.get_stat_by_region(current_region, Stats.VDAT_COMMON)
        if chrom_name not in vdat_common_per_chrom:
            vdat_common_per_chrom[chrom_name] = []
        vdat_common_per_chrom[chrom_name].append((variant.chr_start, variant.REF, variant.REF, variant.ALT, pop_freq))


def save_stats_to_file(out_pickle, skip_common, regions_stats):
    out_dict = {}
    for region_name, region_stats in regions_stats.get_all_stats().items():
        out_dict[f'{region_name}.{Stats.AVG_MUT_RATE.value}'] = region_stats[Stats.AVG_MUT_RATE.value]
        out_dict[f'{region_name}.{Stats.TOTAL_REFLEN.value}'] = region_stats[Stats.TOTAL_REFLEN.value][0]
        out_dict[f'{region_name}.{Stats.SNP_FREQ.value}'] = region_stats[Stats.SNP_FREQ.value]
        out_dict[f'{region_name}.{Stats.SNP_TRANS_FREQ.value}'] = region_stats[Stats.SNP_TRANS_FREQ.value]
        out_dict[f'{region_name}.{Stats.INDEL_FREQ.value}'] = region_stats[Stats.INDEL_FREQ.value]
        out_dict[f'{region_name}.{Stats.TRINUC_MUT_PROB.value}'] = region_stats[Stats.TRINUC_MUT_PROB.value]
        out_dict[f'{region_name}.{Stats.TRINUC_TRANS_PROBS.value}'] = region_stats[Stats.TRINUC_TRANS_PROBS.value]
        if not skip_common:
            out_dict[f'{region_name}.{Stats.COMMON_VARIANTS.value}'] = region_stats[Stats.COMMON_VARIANTS.value]
            out_dict[f'{region_name}.{Stats.HIGH_MUT_REGIONS.value}'] = region_stats[Stats.HIGH_MUT_REGIONS.value]
    pickle.dump(out_dict, open(out_pickle, "wb"))


if __name__ == "__main__":
    main()
