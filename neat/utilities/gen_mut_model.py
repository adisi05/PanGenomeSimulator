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
from neat.utilities.compute_annotation_regions import seperate_exons_genes_intergenics

# Some constants we'll need later
VALID_NUCL = ['A', 'C', 'G', 'T']
VALID_TRINUC = [VALID_NUCL[i] + VALID_NUCL[j] + VALID_NUCL[k] for i in range(len(VALID_NUCL)) for j in
                range(len(VALID_NUCL)) for k in range(len(VALID_NUCL))]
# if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
VCF_DEFAULT_POP_FREQ = 0.00001

#########################################################
#				VARIOUS HELPER FUNCTIONS				#
#########################################################

def cluster_list(list_to_cluster: list, delta: float) -> list:
    """
    Clusters a sorted list
    :param list_to_cluster: a sorted list
    :param delta: the value to compare list items to
    :return: a clustered list of values
    """
    out_list = [[list_to_cluster[0]]]
    previous_value = list_to_cluster[0]
    current_index = 0
    for item in list_to_cluster[1:]:
        if item - previous_value <= delta:
            out_list[current_index].append(item)
        else:
            current_index += 1
            out_list.append([])
            out_list[current_index].append(item)
        previous_value = item
    return out_list


#####################################
#			Data Structures			#
#####################################

class Regions(Enum):
    EXON = 'exon'
    INTRON = 'intron'
    INTERGENIC = 'intergenic'
    ALL = 'all'

class BasicStats(Enum):
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

class AdditionalStats(Enum):
    SNP_FREQ = 'SNP_FREQ'
    AVG_INDEL_FREQ = 'AVG_INDEL_FREQ'
    INDEL_FREQ = 'INDEL_FREQ'
    AVG_MUT_RATE = 'AVG_MUT_RATE'

class AnnotatedSeqence:
    _chrom_seqeunces = {}
    _code_to_annotation = {0:Regions.EXON, 1:Regions.INTRON, 2:Regions.INTERGENIC}
    _annotation_to_code = {Regions.EXON:0, Regions.INTRON:1, Regions.INTERGENIC:2}

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
        self._regions_stats = {Regions.ALL: Stats()}
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
            BasicStats.TRINUC_REF_COUNT: {},
            BasicStats.TRINUC_TRANSITION_COUNT: {},
            BasicStats.SNP_COUNT: [0],
            BasicStats.SNP_TRANSITION_COUNT: {},
            BasicStats.INDEL_COUNT: {},
            BasicStats.TOTAL_REFLEN: [0],
            BasicStats.COMMON_VARIANTS: [],
            BasicStats.HIGH_MUT_REGIONS: []
        }

#####################################
#				main()				#
#####################################

def main(working_dir):
    os.chdir(working_dir)

    args = parse_arguments()
    (ref, vcf, out_pickle, save_trinuc, skip_common) = (
        args.r, args.m, args.o, args.save_trinuc, args.skip_common)

    # Process bed file
    annotations_df = None
    if args.b:
        print('Processing bed file...')
        try:
            annotations_df = seperate_exons_genes_intergenics(args.b, working_dir)
            annotations_df[['chrom', 'feature']] = annotations_df[['chrom', 'feature']].astype(str)
            # TODO validate chromosome length are same as in reference
        except ValueError:
            print('Problem parsing bed file. Ensure bed file is tab separated, standard bed format')

    regions_stats = RegionStats(annotations_df)

    # is_bed = False
    # annotations_df = None
    # if args.b:
    #     print('Processing bed file...')
    #     try:
    #         annotations_df = seperate_exons_genes_intergenics(args.b, working_dir)
    #         annotations_df[['chrom', 'feature']] = annotations_df[['chrom', 'feature']].astype(str)
    #         regions_stats = {
    #             'exon' : Stats(),
    #             'intron' : Stats(),
    #             'intergenic' : Stats()
    #         }
    #         annotated_sequence = AnnotatedSeqence(annotations_df)
    #         #TODO validate chromosome length are same as in reference
    #         is_bed = True
    #     except ValueError:
    #         print('Problem parsing bed file. Ensure bed file is tab separated, standard bed format')
    # else:
    #     regions_stats = {'all': Stats()}
    #     annotated_sequence = AnnotatedSeqence(None)
    #
    #     # my_bed = my_bed.rename(columns={0: 'chrom', 1: 'start', 2: 'end'})
    #     # Adding a couple of columns we'll need for later calculations
    #     # my_bed['coords'] = list(zip(my_bed.start, my_bed.end))
    #     annotations_df['track_len'] = annotations_df.end - annotations_df.start + 1
    #     # my_bed['chrom'] = my_bed['chrom'].astype(str)

    # Process reference file
    print('Processing reference...')
    try:
        reference = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
    except ValueError:
        print("Problems parsing reference file. Ensure reference is in proper fasta format")

    ref_dict = {}
    for key in reference.keys():
        key_split = key.split("|")
        ref_dict[key_split[0]] = reference[key]

    ref_list = list(map(str, ref_dict.keys()))

    # Process VCF file. First check if it's been entered as a TSV
    if vcf[-3:] == 'tsv':
        print("Warning! TSV file must follow VCF specifications.")

    # Pre-parsing to find all the matching chromosomes between ref and vcf
    print('Processing VCF file...')
    try:
        variants = pd.read_csv(vcf, sep='\t', comment='#', index_col=None, header=None)
        variants[0] = variants[0].map(str)
    except ValueError:
        print("VCF must be in standard VCF format with tab-separated columns")

    # Narrow chromosomes to those matching the reference
    # This is in part to make sure the names match
    variant_chroms = variants[0].to_list()
    variant_chroms = list(set(variant_chroms))
    matching_chromosomes = []
    for ref_name in ref_list:
        if ref_name not in variant_chroms:
            continue
        else:
            matching_chromosomes.append(ref_name)

    # Check to make sure there are some matches
    if not matching_chromosomes:
        print("Found no chromosomes in common between VCF and Fasta. Please fix the chromosome names and try again")
        exit(1)

    # Double check that there are matches
    try:
        matching_variants = variants[variants[0].isin(matching_chromosomes)]
    except ValueError:
        print("Problem matching variants with reference.")

    if matching_variants.empty:
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

    # # Now we check that the bed and vcf have matching regions
    # # This also checks that the vcf and bed have the same naming conventions and cuts out scaffolding.
    # if is_bed:#TODO remove?
    #     annotations_chroms = list(set(annotations_df['chrom'])) #[str(val) for val in set(my_bed['chrom'])]
    #     relevant_chroms = list(set(annotations_chroms) & set(variant_chroms))
    #     try:
    #         matching_annotations = annotations_df[annotations_df['chrom'].isin(relevant_chroms)]
    #     except ValueError:
    #         print('Problem matching bed chromosomes to variant file.')
    #
    #     if matching_annotations.empty:
    #         print("There is no overlap between bed and variant file. "
    #               "This could be a chromosome naming problem")
    #         exit(1)
    #
    # # Count Trinucleotides in reference, based on bed or not
    # print('Counting trinucleotides in reference...')

    if annotations_df:
        print("since you're using a bed input, we have to count trinucs in bed region even if "
              "you already have a trinuc count file for the reference...")
        for ref_name in matching_chromosomes:
            chrom_annotations = annotations_df[annotations_df['chrom'] == ref_name]
            for i, annotation in chrom_annotations.iterrows():
                sub_seq = ref_dict[ref_name][annotation['start']: annotation['end']].seq
                region = regions_stats.get_region(chrom_annotations, annotation['start'])
                update_trinuc_ref_count(sub_seq, regions_stats, region)

    elif not os.path.isfile(ref + '.trinucCounts'):
        for ref_name in matching_chromosomes:
            sub_seq = ref_dict[ref_name].seq
            update_trinuc_ref_count(sub_seq, regions_stats)
    else:
        print('Found trinucCounts file, using that.')

    # Load and process variants in each reference sequence individually, for memory reasons...
    print('Creating mutational model...')
    for ref_name in matching_chromosomes:
        update_total_reflen(ref_dict, ref_name, regions_stats, annotations_df)

        # list to be used for counting variants that occur multiple times in file (i.e. in multiple samples)
        VDAT_COMMON = []

        # Create a view that narrows variants list to current ref
        variants_to_process = matching_variants[matching_variants['CHROM'] == ref_name].copy()
        ref_sequence = str(ref_dict[ref_name].seq)

        # we want only snps
        # so, no '-' characters allowed, and chrStart must be same as chrEnd
        snp_df = variants_to_process[~variants_to_process.index.isin(indices_to_indels)]
        snp_df = snp_df.loc[snp_df['chr_start'] == snp_df['chr_end']]
        # print(snp_df[pd.isnull(snp_df['chr_start'])])
        # if is_bed:
        #     annotations_to_process = matching_annotations[matching_annotations['chrom'] == ref_name].copy()
            # # TODO fix this line (need the intersection of these two, I think)
            # snp_df = annotations_to_process.join(snp_df) #TODO change here
        # print(snp_df[pd.isnull(snp_df['chr_start'])])

        if not snp_df.empty:
            # snp_df['chr_start'] = snp_df['chr_start'].astype("Int64")
            # only consider positions where ref allele in vcf matches the nucleotide in our reference
            for index, row in snp_df.iterrows():
                trinuc_to_analyze = str(ref_sequence[row.chr_start - 1: row.chr_start + 2])
                if trinuc_to_analyze not in VALID_TRINUC:
                    continue
                if row.REF == trinuc_to_analyze[1]:
                    trinuc_ref = trinuc_to_analyze
                    trinuc_alt = trinuc_to_analyze[0] + snp_df.loc[index, 'ALT'] + trinuc_to_analyze[2]
                    if trinuc_alt not in VALID_TRINUC:
                        continue
                    region = regions_stats.get_annotation(ref_name, index)
                    update_trinuc_transition_count(trinuc_alt, trinuc_ref, regions_stats, region)
                    update_snp_count(regions_stats, region)
                    update_snp_transition_count(str(row.REF), str(row.ALT), regions_stats, region)

                    my_pop_freq = VCF_DEFAULT_POP_FREQ
                    if ';CAF=' in snp_df.loc[index, 'INFO']:
                        caf_str = re.findall(r";CAF=.*?(?=;)", row.INFO)[0]
                        if ',' in caf_str:
                            my_pop_freq = float(caf_str[5:].split(',')[1])
                    VDAT_COMMON.append(
                        (row.chr_start, row.REF, row.REF, row.ALT, my_pop_freq))
                else:
                    print('\nError: ref allele in variant call does not match reference.\n')
                    exit(1)

        # now let's look for indels...
        indel_df = variants_to_process[variants_to_process.index.isin(indices_to_indels)]
        if not indel_df.empty:
            for index, row in indel_df.iterrows():
                if "-" in row.REF:
                    len_ref = 0
                else:
                    len_ref = len(row.REF)
                if "-" in row.ALT:
                    len_alt = 0
                else:
                    len_alt = len(row.ALT)
                if len_ref != len_alt:
                    indel_len = len_alt - len_ref
                    region = regions_stats.get_annotation(ref_name, index)
                    update_indel_count(indel_len, regions_stats, region)

                    my_pop_freq = VCF_DEFAULT_POP_FREQ
                    if ';CAF=' in row.INFO:
                        caf_str = re.findall(r";CAF=.*?(?=;)", row.INFO)[0]
                        if ',' in caf_str:
                            my_pop_freq = float(caf_str[5:].split(',')[1])
                    VDAT_COMMON.append((row.chr_start, row.REF, row.REF, row.ALT, my_pop_freq))

        # if we didn't find anything, skip ahead along to the next reference sequence
        if not len(VDAT_COMMON):
            print('Found no variants for this reference.')
            continue

        # identify common mutations
        percentile_var = 95
        min_value = np.percentile([n[4] for n in VDAT_COMMON], percentile_var)
        for k in sorted(VDAT_COMMON):
            if k[4] >= min_value:
                COMMON_VARIANTS.append((ref_name, k[0], k[1], k[3], k[4]))
        VDAT_COMMON = {(n[0], n[1], n[2], n[3]): n[4] for n in VDAT_COMMON}

        # identify areas that have contained significantly higher random mutation rates
        dist_thresh = 2000
        percentile_clust = 97
        scaler = 1000
        # identify regions with disproportionately more variants in them
        VARIANT_POS = sorted([n[0] for n in VDAT_COMMON.keys()])
        clustered_pos = cluster_list(VARIANT_POS, dist_thresh)
        by_len = [(len(clustered_pos[i]), min(clustered_pos[i]), max(clustered_pos[i]), i) for i in
                 range(len(clustered_pos))]
        # Not sure what this was intended to do or why it is commented out. Leaving it here for now.
        # by_len  = sorted(by_len,reverse=True)
        # minLen = int(np.percentile([n[0] for n in by_len],percentile_clust))
        # by_len  = [n for n in by_len if n[0] >= minLen]
        candidate_regions = []
        for n in by_len:
            bi = int((n[1] - dist_thresh) / float(scaler)) * scaler
            bf = int((n[2] + dist_thresh) / float(scaler)) * scaler
            candidate_regions.append((n[0] / float(bf - bi), max([0, bi]), min([len(ref_dict[ref_name]), bf])))
        minimum_value = np.percentile([n[0] for n in candidate_regions], percentile_clust)
        for n in candidate_regions:
            if n[0] >= minimum_value:
                HIGH_MUT_REGIONS.append((ref_name, n[1], n[2], n[0]))
        # collapse overlapping regions
        for i in range(len(HIGH_MUT_REGIONS) - 1, 0, -1):
            if HIGH_MUT_REGIONS[i - 1][2] >= HIGH_MUT_REGIONS[i][1] and HIGH_MUT_REGIONS[i - 1][0] == \
                    HIGH_MUT_REGIONS[i][0]:
                # Might need to research a more accurate way to get the mutation rate for this region
                avg_mut_rate = 0.5 * HIGH_MUT_REGIONS[i - 1][3] + 0.5 * HIGH_MUT_REGIONS[i][
                    3]
                HIGH_MUT_REGIONS[i - 1] = (
                    HIGH_MUT_REGIONS[i - 1][0], HIGH_MUT_REGIONS[i - 1][1], HIGH_MUT_REGIONS[i][2], avg_mut_rate)
                del HIGH_MUT_REGIONS[i]

    # if we didn't count ref trinucs because we found file, read in ref counts from file now
    if os.path.isfile(ref + '.trinucCounts'):
        print('reading pre-computed trinuc counts...')
        f = open(ref + '.trinucCounts', 'r')
        for line in f:
            splt = line.strip().split('\t')
            TRINUC_REF_COUNT[splt[0]] = int(splt[1])
        f.close()
    # otherwise, save trinuc counts to file, if desired
    elif save_trinuc:
        if annotations_df:
            print('unable to save trinuc counts to file because using input bed region...')
        else:
            print('saving trinuc counts to file...')
            f = open(ref + '.trinucCounts', 'w')
            for trinuc in sorted(TRINUC_REF_COUNT.keys()):
                f.write(trinuc + '\t' + str(TRINUC_REF_COUNT[trinuc]) + '\n')
            f.close()


    compute_probabilities(regions_stats)

    #
    # save variables to file
    #
    if skip_common:
        OUT_DICT = {'AVG_MUT_RATE': AVG_MUT_RATE,
                    'SNP_FREQ': SNP_FREQ,
                    'SNP_TRANS_FREQ': SNP_TRANS_FREQ,
                    'INDEL_FREQ': INDEL_FREQ,
                    'TRINUC_MUT_PROB': TRINUC_MUT_PROB,
                    'TRINUC_TRANS_PROBS': TRINUC_TRANS_PROBS}
    else:
        OUT_DICT = {'AVG_MUT_RATE': AVG_MUT_RATE,
                    'SNP_FREQ': SNP_FREQ,
                    'SNP_TRANS_FREQ': SNP_TRANS_FREQ,
                    'INDEL_FREQ': INDEL_FREQ,
                    'TRINUC_MUT_PROB': TRINUC_MUT_PROB,
                    'TRINUC_TRANS_PROBS': TRINUC_TRANS_PROBS,
                    'COMMON_VARIANTS': COMMON_VARIANTS,
                    'HIGH_MUT_REGIONS': HIGH_MUT_REGIONS}
    pickle.dump(OUT_DICT, open(out_pickle, "wb"))


def compute_probabilities(regions_stats):
    # if for some reason we didn't find any valid input variants AT ALL, exit gracefully...
    total_var = regions_stats.get_stat_by_region(Regions.ALL, BasicStats.SNP_COUNT) + \
                sum(regions_stats.get_stat_by_region(Regions.ALL, BasicStats.INDEL_COUNT).values())
    if total_var == 0:
        print(
            '\nError: No valid variants were found, model could not be created. (Are you using the correct reference?)\n')
        exit(1)

    for region_name, region_stats in regions_stats.get_all_stats().items():
        ###	COMPUTE PROBABILITIES

        # frequency that each trinuc mutated into anything else
        region_stats[AdditionalStats.TRINUC_MUT_PROB] = {}
        # frequency that a trinuc mutates into another trinuc, given that it mutated
        region_stats[AdditionalStats.TRINUC_TRANS_PROBS] = {}
        # frequency of snp transitions, given a snp occurs.
        region_stats[AdditionalStats.SNP_TRANS_FREQ] = {}

        for trinuc in sorted(region_stats[BasicStats.TRINUC_REF_COUNT].keys()):
            my_count = 0
            for k in sorted(region_stats[BasicStats.TRINUC_TRANSITION_COUNT].keys()):
                if k[0] == trinuc:
                    my_count += region_stats[BasicStats.TRINUC_TRANSITION_COUNT][k]
            region_stats[AdditionalStats.TRINUC_MUT_PROB][trinuc] = \
                my_count / float(region_stats[BasicStats.TRINUC_REF_COUNT][trinuc])
            for k in sorted(region_stats[BasicStats.TRINUC_TRANSITION_COUNT].keys()):
                if k[0] == trinuc:
                    region_stats[AdditionalStats.TRINUC_TRANS_PROBS][k] = \
                        region_stats[BasicStats.TRINUC_TRANSITION_COUNT][k] / float(my_count)

        for n1 in VALID_NUCL:
            rolling_tot = sum([region_stats[BasicStats.SNP_TRANSITION_COUNT][(n1, n2)] for n2 in VALID_NUCL if (n1, n2)
                               in region_stats[BasicStats.SNP_TRANSITION_COUNT]])
            for n2 in VALID_NUCL:
                key2 = (n1, n2)
                if key2 in region_stats[BasicStats.SNP_TRANSITION_COUNT]:
                    region_stats[AdditionalStats.SNP_TRANS_FREQ][key2] = \
                        region_stats[BasicStats.SNP_TRANSITION_COUNT][key2] / float(rolling_tot)

        # compute average snp and indel frequencies
        total_var = region_stats[BasicStats.SNP_COUNT] + sum(region_stats[BasicStats.INDEL_COUNT].values())
        region_stats[AdditionalStats.SNP_FREQ] = region_stats[BasicStats.SNP_COUNT] / float(total_var)
        region_stats[AdditionalStats.AVG_INDEL_FREQ] = 1. - region_stats[BasicStats.SNP_FREQ]
        region_stats[AdditionalStats.INDEL_FREQ] = \
            {k: (region_stats[BasicStats.INDEL_COUNT][k] / float(total_var)) / region_stats[BasicStats.AVG_INDEL_FREQ]
             for k in region_stats[BasicStats.INDEL_COUNT].keys()}
        region_stats[AdditionalStats.AVG_MUT_RATE] = total_var / region_stats[BasicStats.TOTAL_REFLEN]

        #	if values weren't found in data, appropriately append null entries
        print_trinuc_warning = False
        for trinuc in VALID_TRINUC:
            trinuc_mut = [trinuc[0] + n + trinuc[2] for n in VALID_NUCL if n != trinuc[1]]
            if trinuc not in region_stats[AdditionalStats.TRINUC_MUT_PROB]:
                region_stats[AdditionalStats.TRINUC_MUT_PROB][trinuc] = 0.
                print_trinuc_warning = True
            for trinuc2 in trinuc_mut:
                if (trinuc, trinuc2) not in region_stats[AdditionalStats.TRINUC_TRANS_PROBS]:
                    region_stats[AdditionalStats.TRINUC_TRANS_PROBS][(trinuc, trinuc2)] = 0.
                    print_trinuc_warning = True
        if print_trinuc_warning:
            print(
                'Warning: Some trinucleotides transitions were not encountered in the input dataset, '
                'probabilities of 0.0 have been assigned to these events.')

        #
        #	print some stuff
        #
        print(f'Probabilities for region {region_name}:')

        for k in sorted(region_stats[AdditionalStats.TRINUC_MUT_PROB].keys()):
            print('p(' + k + ' mutates) =', region_stats[AdditionalStats.TRINUC_MUT_PROB][k])

        for k in sorted(region_stats[AdditionalStats.TRINUC_TRANS_PROBS].keys()):
            print('p(' + k[0] + ' --> ' + k[1] + ' | ' + k[0] + ' mutates) =',
                  region_stats[AdditionalStats.TRINUC_TRANS_PROBS][k])

        for k in sorted(region_stats[AdditionalStats.INDEL_FREQ].keys()):
            if k > 0:
                print('p(ins length = ' + str(abs(k)) + ' | indel occurs) =',
                      region_stats[AdditionalStats.INDEL_FREQ][k])
            else:
                print('p(del length = ' + str(abs(k)) + ' | indel occurs) =',
                      region_stats[AdditionalStats.INDEL_FREQ][k])

        for k in sorted(region_stats[AdditionalStats.SNP_TRANS_FREQ].keys()):
            print('p(' + k[0] + ' --> ' + k[1] + ' | SNP occurs) =', region_stats[AdditionalStats.SNP_TRANS_FREQ][k])


        print('p(snp)   =', region_stats[AdditionalStats.SNP_FREQ])
        print('p(indel) =', region_stats[AdditionalStats.AVG_INDEL_FREQ])
        print('overall average mut rate:', region_stats[AdditionalStats.AVG_MUT_RATE])
        print('total variants processed:', total_var)


def update_indel_count(indel_len, regions_stats, region = Regions.ALL):
    regions_to_update = {Regions.ALL, region}
    for current_region in regions_to_update:
        if indel_len not in regions_stats.get_stat_by_region(current_region, BasicStats.INDEL_COUNT):
            regions_stats.get_stat_by_region(current_region, BasicStats.INDEL_COUNT)[indel_len] = 0
        regions_stats.get_stat_by_region(current_region, BasicStats.INDEL_COUNT)[indel_len] += 1


def update_total_reflen(ref_dict, ref_name, regions_stats, annotations_df):
    # Count the number of non-N nucleotides for the reference
    total_reflen_all_chrom = len(ref_dict[ref_name].seq) - ref_dict[ref_name].seq.count('N')
    update_total_reflen_for_region(total_reflen_all_chrom, regions_stats)
    chrom_annotations = annotations_df[annotations_df['chrom'] == ref_name]
    if annotations_df:
        for i, annotation in chrom_annotations.iterrows():
            sub_seq = ref_dict[ref_name][annotation['start']: annotation['end']].seq
            reflen_annotation = len(sub_seq) - sub_seq.count('N')
            update_total_reflen_for_region(reflen_annotation, regions_stats, annotations_df['feature'])


def update_total_reflen_for_region(total_reflen_region_chrom, regions_stats, region = Regions.ALL):
    total_reflen_region = regions_stats.get_stat_by_region(region, BasicStats.TOTAL_REFLEN)
    total_reflen_region[0] += total_reflen_region_chrom


def update_snp_transition_count(row_ref, row_alt, regions_stats, region = Regions.ALL):
    key2 = (row_ref, row_alt)
    regions_to_update = {Regions.ALL, region}
    for current_region in regions_to_update:
        if key2 not in regions_stats.get_stat_by_region(current_region, BasicStats.SNP_TRANSITION_COUNT):
            regions_stats.get_stat_by_region(current_region, BasicStats.SNP_TRANSITION_COUNT)[key2] = 0
        regions_stats.get_stat_by_region(current_region, BasicStats.SNP_TRANSITION_COUNT)[key2] += 1


def update_snp_count(regions_stats, region = Regions.ALL):
    regions_to_update = {Regions.ALL, region}
    for current_region in regions_to_update:
        regions_stats.get_stat_by_region(current_region, BasicStats.SNP_COUNT)[0] += 1


def update_trinuc_transition_count(trinuc_alt, trinuc_ref, regions_stats, region = Regions.ALL):
    regions_to_update = {Regions.ALL, region}
    key = (trinuc_ref, trinuc_alt)
    for current_region in regions_to_update:
        if key not in regions_stats.get_stat_by_region(current_region, BasicStats.TRINUC_TRANSITION_COUNT):
            regions_stats.get_stat_by_region(current_region, BasicStats.TRINUC_TRANSITION_COUNT)[key] = 0
        regions_stats.get_stat_by_region(current_region, BasicStats.TRINUC_TRANSITION_COUNT)[key] += 1


def update_trinuc_ref_count(sub_seq, regions_stats, region = Regions.ALL):
    regions_to_update = {Regions.ALL, region}
    for trinuc in VALID_TRINUC:
        for current_region in regions_to_update:
            if trinuc not in regions_stats.get_stat_by_region(current_region, BasicStats.TRINUC_REF_COUNT):
                regions_stats.get_stat_by_region(current_region, BasicStats.TRINUC_REF_COUNT)[trinuc] = 0
            regions_stats.get_stat_by_region(current_region, BasicStats.TRINUC_REF_COUNT)[trinuc] += sub_seq.count_overlap(
                trinuc)


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
    parser.add_argument('--save-trinuc', required=False, action='store_true', default=False,
                        help='save trinucleotide counts for reference')
    parser.add_argument('--skip-common', required=False, action='store_true', default=False,
                        help='Do not save common snps + high mut regions')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    working_dir = '/groups/itay_mayrose/adisivan/PanGenomeSimulator/test_arabidopsis'
    main(working_dir)
