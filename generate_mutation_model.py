import argparse
import pickle

import pandas as pd
from Bio import SeqIO

from utilities.annotated_sequence import AnnotatedSequence
from utilities.common_data_structues import Region, VALID_TRINUC, VALID_NUCL, ModelStats
from utilities.io.genome_annotations_reader import read_annotations_csv

# if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
VCF_DEFAULT_POP_FREQ = 0.00001
VCF_CHROM_COL = 0


class RegionsStats:
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
        if not self._annotated_sequence_per_chrom:
            return Region.ALL
        region, _ = self._annotated_sequence_per_chrom[chrom].get_region_by_position(pos)
        return region

    def get_stat_by_region(self, region: Region, stat: ModelStats):
        return self._regions_stats[region.value][stat.value]

    def get_stat_by_location(self, chrom: str, pos: int, stat: ModelStats):
        region = self.get_region_by_chrom_and_pos(chrom, pos)
        return self.get_stat_by_region(region, stat)

    def get_all_stats(self):
        return self._regions_stats

    @staticmethod
    def create_stats_dict():
        return {
            ModelStats.TRINUC_REF_COUNT.value: {},
            ModelStats.TRINUC_TRANSITION_COUNT.value: {},
            ModelStats.SNP_COUNT.value: [0],
            ModelStats.SNP_TRANSITION_COUNT.value: {},
            ModelStats.INDEL_COUNT.value: {},
            ModelStats.TOTAL_REFLEN.value: [0],
            ModelStats.TRINUC_MUT_PROB.value: {},
            ModelStats.TRINUC_TRANS_PROBS.value: {},
            ModelStats.SNP_TRANS_FREQ.value: {}
        }


#####################################
#				main()				#
#####################################


def main():
    args = parse_arguments()
    (ref, vcf, out_pickle, annotations_file) = (
        args.r, args.m, args.o, args.a)

    chrom_sequences_dict, chrom_names = process_reference(ref)

    indices_to_indels, matching_chromosomes, matching_variants = process_vcf(chrom_names, vcf)

    annotations_df = read_annotations_csv(annotations_file)

    regions_stats = RegionsStats(annotations_df)

    process_trinuc_counts(chrom_sequences_dict, matching_chromosomes, annotations_df, regions_stats)

    # Load and process variants in each reference sequence individually, for memory reasons...
    print('Creating mutational model...')
    for chrom_name in matching_chromosomes:
        collect_basic_stats_for_chrom(chrom_sequences_dict, chrom_name, indices_to_indels, matching_variants,
                                      annotations_df, regions_stats)

    compute_probabilities(regions_stats)

    save_stats_to_file(out_pickle, regions_stats)


#########################################################
#				VARIOUS HELPER FUNCTIONS				#
#########################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description='mutation_model.source',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-r', type=str, required=True, metavar='/path/to/reference.fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-m', type=str, required=True, metavar='/path/to/variants.vcf',
                        help="Mutation file for organism in VCF format")
    parser.add_argument('-o', type=str, required=True, metavar='/path/to/output/model.p',
                        help="Path of output model")
    parser.add_argument('-a', type=str, required=False, metavar='/path/to/annotations.csv',
                        help='Annotations CSV file of regions to consider', default=None)
    args = parser.parse_args()
    return args


def process_reference(ref_filename: str) -> (dict, list):
    # Process reference file
    print('Processing reference...')
    reference = None
    try:
        reference = SeqIO.to_dict(SeqIO.parse(ref_filename, "fasta"))
    except ValueError:
        print("Problems parsing reference file. Ensure reference is in proper fasta format")
        exit(1)
    chrom_sequences_dict = {}
    for chrom in reference.keys():
        chrom_split = chrom.split("|")
        chrom_sequences_dict[chrom_split[0]] = reference[chrom]
    chrom_names = list(map(str, chrom_sequences_dict.keys()))
    return chrom_sequences_dict, chrom_names


def process_vcf(chrom_list: list, vcf_filename: str) -> (list, list, pd.DataFrame):
    # Process VCF file. First check if it's been entered as a TSV
    if vcf_filename[-3:] == 'tsv':
        print("Warning! TSV file must follow VCF specifications.")
    # Pre-parsing to find all the matching chromosomes between ref and vcf
    print('Processing VCF file...')
    variants = []
    try:
        variants = pd.read_csv(vcf_filename, sep='\t', comment='#', index_col=None, header=None)
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
        print(
            "Found no chromosomes in common between VCF and Fasta. Please fix the chromosome names and try again")
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
    return list(indices_to_indels), matching_chromosomes, matching_variants


def process_trinuc_counts(chrom_sequences_dict: dict, matching_chromosomes: list, annotations_df: pd.DataFrame,
                          regions_stats: RegionsStats):
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


def collect_basic_stats_for_chrom(chrom_sequences_dict: dict, chrom_name: str, indices_to_indels, matching_variants,
                                  annotations_df: pd.DataFrame, regions_stats: RegionsStats):
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


def process_snps(snp_df: pd.DataFrame, chrom_name: str, chrom_sequence: str, regions_stats: RegionsStats):
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

            else:
                print('\nError: ref allele in variant call does not match reference.\n')
                exit(1)


def process_indels(indel_df: pd.DataFrame, chrom_name: str, regions_stats: RegionsStats):
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


def compute_probabilities(regions_stats: RegionsStats):
    # if for some reason we didn't find any valid input variants AT ALL, exit gracefully...
    total_var_all = regions_stats.get_stat_by_region(Region.ALL, ModelStats.SNP_COUNT)[0] + \
                    sum(regions_stats.get_stat_by_region(Region.ALL, ModelStats.INDEL_COUNT).values())
    if total_var_all == 0:
        print(
            '\nError: No valid variants were found, model could not be created. '
            '(Are you using the correct reference?)\n')
        exit(1)

    for region_name, curr_region_stats in regions_stats.get_all_stats().items():
        for trinuc in sorted(curr_region_stats[ModelStats.TRINUC_REF_COUNT.value].keys()):
            count = 0

            # TRINUC_MUT_PROB
            for k in sorted(curr_region_stats[ModelStats.TRINUC_TRANSITION_COUNT.value].keys()):
                if k[0] == trinuc:
                    count += curr_region_stats[ModelStats.TRINUC_TRANSITION_COUNT.value][k]
            curr_region_stats[ModelStats.TRINUC_MUT_PROB.value][trinuc] = \
                count / float(curr_region_stats[ModelStats.TRINUC_REF_COUNT.value][trinuc])

            # TRINUC_TRANS_PROBS
            for k in sorted(curr_region_stats[ModelStats.TRINUC_TRANSITION_COUNT.value].keys()):
                if k[0] == trinuc:
                    curr_region_stats[ModelStats.TRINUC_TRANS_PROBS.value][k] = \
                        curr_region_stats[ModelStats.TRINUC_TRANSITION_COUNT.value][k] / float(count)

        # SNP_TRANS_FREQ
        for n1 in VALID_NUCL:
            rolling_tot = sum([curr_region_stats[ModelStats.SNP_TRANSITION_COUNT.value][(n1, n2)]
                               for n2 in VALID_NUCL if (n1, n2)
                               in curr_region_stats[ModelStats.SNP_TRANSITION_COUNT.value]])
            for n2 in VALID_NUCL:
                key2 = (n1, n2)
                if key2 in curr_region_stats[ModelStats.SNP_TRANSITION_COUNT.value]:
                    curr_region_stats[ModelStats.SNP_TRANS_FREQ.value][key2] = \
                        curr_region_stats[ModelStats.SNP_TRANSITION_COUNT.value][key2] / float(rolling_tot)

        # compute average snp and indel frequencies
        total_indels_region = sum(curr_region_stats[ModelStats.INDEL_COUNT.value].values())
        total_var_region = curr_region_stats[ModelStats.SNP_COUNT.value][0] + total_indels_region
        if total_var_region != 0:
            curr_region_stats[ModelStats.SNP_FREQ.value] = \
                curr_region_stats[ModelStats.SNP_COUNT.value][0] / float(total_var_region)
            curr_region_stats[ModelStats.AVG_INDEL_FREQ.value] = 1. - curr_region_stats[ModelStats.SNP_FREQ.value]
            curr_region_stats[ModelStats.INDEL_FREQ.value] = \
                {k: curr_region_stats[ModelStats.INDEL_COUNT.value][k] / total_indels_region
                 for k in curr_region_stats[ModelStats.INDEL_COUNT.value].keys()}
        else:
            curr_region_stats[ModelStats.SNP_FREQ.value] = 0.
            curr_region_stats[ModelStats.AVG_INDEL_FREQ.value] = 1.
            curr_region_stats[ModelStats.INDEL_FREQ.value] = {}
        curr_region_stats[ModelStats.AVG_MUT_RATE.value] = total_var_region / curr_region_stats[
            ModelStats.TOTAL_REFLEN.value][0]

        # if values weren't found in data, appropriately append null entries
        print_trinuc_warning = False
        for trinuc in VALID_TRINUC:
            trinuc_mut = [trinuc[0] + n + trinuc[2] for n in VALID_NUCL if n != trinuc[1]]
            if trinuc not in curr_region_stats[ModelStats.TRINUC_MUT_PROB.value]:
                curr_region_stats[ModelStats.TRINUC_MUT_PROB.value][trinuc] = 0.
                print_trinuc_warning = True
            for trinuc2 in trinuc_mut:
                if (trinuc, trinuc2) not in curr_region_stats[ModelStats.TRINUC_TRANS_PROBS.value]:
                    curr_region_stats[ModelStats.TRINUC_TRANS_PROBS.value][(trinuc, trinuc2)] = 0.
                    print_trinuc_warning = True
        if print_trinuc_warning:
            print(
                'Warning: Some trinucleotides transitions were not encountered in the input dataset, '
                'probabilities of 0.0 have been assigned to these events.')

        #
        # print some stuff
        #
        print(f'Probabilities for region {region_name}:')

        for k in sorted(curr_region_stats[ModelStats.TRINUC_MUT_PROB.value].keys()):
            print(f'p({k} mutates) = {curr_region_stats[ModelStats.TRINUC_MUT_PROB.value][k]}')

        for k in sorted(curr_region_stats[ModelStats.TRINUC_TRANS_PROBS.value].keys()):
            print(
                f'p({k[0]} --> {k[1]} | {k[0]} mutates) = {curr_region_stats[ModelStats.TRINUC_TRANS_PROBS.value][k]}')

        for k in sorted(curr_region_stats[ModelStats.INDEL_FREQ.value].keys()):
            if k > 0:
                print(
                    f'p(ins length = {str(abs(k))} | indel occurs) = {curr_region_stats[ModelStats.INDEL_FREQ.value][k]}')
            else:
                print(
                    f'p(del length = {str(abs(k))} | indel occurs) = {curr_region_stats[ModelStats.INDEL_FREQ.value][k]}')

        for k in sorted(curr_region_stats[ModelStats.SNP_TRANS_FREQ.value].keys()):
            print(f'p({k[0]} --> {k[1]} | SNP occurs) = {curr_region_stats[ModelStats.SNP_TRANS_FREQ.value][k]}')

        print(f'p(snp) = {curr_region_stats[ModelStats.SNP_FREQ.value]}')
        print(f'p(indel) = {curr_region_stats[ModelStats.AVG_INDEL_FREQ.value]}')
        print(f'overall average mut rate: {curr_region_stats[ModelStats.AVG_MUT_RATE.value]}')
        print(f'total variants processed: {total_var_region}')


def update_trinuc_ref_count(sub_seq, regions_stats, region: Region = Region.ALL, also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    for trinuc in VALID_TRINUC:
        for current_region in regions_to_update:
            if trinuc not in regions_stats.get_stat_by_region(current_region, ModelStats.TRINUC_REF_COUNT):
                regions_stats.get_stat_by_region(current_region, ModelStats.TRINUC_REF_COUNT)[trinuc] = 0
            regions_stats.get_stat_by_region(current_region, ModelStats.TRINUC_REF_COUNT)[trinuc] += sub_seq.count_overlap(
                trinuc)


def update_total_reflen(chrom_sequences_dict, chrom_name, regions_stats, annotations_df):
    # Count the number of non-N nucleotides for the reference
    total_reflen_all_chrom = len(chrom_sequences_dict[chrom_name].seq) - chrom_sequences_dict[chrom_name].seq.count('N')
    update_total_reflen_for_region(total_reflen_all_chrom, regions_stats, Region.ALL)
    if annotations_df is not None and not annotations_df.empty:
        chrom_annotations = annotations_df[annotations_df['chrom'] == chrom_name]
        for i, annotation in chrom_annotations.iterrows():
            sub_seq = chrom_sequences_dict[chrom_name][annotation['start']: annotation['end']].seq
            reflen_annotation = len(sub_seq) - sub_seq.count('N')
            update_total_reflen_for_region(reflen_annotation, regions_stats, Region(annotation['region']))


def update_total_reflen_for_region(total_reflen_region_chrom, regions_stats, region: Region = Region.ALL):
    total_reflen_region = regions_stats.get_stat_by_region(region, ModelStats.TOTAL_REFLEN)
    total_reflen_region[0] += total_reflen_region_chrom


def update_trinuc_transition_count(trinuc_ref, trinuc_alt, regions_stats, region: Region = Region.ALL,
                                   also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    key = (trinuc_ref, trinuc_alt)
    for current_region in regions_to_update:
        if key not in regions_stats.get_stat_by_region(current_region, ModelStats.TRINUC_TRANSITION_COUNT):
            regions_stats.get_stat_by_region(current_region, ModelStats.TRINUC_TRANSITION_COUNT)[key] = 0
        regions_stats.get_stat_by_region(current_region, ModelStats.TRINUC_TRANSITION_COUNT)[key] += 1


def update_snp_transition_count(ref, alt, regions_stats, region: Region = Region.ALL, also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    key2 = (ref, alt)
    for current_region in regions_to_update:
        if key2 not in regions_stats.get_stat_by_region(current_region, ModelStats.SNP_TRANSITION_COUNT):
            regions_stats.get_stat_by_region(current_region, ModelStats.SNP_TRANSITION_COUNT)[key2] = 0
        regions_stats.get_stat_by_region(current_region, ModelStats.SNP_TRANSITION_COUNT)[key2] += 1


def update_snp_count(regions_stats, region: Region = Region.ALL, also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    for current_region in regions_to_update:
        regions_stats.get_stat_by_region(current_region, ModelStats.SNP_COUNT)[0] += 1


def update_indel_count(indel_len, regions_stats, region: Region = Region.ALL, also_update_all: bool = True):
    regions_to_update = {region, Region.ALL} if also_update_all else {region}
    for current_region in regions_to_update:
        if indel_len not in regions_stats.get_stat_by_region(current_region, ModelStats.INDEL_COUNT):
            regions_stats.get_stat_by_region(current_region, ModelStats.INDEL_COUNT)[indel_len] = 0
        regions_stats.get_stat_by_region(current_region, ModelStats.INDEL_COUNT)[indel_len] += 1


def save_stats_to_file(out_pickle: str, regions_stats: RegionsStats):
    out_dict = {}
    for region_name, region_stats in regions_stats.get_all_stats().items():
        out_dict[f'{region_name}.{ModelStats.TOTAL_REFLEN.value}'] = region_stats[ModelStats.TOTAL_REFLEN.value][0]
        out_dict[f'{region_name}.{ModelStats.AVG_MUT_RATE.value}'] = region_stats[ModelStats.AVG_MUT_RATE.value]
        out_dict[f'{region_name}.{ModelStats.SNP_FREQ.value}'] = region_stats[ModelStats.SNP_FREQ.value]
        out_dict[f'{region_name}.{ModelStats.SNP_TRANS_FREQ.value}'] = region_stats[ModelStats.SNP_TRANS_FREQ.value]
        out_dict[f'{region_name}.{ModelStats.INDEL_FREQ.value}'] = region_stats[ModelStats.INDEL_FREQ.value]
        out_dict[f'{region_name}.{ModelStats.TRINUC_MUT_PROB.value}'] = region_stats[ModelStats.TRINUC_MUT_PROB.value]
        out_dict[f'{region_name}.{ModelStats.TRINUC_TRANS_PROBS.value}'] = region_stats[
            ModelStats.TRINUC_TRANS_PROBS.value]
    pickle.dump(out_dict, open(out_pickle, "wb"))


if __name__ == "__main__":
    main()
