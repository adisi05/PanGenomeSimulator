#!/usr/bin/env source
# encoding: utf-8
""" ////////////////////////////////////////////////////////////////////////////////
   ///                                                                          ///
  ///       gen_reads.py                                                       ///
 ///        VERSION 3.0: HARDER, BETTER, FASTER, STRONGER!                    ///
///////                                                                      //////
   ///      Variant and read simulator for benchmarking NGS workflows          ///
  ///                                                                         ///
 ///        Written by:     Zach Stephens                                    ///
///////     For:            DEPEND Research Group, UIUC                     ///////
   ///      Date:           May 29, 2015                                       ///
  ///       Contact:        zstephe2@illinois.edu                             ///
 ///                                                                         ///
/////////////////////////////////////////////////////////////////////////////// """

import sys
import copy
import random
import re
import time
import bisect
import pickle
import numpy as np
import pathlib
import os

from Bio.SeqUtils import seq1

from source.bam_file_writer import BamFileWriter
from source.fasta_file_writer import FastaFileWriter
# from source.fastq_file_writer import FastqFileWriter
from source.input_checking import check_file_open, is_in_range
from source.ref_func import index_ref, read_ref
from source.vcf_file_writer import VcfFileWriter
from source.file_writer_utils import reverse_complement
from source.bam_file_writer import sam_flag
from source.probability import DiscreteDistribution, mean_ind_of_weighted_list
from source.SequenceContainer import SequenceContainer, ReadContainer, parse_input_mutation_model

"""
Some constants needed for analysis
"""

# target window size for read sampling. How many times bigger than read/frag length
WINDOW_TARGET_SCALE = 100

# allowed nucleotides
ALLOWED_NUCL = ['A', 'C', 'G', 'T']

def simulate(args):

    general_params, input_params, output_params, mutation_params, sequencing_params = parse_args(args)
    index_params = process_input_params(input_params) #, output_params["ploids"])
    load_sequencing_model(sequencing_params)
    load_mutation_model(mutation_params)

    # initialize output writers
    bam_file_writer, fasta_file_writer, vcf_file_writer = intialize_reads_writers(index_params,\
                                                                    input_params, output_params, sequencing_params)

    # Using pathlib to make this more machine agnostic
    output_params["out_prefix_name"] = pathlib.Path(output_params["out_prefix"]).name
    # keep track of the number of reads we've sampled, for read-names
    read_name_count = 1
    unmapped_records = []

    for chrom in range(len(index_params["ref_index"])):

        simulate_chrom(general_params, input_params, output_params, mutation_params, sequencing_params, index_params,
                       bam_file_writer, fasta_file_writer, vcf_file_writer, chrom, read_name_count, unmapped_records)

    #TODO translocation feature

    # write unmapped reads to bam file
    write_unmapped_to_bam(bam_file_writer, sequencing_params["paired_end"], output_params["save_bam"], unmapped_records)

    # close output files
    fasta_file_writer.close_file()
    if output_params["save_bam"]:
        bam_file_writer.close_file()
    if not output_params["no_fastq"]:
        generate_reads(fasta_file_writer.get_file(), sequencing_params)
    if output_params["save_vcf"]:
        vcf_file_writer.close_file(add_parent_variants=True)


def generate_reads(fasta_files, sequencing_params):
    print(fasta_files)
    print(sequencing_params)
    # TODO test didn't break anything...
    paired = "-p" if sequencing_params['paired_end'] else ""
    fastq_files = [filename.removesuffix('.fasta') +"_read" for filename in fasta_files]
    read_length = sequencing_params['read_len']
    coverage = sequencing_params['coverage']
    insert_size = sequencing_params['fragment_size']
    insert_std = sequencing_params['fragment_std']
    for fasta, fastq in zip(fasta_files,fastq_files):
        art_command = "ART/art_bin_MountRainier/art_illumina {} -i {} -l {} -f {} -o {} -m {} -s {}".format(paired,fasta,read_length,coverage,fastq,insert_size,insert_std)
        start = time.time()
        os.system(art_command)
        end = time.time()
        print("ART reads simulation took {} seconds.".format(int(end - start)))


def parse_args(args):
    general_params, input_params, output_params, mutation_params, sequencing_params = extract_params(args)
    params_sanity_check(input_params, output_params, general_params, sequencing_params)
    # if coverage val for a given window/position is below this value, consider it effectively zero.
    sequencing_params["low_cov_thresh"] = 50
    return general_params, input_params, output_params, mutation_params, sequencing_params

def extract_params(args):
    input_params = {
        "reference": args.r,
        # "input_vcf": args.v,
        "input_variants": args.input_variants,
        "input_bed": args.tr,
        "discard_bed": args.dr
    }
    mutation_params = {
        "mut_bed": args.Mb,
        "mut_model": args.m,
        "mut_rate": args.M,
        "dist": args.dist
    }
    output_params = {
        "accession": args.name,
        "out_prefix": args.o+'_'+args.name,
        "parent_prefix": args.o+'_'+args.parent_name if args.parent_name else None,
        "save_bam": args.bam,
        "save_vcf": args.vcf,
        "no_fastq": args.no_fastq or args.internal,  # TODO maybe convert to save-fastq?
        "ploids": args.p  # TODO validate it's an output param
    }
    general_params = {
        "rng_seed": args.rng,
        "debug": args.d
    }
    (fragment_size, fragment_std) = args.pe
    sequencing_params = {
        "read_len": args.R,
        "coverage": args.c,
        "se_model": args.e,
        "se_rate": None if args.E == -1 else args.E,
        "off_target_scalar": args.to,
        "off_target_discard": args.discard_offtarget,
        "force_coverage": args.force_coverage,
        "rescale_qual": args.rescale_qual,
        "n_max_qual": args.N,
        "fragment_size": fragment_size,
        "fragment_std": fragment_std,
        "fraglen_model": args.pe_model,
        "gc_bias_model": args.gc_model
    }
    return general_params, input_params, output_params, mutation_params, sequencing_params

def params_sanity_check(input_params, output_params, general_params, sequencing_params):
    # Check that files are real, if provided
    check_file_open(input_params["reference"], 'ERROR: could not open reference, {}'.format(input_params["reference"]), required=True)
    # check_file_open(input_params["input_vcf"], 'ERROR: could not open input VCF, {}'.format(input_params["input_vcf"]), required=False)
    check_file_open(input_params["input_bed"], 'ERROR: could not open input BED, {}'.format(input_params["input_bed"]), required=False)
    # if user specified no fastq, not fasta only, and no bam and no vcf, then print error and exit.
    if output_params["no_fastq"] and not output_params["save_bam"] and not output_params["save_vcf"]:
        print('\nERROR: No files would be written.\n')
        sys.exit(1)
    if (sequencing_params["fragment_size"] is None and sequencing_params["fragment_std"] is not None) or (sequencing_params["fragment_size"] is not None and sequencing_params["fragment_std"] is None):
        print('\nERROR: --pe argument takes 2 space-separated arguments.\n')
        sys.exit(1)
    if general_params["rng_seed"] == -1:
        general_params["rng_seed"] = random.randint(1, 99999999)
    random.seed(general_params["rng_seed"])
    is_in_range(sequencing_params["read_len"], 10, 1000000, 'Error: -R must be between 10 and 1,000,000')
    is_in_range(sequencing_params["coverage"], 0, 1000000, 'Error: -c must be between 0 and 1,000,000')
    is_in_range(output_params["ploids"], 1, 100, 'Error: -p must be between 1 and 100')
    is_in_range(sequencing_params["off_target_scalar"], 0, 1, 'Error: -to must be between 0 and 1')
    if sequencing_params["se_rate"] != None:
        is_in_range(sequencing_params["se_rate"], 0, 0.3, 'Error: -E must be between 0 and 0.3')
    if sequencing_params["n_max_qual"] != -1:
        is_in_range(sequencing_params["n_max_qual"], 1, 40, 'Error: -N must be between 1 and 40')

def process_input_params(input_params): #, ploids):
    index_params = index_reference(input_params)
    # # parse input variants, if present
    # load_input_variants(input_params, ploids)
    # parse input targeted regions, if present
    load_input_regions(input_params, index_params["ref_list"])
    # parse discard bed similarly
    load_discard_regions(input_params)
    return index_params

def index_reference(input_params):
    print("Indexing reference started:", input_params["reference"])
    start = time.time()
    # index reference: [(0: chromosome name, 1: byte index where the contig seq begins,
    #                    2: byte index where the next contig begins, 3: contig seq length),
    #                    (repeat for every chrom)]
    # TODO check to see if this might work better as a dataframe or biopython object
    ref_index, line_width = index_ref(input_params["reference"])
    # print("ANOTHER TEST")
    # print("ref_index =")
    # print(ref_index)
    # TODO check if this index can work, maybe it's faster
    # ref_index2 = SeqIO.index(reference, 'fasta')
    index_params = {
        "ref_index": ref_index,
        "indices_by_ref_name": {ref_index[n][0]: n for n in range(len(ref_index))},
        "ref_list": [n[0] for n in ref_index],
        "line_width": line_width
    }
    end = time.time()
    print("Indexing reference took {} seconds.".format(int(end - start)))
    return index_params

def load_input_regions(input_params, ref_list):
    print("Loading input regions from BED file:", input_params["input_bed"])
    start = time.time()
    # TODO convert bed to pandas dataframe
    input_params["input_regions"] = {}
    if input_params["input_bed"] is not None:
        try:
            with open(input_params["input_bed"], 'r') as f:
                for line in f:
                    [my_chr, pos1, pos2] = line.strip().split('\t')[:3]
                    if my_chr not in input_params["input_regions"]:
                        input_params["input_regions"][my_chr] = [-1]
                    input_params["input_regions"][my_chr].extend([int(pos1), int(pos2)])
        except IOError:
            print("\nProblem reading input target BED file.\n")
            sys.exit(1)

        # some validation
        n_in_bed_only = 0
        n_in_ref_only = 0
        for k in ref_list:
            if k not in input_params["input_regions"]:
                n_in_ref_only += 1
        for k in input_params["input_regions"].keys():
            if k not in ref_list:
                n_in_bed_only += 1
                del input_params["input_regions"][k]
        if n_in_ref_only > 0:
            print('Warning: Reference contains sequences not found in targeted regions BED file.')
        if n_in_bed_only > 0:
            print(
                'Warning: Targeted regions BED file contains sequence names not found in reference (regions ignored).')
    end = time.time()
    print("Loading input regions took {} seconds.".format(int(end - start)))

def load_discard_regions(input_params):
    print("Loading discard regions from BED file:", input_params["discard_bed"])
    start = time.time()
    # TODO convert to pandas dataframe
    input_params["discard_regions"] = {}
    if input_params["discard_bed"] is not None:
        try:
            with open(input_params["discard_bed"], 'r') as f:
                for line in f:
                    [my_chr, pos1, pos2] = line.strip().split('\t')[:3]
                    if my_chr not in input_params["discard_bed"]:
                        input_params["discard_bed"][my_chr] = [-1]
                    input_params["discard_bed"][my_chr].extend([int(pos1), int(pos2)])
        except IOError:
            print("\nProblem reading discard BED file.\n")
            sys.exit(1)
    end = time.time()
    print("Loading discard regions took {} seconds.".format(int(end - start)))

def load_sequencing_model(sequencing_params):
    # absolute path to this script
    sim_path = pathlib.Path(__file__).resolve().parent
    # If user specified mean/std, or specified an empirical model, then the reads will be paired_ended
    # If not, then we're doing single-end reads.
    if (sequencing_params["fragment_size"] is not None and sequencing_params["fragment_std"] is not None) or (
            sequencing_params["fraglen_model"] is not None):
        sequencing_params["paired_end"] = True
    else:
        sequencing_params["paired_end"] = False
    # sequencing error model
    if sequencing_params["se_model"] is None:
        print('Using default sequencing error model.')
        sequencing_params["se_model"] = sim_path / 'models/errorModel_toy.p'
        sequencing_params["se_class"] = ReadContainer(sequencing_params["read_len"], sequencing_params["se_model"],
                                                      sequencing_params["se_rate"], sequencing_params["rescale_qual"])
    else:
        # probably need to do some sanity checking
        sequencing_params["se_class"] = ReadContainer(sequencing_params["read_len"], sequencing_params["se_model"],
                                                      sequencing_params["se_rate"], sequencing_params["rescale_qual"])
    # GC-bias model
    load_gc_model(sequencing_params, sim_path)
    # Assign appropriate values to the needed variables if we're dealing with paired-ended data
    if sequencing_params["paired_end"]:
        load_fraglen_model_PE(sequencing_params)
    if sequencing_params["paired_end"]:
        sequencing_params["n_handling"] = ('random', sequencing_params["fragment_size"])
    else:
        sequencing_params["n_handling"] = ('ignore', sequencing_params["read_len"])

def load_gc_model(sequencing_params, sim_path):
    if sequencing_params["gc_bias_model"] is None:
        print('Using default gc-bias model.')
        sequencing_params["gc_bias_model"] = sim_path / 'models/gcBias_toy.p'
        try:
            [gc_scale_count, gc_scale_val] = pickle.load(open(sequencing_params["gc_bias_model"], 'rb'))
        except IOError:
            print("\nProblem reading the default gc-bias model.\n")
            sys.exit(1)
        gc_window_size = gc_scale_count[-1]
    else:
        try:
            [gc_scale_count, gc_scale_val] = pickle.load(open(sequencing_params["gc_bias_model"], 'rb'))
        except IOError:
            print("\nProblem reading the gc-bias model.\n")
            sys.exit(1)
        gc_window_size = gc_scale_count[-1]
    sequencing_params["gc_scale_val"] = gc_scale_val
    sequencing_params["gc_window_size"] = gc_window_size

def load_fraglen_model_PE(sequencing_params):
    # Empirical fragment length distribution, if input model is specified
    if sequencing_params["fraglen_model"] is not None:
        print('Using empirical fragment length distribution.')
        try:
            [potential_values, potential_prob] = pickle.load(open(sequencing_params["fraglen_model"], 'rb'))
        except IOError:
            print('\nProblem loading the empirical fragment length model.\n')
            sys.exit(1)

        fraglen_values = []
        fraglen_probability = []
        for i in range(len(potential_values)):
            if potential_values[i] > sequencing_params["read_len"]:
                fraglen_values.append(potential_values[i])
                fraglen_probability.append(potential_prob[i])

        # TODO add some validation and sanity-checking code here...
        sequencing_params["fraglen_distribution"] = DiscreteDistribution(fraglen_probability, fraglen_values)
        sequencing_params["fragment_size"] = fraglen_values[mean_ind_of_weighted_list(fraglen_probability)]

    # Using artificial fragment length distribution, if the parameters were specified
    # fragment length distribution: normal distribution that goes out to +- 6 standard deviations
    elif sequencing_params["fragment_size"] is not None and sequencing_params["fragment_std"] is not None:
        print(
            'Using artificial fragment length distribution. mean=' + str(sequencing_params["fragment_size"]) + ', std=' + str(
                sequencing_params["fragment_std"]))
        if sequencing_params["fragment_std"] == 0:
            sequencing_params["fraglen_distribution"] = DiscreteDistribution([1], [sequencing_params["fragment_size"]], degenerate_val=sequencing_params["fragment_size"])
        else:
            potential_values = range(sequencing_params["fragment_size"] - 6 * sequencing_params["fragment_std"], sequencing_params["fragment_size"] + 6 * sequencing_params["fragment_std"] + 1)
            fraglen_values = []
            for i in range(len(potential_values)):
                if potential_values[i] > sequencing_params["read_len"]:
                    fraglen_values.append(potential_values[i])
            fraglen_probability = [np.exp(-(((n - float(sequencing_params["fragment_size"])) ** 2) / (2 * (sequencing_params["fragment_std"] ** 2)))) for n in
                                   fraglen_values]
            sequencing_params["fraglen_distribution"] = DiscreteDistribution(fraglen_probability, fraglen_values)

def load_mutation_model(mutation_params):
    mutation_params["mut_model"] = parse_input_mutation_model(mutation_params["mut_model"], 1)
    if mutation_params["mut_rate"] < 0.:
        mutation_params["mut_rate"] = None
    if mutation_params["mut_rate"] != -1 and mutation_params["mut_rate"] is not None:
        is_in_range(mutation_params["mut_rate"], 0.0, 1.0, 'Error: -M must be between 0 and 0.3')
    load_mutation_regions(mutation_params)

# parse input mutation rate rescaling regions, if present
def load_mutation_regions(mutation_params):
    # TODO convert to pandas dataframe
    mutation_params["mut_rate_regions"] = {}
    mutation_params["mut_rate_values"] = {}
    if mutation_params["mut_bed"] is not None:
        try:
            with open(mutation_params["mut_bed"], 'r') as f:
                for line in f:
                    [my_chr, pos1, pos2, meta_data] = line.strip().split('\t')[:4]
                    mut_str = re.findall(r"mut_rate=.*?(?=;)", meta_data + ';')
                    (pos1, pos2) = (int(pos1), int(pos2))
                    if len(mut_str) and (pos2 - pos1) > 1:
                        # mut_rate = #_mutations / length_of_region, let's bound it by a reasonable amount
                        mutation_params["mut_rate"] = max([0.0, min([float(mut_str[0][9:]), 0.3])])
                        if my_chr not in mutation_params["mut_rate_regions"]:
                            mutation_params["mut_rate_regions"][my_chr] = [-1]
                            mutation_params["mut_rate_values"][my_chr] = [0.0]
                        mutation_params["mut_rate_regions"][my_chr].extend([pos1, pos2])
                        # TODO figure out what the next line is supposed to do and fix
                        mutation_params["mut_rate_values"].extend([mutation_params["mut_rate"] * (pos2 - pos1)] * 2)
        except IOError:
            print("\nProblem reading mutational BED file.\n")
            sys.exit(1)

def intialize_reads_writers(index_params, input_params, output_params, sequencing_params):
    bam_file_writer = None
    if output_params["save_bam"]:
        bam_header = [
            copy.deepcopy(index_params["ref_index"])]  # TODO wondering if this is actually needed in the bam_header
        bam_file_writer = BamFileWriter(output_params["out_prefix"], bam_header)
    # fastq_file_writer = None
    # if not output_params["no_fastq"]:
    #     fastq_file_writer = FastqFileWriter(output_params["out_prefix"], paired=sequencing_params["paired_end"],
    #                                         no_fastq=output_params["no_fastq"])
    # else:
    #     print('Bypassing FASTQ generation...')
    fasta_file_writer = FastaFileWriter(output_params["out_prefix"], output_params["ploids"], index_params["line_width"])
    output_params["no_reads"] = output_params["no_fastq"] and not output_params["save_bam"]

    vcf_file_writer = None
    if output_params["save_vcf"]:
        vcf_header = [input_params["reference"]]
        vcf_file_writer = VcfFileWriter(output_params["out_prefix"],output_params["parent_prefix"],output_params["accession"],vcf_header)

    return bam_file_writer, fasta_file_writer, vcf_file_writer

def simulate_chrom(general_params, input_params, output_params, mutation_params, sequencing_params, index_params,
                       bam_file_writer, fasta_file_writer, vcf_file_writer, chrom, read_name_count, unmapped_records):

    # read in reference sequence and notate blocks of Ns
    index_params["ref_sequence"], index_params["n_regions"] = \
        read_ref(input_params["reference"], index_params["ref_index"][chrom], sequencing_params["n_handling"])
    progress_params = intialize_progress_bar_params(index_params["n_regions"])

    valid_variants_from_vcf = prune_invalid_variants(chrom, input_params["input_variants"], index_params["ref_index"], index_params["ref_sequence"])
    all_variants_out = {}
    sequences = None

    # TODO add large random structural variants

    load_sampling_window_params(sequencing_params)
    print('--------------------------------')
    print('Simulating chromosome {} of sequence started...'.format(chrom))
    t_start = time.time()
    # start the progress bar
    print('[', end='', flush=True)
    # Applying variants to non-N regions
    for i in range(len(index_params["n_regions"]['non_N'])):
        last_non_n_index = index_params["n_regions"]['non_N'][i-1][1] if i>0 else 0
        preliminary_Ns = index_params["n_regions"]['non_N'][i][0] - last_non_n_index
        write_fasta(fasta_file_writer, chrom, index_params["ref_index"], 0, None ,preliminary_Ns)
        sequences = apply_variants_to_region(general_params, input_params, output_params, mutation_params,
                                             sequencing_params, index_params, bam_file_writer, fasta_file_writer,
                                             chrom, read_name_count, unmapped_records, progress_params,
                                             valid_variants_from_vcf, all_variants_out, sequences, i)
    closing_Ns = index_params["ref_index"][chrom][3]-index_params["n_regions"]['non_N'][-1][1]
    write_fasta(fasta_file_writer, chrom, index_params["ref_index"], 0, None, closing_Ns)
    print(']', flush=True)

    # write all output variants for this reference
    if output_params["save_vcf"]:
        write_vcf(vcf_file_writer, all_variants_out, chrom, index_params["ref_index"])

    print("Simulating chromosome {} took {} seconds.".format(index_params["ref_index"][chrom][0], int(time.time() - t_start)))

def intialize_progress_bar_params(n_regions):
    # count total bp we'll be spanning so we can get an idea of how far along we are
    # (for printing progress indicators)
    progress_params = {
        "total_bp_span": sum([n[1] - n[0] for n in n_regions['non_N']]),
        "current_progress": 0,
        "current_percent": 0,
        "have_printed100": False
    }
    return progress_params

def prune_invalid_variants(chrom, input_variants, ref_index, ref_sequence):
    print("Pruning relevant variants from input VCF...")
    start = time.time()
    """Prune invalid input variants, e.g variants that:
                    - try to delete or alter any N characters
                    - don't match the reference base at their specified position
                    - any alt allele contains anything other than allowed characters"""
    valid_variants_from_vcf = []
    n_skipped = [0, 0, 0]
    if input_variants and ref_index[chrom][0] in input_variants:
        for n in input_variants[ref_index[chrom][0]]:
            span = (n[0], n[0] + len(n[1]))
            r_seq = str(ref_sequence[span[0] - 1:span[1] - 1])  # -1 because going from VCF coords to array coords
            # Checks if there are any invalid nucleotides in the vcf items
            any_bad_nucl = any((nn not in ALLOWED_NUCL) for nn in [item for sublist in n[2] for item in sublist])
            # Ensure reference sequence matches the nucleotide in the vcf
            if r_seq != n[1]:
                n_skipped[0] += 1
                continue
            # Ensure that we aren't trying to insert into an N region
            elif 'N' in r_seq:
                n_skipped[1] += 1
                continue
            # Ensure that we don't insert any disallowed characters
            elif any_bad_nucl:
                n_skipped[2] += 1
                continue
            # If it passes the above tests, append to valid variants list
            valid_variants_from_vcf.append(n)

        print('found', len(valid_variants_from_vcf), 'valid variants for ' +
              ref_index[chrom][0] + ' in input VCF...')
        if any(n_skipped):
            print(sum(n_skipped), 'variants skipped...')
            print(' - [' + str(n_skipped[0]) + '] ref allele does not match reference')
            print(' - [' + str(n_skipped[1]) + '] attempting to insert into N-region')
            print(' - [' + str(n_skipped[2]) + '] alt allele contains non-ACGT characters')
    end = time.time()
    print("Done. Pruning took {} seconds.".format(int(end - start)))
    return valid_variants_from_vcf

def load_sampling_window_params(sequencing_params):
    # determine sampling windows based on read length, large N regions, and structural mutations.
    # in order to obtain uniform coverage, windows should overlap by:
    # - read_len, if single-end reads
    # - fragment_size (mean), if paired-end reads
    # ploidy is fixed per large sampling window,
    # coverage distributions due to GC% and targeted regions are specified within these windows
    if sequencing_params["paired_end"]:
        sequencing_params["target_size"] = WINDOW_TARGET_SCALE * sequencing_params["fragment_size"]
        sequencing_params["overlap"] = sequencing_params["fragment_size"]
        sequencing_params["overlap_min_window_size"] = max(sequencing_params["fraglen_distribution"].values) + 10
    else:
        sequencing_params["target_size"] = WINDOW_TARGET_SCALE * sequencing_params["read_len"]
        sequencing_params["overlap"] = sequencing_params["read_len"]
        sequencing_params["overlap_min_window_size"] = sequencing_params["read_len"] + 10

def apply_variants_to_region(general_params, input_params, output_params, mutation_params,
                                             sequencing_params, index_params, bam_file_writer, fasta_file_writer,
                                             chrom, read_name_count, unmapped_records, progress_params,
                                             valid_variants_from_vcf, all_variants_out, sequences, i):
    (initial_position, final_position) = index_params["n_regions"]['non_N'][i]
    number_target_windows = max([1, (final_position - initial_position) // sequencing_params["target_size"]])
    base_pair_distance = int((final_position - initial_position) / float(number_target_windows))
    # if for some reason our region is too small to process, skip it! (sorry)
    if number_target_windows == 1 and (final_position - initial_position) < sequencing_params["overlap_min_window_size"]:
        return sequences
    start = initial_position
    end = min([start + base_pair_distance, final_position])
    vars_from_prev_overlap = []
    v_index_from_prev = 0
    is_last_time = False
    while True:
        # which inserted variants are in this window?
        buffer_added, vars_in_window, v_index_from_prev = get_vars_in_window(end, sequencing_params["overlap"], start,
                                                                             v_index_from_prev, valid_variants_from_vcf)

        buffer_added, end, is_last_time, next_end, next_start = adjust_window_to_vars(base_pair_distance, buffer_added,
                                                                                      end, final_position, is_last_time,
                                                                                      sequencing_params["overlap"],
                                                                                      vars_in_window)

        print_progress_indicator(buffer_added, progress_params, general_params["debug"], start, end, next_start,
                                 next_end, is_last_time)

        coverage_dat, target_hits = compute_coverage_modifiers(chrom, start, end, sequencing_params, input_params,
                                                               index_params["ref_index"])

        skip_this_window = should_skip_this_window(coverage_dat, end, sequencing_params, start, target_hits)

        if skip_this_window:
            # fill the window with Ns - a good practice?
            overlap_with_next = 0 if is_last_time else sequencing_params["overlap"]
            write_fasta(fasta_file_writer, chrom, index_params["ref_index"], overlap_with_next, sequences)
            # skip window, save cpu time
            if is_last_time:
                break
            start, end, is_last_time = prepare_next_window(next_start, next_end, final_position, is_last_time)
            vars_from_prev_overlap = []
            continue

        # construct sequence data that we will sample reads from
        all_inserted_variants, coverage_avg, sequences =\
            update_sequences(coverage_dat, end, output_params, mutation_params, sequencing_params,
                             index_params, sequences, start, vars_from_prev_overlap, vars_in_window)

        # which variants do we need to keep for next time (because of window overlap)?
        vars_from_prev_overlap = get_vars_for_next_window(all_inserted_variants, end, sequencing_params)

        overlap_with_next = 0 if is_last_time else sequencing_params["overlap"]
        write_fasta(fasta_file_writer, chrom, index_params["ref_index"], overlap_with_next, sequences)
        # if we're only producing VCF, no need to go through the hassle of generating reads
        if not output_params["no_reads"]:
            sample_reads(sequencing_params, input_params, index_params, output_params, bam_file_writer, chrom,
                         coverage_avg, coverage_dat, start, end, is_last_time, next_start,
                         read_name_count, sequences, unmapped_records)

        # tally up all the variants that got successfully introduced
        for n in all_inserted_variants:
            all_variants_out[n] = True

        if is_last_time:
            break
        # prepare indices of next window
        start, end, is_last_time = prepare_next_window(next_start, next_end, final_position, is_last_time)
    return sequences

def get_vars_in_window(end, overlap, start, v_index_from_prev, valid_variants_from_vcf):
    vars_in_window = []
    updated = False
    buffer_added = 0
    for j in range(v_index_from_prev, len(valid_variants_from_vcf)):
        variants_position = valid_variants_from_vcf[j][0]
        # update: changed <= to <, so variant cannot be inserted in first position
        if start < variants_position < end:
            # vcf --> array coords
            vars_in_window.append(tuple([variants_position - 1] + list(valid_variants_from_vcf[j][1:])))
        if variants_position >= end - overlap - 1 and updated is False:
            updated = True
            v_index_from_prev = j
        if variants_position >= end:
            break
    return buffer_added, vars_in_window, v_index_from_prev

def adjust_window_to_vars(base_pair_distance, buffer_added, end, final_position, is_last_time, overlap, vars_in_window):
    # determine which structural variants will affect our sampling window positions
    structural_vars = []
    for n in vars_in_window:
        # change: added abs() so that insertions are also buffered.
        buffer_needed = max([max([abs(len(n[1]) - len(alt_allele)), 1]) for alt_allele in n[2]])
        # -1 because going from VCF coords to array coords
        structural_vars.append((n[0] - 1, buffer_needed))
    # adjust end-position of window based on inserted structural mutations
    keep_going = True
    while keep_going:
        keep_going = False
        for n in structural_vars:
            # adding "overlap" here to prevent SVs from being introduced in overlap regions
            # (which can cause problems if random mutations from the previous window land on top of them)
            delta = (end - 1) - (n[0] + n[1]) - 2 - overlap
            if delta < 0:
                buffer_added = -delta
                end += buffer_added
                keep_going = True
                break
    next_start = end - overlap
    next_end = min([next_start + base_pair_distance, final_position])
    if next_end - next_start < base_pair_distance:
        end = next_end
        is_last_time = True
    return buffer_added, end, is_last_time, next_end, next_start

def print_progress_indicator(buffer_added, progress_params, debug, start, end, next_start, next_end, is_last_time):
    # print progress indicator
    if debug:
        print(f'PROCESSING WINDOW: {(start, end), [buffer_added]}, '
              f'next: {(next_start, next_end)}, isLastTime: {is_last_time}')
    progress_params["current_progress"] += end - start
    new_percent = int((progress_params["current_progress"] * 100) / float(progress_params["total_bp_span"]))
    if new_percent > progress_params["current_percent"]:
        if new_percent <= 99 or (new_percent == 100 and not progress_params["have_printed100"]):
            # if new_percent % 10 == 1:
            print('-', end='', flush=True)
        progress_params["current_percent"] = new_percent
        if progress_params["current_percent"] == 100:
            progress_params["have_printed100"] = True

def compute_coverage_modifiers(chrom, start, end, sequencing_params, input_params, ref_index):
    # compute coverage modifiers
    coverage_dat = [sequencing_params["gc_window_size"], sequencing_params["gc_scale_val"], []]
    target_hits = 0
    if input_params["input_bed"] is None:
        coverage_dat[2] = [1.0] * (end - start)
    else:
        if ref_index[chrom][0] not in input_params["input_regions"]:
            coverage_dat[2] = [sequencing_params["off_target_scalar"]] * (end - start)
        else:
            for j in range(start, end):
                if not (bisect.bisect(input_params["input_regions"][ref_index[chrom][0]], j) % 2):
                    coverage_dat[2].append(1.0)
                    target_hits += 1
                else:
                    coverage_dat[2].append(sequencing_params["off_target_scalar"])
    return coverage_dat, target_hits

def should_skip_this_window(coverage_dat, end, sequencing_params, start, target_hits):
    # off-target and we're not interested?
    if sequencing_params["off_target_discard"] and target_hits <= sequencing_params["read_len"]:
        return True
    # print len(coverage_dat[2]), sum(coverage_dat[2])
    if sum(coverage_dat[2]) < sequencing_params["low_cov_thresh"]:
        return True
    # check for small window sizes
    if (end - start) < sequencing_params["overlap_min_window_size"]:
        return True
    return False

def prepare_next_window(next_start, next_end, final_position, is_last_time):
    start = next_start
    end = next_end
    if end >= final_position:
        is_last_time = True
    return start, end, is_last_time

def update_sequences(coverage_dat, end, output_params, mutation_params, sequencing_params,
                             index_params, sequences, start, vars_from_prev_overlap, vars_in_window):
    if sequences is None:
        sequences = SequenceContainer(start, index_params["ref_sequence"][start:end], output_params["ploids"],
                                      sequencing_params["overlap"],sequencing_params["read_len"],
                                      [mutation_params["mut_model"]] * output_params["ploids"],
                                      mutation_params["mut_rate"], mutation_params["dist"],
                                      no_reads=output_params["no_reads"])
        # if [cigar for cigar in sequences.all_cigar[0] if len(cigar) != 100] or \
        #         [cig for cig in sequences.all_cigar[1] if len(cig) != 100]:
        if [hap for hap in range(output_params["ploids"]) if [cigar for cigar in sequences.all_cigar[hap] if len(cigar) != 100]]:
            print("There's a cigar that's off.")
            # pdb.set_trace()
            sys.exit(1)
    else:
        sequences.update(start, index_params["ref_sequence"][start:end], output_params["ploids"], sequencing_params["overlap"],
                         sequencing_params["read_len"], [mutation_params["mut_model"]] * output_params["ploids"],
                         mutation_params["mut_rate"], mutation_params["dist"])
        if [hap for hap in range(output_params["ploids"]) if [cigar for cigar in sequences.all_cigar[hap] if len(cigar) != 100]]:
            print("There's a cigar that's off.")
            # pdb.set_trace()
            sys.exit(1)
    # insert variants
    sequences.insert_mutations(vars_from_prev_overlap + vars_in_window)
    all_inserted_variants = sequences.random_mutations()
    # print all_inserted_variants
    # init coverage
    if sum(coverage_dat[2]) >= sequencing_params["low_cov_thresh"]:
        if sequencing_params["paired_end"]:
            coverage_avg = sequences.init_coverage(tuple(coverage_dat), frag_dist=sequencing_params["fraglen_distribution"])
        else:
            coverage_avg = sequences.init_coverage(tuple(coverage_dat))
    return all_inserted_variants, coverage_avg, sequences

def get_vars_for_next_window(all_inserted_variants, end, sequencing_params):
    vars_from_prev_overlap = []
    for n in all_inserted_variants:
        if n[0] >= end - sequencing_params["overlap"] - 1:
            vars_from_prev_overlap.append(n)
    return vars_from_prev_overlap

def sample_reads(sequencing_params, input_params, index_params, output_params, bam_file_writer, chrom, coverage_avg,
                         coverage_dat, start, end, is_last_time, next_start,
                         read_name_count, sequences, unmapped_records):

    reads_to_sample = compute_reads_to_sample(sequencing_params, coverage_avg, coverage_dat, start, end)
    # sample reads
    for k in range(reads_to_sample):

        my_read_data, is_unmapped = get_read_data(sequencing_params, sequences, start)
        if my_read_data is None:
            continue

        # are we discarding offtargets?
        outside_boundaries = get_outside_boundries(chrom, input_params, index_params, sequencing_params, my_read_data)
        if len(outside_boundaries) and any(outside_boundaries):
            continue

        my_read_name = output_params["out_prefix_name"] + '-' + index_params["ref_index"][chrom][0] + '-' + str(read_name_count)
        read_name_count += len(my_read_data)

        # if desired, replace all low-quality bases with Ns
        handle_low_quality(my_read_data, sequencing_params["n_max_qual"], sequencing_params["se_class"])

        # flip a coin, are we forward or reverse strand?
        is_forward = (random.random() < 0.5)

        # if read (or read + mate for PE) are unmapped, put them at end of bam file
        handle_unmapped_reads(is_forward, is_unmapped, my_read_data, my_read_name, sequencing_params["paired_end"], unmapped_records)

        my_ref_index = index_params["indices_by_ref_name"][index_params["ref_index"][chrom][0]]

        output_reads(bam_file_writer, is_forward, is_unmapped, my_read_data, my_read_name,
                     my_ref_index, output_params["no_fastq"], output_params["save_bam"])
    # if not output_params["no_fastq"]:
    #     fastq_file_writer.flush_buffer(is_last_time)
    if output_params["save_bam"]:
        bam_file_writer.flush_buffer(last_time=is_last_time, bam_max=(next_start if not is_last_time else end + 1))

def compute_reads_to_sample(sequencing_params, coverage_avg, coverage_dat, start, end):
    window_span = end - start
    if sequencing_params["paired_end"]:
        if sequencing_params["force_coverage"]:
            reads_to_sample = int((window_span * float(sequencing_params["coverage"])) / (2 * sequencing_params["read_len"])) + 1
        else:
            reads_to_sample = int((window_span * float(sequencing_params["coverage"]) * coverage_avg) / (2 * sequencing_params["read_len"])) + 1
    else:
        if sequencing_params["force_coverage"]:
            reads_to_sample = int((window_span * float(sequencing_params["coverage"])) / sequencing_params["read_len"]) + 1
        else:
            reads_to_sample = int((window_span * float(sequencing_params["coverage"]) * coverage_avg) / sequencing_params["read_len"]) + 1
    # if coverage is so low such that no reads are to be sampled, skip region
    #      (i.e., remove buffer of +1 reads we add to every window)
    if reads_to_sample == 1 and sum(coverage_dat[2]) < sequencing_params["low_cov_thresh"]:
        reads_to_sample = 0
    return reads_to_sample

def get_read_data(sequencing_params, sequences, start):
    is_unmapped = []
    if sequencing_params["paired_end"]:
        my_fraglen = sequencing_params["fraglen_distribution"].sample()
        my_read_data = sequences.sample_read(sequencing_params["se_class"], my_fraglen)
        # skip if we failed to find a valid position to sample read
        if my_read_data is None:
            return my_read_data, is_unmapped
        if my_read_data[0][0] is None:
            is_unmapped.append(True)
        else:
            is_unmapped.append(False)
            # adjust mapping position based on window start
            my_read_data[0][0] += start
        if my_read_data[1][0] is None:
            is_unmapped.append(True)
        else:
            is_unmapped.append(False)
            my_read_data[1][0] += start
    else:
        my_read_data = sequences.sample_read(sequencing_params["se_class"])
        # skip if we failed to find a valid position to sample read
        if my_read_data is None:
            return my_read_data, is_unmapped
        # unmapped read (lives in large insertion)
        if my_read_data[0][0] is None:
            is_unmapped = [True]
        else:
            is_unmapped = [False]
            # adjust mapping position based on window start
            my_read_data[0][0] += start
    return my_read_data, is_unmapped

def get_outside_boundries(chrom, input_params, index_params, sequencing_params, my_read_data):
    outside_boundaries = []
    if sequencing_params["off_target_discard"] and input_params["input_bed"] is not None:
        outside_boundaries += [bisect.bisect(input_params["input_regions"][index_params["ref_index"][chrom][0]], n[0]) % 2 for n
                               in my_read_data]
        outside_boundaries += [
            bisect.bisect(input_params["input_regions"][index_params["ref_index"][chrom][0]], n[0] + len(n[2])) % 2 for n in
            my_read_data]
    if input_params["discard_bed"] is not None:
        outside_boundaries += [bisect.bisect(input_params["discard_bed"][index_params["ref_index"][chrom][0]], n[0]) % 2 for
                               n in my_read_data]
        outside_boundaries += [
            bisect.bisect(input_params["discard_bed"][index_params["ref_index"][chrom][0]], n[0] + len(n[2])) % 2 for n in
            my_read_data]
    return outside_boundaries

def handle_low_quality(my_read_data, n_max_qual, se_class):
    if n_max_qual > -1:
        for j in range(len(my_read_data)):
            my_read_string = [n for n in my_read_data[j][2]]
            for m in range(len(my_read_data[j][3])):
                adjusted_qual = ord(my_read_data[j][3][m]) - se_class.off_q
                if adjusted_qual <= n_max_qual:
                    my_read_string[m] = 'N'
            my_read_data[j][2] = ''.join(my_read_string)

def handle_unmapped_reads(is_forward, is_unmapped, my_read_data, my_read_name, paired_end, unmapped_records):
    if all(is_unmapped):
        if paired_end:
            if is_forward:
                flag1 = sam_flag(['paired', 'unmapped', 'mate_unmapped', 'first', 'mate_reverse'])
                flag2 = sam_flag(['paired', 'unmapped', 'mate_unmapped', 'second', 'reverse'])
            else:
                flag1 = sam_flag(['paired', 'unmapped', 'mate_unmapped', 'second', 'mate_reverse'])
                flag2 = sam_flag(['paired', 'unmapped', 'mate_unmapped', 'first', 'reverse'])
            unmapped_records.append((my_read_name + '/1', my_read_data[0], flag1))
            unmapped_records.append((my_read_name + '/2', my_read_data[1], flag2))
        else:
            flag1 = sam_flag(['unmapped'])
            unmapped_records.append((my_read_name + '/1', my_read_data[0], flag1))

def output_reads(bam_file_writer, is_forward, is_unmapped, my_read_data, my_read_name, my_ref_index,
                 no_fastq, save_bam):
    # write SE output
    if len(my_read_data) == 1:
        write_SE_output(bam_file_writer, is_forward, is_unmapped, my_read_data,
                        my_read_name, my_ref_index, no_fastq, save_bam)
    # write PE output
    elif len(my_read_data) == 2:
        write_PE_output(bam_file_writer, is_forward, is_unmapped, my_read_data,
                        my_read_name, my_ref_index, no_fastq, save_bam)
    else:
        print('\nError: Unexpected number of reads generated...\n')
        sys.exit(1)

def write_vcf(vcf_file_writer, all_variants_out, chrom, ref_index):
    print('Writing output VCF started...')
    start = time.time()
    for k in sorted(all_variants_out.keys()):
        current_ref = ref_index[chrom][0]
        my_id = '.'
        my_quality = '.'
        my_filter = 'PASS'
        # k[0] + 1 because we're going back to 1-based vcf coords
        vcf_file_writer.write_record(current_ref, str(int(k[0]) + 1), my_id, k[1], k[2], my_quality,
                                     my_filter, k[4])
    end = time.time()
    print("Done. Writing output VCF took {} seconds.".format(int(end - start)))

def write_fasta(fasta_file_writer, chrom, ref_index, overlap_with_next, sequences, N_seq_len=0):
    haploid_sequences = sequences.sequences if sequences else None
    fasta_file_writer.write_record(haploid_sequences, ref_index[chrom][0], N_seq_len, overlap_with_next)

def write_SE_output(bam_file_writer, is_forward, is_unmapped, my_read_data, my_read_name,
                    my_ref_index, no_fastq, save_bam):
    # if not no_fastq:
    #     write_fastq_SE(fastq_file_writer, is_forward, my_read_data, my_read_name)
    if save_bam:
        write_SE_to_bam(bam_file_writer, is_forward, is_unmapped, my_read_data, my_read_name, my_ref_index)

# def write_fastq_SE(fastq_file_writer, is_forward, my_read_data, my_read_name):
#     if is_forward:
#         fastq_file_writer.write_record(my_read_name, my_read_data[0][2],
#                                        my_read_data[0][3])
#     else:
#         fastq_file_writer.write_record(my_read_name,
#                                        reverse_complement(my_read_data[0][2]),
#                                        my_read_data[0][3][::-1])

def write_SE_to_bam(bam_file_writer, is_forward, is_unmapped, my_read_data, my_read_name, my_ref_index):
    if is_unmapped[0] is False:
        if is_forward:
            flag1 = 0
            bam_file_writer.write_record(my_ref_index, my_read_name,
                                         my_read_data[0][0],
                                         my_read_data[0][1], my_read_data[0][2],
                                         my_read_data[0][3],
                                         output_sam_flag=flag1)
        else:
            flag1 = sam_flag(['reverse'])
            bam_file_writer.write_record(my_ref_index, my_read_name,
                                         my_read_data[0][0],
                                         my_read_data[0][1], my_read_data[0][2],
                                         my_read_data[0][3],
                                         output_sam_flag=flag1)

def write_PE_output(bam_file_writer, is_forward, is_unmapped, my_read_data, my_read_name,
                    my_ref_index, no_fastq, save_bam):
    # if not no_fastq:
    #     write_fastq_PE(fastq_file_writer, is_forward, my_read_data, my_read_name)
    if save_bam:
        write_PE_to_bam(bam_file_writer, is_forward, is_unmapped, my_read_data, my_read_name, my_ref_index)

# def write_fastq_PE(fastq_file_writer, is_forward, my_read_data, my_read_name):
#     fastq_file_writer.write_record(my_read_name, my_read_data[0][2],
#                                    my_read_data[0][3],
#                                    read2=my_read_data[1][2],
#                                    qual2=my_read_data[1][3],
#                                    orientation=is_forward)

def write_PE_to_bam(bam_file_writer, is_forward, is_unmapped, my_read_data, my_read_name, my_ref_index):
    if is_unmapped[0] is False and is_unmapped[1] is False:
        if is_forward:
            flag1 = sam_flag(['paired', 'proper', 'first', 'mate_reverse'])
            flag2 = sam_flag(['paired', 'proper', 'second', 'reverse'])
        else:
            flag1 = sam_flag(['paired', 'proper', 'second', 'mate_reverse'])
            flag2 = sam_flag(['paired', 'proper', 'first', 'reverse'])
        bam_file_writer.write_record(my_ref_index, my_read_name, my_read_data[0][0],
                                     my_read_data[0][1], my_read_data[0][2],
                                     my_read_data[0][3],
                                     output_sam_flag=flag1,
                                     mate_pos=my_read_data[1][0])
        bam_file_writer.write_record(my_ref_index, my_read_name, my_read_data[1][0],
                                     my_read_data[1][1], my_read_data[1][2],
                                     my_read_data[1][3],
                                     output_sam_flag=flag2, mate_pos=my_read_data[0][0])
    elif is_unmapped[0] is False and is_unmapped[1] is True:
        if is_forward:
            flag1 = sam_flag(['paired', 'first', 'mate_unmapped', 'mate_reverse'])
            flag2 = sam_flag(['paired', 'second', 'unmapped', 'reverse'])
        else:
            flag1 = sam_flag(['paired', 'second', 'mate_unmapped', 'mate_reverse'])
            flag2 = sam_flag(['paired', 'first', 'unmapped', 'reverse'])
        bam_file_writer.write_record(my_ref_index, my_read_name, my_read_data[0][0],
                                     my_read_data[0][1], my_read_data[0][2],
                                     my_read_data[0][3],
                                     output_sam_flag=flag1, mate_pos=my_read_data[0][0])
        bam_file_writer.write_record(my_ref_index, my_read_name, my_read_data[0][0],
                                     my_read_data[1][1], my_read_data[1][2],
                                     my_read_data[1][3],
                                     output_sam_flag=flag2, mate_pos=my_read_data[0][0],
                                     aln_map_quality=0)
    elif is_unmapped[0] is True and is_unmapped[1] is False:
        if is_forward:
            flag1 = sam_flag(['paired', 'first', 'unmapped', 'mate_reverse'])
            flag2 = sam_flag(['paired', 'second', 'mate_unmapped', 'reverse'])
        else:
            flag1 = sam_flag(['paired', 'second', 'unmapped', 'mate_reverse'])
            flag2 = sam_flag(['paired', 'first', 'mate_unmapped', 'reverse'])
        bam_file_writer.write_record(my_ref_index, my_read_name, my_read_data[1][0],
                                     my_read_data[0][1], my_read_data[0][2],
                                     my_read_data[0][3],
                                     output_sam_flag=flag1, mate_pos=my_read_data[1][0],
                                     aln_map_quality=0)
        bam_file_writer.write_record(my_ref_index, my_read_name, my_read_data[1][0],
                                     my_read_data[1][1], my_read_data[1][2],
                                     my_read_data[1][3],
                                     output_sam_flag=flag2, mate_pos=my_read_data[1][0])

def write_unmapped_to_bam(bam_file_writer, paired_end, save_bam, unmapped_records):
    if save_bam and len(unmapped_records):
        print('writing unmapped reads to bam file...')
        for umr in unmapped_records:
            bam_file_writer.write_record(-1, umr[0], 0, umr[1][1], umr[1][2], umr[1][3], output_sam_flag=umr[2],
                                         mate_pos=0 if paired_end else None, aln_map_quality=0)