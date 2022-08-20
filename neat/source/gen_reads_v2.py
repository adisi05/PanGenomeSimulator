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

import bisect
import copy
import pathlib
import pickle
import random
import re
import sys
import time

import numpy as np
import os
import pandas as pd
from pybedtools import BedTool

from neat.source.fastq_file_writer import FastqFileWriter
from source.SequenceContainer import SequenceContainer, ReadContainer, parse_input_mutation_model
from source.bam_file_writer import BamFileWriter
from source.bam_file_writer import sam_flag
from source.fasta_file_writer import FastaFileWriter
# from source.fastq_file_writer import FastqFileWriter
from source.input_checking import check_file_open, is_in_range
from source.probability import DiscreteDistribution, mean_ind_of_weighted_list
from source.ref_func import index_ref, read_ref
from source.vcf_file_writer import VcfFileWriter

"""
Some constants needed for analysis
"""

# target window size for read sampling. How many times bigger than read/frag length
WINDOW_TARGET_SCALE = 100

# allowed nucleotides
ALLOWED_NUCL = ['A', 'C', 'G', 'T']


def simulate(args):
    general_params, input_params, output_params, mutation_params, sequencing_params = parse_args(args)
    index_params = process_input_params(input_params)  # , output_params["ploids"])
    load_mutation_model(mutation_params)

    # initialize output writers
    fasta_file_writer, vcf_file_writer = intialize_reads_writers(index_params,input_params, output_params, sequencing_params)

    # Using pathlib to make this more machine agnostic
    output_params["out_prefix_name"] = pathlib.Path(output_params["out_prefix"]).name

    for chrom in range(len(index_params["ref_index"])):
        simulate_chrom(general_params, input_params, output_params, mutation_params, index_params, vcf_file_writer, chrom)

    # TODO translocation feature

    # close output files
    if not output_params["no_fastq"]:
        FastqFileWriter.generate_reads(fasta_file_writer.get_file(), sequencing_params)
    if output_params["save_vcf"]:
        vcf_file_writer.close_file(add_parent_variants=True)

    # TODO write mut_bed

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
        "input_variants_path": args.input_variants_path,
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
        "out_prefix": args.o + '_' + args.name,
        "parent_prefix": args.o + '_' + args.parent_name if args.parent_name else None,
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
    check_file_open(input_params["reference"], 'ERROR: could not open reference, {}'.format(input_params["reference"]),
                    required=True)
    # check_file_open(input_params["input_vcf"], 'ERROR: could not open input VCF, {}'.format(input_params["input_vcf"]), required=False)
    check_file_open(input_params["input_bed"], 'ERROR: could not open input BED, {}'.format(input_params["input_bed"]),
                    required=False)
    # if user specified no fastq, not fasta only, and no bam and no vcf, then print error and exit.
    if output_params["no_fastq"] and not output_params["save_bam"] and not output_params["save_vcf"]:
        print('\nERROR: No files would be written.\n')
        sys.exit(1)
    if (sequencing_params["fragment_size"] is None and sequencing_params["fragment_std"] is not None) or (
            sequencing_params["fragment_size"] is not None and sequencing_params["fragment_std"] is None):
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


def process_input_params(input_params):  # , ploids):
    index_params = index_reference(input_params)
    # # parse input variants, if present
    # load_input_variants(input_params, ploids)
    # parse input targeted regions, if present
    load_input_regions(input_params, index_params["ref_list"])
    # parse discard bed similarly
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


def load_mutation_model(mutation_params):
    mutation_params["mut_model"] = parse_input_mutation_model(mutation_params["mut_model"], 1)
    if mutation_params["mut_rate"] < 0.:
        mutation_params["mut_rate"] = None
    if mutation_params["mut_rate"] != -1 and mutation_params["mut_rate"] is not None:
        is_in_range(mutation_params["mut_rate"], 0.0, 1.0, 'Error: -M must be between 0 and 0.3')

def intialize_reads_writers(index_params, input_params, output_params):
    fasta_file_writer = FastaFileWriter(output_params["out_prefix"], output_params["ploids"],
                                        index_params["line_width"])
    output_params["no_reads"] = output_params["no_fastq"]

    vcf_file_writer = None
    if output_params["save_vcf"]:
        vcf_header = [input_params["reference"]]
        vcf_file_writer = VcfFileWriter(output_params["out_prefix"], output_params["parent_prefix"],
                                        output_params["accession"], vcf_header)

    return fasta_file_writer, vcf_file_writer


def simulate_chrom(general_params, input_params, output_params, mutation_params, index_params, chrom):

    # read in reference sequence and notate blocks of Ns
    chrom_sequence, index_params["n_regions"] = read_ref(input_params["reference"], index_params["ref_index"][chrom])

    variants_from_vcf = get_input_variants_from_vcf(input_params)
    current_chrom_given_valid_variants = prune_invalid_variants(chrom, variants_from_vcf, index_params["ref_index"],
                                                     chrom_sequence)
    all_variants_out = {}
    sequences = SequenceContainer(chrom_sequence, output_params["ploids"],
                                      [mutation_params["mut_model"]] * output_params["ploids"],
                                      mutation_params["mut_rate"], mutation_params["dist"])

    #TODO add large random structural variants

    print('--------------------------------')
    print('Simulating chromosome {} of sequence started...'.format(chrom))
    t_start = time.time()
    # Applying variants to non-N regions
    for non_n_region in index_params['n_regions']['non_N']:
        start, end = non_n_region
        sequences = apply_variants_to_non_n_region(output_params, mutation_params, chrom, current_chrom_given_valid_variants, all_variants_out, sequences, start, end)

    # write all output variants for this reference
    #TODO write all_variants_out

    print("Simulating chromosome {} took {} seconds.".format(index_params["ref_index"][chrom][0], int(time.time() - t_start)))


def get_input_variants_from_vcf(input_params):
    if input_params["input_variants"] is None:
        variants_from_vcf = pd.read_csv(input_params["input_variants_path"])
        variants_from_vcf[['chrom', 'allele']] = variants_from_vcf[['chrom', 'allele']].astype(str)
        variants_from_vcf['pos'] = variants_from_vcf['pos'].astype(int)
        os.remove(input_params["input_variants_path"])
        input_params["input_variants"] = variants_from_vcf
    return input_params["input_variants"]


def prune_invalid_variants(chrom, input_variants, ref_index, chrom_sequence):
    print("Pruning relevant variants from input VCF...")
    start = time.time()
    """Prune invalid input variants, e.g variants that:
                    - try to delete or alter any N characters
                    - don't match the reference base at their specified position
                    - any alt allele contains anything other than allowed characters"""
    valid_variants_from_vcf_indexes = []
    n_skipped = [0, 0, 0]
    if (not input_variants.empty) and (ref_index[chrom][0] in input_variants.chrom.unique()):
        for index, variant in input_variants[input_variants.chrom == ref_index[chrom][0]].iterrows():
            span = (variant.pos, variant.pos + len(variant.allele))
            r_seq = str(chrom_sequence[span[0] - 1:span[1] - 1])  # -1 because going from VCF coords to array coords
            # Checks if there are any invalid nucleotides in the vcf items
            all_alternatives_nuleotides = [nuc for alt in variant.alternatives for nuc in alt]
            any_bad_nucl = any([(nuc not in ALLOWED_NUCL) for nuc in all_alternatives_nuleotides])
            # Ensure reference sequence matches the nucleotide in the vcf
            if r_seq != variant.allele:
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
            valid_variants_from_vcf_indexes.append(index)

        print('found', len(valid_variants_from_vcf_indexes), 'valid variants for ' +
              ref_index[chrom][0] + ' in input VCF...')
        if any(n_skipped):
            print(sum(n_skipped), 'variants skipped...')
            print(' - [' + str(n_skipped[0]) + '] ref allele does not match reference')
            print(' - [' + str(n_skipped[1]) + '] attempting to insert into N-region')
            print(' - [' + str(n_skipped[2]) + '] alt allele contains non-ACGT characters')
    end = time.time()
    print("Done. Pruning took {} seconds.".format(int(end - start)))
    return input_variants.loc[valid_variants_from_vcf_indexes, 'pos':].reset_index(drop=True)



def apply_variants_to_non_n_region(output_params, mutation_params, chrom,
                                   valid_variants_from_vcf, all_variants_out, sequences, start, end):

    vars_in_current_window = get_vars_in_window(start, end, valid_variants_from_vcf)

    # construct sequence data that we will sample reads from
    inserted_random_variants, sequences = \
        update_sequences(chrom, start, end, sequences, vars_in_current_window, output_params, mutation_params)

    return sequences


def get_vars_in_window(start, end, valid_variants_from_vcf):
    vars_in_window_indexes = []

    for i, variant in valid_variants_from_vcf.iterrows():
        if start < variant['pos'] < end:
            vars_in_window_indexes.append(i)

    if not valid_variants_from_vcf.empty:
        vars_in_window = valid_variants_from_vcf.loc[vars_in_window_indexes, 'pos':].reset_index(drop=True)
    else:
        vars_in_window = valid_variants_from_vcf
    return vars_in_window



def update_sequences(start, end, sequences, vars_in_current_window):
    # sequences.focus_on_window(start, end) done internaly? #TODO decide
    # insert variants
    sequences.insert_given_mutations(vars_in_current_window)
    inserted_random_variants = sequences.insert_random_mutations(start, end)
    return inserted_random_variants, sequences


def get_vars_for_next_window(all_inserted_variants, end, sequencing_params):
    vars_from_prev_overlap = []
    for n in all_inserted_variants:
        if n[0] >= end - sequencing_params["overlap"] - 1:
            vars_from_prev_overlap.append(n)
    return vars_from_prev_overlap



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