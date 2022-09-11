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

import pathlib
import random
import sys
import time

import os
from typing import List, Dict

import pandas as pd

from neat.source.fastq_file_writer import FastqFileWriter
from source.chromosome_processor import ChromosomeProcessor, parse_input_mutation_model
from source.fasta_file_writer import FastaFileWriter
from source.input_checking import check_file_open, is_in_range
from source.ref_func import index_ref, read_ref
from source.vcf_file_writer import VcfFileWriter
from utilities.annotated_sequence import to_annotations_df

"""
Some constants needed for analysis
"""

# target window size for read sampling. How many times bigger than read/frag length
WINDOW_TARGET_SCALE = 100

# allowed nucleotides
ALLOWED_NUCL = ['A', 'C', 'G', 'T']


def simulate(args):
    general_params, input_params, output_params, mutation_params, sequencing_params = parse_args(args)
    index_params = process_input_params(input_params)
    load_mutation_model(mutation_params)
    annotations_df = to_annotations_df(mutation_params["mut_bed"], os.getcwd())

    # Using pathlib to make this more machine agnostic
    output_params["out_prefix_name"] = pathlib.Path(output_params["out_prefix"]).name


    final_chromosomes = {}
    for chrom in index_params["ref_list"]:
        chromosome_processor = simulate_chrom(input_params, output_params, mutation_params, index_params, sequencing_params, chrom, annotations_df)
        final_chromosomes[chromosome_processor.chromosome_name] = str(chromosome_processor.chromosome_sequence)

    write_output(index_params, input_params, output_params, sequencing_params, final_chromosomes)

    # TODO translocation feature

    # TODO close/finalize writer in some way?



def write_output(index_params, input_params, output_params, sequencing_params, final_chromosomes : Dict[str,str]):

    # FASTA
    fasta_file_writer = FastaFileWriter(output_params["out_prefix"], index_params["line_width"])
    for name, sequence in final_chromosomes.items():
        fasta_file_writer.write_record(sequence, name)
    fasta_file_writer.finalize()

    # FASTQ
    if not output_params["no_fastq"]:
        FastqFileWriter.generate_reads([fasta_file_writer.get_file_name()], sequencing_params)

    # VCF
    if output_params["save_vcf"]:
        vcf_header = [input_params["reference"]]
        vcf_file_writer = VcfFileWriter(output_params["out_prefix"], output_params["parent_prefix"],
                                        output_params["accession"], vcf_header)
        vcf_file_writer.close_file(add_parent_variants=True)


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
    mut_model = args.m
    mut_model = mut_model.replace('"', '').replace("'","")
    mutation_params = {
        "mut_bed": args.Mb,
        "mut_model": mut_model,
        "mut_rate": args.M,
        "dist": args.dist
    }
    output_params = {
        "accession": args.name,
        "out_prefix": args.o + '_' + args.name,
        "parent_prefix": args.o + '_' + args.parent_name if args.parent_name else None,
        "save_vcf": args.vcf,
        "no_fastq": args.no_fastq or args.internal
    }
    general_params = {
        "rng_seed": args.rng,
        "debug": args.d
    }
    read_len = args.R
    coverage = args.c
    (fragment_size, fragment_std) = args.pe
    paired_end = (fragment_size and fragment_std) or args.pe_model
    n_handling = ('random', fragment_size) if paired_end else ('ignore', read_len)
    sequencing_params = {
        "read_len": read_len,
        "coverage": coverage,
        "fragment_size": fragment_size,
        "fragment_std": fragment_std,
        "paired_end": paired_end,
        "n_handling": n_handling
    }
    return general_params, input_params, output_params, mutation_params, sequencing_params


def params_sanity_check(input_params, output_params, general_params, sequencing_params):
    # Check that files are real, if provided
    check_file_open(input_params["reference"], 'ERROR: could not open reference, {}'.format(input_params["reference"]),
                    required=True)
    # check_file_open(input_params["input_vcf"], 'ERROR: could not open input VCF, {}'.format(input_params["input_vcf"]), required=False)
    check_file_open(input_params["input_bed"], 'ERROR: could not open input BED, {}'.format(input_params["input_bed"]),
                    required=False)
    # TODO check mut bed file, and vcf?
    if (sequencing_params["fragment_size"] is None and sequencing_params["fragment_std"] is not None) or (
            sequencing_params["fragment_size"] is not None and sequencing_params["fragment_std"] is None):
        print('\nERROR: --pe argument takes 2 space-separated arguments.\n')
        sys.exit(1)
    if general_params["rng_seed"] == -1:
        general_params["rng_seed"] = random.randint(1, 99999999)
    random.seed(general_params["rng_seed"])
    is_in_range(sequencing_params["read_len"], 10, 1000000, 'Error: -R must be between 10 and 1,000,000')
    is_in_range(sequencing_params["coverage"], 0, 1000000, 'Error: -c must be between 0 and 1,000,000')

def process_input_params(input_params):
    index_params = index_reference(input_params)
    # # parse input variants, if present
    # load_input_variants(input_params, ploids)
    # parse input targeted regions, if present
    # load_input_regions(input_params, index_params["ref_list"])
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
        "indices_by_ref_name": {chrom[0]: chrom for chrom in ref_index},  # TODO chrom[1:] ?
        "ref_list": [chrom[0] for chrom in ref_index],
        "line_width": line_width
    }
    end = time.time()
    print("Indexing reference took {} seconds.".format(int(end - start)))
    return index_params


# def load_input_regions(input_params, ref_list):
#     print("Loading input regions from BED file:", input_params["input_bed"])
#     start = time.time()
#     # TODO convert bed to pandas dataframe
#     input_params["input_regions"] = {}
#     if input_params["input_bed"] is not None:
#         try:
#             with open(input_params["input_bed"], 'r') as f:
#                 for line in f:
#                     [my_chr, pos1, pos2] = line.strip().split('\t')[:3]
#                     if my_chr not in input_params["input_regions"]:
#                         input_params["input_regions"][my_chr] = [-1]
#                     input_params["input_regions"][my_chr].extend([int(pos1), int(pos2)])
#         except IOError:
#             print("\nProblem reading input target BED file.\n")
#             sys.exit(1)
#
#         # some validation
#         n_in_bed_only = 0
#         n_in_ref_only = 0
#         for k in ref_list:
#             if k not in input_params["input_regions"]:
#                 n_in_ref_only += 1
#         for k in input_params["input_regions"].keys():
#             if k not in ref_list:
#                 n_in_bed_only += 1
#                 del input_params["input_regions"][k]
#         if n_in_ref_only > 0:
#             print('Warning: Reference contains sequences not found in targeted regions BED file.')
#         if n_in_bed_only > 0:
#             print(
#                 'Warning: Targeted regions BED file contains sequence names not found in reference (regions ignored).')
#     end = time.time()
#     print("Loading input regions took {} seconds.".format(int(end - start)))


def load_mutation_model(mutation_params):
    mutation_params["mut_model"] = parse_input_mutation_model(mutation_params["mut_model"], 1)
    if mutation_params["mut_rate"] < 0.:
        mutation_params["mut_rate"] = None
    if mutation_params["mut_rate"] != -1 and mutation_params["mut_rate"] is not None:
        is_in_range(mutation_params["mut_rate"], 0.0, 1.0, 'Error: -M must be between 0 and 0.3')


def simulate_chrom(input_params, output_params, mutation_params, index_params, sequencing_params, chrom, annotations_df) -> ChromosomeProcessor:

    # read in reference sequence and notate blocks of Ns
    chrom_sequence, index_params["n_regions"] = read_ref(input_params["reference"], index_params['indices_by_ref_name'][chrom], sequencing_params["n_handling"])

    variants_from_vcf = get_input_variants_from_vcf(input_params)
    current_chrom_given_valid_variants = prune_invalid_variants(chrom, variants_from_vcf, index_params["ref_index"],
                                                     chrom_sequence)
    chrom_annotations_df = annotations_df[annotations_df['chrom']==chrom]
    chromosome_processor = ChromosomeProcessor(chrom, chrom_sequence, chrom_annotations_df, annotations_sorted=True, mut_models=mutation_params["mut_model"], mut_rate = mutation_params["mut_rate"], dist = mutation_params["dist"])

    #TODO add large random structural variants

    print('--------------------------------')
    print('Simulating chromosome {} of sequence started...'.format(chrom))
    t_start = time.time()
    # Applying variants to non-N regions
    for non_n_region in index_params['n_regions']['non_N']:
        start, end = non_n_region
        inserted_mutations = apply_variants_to_non_n_region(current_chrom_given_valid_variants, chromosome_processor, start, end)

    # write all output variants for this reference
    #TODO write inserted_mutations - write_vcf?

    print(f"Simulating chromosome {chrom} took {int(time.time() - t_start)} seconds.")

    return chromosome_processor


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
    if (not input_variants.empty) and (chrom in input_variants.chrom.unique()):
        for index, variant in input_variants[input_variants.chrom == chrom].iterrows():
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
              chrom + ' in input VCF...')
        if any(n_skipped):
            print(sum(n_skipped), 'variants skipped...')
            print(' - [' + str(n_skipped[0]) + '] ref allele does not match reference')
            print(' - [' + str(n_skipped[1]) + '] attempting to insert into N-region')
            print(' - [' + str(n_skipped[2]) + '] alt allele contains non-ACGT characters')
    end = time.time()
    print("Done. Pruning took {} seconds.".format(int(end - start)))
    return input_variants.loc[valid_variants_from_vcf_indexes, 'pos':].reset_index(drop=True)



def apply_variants_to_non_n_region(valid_variants_from_vcf, chromosome_processor, start, end):

    vars_in_current_window = get_vars_in_window(start, end, valid_variants_from_vcf)

    # construct sequence data that we will sample reads from
    inserted_mutations =\
        insert_given_and_random_mutations(start, end, chromosome_processor, vars_in_current_window)


    return inserted_mutations


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



def insert_given_and_random_mutations(start, end, chromosome_processor, vars_in_current_window):
    given_mutations_inserted = chromosome_processor.insert_given_mutations(vars_in_current_window, start=start, end=end)
    random_mutations_inserted = chromosome_processor.generate_random_mutations(start, end)
    return given_mutations_inserted + random_mutations_inserted
    # TODO return random_mutations_inserted only?




def write_vcf(vcf_file_writer, inserted_mutations, chrom, ref_index):
    print('Writing output VCF started...')
    start = time.time()
    for k in sorted(inserted_mutations.keys()):
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