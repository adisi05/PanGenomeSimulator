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
from typing import List, Dict, Tuple, Optional

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


class GenomeSimulator:
    def __init__(self, args):
        self._extract_args(args)
        self._params_sanity_check()

        # If coverage val for a given window/position is below this value, consider it effectively zero.
        self._sequencing_low_cov_thresh = 50

        # Use pathlib to make this more machine agnostic
        self._output_prefix_name = pathlib.Path(self._output_prefix).name

        # Load annotations dataframe
        if not self._workdir:
            workdir = os.path.dirname(self._output_prefix)
        self._annotations_df = to_annotations_df(self._annotations_file, workdir)

        self._load_mutation_model()

        self._index_reference()

        self._load_vcf()

    def _extract_args(self, args):
        self._input_reference = args.r
        self._input_variants = args.input_variants
        self._input_variants_path = args.input_variants_path
        self._annotations_file = args.Mb
        self._mutation_model = args.m
        self._mutation_rate = args.M
        self._relative_distance = args.dist
        self._output_accession = args.name
        self._output_prefix = args.o + '_' + args.name
        self._parent_prefix = args.o + '_' + args.parent_name if args.parent_name else None
        self._output_vcf = args.vcf
        self._output_no_fastq = args.no_fastq or args.internal
        self._workdir = args.w if args.w else None
        self._rng_seed = args.rng
        self._debug = args.d
        self._sequencing_read_len = args.R
        self._sequencing_coverage = args.c
        self._sequencing_fragment_size, self._sequencing_fragment_std = args.pe
        # TODO try use args.pe not null? :
        self._sequencing_paired_end = (self._sequencing_fragment_size and self._sequencing_fragment_std) \
                                      or args.pe_model
        self._sequencing_n_handling = ('random', self._sequencing_fragment_size) if self._sequencing_paired_end \
            else ('ignore', self._sequencing_read_len)

    def _params_sanity_check(self):
        # Check that files are real, if provided
        check_file_open(self._input_reference, f'ERROR: could not open reference, {self._input_reference}',
                        required=True)
        if (self._sequencing_fragment_size is None and self._sequencing_fragment_std is not None) \
                or (self._sequencing_fragment_size is not None and self._sequencing_fragment_std is None):
            print('\nERROR: --pe argument takes 2 space-separated arguments.\n')
            sys.exit(1)
        if self._rng_seed == -1:
            self._rng_seed = random.randint(1, 99999999)
        random.seed(self._rng_seed)
        is_in_range(self._sequencing_read_len, 10, 1000000, 'Error: -R must be between 10 and 1,000,000')
        is_in_range(self._sequencing_coverage, 0, 1000000, 'Error: -c must be between 0 and 1,000,000')
        # TODO check sanity workdir (output_params)
        # TODO check mut bed file, and vcf?

    def _load_mutation_model(self):
        self._mutation_model = parse_input_mutation_model(self._mutation_model, 1)
        if self._mutation_rate < 0.:
            self._mutation_rate = None
        if self._mutation_rate != -1 and self._mutation_rate is not None:
            is_in_range(self._mutation_rate, 0.0, 1.0, 'Error: -M must be between 0 and 0.3')

    def _index_reference(self):
        print('Indexing reference started:', self._input_reference)
        start = time.time()
        # index reference: [(0: chromosome name, 1: byte index where the contig seq begins,
        #                    2: byte index where the next contig begins, 3: contig seq length),
        #                    (repeat for every chrom)]
        # TODO check to see if this might work better as a dataframe or biopython object
        ref_index, line_width = index_ref(self._input_reference)
        # TODO check if this index can work, maybe it's faster
        # ref_index2 = SeqIO.index(reference, 'fasta')
        self._indices_by_ref_name = {chrom[0]: chrom for chrom in ref_index},  # TODO chrom[1:] ?
        self._output_line_width = line_width
        end = time.time()
        print('Indexing reference took {} seconds.'.format(int(end - start)))

    def _load_vcf(self):
        if self._input_variants is None:
            variants_from_vcf = pd.read_csv(self._input_variants_path)
            variants_from_vcf[['chrom', 'allele']] = variants_from_vcf[['chrom', 'allele']].astype(str)
            variants_from_vcf['pos'] = variants_from_vcf['pos'].astype(int)
            os.remove(self._input_variants_path)
            self._input_variants = variants_from_vcf

    def simulate(self):
        final_chromosomes = {}
        inserted_mutations = {}  # TODO change to list / df?

        for chrom in self._indices_by_ref_name.keys():
            chrom_processor, chrom_inserted_mutations = self._simulate_chrom(chrom)
            final_chromosomes[chrom_processor.chrom_name] = str(chrom_processor.chrom_sequence)
            inserted_mutations[chrom_processor.chrom_name] = chrom_inserted_mutations

        self._write_output(final_chromosomes, inserted_mutations)

    def _simulate_chrom(self, chrom) -> ChromosomeProcessor:
        # read in reference sequence and notate blocks of Ns
        chrom_sequence, n_regions = read_ref(self._input_reference, self._indices_by_ref_name[chrom],
                                             self._sequencing_n_handling)

        chrom_valid_variants = self._prune_invalid_variants(chrom, chrom_sequence)
        chrom_annotations_df = self._annotations_df[self._annotations_df['chrom'] == chrom]
        chromosome_processor = ChromosomeProcessor(chrom, chrom_sequence, chrom_annotations_df, annotations_sorted=True,
                                                   mut_models=self._mutation_model, mut_rate=self._mutation_rate,
                                                   dist=self._relative_distance, debug=self._debug)

        print('--------------------------------')
        print('Simulating chromosome {} of sequence started...'.format(chrom))
        t_start = time.time()
        for non_n_region in n_regions['non_N']:
            start, end = non_n_region
            inserted_mutations = self._simulate_window(chrom_valid_variants, chromosome_processor,
                                                       start, end, debug=self._debug)
        print(f'Simulating chromosome {chrom} took {int(time.time() - t_start)} seconds.')

        return chromosome_processor, inserted_mutations

    def _simulate_window(self, valid_variants_from_vcf: pd.DataFrame, chromosome_processor: ChromosomeProcessor,
                         start: int, end: int):
        if self._debug:
            print(f"Current window: start={start}, end={end}")

        vars_in_current_window = get_vars_in_window(start, end, valid_variants_from_vcf)
        given_mutations_inserted = chromosome_processor.insert_given_mutations(vars_in_current_window, start=start,
                                                                               end=end)

        random_mutations_inserted = chromosome_processor.generate_random_mutations(start, end)

        # TODO return random_mutations_inserted only?
        return given_mutations_inserted + random_mutations_inserted

    def _write_output(self, final_chromosomes: Dict[str, str], inserted_mutations):
        # FASTA
        fasta_file_writer = FastaFileWriter(self._output_prefix, self._output_line_width)
        for name, sequence in final_chromosomes.items():
            fasta_file_writer.write_record(sequence, name)
        fasta_file_writer.finalize()

        # FASTQ
        if not self._output_no_fastq:
            FastqFileWriter.generate_reads([fasta_file_writer.get_file_name()])

        # VCF
        if self._output_vcf:
            vcf_header = [self._input_reference]
            vcf_file_writer = VcfFileWriter(self._output_prefix, self._parent_prefix,
                                            self._output_accession, vcf_header)
            # TODO write inserted mutations - ?
            vcf_file_writer.close_file(add_parent_variants=True)

    def _prune_invalid_variants(self, chrom, chrom_sequence) -> pd.DataFrame:
        print('Pruning relevant variants from input VCF...')
        start = time.time()
        """Prune invalid input variants, e.g variants that:
                        - try to delete or alter any N characters
                        - don't match the reference base at their specified position
                        - any alt allele contains anything other than allowed characters"""
        valid_variants_from_vcf_indexes = []
        n_skipped = [0, 0, 0]
        if (not self._input_variants.empty) and (chrom in self._input_variants.chrom.unique()):
            for index, variant in self._input_variants[self._input_variants.chrom == chrom].iterrows():
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
        print('Done. Pruning took {} seconds.'.format(int(end - start)))
        return self._input_variants.loc[valid_variants_from_vcf_indexes, 'pos':].reset_index(drop=True)


def get_vars_in_window(start: int, end: int, valid_variants_from_vcf: pd.DataFrame) -> pd.DataFrame:
    vars_in_window_indexes = []

    for i, variant in valid_variants_from_vcf.iterrows():
        if start <= variant['pos'] < end:
            vars_in_window_indexes.append(i)

    if not valid_variants_from_vcf.empty:
        vars_in_window = valid_variants_from_vcf.loc[vars_in_window_indexes, 'pos':].reset_index(drop=True)
    else:
        vars_in_window = valid_variants_from_vcf
    return vars_in_window


# TODO deprecated?
def write_vcf(vcf_file_writer, inserted_mutations, chrom):
    print('Writing output VCF started...')
    start = time.time()
    for k in sorted(inserted_mutations.keys()):
        my_id = '.'
        my_quality = '.'
        my_filter = 'PASS'
        # k[0] + 1 because we're going back to 1-based vcf coords
        vcf_file_writer.write_record(chrom, str(int(k[0]) + 1), my_id, k[1], k[2], my_quality,
                                     my_filter, k[4])
    end = time.time()
    print('Done. Writing output VCF took {} seconds.'.format(int(end - start)))
