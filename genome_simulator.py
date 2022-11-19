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
import pandas as pd
from typing import List, Dict, Tuple

from writers.fastq_file_writer import FastqFileWriter
from chromosome_processor import ChromosomeProcessor
from mutation_model import load_mutation_model_from_file
from writers.fasta_file_writer import FastaFileWriter
from utilities.input_checking import check_file_open, is_in_range
from utilities.ref_func import index_ref, read_ref
from writers.vcf_file_writer import VcfFileWriter
from utilities.genome_annotations import read_annotations_csv

ANNOTATIONS_FILE_FORMAT = '{}_annotations.csv'

"""
Some constants needed for analysis
"""

# Default window size for simulation
DEFAULT_WINDOW_SIZE = 10000
MIN_WINDOW_SIZE = 100


class GenomeSimulator:
    def __init__(self, args):
        self._extract_args(args)
        self._params_sanity_check()

        # Use pathlib to make this more machine agnostic
        self._output_prefix_name = pathlib.Path(self._output_prefix).name

        # Load annotations dataframe
        self._annotations_df = read_annotations_csv(self._annotations_file)
        self._load_mutation_model_data(args)

        self._index_reference()

    def _extract_args(self, args):
        self._input_reference = args.r
        self._annotations_file = args.Mb
        self._relative_distance = args.dist
        self._output_accession = args.name
        self._output_prefix = args.o + '_' + args.name
        self._parent_prefix = args.o + '_' + args.parent_name if args.parent_name else None
        self._output_vcf = args.vcf
        self._output_no_fastq = args.no_fastq or args.internal
        self._rng_seed = args.rng
        self._debug = args.d
        self._sequencing_read_len = args.R
        self._sequencing_coverage = args.c
        self._sequencing_fragment_size, self._sequencing_fragment_std = args.pe
        # TODO try use args.pe not null? :
        self._sequencing_paired_end = (self._sequencing_fragment_size and self._sequencing_fragment_std) \
            or args.pe_model
        self._sequencing_n_handling = {'method': 'random', 'max_threshold': self._sequencing_fragment_size} \
            if self._sequencing_paired_end else {'method': 'ignore', 'max_threshold': None}
        self._window_size = args.ws if args.ws > MIN_WINDOW_SIZE else DEFAULT_WINDOW_SIZE

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
        check_file_open(self._annotations_file, f'ERROR: could not open annotations file, {self._annotations_file}',
                        required=False)

    def _load_mutation_model_data(self, args):
        self._mutation_model = str(args.m)
        self._mutation_model = load_mutation_model_from_file(self._mutation_model)
        self._mutation_scalar = args.M
        if self._mutation_scalar <= 0.:
            self._mutation_scalar = None

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
        self._indices_by_ref_name = {chrom[0]: chrom for chrom in ref_index}
        self._output_line_width = line_width
        end = time.time()
        print('Indexing reference took {} seconds.'.format(int(end - start)))

    def simulate(self) -> set:
        inserted_mutations = {}
        snps_indels_counts = {}
        fasta_file_writer = FastaFileWriter(self._output_prefix, self._output_line_width)

        dfs_to_concat = []
        accession_genes = set()
        for chrom in self._indices_by_ref_name.keys():
            chrom_processor, chrom_inserted_mutations, chrom_snps_count, chrom_indels_count = \
                self._simulate_chrom(chrom)
            inserted_mutations[chrom_processor.chrom_name] = chrom_inserted_mutations
            snps_indels_counts[chrom_processor.chrom_name] = (chrom_snps_count, chrom_indels_count)
            dfs_to_concat.append(chrom_processor.get_annotations_df())
            fasta_file_writer.write_record(str(chrom_processor.chrom_sequence), chrom_processor.chrom_name)
            accession_genes = accession_genes.union(chrom_processor.get_genes())

        fasta_file_writer.finalize()
        new_annotations_df = pd.concat(dfs_to_concat).reset_index(drop=True)
        self._write_non_fasta_output(fasta_file_writer.get_file_name(), inserted_mutations, snps_indels_counts,
                                     new_annotations_df)
        return accession_genes

    def _simulate_chrom(self, chrom) -> (ChromosomeProcessor, List[Tuple], int, int):
        # read in reference sequence and notate blocks of Ns
        chrom_sequence, n_regions = read_ref(self._input_reference, self._indices_by_ref_name[chrom],
                                             self._sequencing_n_handling)

        chromosome_processor = ChromosomeProcessor(chrom, chrom_sequence, self._annotations_df, annotations_sorted=True,
                                                   mut_model=self._mutation_model, mut_scalar=self._mutation_scalar,
                                                   dist=self._relative_distance, debug=self._debug)

        print('--------------------------------')
        print('Simulating chromosome {} of sequence started...'.format(chrom))
        t_start = time.time()
        inserted_mutations = []
        total_snps_count = 0
        total_indels_count = 0
        for non_n_region in n_regions['non_N']:
            windows_list = self._break_to_windows(non_n_region[0], non_n_region[1])
            for window in windows_list:
                start, end = window
                window_mutations, snps_count, indels_count = self._simulate_window(chromosome_processor, start, end)
                inserted_mutations.extend(window_mutations)
                total_snps_count += snps_count
                total_indels_count += indels_count
        print(f'Simulating chromosome {chrom} took {int(time.time() - t_start)} seconds.')

        return chromosome_processor, inserted_mutations, total_snps_count, total_indels_count

    def _break_to_windows(self, start: int, end: int) -> List[Tuple[int, int]]:
        remained_length = end - start
        windows_list = [(start, end)]
        while remained_length > 1.01 * self._window_size and remained_length > self._window_size + MIN_WINDOW_SIZE:
            prev_start, prev_end = windows_list.pop()
            next_end = prev_end
            next_start = prev_start + self._window_size
            prev_end = next_start
            windows_list.append((prev_start, prev_end))
            windows_list.append((next_start, next_end))
            remained_length = next_end - next_start
        return windows_list

    def _simulate_window(self, chromosome_processor: ChromosomeProcessor, start: int, end: int) ->\
            (List[Tuple], int, int):
        if self._debug:
            print(f"New simulation window: start={start}, end={end}")
        chromosome_processor.next_window(start=start, end=end)
        return chromosome_processor.generate_random_mutations()

    def _write_non_fasta_output(self, fasta_file_name: str, inserted_mutations: Dict[str, List[Tuple]],
                                snps_indels_counts: Dict[str, Tuple], new_annotations_df: pd.DataFrame = None):

        # FASTQ
        if not self._output_no_fastq:
            FastqFileWriter.generate_reads([fasta_file_name], self._sequencing_paired_end,
                                           self._sequencing_read_len, self._sequencing_coverage,
                                           self._sequencing_fragment_size, self._sequencing_fragment_std)

        # VCF
        if self._output_vcf:
            self._write_vcf(inserted_mutations)

        # annotations CSV
        new_annotations_df.to_csv(ANNOTATIONS_FILE_FORMAT.format(self._output_prefix))

        # SNPs and Indels counts
        overall_snps = 0
        overall_indels = 0
        for chrom_name, counts in snps_indels_counts.items():
            print(f'Chromosome {chrom_name}: {counts[0]} SNPs and {counts[1]} Indels were inserted')
            overall_snps += counts[0]
            overall_indels += counts[1]
        print(f'Overall {overall_snps} SNPs and {overall_indels} were inserted for this genome.')

    def _write_vcf(self, inserted_mutations: Dict[str, List[Tuple]]):
        vcf_header = [self._input_reference]
        vcf_file_writer = VcfFileWriter(self._output_prefix, self._parent_prefix,
                                        self._output_accession, vcf_header)
        print('Writing output VCF started...')
        start = time.time()
        for chrom in inserted_mutations.keys():
            for mutation in inserted_mutations[chrom]:
                my_id = '.'
                my_quality = '.'
                my_filter = 'PASS'
                # mutation[0] + 1 because we're going back to 1-based vcf coords
                vcf_file_writer.write_record(chrom, str(int(mutation[0]) + 1), my_id, mutation[1], mutation[2],
                                             my_quality, my_filter, 'WP=1')

        end = time.time()
        print('Done. Writing output VCF took {} seconds.'.format(int(end - start)))
        vcf_file_writer.close_file(add_parent_variants=True)
