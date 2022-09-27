import argparse
import copy
import multiprocessing
from itertools import repeat

import ete3
import time
import os
import pandas as pd

from genome_simulator import GenomeSimulator
from writers.fasta_file_writer import FastaFileWriter
from writers.fastq_file_writer import FastqFileWriter
from utilities.input_checking import check_file_open
from writers.vcf_file_writer import VcfFileWriter
from utilities.vcf_func import parse_vcf


def parse_args(raw_args=None):
    parser = argparse.ArgumentParser(description='NEAT-genReads V3.0',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-r', type=str, required=True, metavar='reference', help="Path to reference fasta")
    parser.add_argument('-R', type=int, required=True, metavar='read length', help="The desired read length")
    parser.add_argument('-o', type=str, required=True, metavar='output_prefix',
                        help="Prefix for the output files (can be a path)")
    parser.add_argument('-c', type=float, required=False, metavar='coverage', default=10.0,
                        help="Average coverage, default is 10.0")
    parser.add_argument('-e', type=str, required=False, metavar='error_model', default=None,
                        help="Location of the file for the sequencing error model (omit to use the default)")
    parser.add_argument('-E', type=float, required=False, metavar='Error rate', default=-1,
                        help="Rescale avg sequencing error rate to this, must be between 0.0 and 0.3")
    parser.add_argument('-tr', type=str, required=False, metavar='target.bed', default=None,
                        help="Bed file containing targeted regions")
    parser.add_argument('-dr', type=str, required=False, metavar='discard_regions.bed', default=None,
                        help="Bed file with regions to discard")
    parser.add_argument('-to', type=float, required=False, metavar='off-target coverage scalar', default=0.00,
                        help="off-target coverage scalar")
    parser.add_argument('-m', type=str, required=False, metavar='model.p', default=None,
                        help="Mutation model pickle file")
    parser.add_argument('-M', type=float, required=False, metavar='avg mut rate', default=-1,
                        help="Rescale avg mutation rate to this (1/bp), must be between 0 and 0.3")
    parser.add_argument('-Mb', type=str, required=False, metavar='mut_rates.bed', default=None,
                        help="Bed file containing positional mut rates")
    parser.add_argument('-N', type=int, required=False, metavar='min qual score', default=-1,
                        help="below this quality score, replace base-calls with N's")
    parser.add_argument('-v', type=str, required=False, metavar='vcf.file', default=None,
                        help="Input VCF file of variants to include")
    parser.add_argument('-w', type=str, required=False, metavar='/path/to/working/dir/', default=None,
                        help="Name of working directory to process annotations. If not given, the directory of the "
                             "output file wll be taken instead")
    parser.add_argument('-ws', type=int, required=False, metavar='<int>', default=1000,
                        help='Window size of simulation. Choose an integer number larger than 10')
    parser.add_argument('--pe', nargs=2, type=int, required=False, metavar=('<int>', '<int>'), default=(None, None),
                        help='Paired-end fragment length mean and std')
    parser.add_argument('--pe-model', type=str, required=False, metavar='<str>', default=None,
                        help='empirical fragment length distribution')
    parser.add_argument('--gc-model', type=str, required=False, metavar='<str>', default=None,
                        help='empirical GC coverage bias distribution')
    parser.add_argument('--bam', required=False, action='store_true', default=False, help='output golden BAM file')
    parser.add_argument('--vcf', required=False, action='store_true', default=False, help='output golden VCF file')
    parser.add_argument('--rng', type=int, required=False, metavar='<int>', default=-1,
                        help='rng seed value; identical RNG value should produce identical runs of the program, so '
                             'things like read locations, variant positions, error positions, etc, '
                             'should all be the same.')
    parser.add_argument('--no-fastq', required=False, action='store_true', default=False,
                        help='bypass fastq generation')
    parser.add_argument('--discard-offtarget', required=False, action='store_true', default=False,
                        help='discard reads outside of targeted regions')
    parser.add_argument('--force-coverage', required=False, action='store_true', default=False,
                        help='[debug] ignore fancy models, force coverage to be constant')
    parser.add_argument('--rescale-qual', required=False, action='store_true', default=False,
                        help='Rescale quality scores to match -E input')
    # TODO implement a broader debugging scheme for subclasses.
    parser.add_argument('-d', required=False, action='store_true', default=False, help='Activate Debug Mode')
    parser.add_argument('-newick', type=str, required=False, metavar='newick tree', help="Path to reference newick")
    parser.add_argument('-a', type=str, required=False, metavar='leaf.name', default=None, help='reference accession')
    parser.add_argument('--max-threads', type=int, required=False, metavar='maximum threads number', default=1,
                        help='maximum threads number')

    return parser.parse_args(raw_args)


def main(raw_args=None):
    args = parse_args(raw_args)
    print("Using the next args:", args)
    load_default_args(args)

    # parse input variants, if present
    if args.v:
        check_file_open(args.v, 'ERROR: could not open input VCF, {}'.format(args.v), required=False)
        all_input_variants = load_input_variants(args.v)  # input_vcf = args.v
    else:
        all_input_variants = pd.DataFrame(columns=['dummy'])

    if not args.newick:
        print("No phylogenetic tree supplied")
        args.input_variants = all_input_variants
        print("Generating sequence started")
        start = time.time()
        genome_simulator = GenomeSimulator(args)
        genome_simulator.simulate()
        end = time.time()
        print("Done. Generating sequence took {} seconds.".format(int(end - start)))
        return

    t = ete3.Tree(args.newick, format=1)
    internal_count = 0
    for node in t.traverse("preorder"):
        if not node.name:
            internal_count += 1
            node.name = f"internal_{internal_count}"
    print("Using the next phylogenetic tree:\n", t.get_ascii(show_internal=True))
    clear_previous_tree_output(args.o, t)
    args.total_dist = tree_total_dist(t)

    # Relate to one of the accessions (given by -a param) as the reference
    args.root_to_ref_dist = set_ref_as_accession(args.a, t)

    args.root_fasta = args.r
    if not all_input_variants.empty:
        start = time.time()
        all_input_variants = all_input_variants.sample(frac=1).reset_index(drop=True)
        end = time.time()
        print("Suffling input variants was done in {} seconds.".format(int(end - start)))

    task_list = load_task_list(t, args, all_input_variants)
    del all_input_variants

    if args.max_threads > 1:
        generate_concurrently(args.max_threads, task_list)
    else:
        generate_sequentially(task_list)

    print('================================')
    # TODO remove all mut_bed files
    # TODO remove all csv/whatever files now, and not before!


def load_default_args(args):
    args.name = "simulation"
    args.dist = 1
    args.internal = False
    args.input_variants = None
    args.input_variants_path = None
    args.parent_name = None


def generate_sequentially(task_list):
    print("TEST this is a single-process simulation")
    for args in task_list:
        generate_for_node(args)


def generate_concurrently(ncpu, task_list):
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(processes=ncpu)
    cond = manager.Condition()
    pool.imap(process_handler, zip(task_list, repeat(cond)))
    pool.close()
    pool.join()


def process_handler(params_with_cond):
    simulation_params, cond = params_with_cond
    curr_proc = multiprocessing.current_process()
    print('TEST - current process:', curr_proc.name, curr_proc._identity)
    # When ancestor is ready - then start simulating
    ancestor_path = simulation_params.r
    while not os.path.exists(ancestor_path):
        print("TEST: ", curr_proc.name, 'checking if ancestor exists')
        with cond:
            cond.wait()
            print("TEST: ", curr_proc.name, 'checking again if ancestor exists')
    print("TEST: ", curr_proc, ", ancestor_path=", ancestor_path, "is ready")
    generate_for_node(simulation_params)
    with cond:
        print("TEST: ", curr_proc, " has finished simulation for ", simulation_params.name)
        cond.notify_all()

    return "Done"


def generate_for_node(args):
    print("TEST inside generate_for_node, args.r=", args.r)
    print("Generating sequence for taxon (node):", args.name)
    start = time.time()
    genome_simulator = GenomeSimulator(args)
    genome_simulator.simulate()
    end = time.time()
    print("Done. Generating sequence for taxon {} took {} seconds.".format(args.name, int(end - start)))


def get_node_args_for_simulation(node, args, all_input_variants, input_variants_used):
    new_args = copy.copy(args)
    new_args.name = node.name
    new_args.internal = len(node.children) > 0
    if not node.up.up:
        # Root direct descendants
        new_args.dist = (node.dist + new_args.root_to_ref_dist) / new_args.total_dist
        new_args.r = new_args.root_fasta
        new_args.parent_name = None
    else:
        new_args.dist = node.dist / new_args.total_dist
        new_args.r = FastaFileWriter.get_output_filenames(new_args.o, node.up.name)[0]
        new_args.parent_name = node.up.name
    branch_input_variants, input_variants_used = get_branch_input_variants(new_args.dist, all_input_variants,
                                                                           input_variants_used)
    if args.max_threads > 1:
        new_args.input_variants = None
        new_args.input_variants_path = args.o + "_" + node.name + "_input_variants.csv"
        branch_input_variants.to_csv(new_args.input_variants_path, index=False)
    else:
        new_args.input_variants = branch_input_variants
        new_args.input_variants_path = None
    return new_args, input_variants_used


def get_output_filenames(prefix, name):
    res = []
    res.extend(FastaFileWriter.get_output_filenames(prefix, name))
    res.extend(FastqFileWriter.get_output_filenames(prefix, name))
    res.extend(VcfFileWriter.get_output_filenames(prefix, name))
    return res


def clear_previous_tree_output(prefix, t):
    for node in t.traverse("preorder"):
        if node is t:
            continue
        else:
            output_filenames = get_output_filenames(prefix, node.name)
            for filename in output_filenames:
                if os.path.exists(filename):
                    os.remove(filename)


def load_input_variants(input_vcf):
    print("Loading input variants...")
    start = time.time()
    # TODO read this in as a pandas dataframe
    input_variants = None
    if input_vcf:
        (sample_names, input_variants) = parse_vcf(input_vcf, ploidy=1)
    end = time.time()
    print("Loading input variants took {} seconds.".format(int(end - start)))
    return input_variants


def get_branch_input_variants(branch_dist, input_variants, input_variants_used):
    if input_variants.empty:
        return input_variants, 0
    branch_input_variants_amount = round(branch_dist * len(input_variants))
    variants_used_including_this = min(input_variants_used + branch_input_variants_amount, len(input_variants))
    return input_variants.loc[input_variants_used:variants_used_including_this, :].reset_index(
        drop=True), variants_used_including_this


def set_ref_as_accession(accession, t):
    root_to_ref_dist = 0
    if accession:
        ref_node = t & accession
        t.set_outgroup(ref_node)
        root_to_ref_dist = ref_node.dist
        ref_node.delete()
    return root_to_ref_dist


def tree_total_dist(t):
    total_dist = 0
    for node in t.traverse("preorder"):
        total_dist += node.dist
    print("Total branches distance =", total_dist)
    return total_dist


def load_task_list(t, args, all_input_variants):
    task_list = []
    input_variants_used = 0
    for node in t.traverse("levelorder"):  # AKA breadth first search
        if node is t:
            continue  # This is the root - don't simulate for it
        else:
            node_params, input_variants_used = get_node_args_for_simulation(node, args, all_input_variants,
                                                                            input_variants_used)
            print("TEST, inside create_task_queue, node_params.name=", node_params.name)
            task_list.append(node_params)
    return task_list


if __name__ == "__main__":
    tt = time.time()
    print('Simulation begins')
    main()
    print("Simulation run time is : {0:.3f} (sec)".format(time.time() - tt))
