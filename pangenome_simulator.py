import argparse
import copy
import multiprocessing
from itertools import repeat

import ete3
import time
import os

from genome_simulator import GenomeSimulator, ANNOTATIONS_FILE_FORMAT
from writers.fasta_file_writer import FastaFileWriter
from writers.fastq_file_writer import FastqFileWriter
from writers.vcf_file_writer import VcfFileWriter


def main(raw_args=None):
    args = parse_args(raw_args)
    print("Got the next args from the user:", args)
    load_default_args(args)

    # Single simulation
    if not args.newick:
        print("No phylogenetic tree supplied")
        print("Generating sequence started")
        start = time.time()
        genome_simulator = GenomeSimulator(args)
        genome_simulator.simulate()
        end = time.time()
        print("Done. Generating sequence took {} seconds.".format(int(end - start)))
        return

    # Tree simulation
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

    task_list = load_task_list(t, args)
    if args.max_threads > 1:
        generate_concurrently(args.max_threads, task_list)
    else:
        generate_sequentially(task_list)

    print('================================')
    # TODO remove all mut_bed files
    # TODO remove all csv/whatever files now, and not before!


def parse_args(raw_args=None):
    parser = argparse.ArgumentParser(description='PanGenomeSimulator',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-r', type=str, required=True, metavar='reference', help="Path to reference fasta")
    parser.add_argument('-R', type=int, required=True, metavar='read length', help="The desired read length")
    parser.add_argument('-o', type=str, required=True, metavar='output_prefix',
                        help="Prefix for the output files (can be a path)")
    parser.add_argument('-c', type=float, required=False, metavar='coverage', default=10.0,
                        help="Average coverage, default is 10.0")
    parser.add_argument('-m', type=str, required=False, metavar='model.p', default=None,
                        help="Mutation model pickle file")
    parser.add_argument('-M', type=float, required=False, metavar='avg mut rate', default=-1,
                        help="Rescale avg mutation rate to this (1/bp), must be between 0 and 0.3")
    parser.add_argument('-Mb', type=str, required=False, metavar='mut_rates.bed', default=None,
                        help="Bed file containing positional mut rates")
    parser.add_argument('-ws', type=int, required=False, metavar='<int>', default=10000,
                        help='Window size of simulation. Choose an integer number larger than 10')
    parser.add_argument('--pe', nargs=2, type=int, required=False, metavar=('<int>', '<int>'), default=(None, None),
                        help='Paired-end fragment length mean and std')
    parser.add_argument('--vcf', required=False, action='store_true', default=False, help='output golden VCF file')
    parser.add_argument('--rng', type=int, required=False, metavar='<int>', default=-1,
                        help='rng seed value; identical RNG value should produce identical runs of the program, so '
                             'things like read locations, variant positions, error positions, etc, '
                             'should all be the same.')
    parser.add_argument('--no-fastq', required=False, action='store_true', default=False,
                        help='bypass fastq generation')
    parser.add_argument('-d', required=False, action='store_true', default=False, help='Activate Debug Mode')
    parser.add_argument('-newick', type=str, required=False, metavar='newick tree', help="Path to reference newick")
    parser.add_argument('-a', type=str, required=False, metavar='leaf.name', default=None, help='reference accession')
    parser.add_argument('--max-threads', type=int, required=False, metavar='maximum threads number', default=1,
                        help='maximum threads number')

    return parser.parse_args(raw_args)


def load_default_args(args):
    args.name = "simulation"
    args.dist = 1
    args.internal = False
    args.parent_name = None


def clear_previous_tree_output(prefix, t):
    for node in t.traverse("preorder"):
        if node is t:
            continue
        else:
            output_filenames = get_output_filenames(prefix, node.name)
            for filename in output_filenames:
                if os.path.exists(filename):
                    os.remove(filename)


def get_output_filenames(prefix, name):
    res = []
    res.extend(FastaFileWriter.get_output_filenames(prefix, name))
    res.extend(FastqFileWriter.get_output_filenames(prefix, name))
    res.extend(VcfFileWriter.get_output_filenames(prefix, name))
    res.append(ANNOTATIONS_FILE_FORMAT.format(prefix + '_' + name))
    return res


def tree_total_dist(t):
    total_dist = 0
    for node in t.traverse("preorder"):
        total_dist += node.dist
    print("Total branches distance =", total_dist)
    return total_dist


def set_ref_as_accession(accession, t):
    root_to_ref_dist = 0
    if accession:
        ref_node = t & accession
        t.set_outgroup(ref_node)
        root_to_ref_dist = ref_node.dist
        ref_node.delete()
    return root_to_ref_dist


def load_task_list(t, args):
    task_list = []
    for node in t.traverse("levelorder"):  # AKA breadth first search
        if node is t:
            continue  # This is the root - don't simulate for it
        else:
            node_params = get_node_args_for_simulation(node, args)
            task_list.append(node_params)
    return task_list


def get_node_args_for_simulation(node, args):
    new_args = copy.copy(args)
    new_args.name = node.name
    new_args.internal = len(node.children) > 0
    if not node.up.up:
        # Root direct descendants
        new_args.parent_name = None
        new_args.dist = (node.dist + new_args.root_to_ref_dist) / new_args.total_dist
        new_args.r = new_args.root_fasta
    else:
        new_args.parent_name = node.up.name
        new_args.dist = node.dist / new_args.total_dist
        new_args.r = FastaFileWriter.get_output_filenames(new_args.o, new_args.parent_name)[0]
        new_args.Mb = ANNOTATIONS_FILE_FORMAT.format(new_args.o + '_' + new_args.parent_name)
    return new_args


def generate_concurrently(ncpu, task_list):
    print(f"Starting a multi-process simulation, number of processes: {ncpu}")
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(processes=ncpu)
    cond = manager.Condition()
    pool.imap(process_handler, zip(task_list, repeat(cond)))
    pool.close()
    pool.join()


def process_handler(params_with_cond):
    simulation_params, cond = params_with_cond
    curr_proc = multiprocessing.current_process()
    print('MULTIPROCESS LOG - current process:', curr_proc.name, curr_proc._identity)
    # When ancestor is ready - then start simulating
    ancestor_path = simulation_params.r
    while not os.path.exists(ancestor_path):
        print("MULTIPROCESS LOG: ", curr_proc.name, 'checking if ancestor exists')
        with cond:
            cond.wait()
            print("MULTIPROCESS LOG: ", curr_proc.name, 'checking again if ancestor exists')
    print("MULTIPROCESS LOG: ", curr_proc, ", ancestor_path=", ancestor_path, "is ready")
    generate_for_node(simulation_params)
    with cond:
        print("MULTIPROCESS LOG: ", curr_proc, " has finished simulation for ", simulation_params.name)
        cond.notify_all()

    return "Done"


def generate_sequentially(task_list):
    print("Starting a single-process simulation")
    for args in task_list:
        generate_for_node(args)


def generate_for_node(args):
    print("Generating sequence for taxon (node):", args.name)
    start = time.time()
    genome_simulator = GenomeSimulator(args)
    genome_simulator.simulate()
    end = time.time()
    print("Done. Generating sequence for taxon {} took {} seconds.".format(args.name, int(end - start)))


if __name__ == "__main__":
    tt = time.time()
    print('Simulation begins')
    main()
    print("Simulation run time is : {0:.3f} (sec)".format(time.time() - tt))
