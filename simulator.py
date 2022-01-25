import argparse
from neat import gen_reads
import ete3
import time
import vcf
from vcf import utils
import os
import Bio.bgzf as bgzf

def parse_args(raw_args=None):
    parser = argparse.ArgumentParser(description='NEAT-genReads V3.0',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
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
    # TODO add to the documentation that now ploidy default is 1
    # TODO Consider adding feature of working with a multiple ploidy *reference* (doesn't exist now, only multi-ploidy output exists)
    parser.add_argument('-p', type=int, required=False, metavar='ploidy', default=1,
                        help="Desired ploidy, default = 1")
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
    parser.add_argument('--save-fasta', required=False, action='store_true', default=False,
                        help='outputs FASTA')
    parser.add_argument('-a', type=str, required=False, metavar='leaf.name', default=None, help='reference accession')

    return parser.parse_args(raw_args)

def main(raw_args=None):
    args = parse_args(raw_args)

    if not args.newick:
        print("No phylogenetic tree supplied")
        args.name = "simulation"
        args.dist = 1
        print("Using the next args:",args)
        gen_reads.simulate(args)
        print('================================')
        print('Done.')
        return

    t = ete3.Tree(args.newick, format=1)
    print("Using the next phylogenetic tree:\n",t.get_ascii(show_internal=True))
    total_dist = tree_total_dist(t)

    # Relate to one of the accessions (given by -a param) as the reference
    root_to_ref_dist = set_ref_as_accession(args.a, t)

    ancestor_fasta = args.r
    for node in t.traverse("preorder"):
        if node is t:
            continue # This is the root
        print('================================')
        print("Generating sequence for taxon (node):",node.name)
        if node.up is t:
            print("The parent is the root")
            args.r = ancestor_fasta
            args.dist = (node.dist + root_to_ref_dist)/ total_dist
        else:
            print("The parent is:", node.up.name)
            args.r = args.o + "_" + node.up.name + ".fasta"
            args.dist = node.dist / total_dist
        args.name = node.name
        print("Using the next args:",args)
        gen_reads.simulate(args)
    print('================================')
    print('Merging VCF files across the generations...')
    for node in t.traverse("preorder"):
        if node is t:
            continue # This is the root or its direct descendants
        if node.up is t:
            newname = args.o + "_" + node.name + "_golden_final.vcf.gz"
            oldname = args.o + "_" + node.name + "_golden.vcf.gz"
            if os.path.exists(newname):
                os.remove(newname)
            os.rename(oldname, newname)
        else:
            parent_file = args.o + "_" + node.up.name + "_golden_final.vcf.gz"
            child_file = args.o + "_" + node.name + "_golden.vcf.gz"
            out_file = args.o + "_" + node.name + "_golden_final.vcf.gz"
            add_parent_variants(parent_file, child_file, out_file)
    print('Done.')


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

def add_parent_variants(parent_file, child_file, out_file):
    vcf_reader_parent = vcf.Reader(filename=parent_file, strict_whitespace=True)
    vcf_reader_child = vcf.Reader(filename=child_file, strict_whitespace=True)
    out = bgzf.open(out_file, 'wb')
    vcf_writer = vcf.Writer(out, vcf_reader_child)

    iterate_simulatnously = utils.walk_together(vcf_reader_parent,vcf_reader_child)
    for readers in iterate_simulatnously:
        if readers[1]:
            vcf_writer.write_record(readers[1])
        elif (readers[0]):
            vcf_writer.write_record(readers[0])

    out.close()


if __name__ == "__main__":
    tt = time.time()
    print('Simulation begins')
    main()
    print("Simulation run time is : {0:.3f} (sec)".format(time.time() - tt))
