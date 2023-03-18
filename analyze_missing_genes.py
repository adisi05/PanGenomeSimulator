import re
import argparse
import matplotlib.pyplot as plt
import csv

GENE_PREF = r'>?transcript[_:]'
NOVEL_GENE = 'PanGene'
NO_GENE = '0'
UNMAPPED = {
    'AT1G33355',
    'AT1G64633',
    'AT2G07642',
    'AT2G21105',
    'AT4G09355',
    'AT4G12485',
    'AT4G33735',
    'AT5G07545',
    'AT5G23115',
    'AT5G24575',
    'AT5G40315',
    'ATMG00665'
}
original_genes = {}


def plot_histogram(d):
    plt.bar(d.keys(), d.values(), width=0.4)
    plt.xlabel('Key')
    plt.ylabel('Value')
    plt.title('Histogram of a Dictionary')
    plt.show()


parser = argparse.ArgumentParser(description='Plot and compare missing genes',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
parser.add_argument('--original', type=str, required=True, metavar='<str>', nargs='+',
                    help="* proteins1.fasta [proteins2.fasta] [proteins3.fasta] ...")
parser.add_argument('--pav', type=str, required=True, metavar='<str>', nargs='+',
                    help="* pan_PAV_1.tsv [pan_PAV_2.tsv] [pan_PAV_3.tsv] ...")
args = parser.parse_args()
original_gene_files = args.original
panoramic_pav_files = args.pav


curr_gene_name = None
curr_gene_len = 0
pattern = re.compile(GENE_PREF)
for gene_file in original_gene_files:
    with open(gene_file) as file:
        for line in file:
            if pattern.match(line):
                if curr_gene_name:
                    original_genes[curr_gene_name] = curr_gene_len

                curr_gene_name = re.sub(GENE_PREF, '', line)
                curr_gene_name = curr_gene_name.split('.')[0]
                curr_gene_len = 0
            else:
                curr_gene_len += len(line.strip())


for pav_file in panoramic_pav_files:
    gene_stats = {}

    with open(pav_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:

            # Novel genes
            if line[0].startswith(NOVEL_GENE):
                continue

            # Known genes
            gene_name = re.sub(GENE_PREF, '', line[0])
            gene_name = gene_name.split('.')[0]
            if gene_name in UNMAPPED:
                print("Unmapped", gene_name)
            if NO_GENE in line:
                if gene_name not in original_genes:
                    print("Missing", gene_name)
                else:
                    gene_length = original_genes[gene_name]
                    if gene_length < 50:
                        print("Short:", gene_name, "length:", gene_length)
                    if gene_length not in gene_stats:
                        gene_stats[gene_length] = 0
                    gene_stats[gene_length] += 1

    plot_histogram(gene_stats)
    print(pav_file)
    print("Min value:", min(gene_stats))
    print("Max:", max(gene_stats))