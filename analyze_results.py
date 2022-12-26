import argparse

import csv
import sys

REF_LABEL = 'TAIR10'
NO_GENE = '0'
GENE_PREF = 'transcript_'
NOVEL_GENE = 'PanGene'
parser = argparse.ArgumentParser(description='Plot and compare gene stats',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
parser.add_argument('-f', type=str, required=True, metavar='<str>', nargs='+',
                    help="* pan_PAV_1.tsv [pan_PAV_2.tsv] [pan_PAV_3.tsv] ...")
parser.add_argument('-s', type=str, required=True, metavar='<str>', nargs='+',
                    help="* simulator_PAV_1.csv [simulator_PAV_2.csv] [simulator_PAV_3.tsv] ...")
parser.add_argument('-l', type=str, required=True, metavar='<str>', nargs='+',
                    help="* simulator_legend_1.csv [simulator_legend_2.csv] [simulator_legend_3.tsv] ...")

args = parser.parse_args()
panoramic_pav_files = args.f
simulator_pav_files = args.s
simulator_legend_files = args.l
if len(panoramic_pav_files) != len(simulator_pav_files) or len(simulator_pav_files) != len(simulator_legend_files):
    print("Number of parameters incorrect")
    sys.exit(1)

for i in range(len(simulator_pav_files)):
    # TODO read simulator_legend_files[i] and simulator_pav_files[i] and see which genes are really presence in each
    #  accession and which aren't
    pass

for panoramic_pav in panoramic_pav_files:
    missing_genes = []
    with open(panoramic_pav) as file:

        tsv_file = csv.reader(file, delimiter="\t")

        ref_index = None
        # i = 0
        for line in tsv_file:

            if not ref_index:
                ref_index = line.index(REF_LABEL)
                continue

            if line[0].startswith(NOVEL_GENE):
                continue

            line.pop(ref_index)
            # print(line)  # TEST
            if NO_GENE in line:
                missing_genes.append(line[0])
            # i += 1  # TEST
            # print(line)  # TEST
            # if i == 5:  # TEST
            #     break
    print(missing_genes)

    # TODO compute which genes are also absent in the simulator PAVs.
    #  Not sure if to do it accession-level or general-level.
