import argparse
import pandas as pd

import csv
import sys

REF_LABEL = 'TAIR10'
NO_GENE = '0'
GENE_PREF = 'transcript_'
NOVEL_GENE = 'PanGene'
ACCESSION_NAMES = ['A-lyrata', 'Cvi', 'C24', 'An-1', 'Col-0', 'Sha', 'Ler', 'Kyo', 'Eri']
SIMULATOR_PAV_COLUMNS = ['gene_id'] + ACCESSION_NAMES
PANORAMIC_PAV_COLUMNS = [f'{name}_{name}' for name in ACCESSION_NAMES] + [REF_LABEL, 'gene']
parser = argparse.ArgumentParser(description='Plot and compare gene stats]',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
parser.add_argument('-f', type=str, required=True, metavar='<str>', nargs='+',
                    help="* pan_PAV_1.tsv [pan_PAV_2.tsv] [pan_PAV_3.tsv] ...")
parser.add_argument('-s', type=str, required=True, metavar='<str>', nargs='+',
                    help="* simulator_PAV_1.csv [simulator_PAV_2.csv] [simulator_PAV_3.tsv] ...")
parser.add_argument('-l', type=str, required=True, metavar='<str>', nargs='+',
                    help="* simulator_legend_1.csv [simulator_legend_2.csv] [simulator_legend_3.tsv] ...")

args = parser.parse_args()
panoramic_pav_files = args.f
simulator_pav_dfs = []
simulator_raw_pav_files = args.s
simulator_legend_files = args.l
if len(panoramic_pav_files) != len(simulator_raw_pav_files) or len(simulator_raw_pav_files) != len(simulator_legend_files):
    print("Number of parameters incorrect")
    sys.exit(1)

for i in range(len(simulator_raw_pav_files)):
    raw_pav_df = pd.read_csv(simulator_raw_pav_files[i], index_col=0,
                             dtype={'gene_id': 'int', 'chrom': 'str', 'ref_start': 'int', 'ref_end': 'int'})
    legend_df = pd.read_csv(simulator_legend_files[i], index_col=0, dtype={'Name': 'str', 'ID': 'int'})
    legend_df = legend_df.rename(columns={'ID': 'gene_id', 'oldName2': 'newName2'})
    pav_df = pd.merge(raw_pav_df, legend_df, how='inner', on=['gene_id'])
    pav_df = pav_df[SIMULATOR_PAV_COLUMNS]
    simulator_pav_dfs.append(pav_df)

for i in range(len(panoramic_pav_files)):
    missing_genes = []
    with open(panoramic_pav_files[i]) as file:

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
            if NO_GENE in line:
                missing_genes.append(line[0])
                # TODO test per accession against simulator_pav_dfs[i]
    print(missing_genes)

    # TODO compute which genes are also absent in the simulator PAVs.
    #  Not sure if to do it accession-level or general-level.
