import argparse

import os
import pandas as pd

import csv
import sys

REF_LABEL = 'TAIR10'
NO_GENE = '0'
GENE_PREF = 'transcript_'
NOVEL_GENE = 'PanGene'
ACCESSION_NAMES = ['A-lyrata', 'Cvi', 'C24', 'An-1', 'Col-0', 'Sha', 'Ler', 'Kyo', 'Eri']
SIMULATOR_PAV_COLUMNS = ['Name', 'Transcript'] + ACCESSION_NAMES
PANORAMIC_PAV_COLUMNS = [f'{name}_{name}' for name in ACCESSION_NAMES] + [REF_LABEL, 'gene']
parser = argparse.ArgumentParser(description='Plot and compare gene stats]',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
parser.add_argument('-f', type=str, required=True, metavar='<str>', nargs='+',
                    help="* pan_PAV_1.tsv [pan_PAV_2.tsv] [pan_PAV_3.tsv] ...")
parser.add_argument('-s', type=str, required=True, metavar='<str>', nargs='+',
                    help="* simulator_PAV_1.csv [simulator_PAV_2.csv] [simulator_PAV_3.tsv] ...")
parser.add_argument('-l', type=str, required=True, metavar='<str>', nargs='+',
                    help="* simulator_legend_1.csv [simulator_legend_2.csv] [simulator_legend_3.tsv] ...")
parser.add_argument('-t', type=str, required=True, metavar='<str>', nargs='+',
                    help="* 'De-novo 01' ['Map-to-pan 01'] ['De-novo 02'] ...")


args = parser.parse_args()
panoramic_pav_files = args.f
simulator_raw_pav_files = args.s
simulator_legend_files = args.l
file_tags = args.t

simulator_pav_dfs = []
stats_per_file = {}
if len(panoramic_pav_files) != len(simulator_raw_pav_files) or len(simulator_raw_pav_files) != len(simulator_legend_files):
    print("Number of parameters incorrect")
    sys.exit(1)

for i in range(len(simulator_raw_pav_files)):
    raw_pav_df = pd.read_csv(simulator_raw_pav_files[i], index_col=0,
                             dtype={'gene_id': 'int', 'chrom': 'str', 'ref_start': 'int', 'ref_end': 'int'})
    legend_df = pd.read_csv(simulator_legend_files[i], index_col=0, dtype={'Name': 'str', 'ID': 'int'})
    legend_df = legend_df.rename(columns={'ID': 'gene_id'})
    pav_df = pd.merge(raw_pav_df, legend_df, how='inner', on=['gene_id'])
    pav_df[['Name', 'Transcript']] = pav_df.apply(lambda x: x['Name'].split('.'), axis=1, result_type='expand')
    pav_df = pav_df.loc[:, pav_df.columns.isin(SIMULATOR_PAV_COLUMNS)]  # select this columns if exist

    simulator_pav_dfs.append(pav_df)

for i in range(len(panoramic_pav_files)):
    with open(panoramic_pav_files[i]) as file:
        file_stats = {}
        missing_genes = {}
        relevant_accessions = []
        ref_index = None
        accession_indices = {}
        simulator_pav_df = simulator_pav_dfs[i]

        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:

            if not ref_index:
                ref_index = line.index(REF_LABEL)
                for a_name in ACCESSION_NAMES:
                    if f'{a_name}_{a_name}' in line:
                        accession_indices[a_name] = line.index(f'{a_name}_{a_name}')
                relevant_accessions = list(accession_indices.keys())
                continue

            if line[0].startswith(NOVEL_GENE):
                continue

            line.pop(ref_index)
            if NO_GENE in line:
                gene_name = line[0].strip(GENE_PREF).split('.')[0]
                simulator_result = simulator_pav_df[simulator_pav_df['Name'] == gene_name]
                accession_dict = {}
                if not simulator_result.empty:
                    accession_dict['All'] = False in simulator_result[relevant_accessions].values[0]
                    for a_name in relevant_accessions:
                        if NO_GENE == line[accession_indices[a_name]]:
                            accession_dict[a_name] = not simulator_result[a_name].item()
                else:
                    accession_dict['All'] = True
                    for a_name in relevant_accessions:
                        if NO_GENE == line[accession_indices[a_name]]:
                            accession_dict[a_name] = True
                missing_genes[gene_name] = accession_dict

    print("Total missing genes by Panoramic:", len(missing_genes))
    simulator_agree_count = 0
    for _, v in missing_genes.items():
        if v['All']:
            simulator_agree_count += 1
    print("The numbers of missing genes that are defined missing also by the simulator:", simulator_agree_count)
    file_stats['All'] = {'Overall': len(missing_genes), 'Agreed': simulator_agree_count}
    for a_name in relevant_accessions:
        missing_genes_in_accession = 0
        simulator_agree_count_accession = 0
        for _, v in missing_genes.items():
            if a_name in v:
                missing_genes_in_accession += 1
                if v[a_name]:
                    simulator_agree_count_accession += 1
        print("For accession:", a_name)
        print("Total missing genes by Panoramic:", missing_genes_in_accession)
        print("The numbers of missing genes that are defined missing also by the simulator:",
              simulator_agree_count_accession)
        file_stats[a_name] = {'Overall': missing_genes_in_accession, 'Agreed': simulator_agree_count_accession}

    file_name = file_tags[i]
    stats_per_file[file_name] = file_stats

import matplotlib.pyplot as plt

bar_width = 0.35
for file_name, file_stats in stats_per_file.items():
    # Set the number of groups and the values for each group
    N = len(file_stats)
    values1 = [v['Overall'] for _, v in file_stats.items()]
    values2 = [v['Agreed'] for _, v in file_stats.items()]

    # Set the labels for each group and the bar width
    labels = [k for k, _ in file_stats.items()]

    # Set the position of the bars on the x-axis
    x_pos = [i for i in range(N)]

    # Build the first set of bars
    plt.bar(x_pos, values1, bar_width, alpha=0.5, color='b')

    # Build the second set of bars, positioned slightly to the right of the first set
    plt.bar([i+bar_width for i in x_pos], values2, bar_width, alpha=0.5, color='r')

    # Set the x-axis tick marks and labels
    plt.xticks([i+bar_width/2 for i in x_pos], labels)

    # Add a legend and title
    plt.legend(['Overall', 'Agreed'])
    plt.title(f'{file_name}: Missing genes, overall and agreed')

    # Show the plot
    plt.show()




# Set the number of groups and the values for each group
N = len(stats_per_file)
values1 = [v['All']['Overall'] for _, v in stats_per_file.items()]
values2 = [v['All']['Agreed'] for _, v in stats_per_file.items()]

# Set the labels for each group and the bar width
labels = [k for k, _ in stats_per_file.items()]

# Set the position of the bars on the x-axis
x_pos = [i for i in range(N)]

# Build the first set of bars
plt.bar(x_pos, values1, bar_width, alpha=0.5, color='b')

# Build the second set of bars, positioned slightly to the right of the first set
plt.bar([i + bar_width for i in x_pos], values2, bar_width, alpha=0.5, color='r')

# Set the x-axis tick marks and labels
plt.xticks([i + bar_width / 2 for i in x_pos], labels)

# Add a legend and title
plt.legend(['Overall', 'Agreed'])
plt.title('Missing genes, overall and agreed')

# Show the plot
plt.show()
