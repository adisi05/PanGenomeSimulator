import argparse
import pandas as pd
import csv
import sys
import matplotlib.pyplot as plt
import re

REF_LABEL = 'TAIR10'
NO_GENE = '0'
GENE_PREF = r'transcript[_:]'
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
parser.add_argument('-T', type=str, required=True, metavar='<str>', help="Graph Title")



args = parser.parse_args()
panoramic_pav_files = args.f
simulator_raw_pav_files = args.s
simulator_legend_files = args.l
file_tags = args.t
graph_title = args.T

stats_per_file = {tag: {} for tag in file_tags}
new_stats_per_file = {tag: {} for tag in file_tags}

if len(panoramic_pav_files) != len(simulator_raw_pav_files) or len(simulator_raw_pav_files) != len(simulator_legend_files):
    print("Number of parameters incorrect")
    sys.exit(1)

for i in range(len(file_tags)):

    # Simulator "true" values
    sim_pav_df = pd.read_csv(simulator_raw_pav_files[i], index_col=0,
                             dtype={'gene_id': 'int', 'chrom': 'str', 'ref_start': 'int', 'ref_end': 'int'})
    legend_df = pd.read_csv(simulator_legend_files[i], index_col=0, dtype={'Name': 'str', 'ID': 'int'})
    legend_df = legend_df.rename(columns={'ID': 'gene_id'})
    sim_pav_df = pd.merge(sim_pav_df, legend_df, how='inner', on=['gene_id'])
    sim_pav_df[['Name', 'Transcript']] = sim_pav_df.apply(lambda x: x['Name'].split('.'), axis=1, result_type='expand')
    sim_pav_df = sim_pav_df.loc[:, sim_pav_df.columns.isin(SIMULATOR_PAV_COLUMNS)]  # select this columns if exist

    relevant_accessions = [column for column in sim_pav_df if column in ACCESSION_NAMES]

    with open(panoramic_pav_files[i]) as file:
        file_stats = stats_per_file[file_tags[i]]
        pano_pav_index = {}
        missing_genes = {}

        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:

            # Indices
            if not bool(pano_pav_index):
                pano_pav_index[REF_LABEL] = line.index(REF_LABEL)
                for a_name in relevant_accessions:
                    if f'{a_name}_{a_name}' not in line:
                        raise Exception(f"Accession named {a_name} was not found in {panoramic_pav_files[i]}")
                    pano_pav_index[a_name] = line.index(f'{a_name}_{a_name}')
                    file_stats[a_name] = {
                        'Real Absent': 0,
                        'Real Present': 0,
                        'Predicted Absent': 0,
                        'Predicted Present': 0,
                        'Correct Absent': 0,
                        'Correct Present': 0
                    }
                continue

            # Novel genes
            if line[0].startswith(NOVEL_GENE):
                continue

            # Known genes
            line.pop(pano_pav_index[REF_LABEL])
            gene_name = re.sub(GENE_PREF, '', line[0])
            gene_name = gene_name.split('.')[0]
            simulator_result = sim_pav_df[sim_pav_df['Name'] == gene_name]
            if not simulator_result.empty:
                for a_name in relevant_accessions:
                    if NO_GENE == line[pano_pav_index[a_name]]:
                        file_stats[a_name]['Predicted Absent'] += 1
                        if not simulator_result[a_name].item():
                            file_stats[a_name]['Real Absent'] += 1
                            file_stats[a_name]['Correct Absent'] += 1
                        else:
                            file_stats[a_name]['Real Present'] += 1

                    else:
                        file_stats[a_name]['Predicted Present'] += 1
                        if not simulator_result[a_name].item():
                            file_stats[a_name]['Real Absent'] += 1
                        else:
                            file_stats[a_name]['Real Present'] += 1
                            file_stats[a_name]['Correct Present'] += 1

            # else:   TODO?
            #     print(f"Ignoring gene {gene_name}")

    file_stats['All'] = {
        'Real Absent': sum([file_stats[a_name]['Real Absent'] for a_name in relevant_accessions]),
        'Predicted Absent': sum([file_stats[a_name]['Predicted Absent'] for a_name in relevant_accessions]),
        'Correct Absent': sum([file_stats[a_name]['Correct Absent'] for a_name in relevant_accessions]),
        'Real Present': sum([file_stats[a_name]['Real Present'] for a_name in relevant_accessions]),
        'Predicted Present': sum([file_stats[a_name]['Predicted Present'] for a_name in relevant_accessions]),
        'Correct Present': sum([file_stats[a_name]['Correct Present'] for a_name in relevant_accessions])
    }

    print("======================================")
    print(file_tags[i])
    print("======================================")
    print("Num of genes absent by the simulator (real absent):", file_stats['All']['Real Absent'])
    print("Num of genes absent by the simulator (real present):", file_stats['All']['Real Present'])
    print("Num of genes absent by Panoramic (predicted absent):", file_stats['All']['Predicted Absent'])
    print("Num of genes absent by Panoramic (predicted present):", file_stats['All']['Predicted Present'])
    print("Num of genes predicted absent correctly (TP):", file_stats['All']['Correct Absent'])
    print("Num of genes predicted present correctly (TN):", file_stats['All']['Correct Present'])
    file_stats['All']['Accuracy'] = (file_stats['All']['Correct Absent'] + file_stats['All']['Correct Present']) / \
                                    (file_stats['All']['Predicted Absent'] + file_stats['All']['Predicted Present'])
    file_stats['All']['Precision'] = 'N/A' if 0 == file_stats['All']['Predicted Absent'] else \
        file_stats['All']['Correct Absent'] / file_stats['All']['Predicted Absent']
    print("Accuracy:", file_stats['All']['Accuracy'])
    print("Precision:", file_stats['All']['Precision'])

    for a_name in relevant_accessions:
        print("----------")
        print("For accession:", a_name)
        print("Num of genes absent by the simulator (real absent):", file_stats[a_name]['Real Absent'])
        print("Num of genes absent by the simulator (real present):", file_stats[a_name]['Real Present'])
        print("Num of genes absent by Panoramic (predicted absent):", file_stats[a_name]['Predicted Absent'])
        print("Num of genes absent by Panoramic (predicted present):", file_stats[a_name]['Predicted Present'])
        print("Num of genes predicted absent correctly (TP):", file_stats[a_name]['Correct Absent'])
        print("Num of genes predicted present correctly (TN):", file_stats[a_name]['Correct Present'])
        file_stats[a_name]['Accuracy'] = \
            (file_stats[a_name]['Correct Absent'] + file_stats[a_name]['Correct Present']) / \
            (file_stats[a_name]['Predicted Absent'] + file_stats[a_name]['Predicted Present'])
        file_stats[a_name]['Precision'] = 'N/A' if 0 == file_stats[a_name]['Predicted Absent'] else \
            file_stats[a_name]['Correct Absent'] / file_stats[a_name]['Predicted Absent']
        print("Accuracy:", file_stats[a_name]['Accuracy'])
        print("Precision:", file_stats[a_name]['Precision'])

bar_width = 0.25
for file_name, file_stats in stats_per_file.items():
    # Set the number of groups and the values for each group
    N = len(file_stats)
    values1 = [v['Real Absent'] for _, v in file_stats.items()]
    values2 = [v['Predicted Absent'] for _, v in file_stats.items()]
    values3 = [v['Correct Absent'] for _, v in file_stats.items()]

    # Set the labels for each group and the bar width
    labels = [k for k, _ in file_stats.items()]

    # Set the position of the bars on the x-axis
    x_pos = [i for i in range(N)]

    # Build the first set of bars
    plt.bar(x_pos, values1, bar_width, alpha=0.5, color='#4682B4')

    # Build the second set of bars, positioned slightly to the right of the first set
    plt.bar([i + bar_width for i in x_pos], values2, bar_width, alpha=0.5, color='#FA8072')

    # Build the third set of bars, positioned slightly to the right of the second set
    plt.bar([i + bar_width * 2 for i in x_pos], values3, bar_width, alpha=0.5, color='#808000')

    # Set the x-axis tick marks and labels
    plt.xticks([i + bar_width for i in x_pos], labels)

    # Add a legend and title
    plt.legend(['Real Absent', 'Predicted Absent', 'Correct Absent'])
    plt.title(graph_title + f'\n{file_name}')

    # Show the plot
    plt.show()


# Set the number of groups and the values for each group
N = len(stats_per_file)
values1 = [v['All']['Real Absent'] for _, v in stats_per_file.items()]
values2 = [v['All']['Predicted Absent'] for _, v in stats_per_file.items()]
values3 = [v['All']['Correct Absent'] for _, v in stats_per_file.items()]

# Set the labels for each group and the bar width
labels = [k for k, _ in stats_per_file.items()]

# Set the position of the bars on the x-axis
x_pos = [i for i in range(N)]

# Build the first set of bars
plt.bar(x_pos, values1, bar_width, alpha=0.5, color='#4682B4')

# Build the second set of bars, positioned slightly to the right of the first set
plt.bar([i + bar_width for i in x_pos], values2, bar_width, alpha=0.5, color='#FA8072')

# Build the third set of bars, positioned slightly to the right of the second set
plt.bar([i + bar_width * 2 for i in x_pos], values3, bar_width, alpha=0.5, color='#808000')

# Set the x-axis tick marks and labels
plt.xticks([i + bar_width for i in x_pos], labels)

# Add a legend and title
plt.legend(['Real Absent', 'Predicted Absent', 'Correct Absent'])
plt.title(graph_title+'\nComparison')

# Show the plot
plt.show()
