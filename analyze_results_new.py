import argparse
import pandas as pd
import csv
import sys
import matplotlib.pyplot as plt
import re

# TODO remove later
spreadsheet_input = ''

GENES_START_NUMBER = 27655
REF_LABEL = 'TAIR10'
NO_GENE = '0'
GENE_PREF = r'transcript[_:]'
NOVEL_GENE = 'PanGene'
ACCESSION_NAMES = ['A-lyrata', 'Cvi', 'C24', 'An-1', 'Col-0', 'Sha', 'Ler', 'Kyo', 'Eri']
SIMULATOR_PAV_COLUMNS = ['Name', 'Transcript'] + ACCESSION_NAMES
PANORAMIC_PAV_COLUMNS = [f'{name}_{name}' for name in ACCESSION_NAMES] + [REF_LABEL, 'gene']
parser = argparse.ArgumentParser(description='Plot and compare gene stats',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
parser.add_argument('-f', type=str, required=True, metavar='<str>', nargs='+',
                    help="* pan_PAV_1.tsv [pan_PAV_2.tsv] [pan_PAV_3.tsv] ...")
parser.add_argument('-s', type=str, required=True, metavar='<str>', nargs='+',
                    help="* simulator_PAV_1.csv [simulator_PAV_2.csv] [simulator_PAV_3.tsv] ...")
parser.add_argument('-l', type=str, required=True, metavar='<str>', nargs='+',
                    help="* simulator_legend_1.csv [simulator_legend_2.csv] [simulator_legend_3.tsv] ...")
parser.add_argument('-t', type=str, required=True, metavar='<str>', nargs='+',
                    help="* 'De-novo 01' ['Map-to-pan 01'] ['De-novo 02'] ...")
parser.add_argument('-T', type=str, required=True, metavar='<str>', help="Graph Title")
parser.add_argument('-g', type=str, required=True, metavar='<str>', help="Graph Mode (P/N")

args = parser.parse_args()
panoramic_pav_files = args.f
simulator_raw_pav_files = args.s
simulator_legend_files = args.l
file_tags = args.t
graph_title = args.T
graph_mode = args.g

if graph_mode.upper() == 'P':
    mode = 'Positive'
    mode_key = 'Present'
elif graph_mode.upper() == 'N':
    mode = 'Negative'
    mode_key = 'Absent'
else:
    raise Exception("Choose graph mode")
graph_title = mode_key+' '+graph_title

stats_per_file = {tag: {} for tag in file_tags}
new_stats_per_file = {tag: {} for tag in file_tags}

if len(panoramic_pav_files) != len(simulator_raw_pav_files) or len(simulator_raw_pav_files) != len(
        simulator_legend_files):
    print("Number of parameters incorrect")
    sys.exit(1)

for i in range(len(file_tags)):

    ###############################
    ### Simulator real values ###
    ###############################
    sim_pav_df = pd.read_csv(simulator_raw_pav_files[i], index_col=0,
                             dtype={'gene_id': 'int', 'chrom': 'str', 'ref_start': 'int', 'ref_end': 'int'})
    legend_df = pd.read_csv(simulator_legend_files[i], index_col=0, dtype={'Name': 'str', 'ID': 'int'})
    legend_df = legend_df.rename(columns={'ID': 'gene_id'})
    sim_pav_df = pd.merge(sim_pav_df, legend_df, how='inner', on=['gene_id'])
    sim_pav_df[['Name', 'Transcript']] = sim_pav_df.apply(lambda x: x['Name'].split('.'), axis=1, result_type='expand')
    sim_pav_df = sim_pav_df.loc[:, sim_pav_df.columns.isin(SIMULATOR_PAV_COLUMNS)]  # select this columns if exist
    sim_known_genes = len(sim_pav_df)

    relevant_accessions = [column for column in sim_pav_df if column in ACCESSION_NAMES]


    ###################################
    ### Panoramic predicted values ###
    ###################################
    with open(panoramic_pav_files[i]) as file:
        file_stats = stats_per_file[file_tags[i]]
        pano_pav_index = {}
        common_known_genes = 0
        pano_known_genes = 0
        pano_new_genes = 0

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
                        'Correct Present': 0,
                        'Novel': 0
                    }
                continue

            if line[0].startswith(NOVEL_GENE):  # Novel genes
                pano_new_genes += 1
                for a_name in relevant_accessions:
                    if NO_GENE == line[pano_pav_index[a_name]]:
                        pass
                        # 15.03.2023 - Itay decided not to count TN in this case
                        # file_stats[a_name]['Real Absent'] += 1
                        # file_stats[a_name]['Novel'] += 1
                        # file_stats[a_name]['Predicted Absent'] += 1
                        # file_stats[a_name]['Correct Absent'] += 1
                    else:
                        file_stats[a_name]['Real Absent'] += 1
                        file_stats[a_name]['Novel'] += 1
                        file_stats[a_name]['Predicted Present'] += 1

            else:  # Known genes
                pano_known_genes += 1
                line.pop(pano_pav_index[REF_LABEL])
                gene_name = re.sub(GENE_PREF, '', line[0])
                gene_name = gene_name.split('.')[0]
                simulator_result = sim_pav_df[sim_pav_df['Name'] == gene_name]
                if not simulator_result.empty:
                    common_known_genes += 1
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
        'Novel': sum([file_stats[a_name]['Novel'] for a_name in relevant_accessions]),
        'Predicted Absent': sum([file_stats[a_name]['Predicted Absent'] for a_name in relevant_accessions]),
        'Correct Absent': sum([file_stats[a_name]['Correct Absent'] for a_name in relevant_accessions]),
        'Real Present': sum([file_stats[a_name]['Real Present'] for a_name in relevant_accessions]),
        'Predicted Present': sum([file_stats[a_name]['Predicted Present'] for a_name in relevant_accessions]),
        'Correct Present': sum([file_stats[a_name]['Correct Present'] for a_name in relevant_accessions]),
        'Simulator Known': sim_known_genes,
        'Panoramic Known': pano_known_genes,
        'Common Known': common_known_genes,
        'Panoramic New': pano_new_genes
    }

    print("======================================")
    print(file_tags[i])
    print("======================================")

    accuracy_sum = 0
    accuracy_count = 0
    precision_sum = 0
    precision_count = 0
    f1_sum = 0
    f1_count = 0
    for a_name in relevant_accessions:
        print("For accession:", a_name)
        print("Num of genes absent by the simulator (actual N):", file_stats[a_name]['Real Absent'])
        print("Num of reference genes out of them:", file_stats[a_name]['Real Absent'] - file_stats[a_name]['Novel'])
        print("Num of novel genes out of them:", file_stats[a_name]['Novel'])
        print("Num of genes present by the simulator (actual P):", file_stats[a_name]['Real Present'])
        print("Num of genes absent by Panoramic (predicted N):", file_stats[a_name]['Predicted Absent'])
        print("Num of genes present by Panoramic (predicted P):", file_stats[a_name]['Predicted Present'])
        print("Num of reference genes out of them:", file_stats[a_name]['Predicted Present'] - file_stats[a_name]['Novel'])
        print("Num of genes predicted absent correctly (TN):", file_stats[a_name]['Correct Absent'])
        print("Num of genes predicted absent wrongly (FN):", file_stats[a_name]['Predicted Absent'] - file_stats[a_name]['Correct Absent'])
        print("Num of genes predicted present correctly (TP):", file_stats[a_name]['Correct Present'])
        print("Num of genes predicted present wrongly (FP):", file_stats[a_name]['Predicted Present'] - file_stats[a_name]['Correct Present'])
        print("Among them reference genes:", file_stats[a_name]['Predicted Present'] - file_stats[a_name]['Correct Present'] - file_stats[a_name]['Novel'])

        # Accuracy = (TP + TN) / (TP + FP + TN + FN)
        # note that Predicted Absent is (TN + FN) while as Predicted Present is TP + FP
        file_stats[a_name]['Accuracy'] = \
            (file_stats[a_name]['Correct Absent'] + file_stats[a_name]['Correct Present']) / \
            (file_stats[a_name]['Predicted Absent'] + file_stats[a_name]['Predicted Present'])
        accuracy_sum += file_stats[a_name]['Accuracy']
        accuracy_count += 1

        # Precision = TP / (TP + FP)
        # note that Predicted Present is TP + FP
        if file_stats[a_name]['Predicted Present'] != 0:
            file_stats[a_name]['Precision'] = \
                (file_stats[a_name]['Correct Present'] / file_stats[a_name]['Predicted Present'])
            precision_sum += file_stats[a_name]['Precision']
            precision_count += 1
        else:
            file_stats[a_name]['Precision'] = 'N/A'

        # Recall = TP / (TP + FN)
        # note that (Predicted Absent - Correct Absent) is FN
        if file_stats[a_name]['Predicted Present'] != 0:
            file_stats[a_name]['Recall'] =\
                (file_stats[a_name]['Correct Present'] / (file_stats[a_name]['Correct Present']
                + file_stats[a_name]['Predicted Absent'] - file_stats[a_name]['Correct Absent']))
        else:
            file_stats[a_name]['Recall'] = 'N/A'

        # F1 = 2 * (Precision * Recall) / (Precision + Recall))
        if file_stats[a_name]['Correct Present'] != 0 and file_stats[a_name]['Recall'] != 'N/A' and \
                file_stats[a_name]['Precision'] != 'N/A':
            file_stats[a_name]['F1'] = (2 * file_stats[a_name]['Recall'] * file_stats[a_name]['Precision'] /
                  (file_stats[a_name]['Recall'] + file_stats[a_name]['Precision']))
            f1_sum += file_stats[a_name]['F1']
            f1_count += 1
        else:
            file_stats[a_name]['F1'] = 'N/A'

        print("Accuracy:", file_stats[a_name]['Accuracy'])
        print("Precision:", file_stats[a_name]['Precision'])
        print("Recall:", file_stats[a_name]['Recall'])
        print("F1:", file_stats[a_name]['F1'])
        print("----------")

    file_stats['All']['Accuracy'] = accuracy_sum / accuracy_count
    file_stats['All']['Precision'] = precision_sum / precision_count
    file_stats['All']['F1'] = f1_sum / f1_count

    print("Simulator\'s known genes:", file_stats['All']['Simulator Known'], "which are the 'true' pan-genes")
    simulation_discarded = GENES_START_NUMBER - file_stats['All']['Simulator Known']
    print(f"{simulation_discarded} genes were discarded before simulation even started, "
          f"which are {simulation_discarded/GENES_START_NUMBER} of the genes.")
    print("Panoramic\'s known genes:", file_stats['All']['Panoramic Known'])
    print("Common known genes:", file_stats['All']['Common Known'])
    common_percentage = file_stats['All']['Common Known'] / file_stats['All']['Panoramic Known']
    print(f"Out of Panoramic known genes, only {common_percentage} are relevant")
    print("Panoramic\'s new genes:", file_stats['All']['Panoramic New'])
    print("--------------------------------------")
    print("Num of genes absent by the simulator (actual N):", file_stats['All']['Real Absent'])
    print("Num of novel genes out of them:", file_stats['All']['Novel'])
    print("Num of reference genes out of them:", file_stats['All']['Real Absent'] - file_stats['All']['Novel'])
    print("Num of genes present by the simulator (actual P):", file_stats['All']['Real Present'])
    print("Num of genes absent by Panoramic (predicted N):", file_stats['All']['Predicted Absent'])
    print("Num of genes present by Panoramic (predicted P):", file_stats['All']['Predicted Present'])
    print("Num of reference genes out of them:", file_stats['All']['Predicted Present'] - file_stats['All']['Novel'])
    print("Num of genes predicted absent correctly (TN):", file_stats['All']['Correct Absent'])
    print("Num of genes predicted absent wrongly (FN):", file_stats['All']['Predicted Absent'] - file_stats['All']['Correct Absent'])
    print("Num of genes predicted present correctly (TP):", file_stats['All']['Correct Present'])
    print("Num of genes predicted present wrongly (FP):", file_stats['All']['Predicted Present'] - file_stats['All']['Correct Present'])
    print("Among them reference genes:", file_stats['All']['Predicted Present'] - file_stats['All']['Correct Present'] - file_stats['All']['Novel'])

    # file_stats['All']['Accuracy'] = (file_stats['All']['Correct Absent'] + file_stats['All']['Correct Present']) / \
    #                                 (file_stats['All']['Predicted Absent'] + file_stats['All']['Predicted Present'])
    # file_stats['All']['Precision'] = 'N/A' if 0 == file_stats['All']['Predicted Absent'] else \
    #     (file_stats['All']['Correct Absent'] / file_stats['All']['Predicted Absent'])
    # file_stats['All']['Recall'] = 'N/A' if 0 == file_stats['All']['Predicted Absent'] else \
    #     (file_stats['All']['Correct Absent'] /
    #      (file_stats['All']['Correct Absent'] + file_stats['All']['Predicted Present']
    #       - file_stats['All']['Correct Present']))
    # file_stats['All']['F1'] = 'N/A'\
    #     if file_stats['All']['Correct Absent'] == 0 or 'N/A' == file_stats['All']['Recall'] \
    #        or 'N/A' == file_stats['All']['Precision'] \
    #     else (2 * file_stats['All']['Recall'] * file_stats['All']['Precision'] /
    #           (file_stats['All']['Recall'] + file_stats['All']['Precision']))
    print("Accuracy:", file_stats['All']['Accuracy'])
    print("Precision:", file_stats['All']['Precision'])
    print("F1:", file_stats['All']['F1'])

    # TODO remove later
    # if not spreadsheet_input:
    #     spreadsheet_input = spreadsheet_input \
    #                     + str(file_stats['All']['Real Present']) + '\n' \
    #                     + str(file_stats['All']['Real Absent'] - file_stats['All']['Novel']) + '\n' #\
    #                     #+ str(file_stats['All']['Novel']) + '\n'
    #     for a_name in relevant_accessions:
    #         spreadsheet_input = spreadsheet_input \
    #                         + str(file_stats[a_name]['Real Present']) + '\n' \
    #                         + str(file_stats[a_name]['Real Absent'] - file_stats[a_name]['Novel']) + '\n' #\
    #                         #+ str(file_stats[a_name]['Novel']) + '\n'

    # numbers
    # spreadsheet_input = spreadsheet_input \
    #                     + str(file_stats['All']['Correct Present']) + '\n' \
    #                     + str(file_stats['All']['Correct Absent']) + '\n'
    for a_name in relevant_accessions:
        spreadsheet_input = spreadsheet_input \
                        + str(file_stats[a_name]['Correct Present']) + '\n' \
                        + str(file_stats[a_name]['Correct Absent']) + '\n'
    for a_name in relevant_accessions:
        spreadsheet_input = spreadsheet_input \
                        + str(file_stats[a_name]['Predicted Present'] - file_stats[a_name]['Correct Present']) + '\n' \
                        + str(file_stats[a_name]['Predicted Absent'] - file_stats[a_name]['Correct Absent']) + '\n'
    # # percentages
    # spreadsheet_input = spreadsheet_input \
    #                     + str(file_stats['All']['Correct Present']/file_stats['All']['Real Present']) + '\n' \
    #                     + str(file_stats['All']['Correct Absent']/file_stats['All']['Real Absent']) + '\n'
    # for a_name in relevant_accessions:
    #     spreadsheet_input = spreadsheet_input \
    #                     + str(file_stats[a_name]['Correct Present']/file_stats[a_name]['Real Present']) + '\n' \
    #                     + str(file_stats[a_name]['Correct Absent']/file_stats[a_name]['Real Absent']) + '\n'


# TODO remove later
print("----------")
print('Spreadsheet input:')
print(spreadsheet_input)

bar_width = 0.35
for file_name, file_stats in stats_per_file.items():
    # Set the number of groups and the values for each group
    N = len(file_stats) - 1
    values1 = [v[f'Real {mode_key}'] for k, v in file_stats.items() if k != 'All']
    if mode_key == 'Present':
        values1 = [v[f'Real Present'] for k, v in file_stats.items() if k != 'All']
    else:
        values1 = [v[f'Real Absent'] - v['Novel'] for k, v in file_stats.items() if k != 'All']

    values2 = [v[f'Predicted {mode_key}'] for k, v in file_stats.items() if k != 'All']
    values3 = [v[f'Correct {mode_key}'] for k, v in file_stats.items() if k != 'All']

    # Set the labels for each group and the bar width
    labels = [k for k, _ in file_stats.items() if k != 'All']

    # Set the position of the bars on the x-axis
    x_pos = [i for i in range(N)]

    # Build the first set of bars
    plt.bar(x_pos, values1, bar_width, alpha=0.5, color='#4682B4')

    # Build the second set of bars, positioned slightly to the right of the first set
    plt.bar([i + bar_width for i in x_pos], values2, bar_width, alpha=0.5, color='#FA8072')

    # Build the third set of bars, positioned at the same x-position as the second set of bars
    plt.bar([i + bar_width for i in x_pos], values3, bar_width, alpha=0.5, color='#808000')

    # Set the x-axis tick marks and labels
    plt.xticks([i + bar_width/2 for i in x_pos], labels)

    # Add a legend and title
    plt.legend(['Actual Negative\n(Simulator)', 'Predicted Negative\n(Panoramic pipeline)', 'True Negative'],
           bbox_to_anchor=(0.5, -0.15), loc='center', ncol=3)
    plt.title(graph_title + f'\n{file_name}')

    plt.tight_layout()

    # Show the plot
    plt.show()

# Set the number of groups and the values for each group
N = len(stats_per_file)
if mode_key == 'Present':
    values1 = [v['All']['Real Present'] for _, v in stats_per_file.items()]
else:
    values1 = [v['All']['Real Absent'] - v['All']['Novel'] for _, v in stats_per_file.items()]
values2 = [v['All'][f'Predicted {mode_key}'] for _, v in stats_per_file.items()]
values3 = [v['All'][f'Correct {mode_key}'] for _, v in stats_per_file.items()]

# Set the labels for each group and the bar width
labels = [k for k, _ in stats_per_file.items()]

# Set the position of the bars on the x-axis
x_pos = [i for i in range(N)]

# Build the first set of bars
plt.bar(x_pos, values1, bar_width, alpha=0.5, color='#4682B4')

# Build the second set of bars, positioned slightly to the right of the first set
plt.bar([i + bar_width for i in x_pos], values2, bar_width, alpha=0.5, color='#FA8072')

# Build the third set of bars, positioned at the same x-position as the second set of bars
plt.bar([i + bar_width for i in x_pos], values3, bar_width, alpha=0.5, color='#808000')

# Set the x-axis tick marks and labels
plt.xticks([i + bar_width/2 for i in x_pos], labels)

# Add a legend and title
plt.legend([f'Actual {mode}\n(Simulator)', f'Predicted {mode}\n(Panoramic pipeline)', f'True {mode}'],
           bbox_to_anchor=(0.5, -0.15), loc='center', ncol=3)
plt.title(graph_title + '\nComparison')

plt.tight_layout()

# Show the plot
plt.show()
