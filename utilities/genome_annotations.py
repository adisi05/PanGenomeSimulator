import argparse
import os
import time
import re
import numpy as np
from os.path import exists
from typing import Dict, List, Tuple

import operator
import pandas as pd
import pybedtools

from common_data_structues import Strand, Region

DEFAULT_INPUT_FILE = '/groups/itay_mayrose/adisivan/arabidopsis/ensemblgenomes/gff3/Arabidopsis_thaliana.TAIR10.53.gff3'
DEFAULT_OUTPUT_FILE = 'all_chroms_annotations.csv'


def workflow(input_path: str, output_path: str, legend_output_path: str):
    genes_df, cds_elements_df, chroms_df = extract_data_frames(input_path)

    genes_df = get_genes_without_overlaps(genes_df)
    relevant_genes_df, sub_gene_annotations_df = get_sub_gene_annotations(genes_df, cds_elements_df)
    intergenics_df = get_intergenic_annotations(chroms_df, relevant_genes_df)

    all_annotations_df = pd.concat([sub_gene_annotations_df, intergenics_df])  # ignore_index=True ?
    all_annotations_df = all_annotations_df.sort_values(by=['chrom', 'start', 'end']).reset_index(drop=True)
    all_annotations_df, genes_legend_df = gene_ids_to_numbers(all_annotations_df)
    all_annotations_df[['start', 'end', 'gene_id']] = all_annotations_df[['start', 'end', 'gene_id']].astype(int)
    all_annotations_df = assign_frame_to_cds(all_annotations_df)
    sanity_check(all_annotations_df)
    all_annotations_df.to_csv(output_path)
    genes_legend_df.to_csv(legend_output_path)


def extract_data_frames(file_path: str) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    if file_path is not None:
        try:
            using_bed = file_path.lower().endswith('.bed')
            annotations = pybedtools.example_bedtool(file_path)

            # genes
            gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene')
            gene_elements = gene_elements.sort()
            gene_elements_df = gene_elements.to_dataframe()
            gene_elements_df = gene_elements_df[['seqname', 'start', 'end', 'attributes', 'strand']]
            gene_elements_df.rename(columns={'seqname': 'chrom'}, inplace=True)
            gene_elements_df['gene_id'] = gene_elements_df.apply(
                lambda x: extract_gene_id_from_gene(x['attributes']), axis=1, result_type='expand')
            del gene_elements_df['attributes']
            if not using_bed:
                gene_elements_df['start'] -= 1

            # CDS elements
            cds_elements = annotations.filter(lambda x: x.fields[2] == Region.CDS.value)
            cds_elements = cds_elements.sort()
            cds_elements_df = cds_elements.to_dataframe()
            cds_elements_df = cds_elements_df[['seqname', 'start', 'end', 'strand', 'frame', 'attributes']]
            cds_elements_df.rename(columns={'seqname': 'chrom'}, inplace=True)
            cds_elements_df[['gene_id', 'variant_id']] = cds_elements_df.apply(
                lambda x: extract_gene_id_variant_id_from_cds(x['attributes']), axis=1, result_type='expand')
            del cds_elements_df['attributes']
            if not using_bed:
                cds_elements_df['start'] -= 1

            # chromosomes
            chrom_elements = annotations.filter(lambda x: x.fields[2] == 'chromosome')
            chrom_elements = chrom_elements.sort()
            chrom_elements_df = chrom_elements.to_dataframe()
            chrom_elements_df = chrom_elements_df[['seqname', 'start', 'end']]
            chrom_elements_df.rename(columns={'seqname': 'chrom'}, inplace=True)
            if not using_bed:
                chrom_elements_df['start'] -= 1

            return gene_elements_df, cds_elements_df, chrom_elements_df

        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")


def extract_gene_id_from_gene(gff_attributes: str) -> str:
    id_attribute = gff_attributes.split(';')[0]
    return id_attribute.split(':')[1] if len(id_attribute.split(':')) == 2 else None


def extract_gene_id_variant_id_from_cds(gff_attributes: str) -> (str, int):
    if gff_attributes:
        if len(gff_attributes.split(';')) > 1:
            id_attribute = gff_attributes.split(';')[0]
            if len(id_attribute.split(':')) == 2:
                gene_with_variant = id_attribute.split(':')[1]
                if len(gene_with_variant.split('.')) == 2:
                    gene_variant = gene_with_variant.split('.')
                    return gene_variant[0], int(gene_variant[1])

    return None


def get_genes_without_overlaps(genes_df: pd.DataFrame) -> pd.DataFrame:
    overlapping_indexes = set()
    prev_gene_end = 0
    curr_chrom = ''
    for i, gene in genes_df.iterrows():
        if gene['chrom'] != curr_chrom:
            curr_chrom = gene['chrom']
            prev_gene_end = gene['end']
            continue
        if gene['start'] < prev_gene_end:
            # this gene is to discard due to overlapping with previous genes
            overlapping_indexes.add(i)
        else:
            # this gene is to keep
            prev_gene_end = gene['end']

    gene_count_before = len(genes_df)
    genes_df = genes_df[~genes_df.index.isin(overlapping_indexes)] \
        .reset_index(drop=True)
    print(f'Found {len(genes_df)} non-overlapping genes out of {gene_count_before}.')
    return genes_df


def get_sub_gene_annotations(genes_df: pd.DataFrame, cds_elements_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    valid_gene_indices = []
    elements_list = []
    for i, gene in genes_df.iterrows():
        gene_id = gene['gene_id']

        # choose the longest variants that its length is dividable by 3
        var_lenths_dict = get_lengths_of_gene_variants(gene_id, cds_elements_df)
        var_lenths_dict = dict(sorted(var_lenths_dict.items(), key=operator.itemgetter(1), reverse=True))
        chosen_variant = -1
        for var_id, length in var_lenths_dict.items():
            if length % 3 == 0:
                chosen_variant = var_id
                break
        if chosen_variant == -1:
            print(f"Did not found an appropriate variant for gene {gene_id} on chromosome {gene['chrom']}."
                  f" Ignoring this gene")
            continue

        valid_gene_indices.append(i)
        current_gene_elements_list = get_sub_gene_elements(gene_id, genes_df, chosen_variant, cds_elements_df)
        elements_list.extend(current_gene_elements_list)

    valid_genes_df = genes_df[genes_df.index.isin(valid_gene_indices)]
    sub_gene_annotations_df = pd.DataFrame(elements_list)

    return valid_genes_df, sub_gene_annotations_df


def get_lengths_of_gene_variants(gene_id: str, cds_elements_df: pd.DataFrame) -> Dict[int, int]:
    var_lenth_dict = {}
    variant_ids = list(cds_elements_df[cds_elements_df['gene_id'] == gene_id]['variant_id'].unique())
    for var_id in variant_ids:
        variant_elements = cds_elements_df[(cds_elements_df['gene_id'] == gene_id) &
                                           (cds_elements_df['variant_id'] == var_id)]
        var_lenth_dict[int(var_id)] = int(variant_elements.apply(
            lambda x: x['end'] - x['start'], axis=1, result_type='expand').sum())
    return var_lenth_dict


def get_sub_gene_elements(gene_id: str, genes_df: pd.DataFrame, var_id: int, cds_elements_df: pd.DataFrame)\
        -> List[dict]:
    gene = genes_df[genes_df['gene_id'] == gene_id]
    gene_start = gene['start'].item()
    gene_end = gene['end'].item()
    gene_strand = gene['strand'].item()
    gene_chrom = gene['chrom'].item()
    elements_list = []

    variants_cds_elements_df = cds_elements_df[(cds_elements_df['gene_id'] == gene_id) &
                                               (cds_elements_df['variant_id'] == var_id)]
    cursor = gene_start
    for i, cds_element in variants_cds_elements_df.iterrows():
        cds_start = cds_element['start']
        cds_end = cds_element['end']

        # NON_CODING_GENE element before current CDS
        if cursor < cds_start:
            non_coding_gene_dict = {
                'chrom': gene_chrom,
                'gene_id': gene_id + '.' + str(var_id),
                'region': Region.NON_CODING_GENE.value,
                'start': cursor,
                'end': cds_start,
                'strand': gene_strand,
                'frame': np.nan
            }
            elements_list.append(non_coding_gene_dict)

        # CDS element
        cds_dict = {
            'chrom': gene_chrom,
            'gene_id': gene_id + '.' + str(var_id),
            'region': Region.CDS.value,
            'start': cds_start,
            'end': cds_end,
            'strand': gene_strand,
            'frame': cds_element['frame']
        }
        elements_list.append(cds_dict)
        cursor = cds_end

    # NON_CODING_GENE element after last CDS
    if cursor != gene_end:
        non_coding_gene_dict = {
            'chrom': gene_chrom,
            'gene_id': gene_id + '.' + str(var_id),
            'region': Region.NON_CODING_GENE.value,
            'start': cursor,
            'end': gene_end,
            'strand': gene_strand,
            'frame': np.nan
        }
        elements_list.append(non_coding_gene_dict)

    return elements_list


def get_intergenic_annotations(chroms_df: pd.DataFrame, genes_df: pd.DataFrame) -> pd.DataFrame:
    elements_list = []

    for _, chrom in chroms_df.iterrows():
        chrom_name = chrom['chrom']
        chrom_end = chrom['end']
        cursor = 0  # chrom start at 0 always
        for _, gene in genes_df[genes_df['chrom'] == chrom_name].iterrows():

            # INTERGENIC element before current gene
            if cursor < gene['start']:
                intergenic_element = {
                    'chrom': chrom_name,
                    'gene_id': '',
                    'region': Region.INTERGENIC.value,
                    'start': cursor,
                    'end': gene['start'],
                    'strand': Strand.UNKNOWN.value
                }
                elements_list.append(intergenic_element)

            cursor = gene['end']

        # INTERGENIC element after last gene
        if cursor != chrom_end:
            intergenic_element = {
                'chrom': chrom_name,
                'gene_id': '',
                'region': Region.INTERGENIC.value,
                'start': cursor,
                'end': chrom_end,
                'strand': Strand.UNKNOWN.value
            }
            elements_list.append(intergenic_element)

    intergenics_df = pd.DataFrame(elements_list)
    return intergenics_df


def gene_ids_to_numbers(all_annotations_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    name_to_number = {'': 0}
    id_counter = 1
    for i, annotation in all_annotations_df.iterrows():
        if annotation['gene_id'] in name_to_number:
            all_annotations_df.loc[i, 'gene_id'] = name_to_number[annotation['gene_id']]
        else:
            name_to_number[annotation['gene_id']] = id_counter
            all_annotations_df.loc[i, 'gene_id'] = id_counter
            id_counter += 1
    name_number_df = pd.DataFrame([(name, number) for name, number in name_to_number.items()], columns=['Name', 'ID'])

    return all_annotations_df, name_number_df


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='Genome annotations to CSV format',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-a', type=str, required=True, metavar='/path/to/annotations.gff',
                        help="Annotations file in one of the next formats: bed, gff, gtf, gff3")
    parser.add_argument('-o', type=str, required=False, metavar='/path/to/output.csv',
                        help="Path to output csv file")
    parser.add_argument('--test', required=False, default=False, help="Output overlapping genes",
                        action='store_true')
    args = parser.parse_args(raw_args)

    input_file = args.a if args.a else DEFAULT_INPUT_FILE
    output_file = args.o if args.o else DEFAULT_OUTPUT_FILE
    output_file = output_file if output_file.lower().endswith('.csv') else DEFAULT_OUTPUT_FILE
    legend_file = re.split('.csv', output_file, flags=re.IGNORECASE)[0] + '_legend.csv'
    test_mode = args.test

    if not test_mode:
        workflow(input_file, output_file, legend_file)
    else:
        test_overlapping_genes(input_file)


def assign_frame_to_cds(annotations_df: pd.DataFrame) -> pd.DataFrame:
    annotations_df = assign_frame_to_cds_on_strand(annotations_df, Strand.FORWARD)
    annotations_df = assign_frame_to_cds_on_strand(annotations_df, Strand.REVERSE)
    return annotations_df


def assign_frame_to_cds_on_strand(annotations_df: pd.DataFrame, strand: Strand) -> pd.DataFrame:
    print(f"Assigning frames for CDS elements in strand {strand.value}...")
    t_start = time.time()
    index = annotations_df.index
    if strand.value == Strand.REVERSE.value:
        index = reversed(index)

    current_gene_id = 0
    previous_offset = 0
    for i in index:
        annotation = annotations_df.iloc[i, :]
        if annotation['strand'] != strand.value:
            continue

        # encountered a new gene
        if annotation['gene_id'] != current_gene_id and annotation['gene_id'] != 0:
            current_gene_id = annotation['gene_id']
            previous_offset = 0

        if annotation['region'] == Region.CDS.value:
            if annotation['frame'] != previous_offset:
                annotations_df.iloc[i, annotations_df.columns.get_loc('frame')] = previous_offset

            previous_offset = (previous_offset + annotation['end'] - annotation['start']) % 3
    print(f"Completed frames assignment, took {int(time.time() - t_start)} seconds.")
    return annotations_df


#####################
### Sanity Checks ###
#####################


def sanity_check(annotations_df: pd.DataFrame):
    start_t = time.time()
    print('Starting sanity check...')
    sanity_check_elements_positions(annotations_df)
    sanity_check_gene_strand(annotations_df)
    sanity_check_cds_frames(annotations_df)
    print('Sanity check completed. Took {0:.3f} (sec)'.format(time.time() - start_t))


def sanity_check_elements_positions(annotations_df):
    print("Starting positions sanity check")
    count = 0
    current_chrom = ''
    previous_end = -1
    for i, annotation in annotations_df.iterrows():
        count += 1
        # first annotation in chromosome
        if current_chrom != annotation['chrom']:
            current_chrom = annotation['chrom']
            previous_end = annotation['end']
            continue

        if annotation['start'] != previous_end:
            print(f'Wrong start for annotation index {i}')

        previous_end = annotation['end']
    print(f"Positions sanity check completed, {count} annotations were checked")


def sanity_check_gene_strand(annotations_df):
    print("Starting strand sanity check")
    count = 0
    for i, annotation in annotations_df.iterrows():
        count += 1
        if annotation['region'] in [Region.CDS.value, Region.NON_CODING_GENE.value] and \
                annotation['strand'] not in [Strand.FORWARD.value, Strand.REVERSE.value]:
            print(f'Wrong strand for annotation index {i}')
    print(f"Strand sanity check completed, {count} annotations were checked")


def sanity_check_cds_frames(annotations_df):

    sanity_check_cds_frames_for_strand(annotations_df, Strand.FORWARD)

    reversed_df = annotations_df.reindex(index=annotations_df.index[::-1])
    sanity_check_cds_frames_for_strand(reversed_df, Strand.REVERSE)


def sanity_check_cds_frames_for_strand(annotations_df: pd.DataFrame, strand: Strand):
    print(f"Sanity check for CDS frames in strand {strand.value}...")
    current_gene_id = 0
    previous_offset = 0
    count = 0
    invalid_count = 0
    for i, annotation in annotations_df.iterrows():
        if annotation['strand'] != strand.value:
            continue

        # encountered a new gene
        if annotation['gene_id'] != current_gene_id and annotation['gene_id'] != 0:
            current_gene_id = annotation['gene_id']
            previous_offset = 0

        if annotation['region'] == Region.CDS.value:
            count += 1
            if annotation['frame'] != previous_offset:
                # print(f'Wrong frame for annotation index {i}')
                invalid_count += 1
            previous_offset = (previous_offset + annotation['end'] - annotation['start']) % 3
    print(f"Out of {count} CDS elements, {invalid_count} have invalid frames")


def external_sanity_check(raw_args=None):
    parser = argparse.ArgumentParser(description='Genome annotations sanity check',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-csv', type=str, required=False, metavar='/path/to/output.csv',
                        help="Path to output csv file")
    args = parser.parse_args(raw_args)

    output_file = args.csv if args.csv else DEFAULT_OUTPUT_FILE
    annotations_df = read_annotations_csv(output_file)

    sanity_check(annotations_df)


############
### Test ###
############


def test_overlapping_genes(in_file: str = DEFAULT_INPUT_FILE):
    if in_file is not None:
        try:
            using_bed = in_file.lower().endswith('.bed')
            annotations = pybedtools.example_bedtool(in_file)
            gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene')
            gene_elements = gene_elements.sort()
            gene_elements_df = gene_elements.to_dataframe()
            gene_elements_df['gene_id'] = gene_elements_df.apply(
                lambda x: extract_gene_id_from_gene(x['attributes']), axis=1, result_type='expand')
            del gene_elements_df['attributes']
            if not using_bed:
                gene_elements_df['start'] -= 1

            overlap_count = 0
            overlap_count_same_strand = 0
            chrom_overlap_count = 0
            chrom_overlap_count_same_strand = 0
            chrom_gene_count = 0
            chrom = ''
            overlapping_indexes = set()
            prev_gene_end = 0
            for i, gene in gene_elements_df.iterrows():
                if i == 0:
                    chrom = gene['chrom']
                    chrom_gene_count += 1
                    continue
                if chrom != gene['chrom']:
                    print(f'Chromosome: {chrom}, overlapping occurrences: {chrom_overlap_count}')
                    print(f'On the same strand: {chrom_overlap_count_same_strand}')
                    print(f'Chromosome: {chrom}, all genes count: {chrom_gene_count}')
                    print('----------')

                    chrom = gene['chrom']
                    chrom_overlap_count = 0
                    chrom_overlap_count_same_strand = 0
                    chrom_gene_count = 0
                    prev_gene_end = 0

                chrom_gene_count += 1
                prev_gene = gene_elements_df.iloc[i - 1, :]
                if gene['chrom'] == prev_gene['chrom'] and gene['start'] < prev_gene_end:
                    overlapping_indexes.add(i)
                    overlapping_indexes.add(i - 1)
                    overlap_count += 1
                    chrom_overlap_count += 1
                    if gene['strand'] == prev_gene['strand']:
                        overlap_count_same_strand += 1
                        chrom_overlap_count_same_strand += 1

            print(f'Chromosome: {chrom}, overlapping occurrences: {chrom_overlap_count}')
            print(f'On the same strand: {chrom_overlap_count_same_strand}')
            print(f'Chromosome: {chrom}, all genes count: {chrom_gene_count}')
            print('----------')
            print(f'Overall overlapping occurrences: {overlap_count}')
            print(f'On the same strand: {overlap_count_same_strand}')
            print(f'All genes count: {len(gene_elements_df)}')
            gene_elements_df = gene_elements_df.iloc[list(overlapping_indexes), :].reset_index(drop=True)
            gene_elements_df.to_csv('overlapping_genes.csv')

        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")


if __name__ == "__main__":
    main()


def read_annotations_csv(file_path: str):
    annotations_df = None
    if file_path:
        print('Loading csv...')
        try:
            _, extension = os.path.splitext(file_path)
            if extension != '.csv':
                print(f'Annotations file must be a csv file. Got extension: {extension}')
                raise Exception
            if exists(file_path):
                annotations_df = pd.read_csv(file_path, index_col=0)
                annotations_df[['chrom', 'region', 'strand']] = \
                    annotations_df[['chrom', 'region', 'strand']].astype(str)
                annotations_df[['start', 'end', 'gene_id']] = \
                    annotations_df[['start', 'end', 'gene_id']].astype(int)
            else:
                print(f'Annotations file does not exist. File path = {file_path}')
                raise Exception
        except Exception:
            print('Problem parsing annotations file')
    return annotations_df


def get_all_genes(file_path: str) -> Dict[str, Tuple[str, int, int]]:
    annotations_df = read_annotations_csv(file_path)
    relevant_genes = {}
    for gene_id in annotations_df.gene_id.unique():
        if gene_id == 0:
            continue
        gene_annotations = annotations_df[annotations_df['gene_id'] == gene_id]
        chrom = gene_annotations.iloc[0]['chrom']
        start = gene_annotations.iloc[0]['start']
        end = gene_annotations.iloc[-1]['end']
        relevant_genes[gene_id] = (chrom, start, end)
    return relevant_genes

# if __name__ == "__main__":
#     external_sanity_check()
