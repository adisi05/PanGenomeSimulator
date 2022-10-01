import argparse
import os
from os.path import exists
from typing import Dict, List

import operator
import pandas as pd
import pybedtools

from utilities.common_data_structues import Strand, Region

DEFAULT_INPUT_FILE = '/groups/itay_mayrose/adisivan/arabidopsis/ensemblgenomes/gff3/Arabidopsis_thaliana.TAIR10.53.gff3'
DEFAULT_OUTPUT_FILE = 'all_chroms_annotations.csv'


def test_overlapping_genes(in_file: str = DEFAULT_INPUT_FILE, using_bed: bool = False):
    if in_file is not None:
        try:
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
            for i, gene in gene_elements_df.iterrows():
                if i == 0:
                    chrom = gene['seqname']
                    chrom_gene_count += 1
                    continue
                if chrom != gene['seqname']:
                    print(f'Chromosome: {chrom}, overlapping occurrences: {chrom_overlap_count}')
                    print(f'On the same strand: {chrom_overlap_count_same_strand}')
                    print(f'Chromosome: {chrom}, all genes count: {chrom_gene_count}')
                    print('----------')

                    chrom = gene['seqname']
                    chrom_overlap_count = 0
                    chrom_overlap_count_same_strand = 0
                    chrom_gene_count = 0

                chrom_gene_count += 1
                prev_gene = gene_elements_df.iloc[i - 1, :]
                if genes_overlap(prev_gene['seqname'], prev_gene['start'], prev_gene['end'],
                                 gene['seqname'], gene['start'], gene['end']):
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


def get_gene_ids_without_overlaps(genes_df: pd.DataFrame) -> pd.DataFrame:
    overlapping_indexes = set()
    for i, gene in genes_df.iterrows():
        if i == 0:
            continue
        prev_gene = genes_df.iloc[i - 1, :]
        if genes_overlap(prev_gene['seqname'], prev_gene['start'], prev_gene['end'],
                         gene['seqname'], gene['start'], gene['end']):
            overlapping_indexes.add(i)

    genes_df = genes_df[~genes_df.index.isin(overlapping_indexes)] \
        .reset_index(drop=True)
    genes_df = genes_df[['seqname', 'start', 'end', 'gene_id', 'strand']]
    return genes_df


def extract_gene_id_from_gene(gff_attributes: str) -> str:
    id_attribute = gff_attributes.split(';')[0]
    return id_attribute.split(':')[1] if len(id_attribute.split(':')) == 2 else None


def genes_overlap(chrom1: str, start1: int, end1: int, chrom2: str, start2: int, end2: int) -> bool:
    if chrom1 != chrom2:
        return False

    return start1 <= start2 < end1 or start2 <= start1 < end2


def gff_to_bed_coordinates(start: int, end: int) -> (int, int):
    return start - 1, end


def get_sub_gene_annotations(genes_df: pd.DataFrame, cds_elements_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    valid_gene_indices = []
    elements_list = []
    for i, gene in genes_df.iterrows():
        gene_id = gene['gene_id']
        print(gene_id)

        # choose the longest variants that its length is dividable by 3
        var_lenths_dict = get_lengths_of_gene_variants(gene_id, cds_elements_df)
        var_lenths_dict = dict(sorted(var_lenths_dict.items(), key=operator.itemgetter(1), reverse=True))
        chosen_variant = -1
        for var_id, length in var_lenths_dict.items():
            if length % 3 == 0:  # TODO maybe ignore this check?
                chosen_variant = var_id
                break
        if chosen_variant == -1:
            print(f'Did not found an appropriate variant for gene {gene_id}. Ignoring this gene')
            continue

        valid_gene_indices.append(i)
        current_gene_elements_list = get_sub_gene_elements(gene_id, genes_df, var_id, cds_elements_df)
        elements_list.extend(current_gene_elements_list)

    valid_genes_df = genes_df[genes_df.index.isin(valid_gene_indices)]
    sub_gene_annotations_df = pd.DataFrame(elements_list)

    return valid_genes_df, sub_gene_annotations_df


def extract_gene_id_variant_id_from_cds(gff_attributes: str) -> str:
    if gff_attributes:
        if len(gff_attributes.split(';')) > 1:
            id_attribute = gff_attributes.split(';')[0]
            if len(id_attribute.split(':')) == 2:
                gene_with_variant = id_attribute.split(':')[1]
                if len(gene_with_variant.split('.')) == 2:
                    gene_variant = gene_with_variant.split('.')
                    return gene_variant[0], gene_variant[1]

    return None


def get_lengths_of_gene_variants(gene_id: str, cds_elements_df: pd.DataFrame) -> Dict[str, int]:
    var_lenth_dict = {}
    variant_ids = list(cds_elements_df[cds_elements_df['gene_id'] == gene_id]['variant_id'].unique())
    for var_id in variant_ids:
        variant_elements = cds_elements_df[(cds_elements_df['gene_id'] == gene_id) &
                                           (cds_elements_df['variant_id'] == var_id)]
        var_lenth_dict[var_id] = int(variant_elements.apply(
            lambda x: x['end'] - x['start'], axis=1, result_type='expand').sum())
    return var_lenth_dict


def get_sub_gene_elements(gene_id: str, genes_df: pd.DataFrame, var_id: int, cds_elements_df: pd.DataFrame) -> List[dict]:
    gene = genes_df[genes_df['gene_id'] == gene_id]
    gene_start = gene['start'].item()
    gene_end = gene['end'].item()
    strand = gene['strand'].item()
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
                'gene_id': gene_id,
                'region': Region.NON_CODING_GENE.value,
                'start': cursor,
                'end': cds_start,
                'strand': strand
            }
            elements_list.append(non_coding_gene_dict)

        # CDS element
        cds_dict = {
            'gene_id': gene_id,
            'region': Region.CDS.value,
            'start': cds_start,
            'end': cds_end,
            'strand': strand
        }
        elements_list.append(cds_dict)
        cursor = cds_end

    # NON_CODING_GENE element after last CDS
    if cursor != gene_end:
        non_coding_gene_dict = {
            'gene_id': gene_id,
            'region': Region.NON_CODING_GENE.value,
            'start': cursor,
            'end': gene_end,
            'strand': strand
        }
        elements_list.append(non_coding_gene_dict)

    return elements_list


def get_intergenic_annotations(chroms_df: pd.DataFrame, genes_df: pd.DataFrame) -> pd.DataFrame:
    elements_list = []

    for _, chrom in chroms_df.iterrows():
        chrom_name = chrom['seqname']
        chrom_end = chrom['end']
        cursor = 0  # chrom start at 0 always
        for _, gene in genes_df[genes_df['seqname'] == chrom_name].iterrows():

            # INTERGENIC element before current gene
            if cursor < gene['start']:
                intergenic_element = {
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
                'gene_id': '',
                'region': Region.INTERGENIC.value,
                'start': cursor,
                'end': chrom_end,
                'strand': Strand.UNKNOWN.value
            }
            elements_list.append(intergenic_element)

    intergenics_df = pd.DataFrame(elements_list)
    return intergenics_df


def test_cds_frames(in_file: str = DEFAULT_INPUT_FILE):
    if in_file is not None:
        try:
            annotations = pybedtools.example_bedtool(in_file)
            cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS')
            cds_elements = cds_elements.sort()
            cds_elements_df = cds_elements.to_dataframe()
            cds_elements_df = cds_elements_df[['seqname', 'strand', 'start', 'end', 'frame']]
            # gene_ids = cds_elements_df.seqname.unique()
            # for gene_id in gene_ids:
            #     if gene_id == 0:
            #         continue  # intergenic region
            #
            #     cds_annotations = self._annotations_df[(self._annotations_df['gene'] == gene_id) &
            #                                            (self._annotations_df['region'] == Region.CDS.value)]
            #     cds_total_len = cds_annotations['end'].sum() - cds_annotations['start'].sum()
            #
            #     if cds_total_len % 3 == 0:
            #         reading_offset = 0
            #         for idx, annotation in cds_annotations.iterrows():
            #             self._annotations_df.loc[idx, 'reading_offset'] = reading_offset
            #             reading_offset += annotation['end'] - annotation['start']
            #             reading_offset = reading_offset % 3
            #
            #     else:
            # for i, cds in cds_elements_df.iterrows():
            #     if i == 0:
            #         chrom = cds['seqname']
            #         chrom_gene_count += 1
            #         continue
            #     if chrom != cds['seqname']:
            #         print(f'Chromosome: {chrom}, overlapping occurrences: {chrom_overlap_count}')
            #         print(f'On the same strand: {chrom_overlap_count_same_strand}')
            #         print(f'Chromosome: {chrom}, all genes count: {chrom_gene_count}')
            #
            #         chrom = cds['seqname']
            #         chrom_overlap_count = 0
            #         chrom_overlap_count_same_strand = 0
            #         chrom_gene_count = 0
            #
            #     chrom_gene_count += 1
            #     prev_gene = gene_elements_df.iloc[i-1, :]
            #     if cds['start'] < prev_gene['end'] and cds['seqname'] == prev_gene['seqname']:
            #         # print('overlapping genes:')
            #         # print(gene)
            #         # print(prev_gene)
            #         overlapping_indexes.add(i)
            #         overlapping_indexes.add(i-1)
            #         overlap_count += 1
            #         chrom_overlap_count += 1
            #         if cds['strand'] == prev_gene['strand']:
            #             overlap_count_same_strand += 1
            #             chrom_overlap_count_same_strand += 1
            #

        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")


def separate_cds_genes_intergenics(in_file: str, out_file: str = 'all_chroms_annotations.csv'):
    # GENE_ID = 'gene_id'
    # def extract_gene_id(attributes: str):
    #     return next(filter(lambda x: x.startswith(GENE_ID), attributes.split(';')), None).split('=')[1]

    if in_file is not None:
        try:
            annotations = pybedtools.example_bedtool(in_file)
            all_chroms_annotations = []
            for chrom in annotations.filter(lambda x: x.fields[2] == 'chromosome'):
                cds_elements = annotations.filter(lambda x: x.fields[2] == 'CDS' and x.chrom == chrom.chrom)
                cds_elements = cds_elements.sort()
                cds_elements = cds_elements.merge()
                cds_elements_df = cds_elements.to_dataframe()
                cds_elements_df['region'] = Region.CDS.value
                del cds_elements_df['chrom']
                cds_elements_df = cds_elements_df.drop_duplicates()

                gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene' and x.chrom == chrom.chrom)
                gene_elements = gene_elements.sort()
                non_coding_gene_elements = gene_elements.subtract(cds_elements)
                non_coding_gene_elements = non_coding_gene_elements.sort()
                non_coding_gene_elements = non_coding_gene_elements.merge()
                non_coding_gene_elements_df = non_coding_gene_elements.to_dataframe()
                non_coding_gene_elements_df['region'] = Region.NON_CODING_GENE.value
                del non_coding_gene_elements_df['chrom']
                non_coding_gene_elements_df = non_coding_gene_elements_df.drop_duplicates()

                all_elements = annotations.filter(lambda x: x.chrom == chrom.chrom)
                all_elements = all_elements.sort()
                intergenic_elements = all_elements.subtract(gene_elements).subtract(cds_elements)
                intergenic_elements = intergenic_elements.sort()
                intergenic_elements = intergenic_elements.merge()
                intergenic_elements_df = intergenic_elements.to_dataframe()
                intergenic_elements_df['region'] = Region.INTERGENIC.value
                del intergenic_elements_df['chrom']
                intergenic_elements_df = intergenic_elements_df.drop_duplicates()

                cds_genes_intergenics = pd.concat([cds_elements_df, non_coding_gene_elements_df,
                                                   intergenic_elements_df], ignore_index=True)
                cds_genes_intergenics.sort_values(by=['start', 'end'], inplace=True)
                cds_genes_intergenics.insert(0, 'chrom', str(chrom.chrom))
                cds_genes_intergenics[['chrom', 'region']] = cds_genes_intergenics[['chrom', 'region']].astype(str)

                # cds_genes_intergenics[['chrom', 'region']] = cds_genes_intergenics[['chrom', 'region']].astype(str)

                gene_elements_df = gene_elements.to_dataframe()
                gene_elements_df = gene_elements_df[['strand', 'start', 'end']]
                cds_genes_intergenics[['gene', 'strand']] = cds_genes_intergenics.apply(
                    lambda x: assign_gene(x['start'], x['end'], gene_elements_df), axis=1, result_type='expand')

                # add_reading_frames_test(cds_genes_intergenics, str(chrom.chrom))

                all_chroms_annotations.append(cds_genes_intergenics)
            all_chroms_annotations = pd.concat(all_chroms_annotations).reset_index(drop=True)
            all_chroms_annotations.to_csv(out_file)
            return all_chroms_annotations
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")


def assign_gene(start: int, end: int, gene_elements_df: pd.DataFrame) -> (int, int):
    gene_ids = gene_elements_df.index[
        (gene_elements_df['start'] - 1 <= start) & (end <= gene_elements_df['end'])].tolist()
    # gene_elements_df['start'] - 1 correction is because of some bug in pybedtools

    if len(gene_ids) < 1:
        return 0, Strand.UNKNOWN.value
    if len(gene_ids) > 1:
        # We assume genes from different strands don't overlap, but if it happens nonetheless - ignore
        return 0, Strand.UNKNOWN.value
    strand = Strand(gene_elements_df.loc[gene_ids[0], 'strand'])
    return gene_ids[0] + 1, strand.value  # indices in pandas are 0-based


def workflow(input_path: str, output_path: str):
    genes_df, cds_elements_df, chroms_df = extract_data_frames(input_path)

    genes_df = get_gene_ids_without_overlaps(genes_df)
    genes_df, sub_gene_annotations_df = get_sub_gene_annotations(genes_df, cds_elements_df)
    intergenics_df = get_intergenic_annotations(chroms_df, genes_df)

    # TODO merge and reorder


def extract_data_frames(file_path: str) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    if file_path is not None:
        try:
            using_bed = file_path.lower().endswith('.bed')
            annotations = pybedtools.example_bedtool(file_path)

            # gene df
            gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene')
            gene_elements = gene_elements.sort()
            gene_elements_df = gene_elements.to_dataframe()
            gene_elements_df['gene_id'] = gene_elements_df.apply(
                lambda x: extract_gene_id_from_gene(x['attributes']), axis=1, result_type='expand')
            del gene_elements_df['attributes']
            if not using_bed:
                gene_elements_df['start'] -= 1

            # CDS df
            cds_elements = annotations.filter(lambda x: x.fields[2] == Region.CDS.value)
            cds_elements = cds_elements.sort()
            cds_elements_df = cds_elements.to_dataframe()  # TODO remove unnessesary columns: score, source
            cds_elements_df[['gene_id', 'variant_id']] = cds_elements_df.apply(
                lambda x: extract_gene_id_variant_id_from_cds(x['attributes']), axis=1, result_type='expand')
            del cds_elements_df['attributes']
            if not using_bed:
                cds_elements_df['start'] -= 1

            chrom_elements = annotations.filter(lambda x: x.fields[2] == 'chromosome')
            chrom_elements = chrom_elements.sort()
            chrom_elements_df = chrom_elements.to_dataframe()
            # TODO some preocessing in chrom_elements_df, fields and such
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if not using_bed:
                chrom_elements_df['start'] -= 1
        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")
    return gene_elements_df, cds_elements_df, chrom_elements_df


def main(raw_args=None):
    print("Testing...")
    parser = argparse.ArgumentParser(description='Genome annotations to CSV format',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-a', type=str, required=True, metavar='/path/to/annotations.gff',
                        help="Annotations file in one of the next formats: bed, gff, gtf, gff3")
    parser.add_argument('-o', type=str, required=True, metavar='/path/to/output.csv',
                        help="The path of the output csv file")

    args = parser.parse_args(raw_args)

    input_file = args.a if args.a else DEFAULT_INPUT_FILE
    output_file = args.o if args.o else DEFAULT_OUTPUT_FILE
    output_file = output_file if output_file.lower().endswith('.csv') else DEFAULT_OUTPUT_FILE

    # test_overlapping_genes(using_bed = False)
    # test_cds_frames()
    workflow(input_file, output_file)

if __name__ == "__main__":
    main()


def read_annotations_csv(file_path: str, output_dir: str = None):
    annotations_df = None
    if file_path:
        print('Loading csv...')
        try:
            _, extension = os.path.splitext(file_path)
            if extension != '.csv':
                print(f'Annotations file must be a csv file. Got extension: {extension}')
                raise Exception
            if output_dir:
                file_path = os.path.join(output_dir, file_path)
            if exists(file_path):
                annotations_df = pd.read_csv(file_path)
                # TODO first/index column - cast?
                annotations_df[['chrom', 'region']] = annotations_df[['chrom', 'region']].astype(str)
                # TODO add in the future: (should be calculated here and not in annotated_seq)
                # annotations_df[['strand', 'gene']] = annotations_df[['strand', 'gene']].astype(str)
                annotations_df[['start', 'end']] = annotations_df[['start', 'end']].astype(int)
            else:
                print(f'Annotations file does not exist. File path = {file_path}')
                raise Exception
        except Exception:
            print('Problem parsing annotations file')
    return annotations_df
