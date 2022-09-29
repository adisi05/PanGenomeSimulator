import os
from os.path import exists

import pandas as pd
import pybedtools

from utilities.common_data_structues import Strand, Region

DEFAULT_FILE_PATH = '/groups/itay_mayrose/adisivan/arabidopsis/ensemblgenomes/gff3/Arabidopsis_thaliana.TAIR10.53.gff3'

def test_overlapping_genes(in_file: str = DEFAULT_FILE_PATH):
    if in_file is not None:
        try:
            annotations = pybedtools.example_bedtool(in_file)
            gene_elements = annotations.filter(lambda x: x.fields[2] == 'gene')
            gene_elements = gene_elements.sort()
            gene_elements_df = gene_elements.to_dataframe()
            gene_elements_df = gene_elements_df[['seqname', 'strand', 'start', 'end']]
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

                    chrom = gene['seqname']
                    chrom_overlap_count = 0
                    chrom_overlap_count_same_strand = 0
                    chrom_gene_count = 0

                chrom_gene_count += 1
                prev_gene = gene_elements_df.iloc[i-1, :]
                if gene['start'] < prev_gene['end'] and gene['seqname'] == prev_gene['seqname']:
                    # print('overlapping genes:')
                    # print(gene)
                    # print(prev_gene)
                    overlapping_indexes.add(i)
                    overlapping_indexes.add(i-1)
                    overlap_count += 1
                    chrom_overlap_count += 1
                    if gene['strand'] == prev_gene['strand']:
                        overlap_count_same_strand += 1
                        chrom_overlap_count_same_strand += 1

            print(f'Chromosome: {chrom}, overlapping occurrences: {chrom_overlap_count}')
            print(f'On the same strand: {chrom_overlap_count_same_strand}')
            print(f'Chromosome: {chrom}, all genes count: {chrom_gene_count}')
            print(f'Overall overlapping occurrences: {overlap_count}')
            print(f'On the same strand: {overlap_count_same_strand}')
            print(f'All genes count: {len(gene_elements_df)}')
            gene_elements_df = gene_elements_df.iloc[list(overlapping_indexes), :].reset_index(drop=True)
            gene_elements_df.to_csv('overlapping_genes.csv')


        except IOError:
            print("\nProblem reading annotation (BED/GFF) file.\n")


def test_cds_frames(in_file: str = DEFAULT_FILE_PATH):
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

if __name__ == "__main__":
    print("Testing...")
    # test_overlapping_genes()
    test_cds_frames()


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