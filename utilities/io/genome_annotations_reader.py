import os
from os.path import exists
from typing import Dict, Tuple

import pandas as pd

from utilities.io.logger import Logger


def read_annotations_csv(file_path: str, logger: Logger = None):
    logger = logger if logger else Logger()
    annotations_df = None
    if file_path:
        logger.message('Loading annotations CSV file...')
        try:
            _, extension = os.path.splitext(file_path)
            if extension != '.csv':
                logger.message(f'Annotations file must be a csv file. Got extension: {extension}')
                raise Exception
            if exists(file_path):
                annotations_df = pd.read_csv(file_path, index_col=0,
                                             dtype={'chrom': 'str', 'gene_id': 'int', 'region': 'str',
                                                    'start': 'int', 'end': 'int', 'strand': 'str', 'frame': 'Int8'})
            else:
                logger.message(f'Annotations file does not exist. File path = {file_path}')
                raise Exception
        except Exception as e:
            logger.message(f'Problem parsing annotations file. Error details: {e}')
    return annotations_df


def get_all_genes(file_path: str) -> Dict[str, Tuple[str, int, int]]:
    relevant_genes = {}
    annotations_df = read_annotations_csv(file_path)
    if annotations_df is None or annotations_df.empty:
        return relevant_genes
    for gene_id in annotations_df.gene_id.unique():
        if gene_id == 0:
            continue
        gene_annotations = annotations_df[annotations_df['gene_id'] == gene_id]
        chrom = gene_annotations.iloc[0]['chrom']
        start = gene_annotations.iloc[0]['start']
        end = gene_annotations.iloc[-1]['end']
        relevant_genes[gene_id] = (chrom, start, end)
    return relevant_genes
