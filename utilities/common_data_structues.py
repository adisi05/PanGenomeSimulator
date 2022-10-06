from enum import Enum


class Strand(Enum):
    FORWARD = '+'
    REVERSE = '-'
    UNKNOWN = '.'


class MutType(Enum):
    SNP = 'snp'
    INDEL = 'indel'
    SV = 'SV'


class Region(Enum):
    CDS = 'CDS'
    NON_CODING_GENE = 'non_coding_gene'
    INTERGENIC = 'intergenic'
    ALL = 'all'
    # TODO see how we pass through evey Region which is not ALL, or only through ALL, when inserting mutations.
    # need some kind of "strategy" solution


VALID_NUCL = ['A', 'C', 'G', 'T']
VALID_TRINUC = [VALID_NUCL[i] + VALID_NUCL[j] + VALID_NUCL[k] for i in range(len(VALID_NUCL)) for j in
                range(len(VALID_NUCL)) for k in range(len(VALID_NUCL))]

STOP_CODONS_FORWARD_STRAND = ['TAG', 'TAA', 'TGA']
STOP_CODONS_REVERSE_STRAND = ['CTA', 'TTA', 'TCA']

NUCL = ['A', 'C', 'G', 'T']
TRI_IND = {'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3, 'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
           'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11, 'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15}
NUC_IND = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
ALL_TRI = [NUCL[i] + NUCL[j] + NUCL[k] for i in range(len(NUCL)) for j in range(len(NUCL)) for k in range(len(NUCL))]
ALL_IND = {ALL_TRI[i]: i for i in range(len(ALL_TRI))}