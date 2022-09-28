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