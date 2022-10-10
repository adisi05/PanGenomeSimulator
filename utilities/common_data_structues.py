from dataclasses import dataclass
from enum import Enum


VALID_NUCL = ['A', 'C', 'G', 'T']
VALID_TRINUC = [VALID_NUCL[i] + VALID_NUCL[j] + VALID_NUCL[k] for i in range(len(VALID_NUCL)) for j in
                range(len(VALID_NUCL)) for k in range(len(VALID_NUCL))]

STOP_CODONS_FORWARD_STRAND = ['TAG', 'TAA', 'TGA']
STOP_CODONS_REVERSE_STRAND = ['CTA', 'TTA', 'TCA']

TRI_IND = {'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3, 'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
           'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11, 'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15}
NUC_IND = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
ALL_TRI = [VALID_NUCL[i] + VALID_NUCL[j] + VALID_NUCL[k] for i in range(len(VALID_NUCL)) for j in range(len(VALID_NUCL))
           for k in range(len(VALID_NUCL))]
ALL_IND = {ALL_TRI[i]: i for i in range(len(ALL_TRI))}


class MutType(Enum):
    SNP = 'snp'
    INDEL = 'indel'
    SV = 'SV'


# see https://stackoverflow.com/questions/68840767/what-is-the-purpose-of-sort-index-of-a-dataclass-in-python
@dataclass(order=True)
@dataclass
class Mutation:
    position: int
    ref_nucl: str
    new_nucl: str
    mut_type: MutType

    def get_offset_change(self):
        return len(self.new_nucl) - len(self.ref_nucl)


class Region(Enum):
    CDS = 'CDS'
    NON_CODING_GENE = 'non_coding_gene'
    INTERGENIC = 'intergenic'
    ALL = 'all'
    # TODO see how we pass through evey Region which is not ALL, or only through ALL, when inserting mutations.
    # need some kind of "strategy" solution


class Strand(Enum):
    FORWARD = '+'
    REVERSE = '-'
    UNKNOWN = '.'
