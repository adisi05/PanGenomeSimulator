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


class ModelStats(Enum):
    # how many times do we observe each trinucleotide in the reference (or genomic region)?
    TRINUC_REF_COUNT = 'TRINUC_REF_COUNT'
    # [(trinuc_ref, trinuc_alt)] = # of times we observed a mutation from trinuc_ref into trinuc_alt
    TRINUC_TRANSITION_COUNT = 'TRINUC_TRANSITION_COUNT'
    # total count of SNPs
    SNP_COUNT = 'SNP_COUNT'
    # overall SNP transition probabilities
    SNP_TRANSITION_COUNT = 'SNP_TRANSITION_COUNT'
    # total count of indels, indexed by length
    INDEL_COUNT = 'INDEL_COUNT'
    # tabulate how much non-N reference sequence we've eaten through
    TOTAL_REFLEN = 'TOTAL_REFLEN'
    # the number of SNPs divided to the overall number of mutations
    SNP_FREQ = 'SNP_FREQ'
    # the number of indels divided to the overall number of mutations (1 - SNP_FREQ)
    AVG_INDEL_FREQ = 'AVG_INDEL_FREQ'
    # number of a specific indel divided to the number of the overall indels
    INDEL_FREQ = 'INDEL_FREQ'
    # the number of mutations divided to the overall length
    AVG_MUT_RATE = 'AVG_MUT_RATE'
    # frequency of snp transitions, given a snp occurs.
    SNP_TRANS_FREQ = 'SNP_TRANS_FREQ'
    # frequency that each trinuc mutated into anything else
    TRINUC_MUT_PROB = 'TRINUC_MUT_PROB'
    # frequency that a trinuc mutates into another trinuc, given that it mutated
    TRINUC_TRANS_PROBS = 'TRINUC_TRANS_PROBS'


class ModelKeys(Enum):
    # average mutation rate
    AVG_MUT_RATE = 'AVG_MUT_RATE'
    # p(mut is indel | mut occurs)
    P_INDEL = 'P_INDEL'
    # p(insertion | indel occurs)
    P_INSERTION = 'P_INSERTION'
    # distribution of insertion lengths
    INS_LEN_DSTRBTN = 'INS_LEN_DSTRBTN'
    # distribution of deletion lengths
    DEL_LEN_DSTRBTN = 'DEL_LEN_DSTRBTN'
    # distribution of trinucleotide SNP transitions
    TRINUC_TRANS_DSTRBTN = 'TRINUC_TRANS_DSTRBTN'
    # p(trinuc mutates)
    P_TRINUC_MUT = 'P_TRINUC_MUT'
