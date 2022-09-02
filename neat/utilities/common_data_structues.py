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
    #TODO see how we pass through evey Region which is not ALL, or only through ALL, when inserting mutations.
    # need some kind of "strategy" solution
