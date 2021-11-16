from struct import pack
import Bio.bgzf as bgzf
import pathlib
import re

from neat.source.neat_cigar import CigarString
from neat.source.output_file_writer import OutputFileWriter


class FastaFileWriter(OutputFileWriter):
    def __init__(self):
        pass

    def write_record(self, sequences, chrom, out_prefix):

        for hapl_idx in range(sequences.ploidy):
            file_name = '{0}.fasta'.format(out_prefix) if sequences.ploidy == 1 else '{0}_{hapl_idx}.fasta'.format(out_prefix,hapl_idx)
            self.fasta_file = open(file_name, 'w')
            strID = '>{0}\n'.format(chrom)
            self.fasta_file.write(strID)

            strDNA = str(sequences.sequences[hapl_idx])
            strDNA_len = len(strDNA)
            i = 0

            while i + 50 <= strDNA_len:
                self.fasta_file.write(strDNA[i:i + 50] + '\n')
                i = i + 50

            if i < strDNA_len:
                self.fasta_file.write(strDNA[i:strDNA_len] + '\n')
            self.fasta_file.close()


    def flush_buffer(self, last_time=False):
        pass

    def close_file(self):
        pass
