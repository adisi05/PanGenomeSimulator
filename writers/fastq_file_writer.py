import time
import os
from typing import List


class FastqFileWriter:

    @staticmethod
    def get_output_filenames(prefix, name) -> List[str]:
        return [prefix + "_" + name + "_read.fq",
                prefix + "_" + name + "_read.aln",
                prefix + "_" + name + "_read1.fq",
                prefix + "_" + name + "_read2.fq",
                prefix + "_" + name + "_read1.aln",
                prefix + "_" + name + "_read2.aln"]

    @staticmethod
    def generate_reads(fasta_files: List[str], paired_end: bool, read_length: int, coverage: float, insert_size: int,
                       insert_std: int):
        print(fasta_files)
        # print(sequencing_params)
        paired = "-p" if paired_end else ""
        fastq_files = [filename[:-len('.fasta')] + "_read" for filename in fasta_files]
        for fasta, fastq in zip(fasta_files, fastq_files):
            art_command = "ART/art_bin_MountRainier/art_illumina {} -i {} -l {} -f {} -o {} -m {} -s {}"\
                .format(paired, fasta, read_length, coverage, fastq, insert_size, insert_std)
            start = time.time()
            os.system(art_command)
            end = time.time()
            print("ART reads simulation took {} seconds.".format(int(end - start)))
