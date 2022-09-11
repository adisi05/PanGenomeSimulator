import os


class FastaFileWriter:

    @staticmethod
    def get_output_filenames(prefix, name):
        return [prefix + "_" + name + ".fasta"]

    def __init__(self, out_prefix, line_width):

        self.prev_chrom = None
        self.line_width = line_width
        self.last_line_len = 0
        self.file_name = '{0}.fasta_temp'.format(out_prefix)
        print("TEST file_name",self.file_name)
        open(self.file_name, 'w').close()

    def write_record(self, sequence, chrom, N_seq_len=0, ignored_ending=0):
        fasta_file = open(self.file_name, 'a')

        if chrom != self.prev_chrom:
            # start new chromosome
            if self.prev_chrom != None:
                fasta_file.write("\n")
                self.last_line_len = 0
            str_id = '>{0}\n'.format(chrom)
            fasta_file.write(str_id)
            self.prev_chrom = chrom

        # determine - nucleotides or 'N' sequence
        if sequence is not None and N_seq_len == 0:
            dna_str = str(sequence)
        else:
            dna_str = "N" * N_seq_len

        dna_str_len = len(dna_str) - ignored_ending
        dna_str = dna_str[:dna_str_len]

        idx = 0
        while idx < dna_str_len:
            if self.last_line_len == self.line_width:
                fasta_file.write('\n')
                self.last_line_len = 0
            chunk_len = min(self.line_width - self.last_line_len, dna_str_len - idx)
            fasta_file.write(dna_str[idx:idx+chunk_len])
            idx += chunk_len
            self.last_line_len += chunk_len
        fasta_file.close()

    def finalize(self):
        final_name = self.file_name[:-len('_temp')]
        os.rename(self.file_name, final_name)
        self.file_name = final_name

    def get_file_name(self):
        return self.file_name