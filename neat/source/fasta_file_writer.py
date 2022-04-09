import os


class FastaFileWriter:

    @staticmethod
    def get_output_filenames(prefix, name):
        return [prefix + "_" + name + ".fasta"]

    def __init__(self, out_prefix, ploidy, line_width):

        self.files = {}
        self.prev_chrom = None
        self.ploidy = ploidy
        self.line_width = line_width
        self.last_line_len = 0
        for hapl_idx in range(self.ploidy):
            file_name = '{0}.fasta_temp'.format(out_prefix) if ploidy == 1 else '{0}_{hapl_idx}.fasta_temp'.format(
                out_prefix, hapl_idx)
            print("TEST file_name",file_name)
            open(file_name, 'w').close()
            self.files[hapl_idx] = file_name

    def write_record(self, haploid_sequences, chrom, N_seq_len=0, ignored_ending=0):
        for hapl_idx in range(self.ploidy):
            fasta_file = open(self.files[hapl_idx], 'a')

            if (not self.prev_chrom) or self.prev_chrom != chrom:
                if self.prev_chrom != None:
                    fasta_file.write("\n")
                    self.last_line_len = 0
                strID = '>{0}\n'.format(chrom)
                fasta_file.write(strID)
                self.prev_chrom = chrom

            if haploid_sequences is None or N_seq_len > 0:
                strDNA = "N" * N_seq_len
                # print("TEST2 - N_seq_len is ",N_seq_len)
            # elif not haploid_sequences:
                # print("TEST2 - haploid_sequences is None, N_seq_len is ",N_seq_len)
            else:
                strDNA = str(haploid_sequences[hapl_idx])
                # print("TEST2 - strDNA start with", strDNA[0:9])
            strDNA_len = len(strDNA)-ignored_ending
            strDNA=strDNA[:strDNA_len]
            i = 0

            if self.last_line_len + strDNA_len < self.line_width:
                pass
                self.last_line_len = self.last_line_len + strDNA_len
            else: # self.last_line_len + strDNA_len >= self.line_width
                fasta_file.write(strDNA[0:self.line_width-self.last_line_len] + '\n')
                i = self.line_width - self.last_line_len
                self.last_line_len = 0
                while i + self.line_width <= strDNA_len:
                    fasta_file.write(strDNA[i:i + self.line_width] + '\n')
                    i = i + self.line_width
                if i < strDNA_len:
                    fasta_file.write(strDNA[i:strDNA_len])
                    self.last_line_len = strDNA_len - i
            fasta_file.close()

    def flush_buffer(self, last_time=False):
        pass

    def close_file(self):
        for hapl_idx in range(self.ploidy):
            final_name = self.files[hapl_idx].removesuffix('_temp')
            os.rename(self.files[hapl_idx], final_name)
            self.files[hapl_idx] = final_name

    def get_file(self):
        return [self.files[hapl_idx] for hapl_idx in range(self.ploidy)]