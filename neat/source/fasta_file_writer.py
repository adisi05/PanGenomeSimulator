class FastaFileWriter:
    def __init__(self, out_prefix, ploidy):
        self.files = {}
        for hapl_idx in range(ploidy):
            file_name = '{0}.fasta'.format(out_prefix) if ploidy == 1 else '{0}_{hapl_idx}.fasta'.format(
                out_prefix, hapl_idx)
            open(file_name, 'w').close()
            self.files[hapl_idx] = file_name

    def write_record(self, sequences, chrom):

        for hapl_idx in range(sequences.ploidy):
            fasta_file = open(self.files[hapl_idx], 'a')
            strID = '>{0}\n'.format(chrom)
            fasta_file.write(strID)

            strDNA = str(sequences.sequences[hapl_idx])
            strDNA_len = len(strDNA)
            i = 0

            while i + 50 <= strDNA_len:
                fasta_file.write(strDNA[i:i + 50] + '\n')
                i = i + 50

            if i < strDNA_len:
                fasta_file.write(strDNA[i:strDNA_len] + '\n')
            fasta_file.close()


    def flush_buffer(self, last_time=False):
        pass

    def close_file(self):
        pass
