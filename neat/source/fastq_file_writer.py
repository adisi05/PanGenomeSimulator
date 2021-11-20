import Bio.bgzf as bgzf
import pathlib

from neat.source.file_writer_utils import BUFFER_BATCH_SIZE


class FastqFileWriter:
    def __init__(self, out_prefix, paired=False, no_fastq=False):
        self.no_fastq = no_fastq
        if self.no_fastq:
            return

        fq1 = pathlib.Path(out_prefix + '_read1.fq.gz')
        fq2 = pathlib.Path(out_prefix + '_read2.fq.gz')
        self._file1 = bgzf.open(fq1, 'w')
        self._file2 = None
        if paired:
            self._file2 = bgzf.open(fq2, 'w')
        self._buffer1 = []
        self._buffer2 = []

    def write_record(self, read_name, read1, qual1, read2=None, qual2=None, orientation=None):
        if self.no_fastq:
            return

        # Since read1 and read2 are Seq objects from Biopython, they have reverse_complement methods built-in
        (read1, quality1) = (read1, qual1)
        if read2 is not None and orientation is True:
            (read2, quality2) = (read2.reverse_complement(), qual2[::-1])
        elif read2 is not None and orientation is False:
            read2_tmp = read2
            qual2_tmp = qual2
            (read2, quality2) = (read1, qual1)
            (read1, quality1) = (read2_tmp.reverse_complement(), qual2_tmp[::-1])

        self._buffer1.append('@' + read_name + '/1\n' + str(read1) + '\n+\n' + quality1 + '\n')
        if read2 is not None:
            self._buffer2.append('@' + read_name + '/2\n' + str(read2) + '\n+\n' + quality2 + '\n')

    def flush_buffer(self, last_time=False):
        if self.no_fastq:
            return

        if (len(self._buffer1) >= BUFFER_BATCH_SIZE) or (len(self._buffer1) and last_time):
            if not self.no_fastq:
                self._file1.write(''.join(self._buffer1))
                if len(self._buffer2):
                    self._file2.write(''.join(self._buffer2))
            self._buffer1 = []
            self._buffer2 = []

    def close_file(self):
        if self.no_fastq:
            return

        self.flush_buffer(last_time=True)
        self._file1.close()
        if self._file2 is not None:
            self._file2.close()