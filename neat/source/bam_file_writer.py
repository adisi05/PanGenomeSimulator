from struct import pack
import Bio.bgzf as bgzf
import pathlib
import re

from neat.source.neat_cigar import CigarString
from neat.source.file_writer_utils import BUFFER_BATCH_SIZE

BAM_COMPRESSION_LEVEL = 6

# SAMtools reg2bin function
def reg2bin(beg: int, end: int):
    """
    Finds the largest superset bin of region. Numeric values taken from hts-specs
    Note: description of this function taken from source code for bamnostic.bai
        (https://bamnostic.readthedocs.io/en/latest/_modules/bamnostic/bai.html)
    :param beg: inclusive beginning position of region
    :param end: exclusive end position of region
    :return: distinct bin ID or largest superset bin of region
    """
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


# takes list of strings, returns numerical flag
def sam_flag(string_list: list) -> int:
    out_val = 0
    string_list = list(set(string_list))
    for n in string_list:
        if n == 'paired':
            out_val += 1
        elif n == 'proper':
            out_val += 2
        elif n == 'unmapped':
            out_val += 4
        elif n == 'mate_unmapped':
            out_val += 8
        elif n == 'reverse':
            out_val += 16
        elif n == 'mate_reverse':
            out_val += 32
        elif n == 'first':
            out_val += 64
        elif n == 'second':
            out_val += 128
        elif n == 'not_primary':
            out_val += 256
        elif n == 'low_quality':
            out_val += 512
        elif n == 'duplicate':
            out_val += 1024
        elif n == 'supplementary':
            out_val += 2048
    return out_val


CIGAR_PACKED = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
SEQ_PACKED = {'=': 0, 'A': 1, 'C': 2, 'M': 3, 'G': 4, 'R': 5, 'S': 6, 'V': 7,
              'T': 8, 'W': 9, 'Y': 10, 'H': 11, 'K': 12, 'D': 13, 'B': 14, 'N': 15}


# TODO find a better way to write output files
class BamFileWriter:
    def __init__(self, out_prefix, header=None):

        bam = pathlib.Path(out_prefix + '_golden.bam')

        # BAM OUTPUT
        self._file = None
        if header is not None:
            self._file = bgzf.BgzfWriter(bam, 'w', compresslevel=BAM_COMPRESSION_LEVEL)

            # WRITE BAM HEADER
            self._file.write("BAM\1")
            header = '@HD\tVN:1.5\tSO:coordinate\n'
            for n in header[0]:
                header += '@SQ\tSN:' + n[0] + '\tLN:' + str(n[3]) + '\n'
            header += '@RG\tID:NEAT\tSM:NEAT\tLB:NEAT\tPL:NEAT\n'
            header_bytes = len(header)
            num_refs = len(header[0])
            self._file.write(pack('<i', header_bytes))
            self._file.write(header)
            self._file.write(pack('<i', num_refs))

            for n in header[0]:
                l_name = len(n[0]) + 1
                self._file.write(pack('<i', l_name))
                self._file.write(n[0] + '\0')
                self._file.write(pack('<i', n[3]))

        # buffers for more efficient writing
        self._buffer = []

    def write_record(self, ref_id, read_name, pos_0, cigar, seq, qual, output_sam_flag,
                         mate_pos=None, aln_map_quality=70):

        my_bin = reg2bin(pos_0, pos_0 + len(seq))
        # my_bin     = 0	# or just use a dummy value, does this actually matter?

        my_map_quality = aln_map_quality
        cigar_string = CigarString.list_to_string(cigar)
        cig_letters = re.split(r"\d+", cigar_string)[1:]
        cig_numbers = [int(n) for n in re.findall(r"\d+", cigar_string)]
        cig_ops = len(cig_letters)
        next_ref_id = ref_id
        if mate_pos is None:
            next_pos = 0
            my_t_len = 0
        else:
            next_pos = mate_pos
            if next_pos > pos_0:
                my_t_len = next_pos - pos_0 + len(seq)
            else:
                my_t_len = next_pos - pos_0 - len(seq)

        encoded_cig = bytearray()
        for i in range(cig_ops):
            encoded_cig.extend(pack('<I', (cig_numbers[i] << 4) + CIGAR_PACKED[cig_letters[i]]))
        encoded_seq = bytearray()
        encoded_len = (len(seq) + 1) // 2
        seq_len = len(seq)
        if seq_len & 1:
            seq += '='
        for i in range(encoded_len):
            # print(seq[2*i], seq[2*i+1])
            encoded_seq.extend(
                pack('<B', (SEQ_PACKED[seq[2 * i].capitalize()] << 4) + SEQ_PACKED[seq[2 * i + 1].capitalize()]))

        # apparently samtools automatically adds 33 to the quality score string...
        encoded_qual = ''.join([chr(ord(n) - 33) for n in qual])

        """
        block_size = 4 +		# refID 		int32
                     4 +		# pos			int32
                     4 +		# bin_mq_nl		uint32
                     4 +		# flag_nc		uint32
                     4 +		# l_seq			int32
                     4 +		# next_ref_id	int32
                     4 +		# next_pos		int32
                     4 +		# tlen			int32
                     len(readName)+1 +
                     4*cig_ops +
                     encoded_len +
                     len(seq)
        """
        # block_size = 32 + len(readName)+1 + 4*cig_ops + encoded_len + len(seq)
        block_size = 32 + len(read_name) + 1 + len(encoded_cig) + len(encoded_seq) + len(encoded_qual)

        """
        Not sure what the point of the following lines are
        # self.bam_file.write(pack('<i',block_size))
        # self.bam_file.write(pack('<i',refID))
        # self.bam_file.write(pack('<i',pos_0))
        # self.bam_file.write(pack('<I',(my_bin<<16) + (my_map_quality<<8) + len(readName)+1))
        # self.bam_file.write(pack('<I',(samFlag<<16) + cig_ops))
        # self.bam_file.write(pack('<i',seq_len))
        # self.bam_file.write(pack('<i',next_ref_id))
        # self.bam_file.write(pack('<i',next_pos))
        # self.bam_file.write(pack('<i',my_tlen))
        # self.bam_file.write(readName+'\0')
        # self.bam_file.write(encoded_cig)
        # self.bam_file.write(encoded_seq)
        # self.bam_file.write(encoded_qual)
        """

        # a horribly compressed line, I'm sorry.
        # (ref_index, position, data)
        self._buffer.append((ref_id, pos_0, pack('<i', block_size) + pack('<i', ref_id) + pack('<i', pos_0) +
                             pack('<I', (my_bin << 16) + (my_map_quality << 8) + len(read_name) + 1) +
                             pack('<I', (output_sam_flag << 16) + cig_ops) + pack('<i', seq_len) + pack('<i', next_ref_id) +
                             pack('<i', next_pos) + pack('<i', my_t_len) + read_name.encode('utf-8') +
                                b'\0' + encoded_cig + encoded_seq + encoded_qual.encode('utf-8')))

    def flush_buffer(self, bam_max=None, last_time=False):
        if (len(self._buffer) >= BUFFER_BATCH_SIZE) or (len(self._buffer) and last_time):
            # bam
            bam_data = sorted(self._buffer)
            if last_time:
                self._file.write(b''.join([n[2] for n in bam_data]))
                self._buffer = []
            else:
                ind_to_stop_at = 0
                for i in range(0, len(bam_data)):
                    # if we are from previous reference, or have coordinates lower
                    # than next window position, it's safe to write out to file
                    if bam_data[i][0] != bam_data[-1][0] or bam_data[i][1] < bam_max:
                        ind_to_stop_at = i + 1
                    else:
                        break
                self._file.write(b''.join([n[2] for n in bam_data[:ind_to_stop_at]]))
                # Debug statement
                # print(f'BAM WRITING: {ind_to_stop_at}/{len(bam_data)}')
                if ind_to_stop_at >= len(bam_data):
                    self._buffer = []
                else:
                    self._buffer = bam_data[ind_to_stop_at:]

    def close_file(self):
        self.flush_buffer(last_time=True)
        if self._file is not None:
            self._file.close()
