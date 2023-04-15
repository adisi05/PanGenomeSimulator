import sys
import time
import gzip
import pathlib
import random
from Bio.Seq import Seq, MutableSeq

from utilities.common_data_structues import VALID_NUCL
from utilities.io.logger import Logger

OK_CHR_ORD = {'A': True, 'C': True, 'G': True, 'T': True, 'U': True}


def index_ref(reference_path: str, logger: Logger = None) -> [list, int]:
    """
    Index reference fasta
    :param logger:
    :param reference_path: string path to the reference
    :return: reference index in list from
    """
    logger = logger if logger else Logger()
    tt = time.time()

    absolute_reference_location = pathlib.Path(reference_path)

    # sanity check
    if not absolute_reference_location.is_file():
        logger.message("\nProblem reading the reference fasta file.\n")
        sys.exit(1)

    index_filename = None

    # check if the reference file already exists
    if absolute_reference_location.with_suffix('.fai').is_file():
        logger.message('found index ' + str(absolute_reference_location.with_suffix('.fai')))
        index_filename = absolute_reference_location.with_suffix('.fai')
    elif absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai').is_file():
        logger.message('found index ' +
                       str(absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai')))
        index_filename = absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai')
    else:
        pass

    ref_indices = []
    if index_filename is not None:
        fai = open(index_filename, 'r')
        for line in fai:
            splt = line[:-1].split('\t')
            # Defined as the number of bases in the contig
            seq_len = int(splt[1])
            # Defined as the byte index where the contig sequence begins
            offset = int(splt[2])
            # Defined as bases per line in the Fasta file
            line_ln = int(splt[3])
            n_lines = seq_len // line_ln
            if seq_len % line_ln != 0:
                n_lines += 1
            # Item 3 in this gives you the byte position of the next contig, I believe
            ref_indices.append((splt[0], offset, offset + seq_len + n_lines, seq_len))
        fai.close()
        return ref_indices, line_ln

    logger.message('Index not found, creating one... ')
    if absolute_reference_location.suffix == ".gz":
        ref_file = gzip.open(absolute_reference_location, 'rt')
    else:
        ref_file = open(absolute_reference_location, 'r')
    prev_r = None
    prev_p = None
    seq_len = 0

    line_width = -1
    while True:
        data = ref_file.readline()
        if not data:
            ref_indices.append((prev_r, prev_p, ref_file.tell() - len(data), seq_len))
            break
        elif data[0] == '>':
            if prev_p is not None:
                ref_indices.append((prev_r, prev_p, ref_file.tell() - len(data), seq_len))
            seq_len = 0
            prev_p = ref_file.tell()  # including '\n' characters
            prev_r = data[1:-1]
        else:
            # look at the first line of the first chromosome
            line_width = len(data) - 1 if line_width == -1 else line_width
            seq_len += len(data) - 1  # -1 is for ignoring the '\n' characters
    ref_file.close()

    logger.message('{0:.3f} (sec)'.format(time.time() - tt))
    return ref_indices, line_width


def read_ref(ref_path, ref_inds_i, n_handling, n_unknowns=True, logger: Logger = None):
    logger = logger if logger else Logger()
    tt = time.time()
    logger.debug_message(f'Reading & indexing reference, chromosome: {ref_inds_i[0]}... ')

    absolute_reference_path = pathlib.Path(ref_path)
    if absolute_reference_path.suffix == '.gz':
        ref_file = gzip.open(absolute_reference_path, 'rt')
    else:
        ref_file = open(absolute_reference_path, 'r')

    # TODO convert to SeqIO containers?
    #  for seq_record in SeqIO.parse(ref_file, "fasta"):
    #     ...

    ref_file.seek(ref_inds_i[1])
    my_dat = ''.join(ref_file.read(ref_inds_i[2] - ref_inds_i[1]).split('\n'))
    my_dat = Seq(my_dat.upper())
    # TODO Mutable seqs have a number of disadvantages. Convert to immutable and see if that helps?
    #  my_dat = MutableSeq(my_dat)

    # find N regions
    # data explanation: my_dat[n_atlas[0][0]:n_atlas[0][1]] = solid block of Ns
    prev_ni = 0
    n_count = 0
    n_atlas = []
    for i in range(len(my_dat)):
        # if my_dat[i] == '>':
        #     break
        if my_dat[i] == 'N' or (n_unknowns and my_dat[i] not in OK_CHR_ORD):
            if n_count == 0:
                prev_ni = i
            n_count += 1
            if i == len(my_dat) - 1:
                n_atlas.append((prev_ni, prev_ni + n_count))
        else:
            if n_count > 0:
                n_atlas.append((prev_ni, prev_ni + n_count))
            n_count = 0

    # handle N base-calls as desired
    # TODO this seems to randomly replace an N with a base. Is this necessary? How to do this in an immutable seq?
    n_info = {'all': [], 'big': [], 'non_N': []}
    if n_handling['method'] == 'random':
        for n_region in n_atlas:
            n_info['all'].extend(n_region)
            if n_region[1] - n_region[0] <= n_handling['max_threshold']:
                for i in range(n_region[0], n_region[1]):
                    temp = MutableSeq(my_dat)
                    temp[i] = random.choice(VALID_NUCL)
                    my_dat = Seq(temp)
            else:
                n_info['big'].extend(n_region)
    elif n_handling['method'] == 'ignore':
        for n_region in n_atlas:
            n_info['all'].extend(n_region)
            n_info['big'].extend(n_region)
    else:
        logger.message('\nERROR: UNKNOWN N_HANDLING MODE\n')
        sys.exit(1)

    habitable_regions = []
    if not n_info['big']:
        n_info['non_N'] = [(0, len(my_dat))]
    else:
        for i in range(0, len(n_info['big']), 2):
            if i == 0:
                habitable_regions.append((0, n_info['big'][0]))
            else:
                habitable_regions.append((n_info['big'][i - 1], n_info['big'][i]))
        habitable_regions.append((n_info['big'][-1], len(my_dat)))
    for n in habitable_regions:
        if n[0] != n[1]:
            n_info['non_N'].append(n)

    ref_file.close()

    logger.debug_message('Reading reference took {0:.3f} (sec)'.format(time.time() - tt))

    return my_dat, n_info
