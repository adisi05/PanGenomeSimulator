import Bio.bgzf as bgzf
import pathlib
import vcf
from vcf import utils


class VcfFileWriter:

    @staticmethod
    def get_output_filenames(prefix, name):
        return [prefix + "_" + name + "_golden.vcf.gz",
                prefix + "_" + name + "_golden_final.vcf.gz"]

    def __init__(self, out_prefix, parent_prefix, accession, header=None):
        self._out_prefix = out_prefix
        self._parent_prefix = parent_prefix

        path = self._out_prefix + '_golden.vcf.gz' if self._parent_prefix else self._out_prefix + '_golden_final.vcf.gz'
        vcf_file = pathlib.Path(path)

        # VCF OUTPUT
        self._file = None
        if header is not None:
            self._file = bgzf.open(vcf_file, 'wb')

            # WRITE VCF HEADER
            self._file.write('##fileformat=VCFv4.1\n'.encode('utf-8'))
            reference = '##reference=' + header[0] + '\n'
            self._file.write(reference.encode('utf-8'))
            self._file.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'.encode('utf-8'))
            self._file.write(
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'.encode('utf-8'))
            self._file.write(
                '##INFO=<ID=VMX,Number=1,Type=String,Description="SNP is Missense in these Read Frames">\n'.encode(
                    'utf-8'))
            self._file.write(
                '##INFO=<ID=VNX,Number=1,Type=String,Description="SNP is Nonsense in these Read Frames">\n'.encode(
                    'utf-8'))
            self._file.write(
                '##INFO=<ID=VFX,Number=1,Type=String,Description="Indel Causes Frameshift">\n'.encode('utf-8'))
            self._file.write(
                '##INFO=<ID=WP,Number=A,Type=Integer,Description="NEAT-GenReads ploidy indicator">\n'.encode(
                    'utf-8'))
            self._file.write('##ALT=<ID=DEL,Description="Deletion">\n'.encode('utf-8'))
            self._file.write('##ALT=<ID=DUP,Description="Duplication">\n'.encode('utf-8'))
            self._file.write('##ALT=<ID=INS,Description="Insertion of novel sequence">\n'.encode('utf-8'))
            self._file.write('##ALT=<ID=INV,Description="Inversion">\n'.encode('utf-8'))
            self._file.write('##ALT=<ID=CNV,Description="Copy number variable region">\n'.encode('utf-8'))
            self._file.write('##ALT=<ID=TRANS,Description="Translocation">\n'.encode('utf-8'))
            self._file.write('##ALT=<ID=INV-TRANS,Description="Inverted translocation">\n'.encode('utf-8'))
            self._file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'.encode('utf-8'))
            self._file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'
                             .format(accession).encode('utf-8'))

    def write_record(self, chrom, pos, id_str, ref, alt, qual, filt, info):
        self._file.write(
            str(chrom) + '\t' + str(pos) + '\t' + str(id_str) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(
                qual) + '\t' + str(filt) + '\t' + str(info) + '\tGT\t1|1\n')

    def flush_buffer(self, last_time=False):
        pass

    def close_file(self, add_parent_variants=False):
        self.flush_buffer(last_time=True)
        if self._file is not None:
            self._file.close()

        if add_parent_variants:
            self.merge_parent_variants()

    def merge_parent_variants(self):
        if not self._parent_prefix:
            return

        file_path = self._out_prefix + '_golden.vcf.gz'
        parent_path = self._parent_prefix + '_golden_final.vcf.gz'
        out_file = self._out_prefix + '_golden_final.vcf.gz'

        vcf_reader_parent = vcf.Reader(filename=parent_path, strict_whitespace=True)
        vcf_reader_child = vcf.Reader(filename=file_path, strict_whitespace=True)
        out = bgzf.open(out_file, 'wb')
        vcf_writer = vcf.Writer(out, vcf_reader_child)

        iterate_simulatnously = utils.walk_together(vcf_reader_parent, vcf_reader_child)
        for readers in iterate_simulatnously:
            if readers[1]:
                vcf_writer.write_record(readers[1])
            elif readers[0]:
                vcf_writer.write_record(readers[0])
        out.close()
