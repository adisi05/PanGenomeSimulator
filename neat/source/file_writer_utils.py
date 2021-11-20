def reverse_complement(dna_string) -> str:
    """
    Return the reverse complement of a string from a DNA strand. Found this method that is slightly faster than
    biopython. Thanks to this stack exchange post:
    https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
    :param dna_string: string of DNA, either in string or Seq format
    :return: the reverse complement of the above string in either string or MutableSeq format
    """
    if type(dna_string) != str:
        dna_string.reverse_complement()
        return dna_string
    else:
        tab = str.maketrans("ACTGN", "TGACN")

        return dna_string.translate(tab)[::-1]

# TODO figure out an optimum batch size
BUFFER_BATCH_SIZE = 8000  # write out to file after this many reads