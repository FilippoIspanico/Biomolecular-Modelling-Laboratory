import numpy as np

def get_frequencies(filename: str):
    """
    This function parses the input file that contains BLAST sequences and returns two lists. The first one contains
    the most frequent amino acid at every position, the second one contains their frequency at every position.


    For
    :param filename: must be a vaild path to the file that stores the result of the query done at:
                    https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
    :return:
    """
    alignments = open(filename)
    lines = alignments.readlines()
    n_lines = len(lines)

    sequences = []

    for idx, line in enumerate(lines):
        if line.startswith('>'):
            current_sequence = ""
            i = 1
            while (idx + i) < n_lines and not lines[idx + i].startswith('>'):  # continua fino alla prossima '>'
                current_sequence += lines[idx + i].strip()  # rimuove newline
                i += 1
            sequences.append(current_sequence)


    number_of_sequences = len(sequences)
    print(f'Found {number_of_sequences} aligned sequences')


    lengths = []
    for sequence in sequences:
        lengths.append(len(sequence))


    most_common = []
    frequencies = []
    for position in range(129):
        amminoacids = []
        for sequence in sequences:
            # we need to count the number of the same amino acids: if a sequence
            if position < len(sequence):
                amminoacids.append(sequence[position])
            else:
                amminoacids.append(str('no amminoacid here'))

        unique, counts = np.unique(amminoacids, return_counts=True)
        frequencies.append(np.max(counts) / number_of_sequences)
        most_common_idx = np.argmax(counts)
        most_common.append(unique[most_common_idx])


    return  most_common, frequencies,
