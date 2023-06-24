import regex


def compile_regex(split_str: str, error: float):
    """
    Compiles the regex to speed up the code, given the string that should be compiled and the error percentage.
    :param split_str: The substring that the sequence should be split on (str).
    :param error: The allowed error percentage (float).
    :return: Compiled regex pattern (Pattern[str]).
    """

    substitutions = int(len(split_str) / 100 * error)
    return regex.compile(f'({split_str}){{s<={substitutions}}}')


def find_and_split(compiled_pattern, sequence: str, part: int):
    """
    Splits a DNA or peptide sequence on a subsequence (primer/linker sequence) given the compiled regex, the sequence
    and whether the part before or after the split should be returned.

    :param compiled_pattern: The compiled regex pattern of the primer and the allowed error (Pattern).
    :param sequence: The sequence that should be checked (str).
    :param part: If the part before (0) or after (-1) the split should be used (int).
    :return: Coding sequence (str) or None.
    """

    # If the subsequence is present in the sequence, split is a list of more than one item.
    split = compiled_pattern.split(sequence)
    # Return the part before/after the subsequence as a string if there was a split done, else return None.
    return split[part] if len(split) > 1 else None


def translate(coding_seq, codon_table, stop_codons=("TAA", "TAG", "TGA")):
    """
    Function to translate the coding DNA sequence to its corresponding protein sequence using a custom codon table.

    :param coding_seq: The DNA coding sequence (str).
    :param codon_table: The (reprogrammed) codon table (dict).
    :param stop_codons: The stop codons (tuple).
    :return: Amino acid sequence (str).
    """

    peptide = ''
    for i in range(0, len(coding_seq), 3):
        codon = coding_seq[i:i + 3]
        # If codon is of length 3, continue translation, else break.
        if len(codon) != 3:
            break
        if codon in codon_table:
            peptide += codon_table[codon]
            if codon in stop_codons:
                break
        # Unknown amino acid.
        else:
            peptide += 'X'
    return peptide


def translate_backwards(coding_seq, codon_table):
    """
        Function to translate the coding DNA sequence starting at the reverse primers to its corresponding protein
        sequence using a custom codon table. This is to account for frameshifts.

        :param coding_seq: The DNA coding sequence (str).
        :param codon_table: The (reprogrammed) codon table (dict).
        :param stop_codons: The stop codons (tuple).
        :return: Amino acid sequence (str).
    """

    peptide = ''
    # Loop over coding sequence backwards.
    for i in range(len(coding_seq) - 1, -1, -3):
        codon = coding_seq[i - 2:i + 1]
        # If codon is of length 3, continue translation.
        if len(codon) != 3:
            break
        if codon in codon_table:
            peptide += codon_table[codon]
            # if codon in stop_codons:
            #     break
        else:
            peptide += 'X'  # unknown amino acid
    # Reverse resulting peptide.
    return peptide[::-1]


def translate_helm(coding_seq, codon_table, stop_codons=("TAA", "TAG", "TGA")):
    """
        Function to translate the coding DNA sequence to its corresponding protein sequence using a custom codon table.
        HELM sequence is used as output, as this way the codon table can contain unnatural amino acids.

        :param coding_seq: The DNA coding sequence (str).
        :param codon_table: The (reprogrammed) codon table (dict).
        :param stop_codons: The stop codons (tuple).
        :return: HELM amino acid sequence (str).
    """

    peptide = 'PEPTIDE1{'
    for i in range(0, len(coding_seq), 3):
        codon = coding_seq[i:i + 3]
        # If codon is of length 3, continue translation, else break.
        if len(codon) != 3:
            break
        if codon in codon_table:
            if codon in stop_codons:
                break
            peptide += codon_table[codon]
            peptide += '.'
        # Unknown amino acid.
        else:
            peptide += 'X'
    peptide = peptide[:-1]
    peptide += '}$$$$'
    return peptide
