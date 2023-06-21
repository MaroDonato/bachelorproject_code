# Input: fastq sequence files.
# Output:
#  - peptide_counts (csv): file with each DNA sequence, peptide sequence and count
#  - peptide_clean (fasta): file with peptide sequences based on linker, sorted on prevalence
#  - peptide_alignment (?): file with top N peptides aligned locally
# Optional (if I have time):
#  - Peptide sequence -> SMILES (including the non-canonical amino acids)

# =====Imports===== #
import csv
import os
from collections import Counter, defaultdict
from multiprocessing import Pool

import regex
from Bio import SeqIO
from Bio.Seq import Seq


# =====Inputs===== #
# Filename(s) todo change this (see above)
input_file = 'SEQ2/S2RF1LinR5.fastq'
output_name = 'S2RF1LinR5_R1'

counts_file = f'peptide_counts_{output_name}.csv'
clean_file = f'peptide_clean_{output_name}.fasta'
clean_file_rev = f'peptide_clean_rev_{output_name}.fasta'
alignment_file = f'peptide_alignment{output_name}.?'
# Linker and primer sequences
linker_seq = 'GGGGGS'
fw_primer = 'ATACTAATACGACTCACTATAGGATTAAGGAGGTGATATTTATG'
rev_primer = 'TAGGACGGGGGGCGGGAGGCGGG'
# Codon table
codon_table = {
    "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
    "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
    "TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*",
    "TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W",

    "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
    "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
    "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
    "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",

    "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",  #  todo AGT -> Z (pentafluorophenylalanine)
    "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
    "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
    "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",

    "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
    "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
    "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
    "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"
}
start_codon = 'ATG'
stop_codons = ("TAA", "TAG", "TGA")

# Automatically changed.
flipped_fw_primer = str(Seq(fw_primer).reverse_complement())


# =====Regex compiling===== #
def compile_regex(split_str: str, error: float):
    """
    Compiles the regex to speed up the code, given the string that should be compiled and the error percentage.
    :param split_str: The substring that the sequence should be split on (str).
    :param error: The allowed error percentage (float).
    :return: Compiled regex pattern (Pattern[str]).
    """

    substitutions = int(len(split_str) / 100 * error)
    return regex.compile(f'({split_str}){{s<={substitutions}}}')

# Compile the regex for find_and_split(), the allowed error percentage can be changed here.
fw_primer_pattern = compile_regex(fw_primer, error=0)
flipped_fw_primer_pattern = compile_regex(flipped_fw_primer, error=0)
rev_primer_pattern = compile_regex(rev_primer, error=0)
linker_pattern = compile_regex(linker_seq + '\\*', error=0) if linker_seq else None


# =====Functions===== #
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


def translate(coding_seq):
    """
    Function to translate the coding DNA sequence to its corresponding protein sequence using a custom codon table.

    :param coding_seq: The DNA coding sequence (str).
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


def translate_backwards(coding_seq):
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


def process_record(record):
    seq_str = str(record.seq)
    # Check if the sequence is flipped using the flipped forward primer.
    if rev_hit := find_and_split(flipped_fw_primer_pattern, seq_str, 0):
        # Flip the sequence if necessary.
        record = record.reverse_complement()
        seq_str = str(record.seq)
        hit = str(Seq(rev_hit).reverse_complement())
        # Translate the hit (coding part of sequence).
        peptide = translate(hit)
        peptide_rev = None
        # Find the part between fw and rev primers.
        if rev_hit := find_and_split(rev_primer_pattern, hit, 0):
            # Translate the hit again, but backwards.
            peptide_rev = translate_backwards(rev_hit)
        # Optional: only if linker sequence is given, check if the peptide contains the linker.
        if linker_seq:
            correct_peptide_split = find_and_split(linker_pattern, peptide, 0)
            correct_peptide_split_rev = None
            if peptide_rev:
                correct_peptide_split_rev = find_and_split(linker_pattern, peptide_rev, 0)
            return peptide, peptide_rev, seq_str, correct_peptide_split, correct_peptide_split_rev
        return peptide, peptide_rev, seq_str, None, None
    # Check if the primer is in the sequence.
    if hit := find_and_split(fw_primer_pattern, seq_str, -1):
        # Translate the hit (coding part of sequence).
        peptide = translate(hit)
        peptide_rev = None
        # Find the part between fw and rev primers.
        if rev_hit := find_and_split(rev_primer_pattern, hit, 0):
            # Translate the hit again, but backwards.
            peptide_rev = translate_backwards(rev_hit)
        # Optional: only if linker sequence is given, check if the peptide contains the linker.
        if linker_seq:
            correct_peptide_split = find_and_split(linker_pattern, peptide, 0)
            correct_peptide_split_rev = None
            if peptide_rev:
                correct_peptide_split_rev = find_and_split(linker_pattern, peptide_rev, 0)
            return peptide, peptide_rev, seq_str, correct_peptide_split, correct_peptide_split_rev
        return peptide, peptide_rev, seq_str, None, None
    return None


def main():
    # Create empty lists for the sequences.
    DNA_seq = []
    pep_seq = []
    pep_rev_seq = []
    correct_peptides_split = []
    correct_peptides_split_rev = []
    # Create defaultdict to count the number of peptide sequences.
    pep_seq_cnt = defaultdict(int)
    pep_rev_seq_cnt = defaultdict(int)
    # The multiprocessing pool is used to use all CPU cores simultaneously.
    with Pool() as pool:
        print(f'Analyzing sequences for {output_name}...')
        # call the function for each item in parallel, get results as tasks complete
        for result in pool.imap_unordered(process_record, SeqIO.parse(input_file, "fastq")):
            # Save the results in their lists.
            if result:
                pep_seq.append(result[0])
                pep_seq_cnt[result[0]] += 1
                DNA_seq.append(result[2])
                if result[3]:
                    correct_peptides_split.append(result[3])
                if result[1]:
                    pep_rev_seq.append(result[1])
                    pep_rev_seq_cnt[result[1]] += 1
                    if result[4]:
                        correct_peptides_split_rev.append(result[4])
    print(f'Parsed all records for {output_name}, writing files...')

    # Dictionary with all unique peptides and their count.
    combined_lists = [
        (DNA, pep, pep_seq_cnt[pep], pep_rev, pep_rev_seq_cnt[pep_rev])
        for DNA, pep, pep_rev in zip(DNA_seq, pep_seq, pep_rev_seq)]
    unique_combined_lists = list(dict.fromkeys(combined_lists))
    sorted_list = sorted(unique_combined_lists, key=lambda x: x[2], reverse=True)
    with open(counts_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=";")
        writer.writerow(("DNA Sequence", "Peptide Sequence", "Count", "Rev Peptide Sequence", "Rev Count"))
        # Writing all lines at once.
        writer.writerows(sorted_list)
    print(f'DNA and peptide sequences with count are located in {counts_file}.')

    # Optional: only if linker sequence is given.
    if linker_seq:
        pep_spl_cnt = dict(Counter(correct_peptides_split).most_common())
        # Faster to write all lines at once, so first make a list.
        peptide_fasta_lines = [f'>{i + 1}_{n}\n{peptide}\n' for i, (peptide, n) in enumerate(pep_spl_cnt.items())]
        with open(clean_file, 'w') as peptide_fasta:
            peptide_fasta.writelines(peptide_fasta_lines)

        # Same for reverse translated peptides.
        pep_spl_cnt_rev = dict(Counter(correct_peptides_split_rev).most_common())
        # Faster to write all lines at once, so first make a list.
        peptide_fasta_lines_rev = [f'>{i + 1}_{n}\n{peptide}\n' for i, (peptide, n) in enumerate(pep_spl_cnt_rev.items())]
        with open(clean_file_rev, 'w') as peptide_fasta:
            peptide_fasta.writelines(peptide_fasta_lines_rev)
        print(f'Cleaved peptides based on given linker are located in {clean_file} and {clean_file_rev}.')


if __name__ == '__main__':
    # main()
    for filename in os.listdir('SEQ'):
        # Skip any non fastq files.
        if filename.split('.')[-1] != 'fastq':
            continue
        # Set input file path.
        input_file = os.path.join('SEQ', filename)
        # Get sample name.
        sample = filename.split('_')[0]
        # Set reverse or forward
        if '_R1_' in filename:
            output_name = f'{sample}_FW'
        else:
            output_name = f'{sample}_REV'
        # Get sample round and library.
        # sample_round = sample[-2:]
        sample_library = sample[-5:-2]
        if sample_library == 'Cyc':
            codon_table = {
                "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
                "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
                "TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*",
                "TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W",

                "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
                "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
                "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
                "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",

                "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "Z",  # AGT -> Z (pentafluorophenylalanine)
                "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
                "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
                "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",

                "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
                "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
                "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
                "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"
            }
        else:
            codon_table = {
                "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
                "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
                "TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*",
                "TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W",

                "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
                "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
                "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
                "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",

                "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",
                "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
                "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
                "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",

                "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
                "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
                "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
                "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"
            }

        counts_file = f'OUT/peptide_counts_{output_name}.csv'
        clean_file = f'OUT/peptide_clean_{output_name}.fasta'
        clean_file_rev = f'OUT/peptide_clean_rev_{output_name}.fasta'
        alignment_file = f'OUT/peptide_alignment{output_name}.?'
        main()
        # break
