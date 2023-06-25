# Translate DNA to protein, with a reprogrammed dictionary.

# How to run:
# Place the fastq files in a folder in the same directory as the python script. For example:
# .
# ├── library
# │   ├── file_1.fastq
# │   ├── ...
# │   └── file_n.fastq
# ├── nanopore_script.py
# └── reverse_transcription.py
# Then, go into the folder containing the python files using cd. For example:
# cd Documents
# Now you can run the python files:
# python3 nanopore_script.py

# NOTE: requires biopython and regex to be installed.
# See https://pypi.org/project/biopython/ and https://pypi.org/project/regex/ for installation instructions.

# Imports
import csv
import os
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from os import system

import regex
from Bio import SeqIO
from Bio.Seq import Seq

# Variables, feel free to edit.
# The name of the folder in which the fastq files are located.
folder_name = 'library'
# File names.
input_file = 'combined.fastq'
output1_file = 'peptide_counts.csv'
output2_file = 'peptide_clean.fasta'
# Forward primer sequence.
fw_primer = 'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG'
# The linker sequence, can be None or a string. If given, filters out peptides without this sequence into output2_file.
linker_seq = 'GAGAGA'
# Customize the codon table if needed!
# Codon table below: start codon M --> Y
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
start_codon = 'ATG'
stop_codons = ("TAA", "TAG", "TGA")


# Functions, no need to edit beyond this part.
def concatenate_files():
    """
    Checks the current platform (linux/macOS or windows) and uses the appropriate command to concatenate the fastq
    files.
    """

    current_platform = sys.platform
    if current_platform == 'win32':
        print('Running on Windows')
        system(f'type {os.getcwd()}{os.sep}{folder_name}{os.sep}*.fastq > {input_file}')
    elif current_platform == 'linux' or current_platform == 'linux2' or current_platform == 'darwin':
        print('Running on Linux/macOS')
        system(f'cat {os.getcwd()}{os.sep}{folder_name}{os.sep}*.fastq > {input_file}')


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


# Flip the forward primer using biopython.
flipped_fw_primer = str(Seq(fw_primer).reverse_complement())
# Compile the regex for find_and_split(), the allowed error percentage can be changed here.
fw_primer_pattern = compile_regex(fw_primer, error=20)
flipped_fw_primer_pattern = compile_regex(flipped_fw_primer, error=20)
linker_pattern = compile_regex(linker_seq + '\\*', error=30) if linker_seq else None


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
        # Optional: only if linker sequence is given, check if the peptide contains the linker.
        if linker_seq:
            if correct_peptide_split := find_and_split(linker_pattern, peptide, 0):
                return peptide, seq_str, correct_peptide_split
        return peptide, seq_str, None
    # Check if the primer is in the sequence.
    if hit := find_and_split(fw_primer_pattern, seq_str, -1):
        # Translate the hit (coding part of sequence).
        peptide = translate(hit)
        # Optional: only if linker sequence is given, check if the peptide contains the linker.
        if linker_seq:
            if correct_peptide_split := find_and_split(linker_pattern, peptide, 0):
                return peptide, seq_str, correct_peptide_split
        return peptide, seq_str, None
    return None


def main():
    # Concatenate the files in the folder.
    concatenate_files()
    print(f'Combined all files from "{folder_name}" folder into {input_file}.')
    # Create empty lists for the sequences.
    DNA_seq = []
    pep_seq = []
    correct_peptides_split = []
    # Create defaultdict to count the number of peptide sequences.
    pep_seq_cnt = defaultdict(int)
    # The multiprocessing pool is used to use all CPU cores simultaneously.
    with Pool() as pool:
        print(f'Analyzing sequences...')
        # call the function for each item in parallel, get results as tasks complete
        for result in pool.imap_unordered(process_record, SeqIO.parse(input_file, "fastq")):
            # Save the results in their lists.
            if result:
                pep_seq.append(result[0])
                pep_seq_cnt[result[0]] += 1
                DNA_seq.append(result[1])
                if result[2]:
                    correct_peptides_split.append(result[2])
    print('Parsed all records, writing files...')

    # for record in SeqIO.parse(input_file, "fastq"):
    #     seq_str = str(record.seq)
    #     # Check if the sequence is flipped using the flipped forward primer.
    #     if rev_hit := find_and_cleave_peptide(flipped_fw_primer_pattern, seq_str):
    #         # Flip the sequence if necessary.
    #         record = record.reverse_complement()
    #         seq_str = str(record.seq)
    #         hit = str(Seq(rev_hit).reverse_complement())
    #         # Translate the hit (coding part of sequence).
    #         peptide = translate(hit, codon_table)
    #         pep_seq.append(peptide)
    #         DNA_seq.append(seq_str)
    #         # Optional: only if linker sequence is given, check if the peptide contains the linker.
    #         if linker_seq:
    #             if correct_peptide_split := find_and_cleave_peptide(linker_pattern, peptide):
    #                 correct_peptides_split.append(correct_peptide_split)
    #         continue
    #     # Check if the primer is in the sequence.
    #     if hit := find_primer(fw_primer_pattern, seq_str):
    #         # Translate the hit (coding part of sequence).
    #         peptide = translate(hit, codon_table)
    #         pep_seq.append(peptide)
    #         DNA_seq.append(seq_str)
    #         # Optional: only if linker sequence is given, check if the peptide contains the linker.
    #         if linker_seq:
    #             if correct_peptide_split := find_and_cleave_peptide(linker_pattern, peptide):
    #                 correct_peptides_split.append(correct_peptide_split)

    # Dictionary with all unique peptides and their count.
    combined_lists = [(DNA, pep, pep_seq_cnt[pep]) for DNA, pep in zip(DNA_seq, pep_seq)]
    unique_combined_lists = list(dict.fromkeys(combined_lists))
    sorted_list = sorted(unique_combined_lists, key=lambda x: x[2], reverse=True)
    with open(output1_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(("DNA Sequence", "Peptide Sequence", "Count"))
        # Writing all lines at once.
        writer.writerows(sorted_list)
    print(f'DNA and peptide sequences with count are located in {output1_file}.')

    # Optional: only if linker sequence is given.
    if linker_seq:
        pep_spl_cnt = dict(Counter(correct_peptides_split).most_common())
        # Faster to write all lines at once, so first make a list.
        peptide_fasta_lines = [f'>{i + 1}_{n}\n{peptide}\n' for i, (peptide, n) in enumerate(pep_spl_cnt.items())]
        with open(output2_file, 'w') as peptide_fasta:
            peptide_fasta.writelines(peptide_fasta_lines)
        print(f'Cleaved peptides based on given linker are located in {output2_file}.')


if __name__ == "__main__":
    main()
    # 1.7 minutes per million sequences on macbook
