# Translate DNA to protein, with a reprogrammed dictionary.

# How to run
# Place the fastq files in a folder in the same directory as the python script.
# python3 nanopore_script.py
# NOTE: requires biopython and regex to be installed.
# See https://pypi.org/project/biopython/ and https://pypi.org/project/regex/ for installation instructions.


# Imports
import csv
import sys
import os
from collections import Counter, defaultdict
from os import system

import regex
from Bio import SeqIO
from Bio.Seq import Seq


# The name of the folder in which the fastq files are located.
folder_name = 'library'
# File names.
input_file = 'combined.fastq'
output1_file = 'peptide_counts.csv'
output2_file = 'peptide_clean.fasta'
# Forward primer sequence.
fw_primer = 'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG'
# The linker sequence, can be empty or string.
# If given, filters out peptides without this sequence into output2_file.
linker_seq = 'GAGAGA'

# input_fastq = sys.argv[1]
# output_txt = sys.argv[2]

# Check the current platform and use the appropriate command to concatenate the fastq files.
current_platform = sys.platform
if current_platform == 'win32':
    print('Running on Windows')
    system(f'type {os.getcwd()}{os.sep}{folder_name}{os.sep}*.fastq > {input_file}')
elif current_platform == 'linux' or current_platform == 'linux2' or current_platform == 'darwin':
    print('Running on Linux/macOS')
    system(f'cat {os.getcwd()}{os.sep}{folder_name}{os.sep}*.fastq > {input_file}')

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
    "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"}
start_codon = 'ATG'
stop_codons = ["TAA", "TAG", "TGA"]
flipped_fw_primer = str(Seq(fw_primer).reverse_complement())


def find_primer(primer_sequence: str, sequence: str, error: float = 20):
    """
    Finds a sequence containing a primer given the primer sequence and the allowed substitutions (defaults to 20%). Also
    works for a flipped primer.

    :param primer_sequence: The sequence of the primer (str).
    :param sequence: The sequence that should be checken (str).
    :param error: The error percentage (float).
    :return: Coding sequence (str) or None.
    """

    substitutions = int(len(primer_sequence) / 100 * error)
    # If there is a primer in the sequence, split is a list of more than one item.
    split = regex.split(f'({primer_sequence}){{s<={substitutions}}}', sequence)
    # Return the part after the primer as a string if there was a split done, else return None.
    return split[-1] if len(split) > 1 else None


def find_and_cleave_peptide(search_str: str, linker: str, error: float = 30):
    """
    Gives the part before the linker of a peptide given its linker sequence and the allowed substitutions in the linker
    (defaults to 30%).

    :param search_str: The peptide sequence that should be split (str).
    :param linker: The sequence of the linker (str).
    :param error: The allowed error (int, n of amino acids).
    :return: Hit peptide sequence (str) or None.
    """

    substitutions = int(len(linker) / 100 * error)
    split = regex.split(f'({linker}\\*){{s<={str(substitutions)}}}', search_str)
    # Return the part before the linker as a string if there was a split done, else return None.
    return split[0] if len(split) > 1 else None


def translate(coding_seq, codon_table):
    """
    Function to translate the coding DNA sequence to its corresponding protein sequence using a custom codon table.

    :param coding_seq: The DNA coding sequence (str).
    :param codon_table: The (custom) codon table for translation (dict).
    :return: Amino acid sequence (str).
    """

    peptide = ''
    for i in range(0, len(coding_seq), 3):
        codon = coding_seq[i:i + 3]
        # If codon is of length 3, continue translation.
        if len(codon) != 3:
            # print(f"Codon encountered of less then 3 amino acids:{codon}. Translation terminated., frameshift may have occured.")
            break
        if codon in codon_table:
            peptide += codon_table[codon]
            if codon in stop_codons:
                break
        else:
            peptide += 'X'  # unknown amino acid
    return peptide


def main():
    print(f'Analyzing sequences from "{folder_name}" folder...')

    DNA_seq = []
    pep_seq = []
    correct_peptides_split = defaultdict(int)

    for record in SeqIO.parse(input_file, "fastq"):
        seq_str = str(record.seq)
        # Check if the sequence is flipped using the flipped forward primer.
        if find_primer(flipped_fw_primer, seq_str):
            # Flip the sequence if necessary.
            record = record.reverse_complement()
            seq_str = str(record.seq)
        # Check if the primer is in the sequence.
        if hit := find_primer(fw_primer, seq_str):
            # Translate the hit (coding part of sequence).
            peptide = translate(hit, codon_table)
            pep_seq.append(peptide)
            DNA_seq.append(seq_str)
            # Optional: only if linker sequence is given, check if the peptide contains the linker.
            if linker_seq:
                if correct_peptide_split := find_and_cleave_peptide(peptide, linker_seq):
                    correct_peptides_split[correct_peptide_split] += 1

    # Dictionary with all unique peptides and their count.
    seq_cnt = []
    for pep in pep_seq:
        seq_cnt.append(pep_seq.count(pep))
    combined_lists = list(zip(DNA_seq, pep_seq, seq_cnt))
    unique_combined_lists = list(dict.fromkeys(combined_lists))
    sorted_list = sorted(unique_combined_lists, key=lambda x: x[2], reverse=True)
    with open(output1_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["DNA Sequence", "Peptide Sequence", "Count"])
        for row in sorted_list:
            writer.writerow(row)
    print(f'DNA and peptide sequences with count are located in {output1_file}.')

    # Optional: only if linker sequence is given.
    if linker_seq:
        pep_spl_cnt = dict(Counter(correct_peptides_split).most_common())
        with open(output2_file, 'w') as peptide_fasta:
            for i, (peptide, n) in enumerate(pep_spl_cnt.items()):
                peptide_fasta.write(f'>{i + 1}_{n}\n{peptide}\n')
        print(f'Cleaved peptides based on given linker are located in {output2_file}.')


if __name__ == "__main__":
    main()
#     TODO speeding up: maybe use regex.finditer() for better performance and pre compile the regex.
