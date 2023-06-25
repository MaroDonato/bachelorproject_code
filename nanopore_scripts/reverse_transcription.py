import csv
from collections import Counter, defaultdict
from multiprocessing import Pool

from nanopore_script_optimized import compile_regex, find_and_split, fw_primer_pattern

# Input file
input_file = 'peptide_counts.csv'
output_file = 'peptides.fasta'
# Tag sequence
tag_sequences = [
    'GCATCA',
    'GAGTCA',
    'TTGCTG',
    'AGAGAT',
    'CTGAAA',
    'CAGTGG',
    'CGTGTA',
    'TTAAGG',
    'TGAGCC',
    'TTCGGA'
]
fw_primer = 'TAATACGACTCACTATAGGGTTAACTTTAAGAAGGAGATATACATATG'
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


def translate_backwards(coding_seq, codon_table):
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


# Create patterns.
tag_seq_patterns = {}
for tag_seq in tag_sequences:
    tag_seq_patterns[tag_seq] = compile_regex(tag_seq, 20)
spacer_pattern = compile_regex('GGYPYDVPDYGAGAG', 20)


def process_peptide(row):
    DNA_seq = row[0]
    for tag_seq in tag_sequences:
        # Check if DNA sequence contains the tag sequence.
        if seq := find_and_split(tag_seq_patterns[tag_seq], DNA_seq, 0):
            # Remove the primer sequence.
            if coding_seq := find_and_split(fw_primer_pattern, seq, -1):
                peptide = translate_backwards(coding_seq, codon_table)
                if peptide_cleaved := find_and_split(spacer_pattern, peptide, 0):
                    return tag_seq, peptide_cleaved
    return None


def main():
    peptides = defaultdict(list)
    with open(input_file, 'r') as csvfile:
        with Pool() as pool:
            # call the function for each item in parallel, get results as tasks complete
            for result in pool.imap_unordered(process_peptide, csv.reader(csvfile, delimiter=',')):
                # Save the results in a list.
                if result:
                    peptides[result[0]].append(result[1])

    # for tag_seq in tag_sequences:
    #     peptides = []
    #     with open(input_file, 'r') as csvfile:
    #         reader = csv.reader(csvfile, delimiter=',')
    #         for row in reader:
    #             DNA_seq = row[0]
    #             # Check if DNA sequence contains the tag sequence.
    #             if seq := find_and_cleave_peptide(DNA_seq, tag_seq):
    #                 # Remove the primer sequence.
    #                 if coding_seq := find_primer(fw_primer, seq):
    #                     peptide = translate_backwards(coding_seq, codon_table)
    #                     if peptide_cleaved := find_and_cleave_peptide(peptide, 'GGYPYDVPDYGAGAG'):
    #                         peptides.append(peptide_cleaved)

    for tag_seq, pep_list in peptides.items():
        pep_cnt = dict(Counter(pep_list).most_common())
        # pep_cnt = dict(Counter(peptides).most_common())
        peptide_fasta_lines = [f'>{i + 1}_{n}\n{peptide}\n' for i, (peptide, n) in enumerate(pep_cnt.items())]
        with open(tag_seq + output_file, 'w') as peptide_fasta:
            peptide_fasta.writelines(peptide_fasta_lines)
            # for i, (peptide, n) in enumerate(pep_cnt.items()):
            #     peptide_fasta.write(f'>{i + 1}_{n}\n{peptide}\n')


if __name__ == "__main__":
    main()
