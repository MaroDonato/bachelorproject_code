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

if __name__ == '__main__':
    for filename in os.listdir('SEQ'):
        input_file = os.path.join('SEQ', filename)
        sample = filename.split('_')[0]
        if '_R1_' in filename:
            output_name = f'{sample}_R1'
        else:
            output_name = f'{sample}_R2'
        print(output_name, input_file)

        counts_file = f'peptide_counts_{output_name}.csv'
        clean_file = f'peptide_clean_{output_name}.fasta'
        alignment_file = f'peptide_alignment{output_name}.?'
        main(input_file, counts_file, clean_file)
        break