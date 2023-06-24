# Input: fastq sequence files.
# Output:
#  - peptide_counts (csv): file with each DNA sequence, peptide sequence and count
#  - peptide_clean (fasta): file with peptide sequences based on linker, sorted on prevalence
# Optional output (not used for report):
#  - peptide_alignment (?): file with top N peptides aligned locally

# =====Imports===== #
import csv
import os
from collections import Counter, defaultdict
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq
from sequence_analysis import compile_regex, find_and_split, translate, translate_backwards

# =====Inputs===== #
# Folder with sequences.
folder = 'SEQ'
# Number of peptides that should be picked to be aligned (top n from top will be picked).
n = 100
# Minimum length of the peptide.
min_length = 11
library = 'Cyclic'
# Linker and primer sequences
linker_seq = 'GGGGGS'
fw_primer = 'ATACTAATACGACTCACTATAGGATTAAGGAGGTGATATTTATG'
rev_primer = 'TAGGACGGGGGGCGGGAGGCGGG'
# Optional: the linker with a frameshift.
frameshift_linker_seq = 'RWRRR*'
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

    "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "Z",  # AGT -> Z (pentafluorophenylalanine)
    "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
    "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
    "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",

    "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
    "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
    "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
    "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"
}
start_codon = 'ATG'

# Automatically changed.
flipped_fw_primer = str(Seq(fw_primer).reverse_complement())


# =====Regex compiling===== #
# Compile the regex for find_and_split(), the allowed error percentage can be changed here.
fw_primer_pattern = compile_regex(fw_primer, error=0)
flipped_fw_primer_pattern = compile_regex(flipped_fw_primer, error=0)
rev_primer_pattern = compile_regex(rev_primer, error=0)
linker_pattern = compile_regex(linker_seq + '\\*', error=0) if linker_seq else None
frameshift_linker_pattern = compile_regex(frameshift_linker_seq + '\\*', error=0) if frameshift_linker_seq else None


# =====Functions===== #
def process_record(record):
    seq_str = str(record.seq)
    # Check if the sequence is flipped using the flipped forward primer.
    if rev_hit := find_and_split(flipped_fw_primer_pattern, seq_str, 0):
        # Flip the sequence if necessary.
        record = record.reverse_complement()
        seq_str = str(record.seq)
        hit = str(Seq(rev_hit).reverse_complement())
        # Translate the hit (coding part of sequence).
        peptide = translate(hit, codon_table)
        peptide_rev = None
        # Find the part between fw and rev primers.
        if rev_hit := find_and_split(rev_primer_pattern, hit, 0):
            # Translate the hit again, but backwards.
            peptide_rev = translate_backwards(rev_hit, codon_table)
        # Optional: only if linker sequence is given, check if the peptide contains the linker.
        if linker_seq:
            correct_peptide_split = find_and_split(linker_pattern, peptide, 0)
            if not correct_peptide_split:
                correct_peptide_split = find_and_split(frameshift_linker_pattern, peptide, 0)
                if correct_peptide_split:
                    correct_peptide_split += '!'
            correct_peptide_split_rev = None
            if peptide_rev:
                correct_peptide_split_rev = find_and_split(linker_pattern, peptide_rev, 0)
            return peptide, peptide_rev, seq_str, correct_peptide_split, correct_peptide_split_rev
        return peptide, peptide_rev, seq_str, None, None
    # Check if the primer is in the sequence.
    if hit := find_and_split(fw_primer_pattern, seq_str, -1):
        # Translate the hit (coding part of sequence).
        peptide = translate(hit, codon_table)
        peptide_rev = None
        # Find the part between fw and rev primers.
        if rev_hit := find_and_split(rev_primer_pattern, hit, 0):
            # Translate the hit again, but backwards.
            peptide_rev = translate_backwards(rev_hit, codon_table)
        # Optional: only if linker sequence is given, check if the peptide contains the linker.
        if linker_seq:
            correct_peptide_split = find_and_split(linker_pattern, peptide, 0)
            if not correct_peptide_split:
                correct_peptide_split = find_and_split(frameshift_linker_pattern, peptide, 0)
                if correct_peptide_split:
                    correct_peptide_split += '!'

            # TODO the thing above this should give the same effect as below
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
                if correct_peptide_split := result[3]:
                    correct_peptides_split.append(correct_peptide_split)
                if peptide_rev := result[1]:
                    pep_rev_seq.append(peptide_rev)
                    pep_rev_seq_cnt[peptide_rev] += 1
                    if correct_peptide_split_rev := result[4]:
                        correct_peptides_split_rev.append(correct_peptide_split_rev)
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
        # peptide_fasta_lines = [f'>{i + 1}_{n}\n{peptide}\n' for i, (peptide, n) in enumerate(pep_spl_cnt.items())]

        peptide_fasta_lines = []
        for i, (peptide, n) in enumerate(pep_spl_cnt.items()):
            if len(peptide) >= min_length:
                if '!' in peptide:
                    peptide = peptide[:-1]
                    peptide_fasta_lines.append(f'>{i + 1}_{n}_fs\n{peptide}\n')
                else:
                    peptide_fasta_lines.append(f'>{i + 1}_{n}\n{peptide}\n')

        with open(clean_file, 'w') as peptide_fasta:
            peptide_fasta.writelines(peptide_fasta_lines)

        # Same for reverse translated peptides.
        # pep_spl_cnt_rev = dict(Counter(correct_peptides_split_rev).most_common())
        # # Faster to write all lines at once, so first make a list.
        # peptide_fasta_lines_rev = [f'>{i + 1}_{n}\n{peptide}\n' for i, (peptide, n) in enumerate(pep_spl_cnt_rev.items())]
        # with open(clean_file_rev, 'w') as peptide_fasta:
        #     peptide_fasta.writelines(peptide_fasta_lines_rev)
        print(f'Cleaved peptides based on given linker are located in {clean_file}.')


if __name__ == '__main__':
    # Create folders for output files if they don't exist yet.
    os.makedirs(f'OUT/{library}/RF1', exist_ok=True)
    os.makedirs(f'OUT/{library}/RF2', exist_ok=True)
    os.makedirs(f'OUT/Alignments/{library}/RF1', exist_ok=True)
    os.makedirs(f'OUT/Alignments/{library}/RF2', exist_ok=True)
    for filename in os.listdir(folder):
        # Skip any non fastq files.
        if filename.split('.')[-1] != 'fastq':
            continue
        # Set input file path.
        input_file = os.path.join(folder, filename)
        # Get sample name.
        sample = filename.split('_')[0]
        # Set reverse or forward
        if '_R1_' in filename:
            output_name = f'{sample}_FW'
        else:
            output_name = f'{sample}_REV'
        # Get sample round and library.
        # sample_round = sample[-2:]
        target = sample[2:5]
        sample_library = sample[-5:-2]
        if library == 'Linear':
            if sample_library == 'Cyc':
                print('Linear only')
                continue
        elif library == 'Cyclic':
            if sample_library == 'Lin':
                print('Cyclic only')
                continue

        counts_file = f'OUT/Linear/{target}/peptide_counts_{output_name}.csv'
        clean_file = f'OUT/Linear/{target}/peptide_clean_{output_name}.fasta'
        alignment_input = f'OUT/Alignments/{library}/{target}/alignment_input_{output_name}.fasta'
        # alignment_file = f'OUT/Alignments/Linear/{target}/peptide_alignment_{output_name}.aln'
        main()

        # Finally, prepare for alignments.
        with open(clean_file, 'r') as f_full, open(alignment_input, 'w') as f_trunc:
            data = f_full.readlines()
            # Number of sequences * 2, as one sequence is 2 lines in the fasta file.
            if len(data) > n*2:
                data = data[0:n*2]
            f_trunc.writelines(data)
