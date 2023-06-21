# Aligns peptides using https://github.com/Merck/PepSeA, as this allows for non-natural amino acids to be in the
# peptide sequences.

from Bio import SeqIO
# All files in the PepSeA directory are from https://github.com/Merck/PepSeA.
from PepSeA.AlignSubPeptides import align_sub_peptides
from PepSeA.ApiUtils import extract_helm_from_json, json_output
from PepSeA.models import MafftMethods

input_file = 'test.fasta'

# Mafft installation path
mafft_path = '/opt/homebrew/bin/'

# Variables (defaults work fine)
mafft_binary = mafft_path + MafftMethods.ginsi
rocs_path = 'PepSeA/ROCS'
monomers_map_file = 'PepSeA/monomers_map.txt'

helm_sequences = []
for record in SeqIO.parse(input_file, 'fasta'):
    # TODO make smiles out of HELM if necessary (using rdmolfiles).
    # Making HELM files in the first place is also easy: change the reprogrammed AA with [NAME] in the dictionary,
    # then add a dot after every amino acid during reprogramming. (A more advanced script could even make them cyclic)
    helm_seq = f'PEPTIDE1{{{".".join([amino_acid for amino_acid in str(record.seq)])}}}$$$$'
    helm_sequences.append({'ID': str(record.id), 'HELM': helm_seq})

align_output, _, alignment_score = align_sub_peptides(
    extract_helm_from_json(helm_sequences),
    gap_opening_penalty=1.53,
    gap_extension_penalty=0.0,
    polymer_to_align=None,
    mafft_options='',
    path_to_mafft=mafft_binary,
    path_to_subst_matrix=rocs_path,
    path_to_monomer_table=monomers_map_file)

# Save alignment result in a clustal alignment file.
alignment_result = {"Alignment": json_output(align_output, helm_sequences), "AlignmentScore": alignment_score}
with open('clustal_alignment.aln', 'w') as f:
    for i in alignment_result['Alignment']:
        f.write(f'{i["ID"]}      {i["AlignedSeq"]}\n')
