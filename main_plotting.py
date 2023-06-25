import pandas as pd
from matplotlib import rc
import os
from plotting import plot_recovery, plot_recovery_split, plot_diversity_horizontal, plot_diversity

# Font settings for all figures.
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text')


def recovery_dict_from_excel(doc, sheet):
    # Read in excel data.
    df = pd.read_excel(doc, sheet_name=sheet)
    # Get the column names
    cols = df.columns
    recoveries = {}
    for i in range(0, len(cols), 2):
        # Add recoveries and corresponding values to dict, omitting NaN values.
        recoveries[cols[i]] = dict(zip(
            df[cols[i]].dropna(),
            # Multiply the values by 100 to get the percentages.
            df[cols[i + 1]] * 100))
    return recoveries


def diversity_dict_from_folder(folder_1, folder_2):
    diversities = {
        'RF1': {
            'Lin': {},
            'Cyc': {}
        },
        'RF2': {
            'Lin': {},
            'Cyc': {}
        }
    }
    for filename in os.listdir(folder_1):
        # Skip any non fastq files.
        if filename.split('.')[-1] != 'fastq':
            continue
        # Get sample name.
        sample = filename.split('_')[0]
        # Skip first selection.
        if sample[:2] != 'S2':
            continue
        # Skip reverse.
        if '_R1_' not in filename:
            continue
        sample_name = f'{sample}_FW'
        # Get sample round and library.
        target = sample[2:5]
        sample_library = sample[-5:-2]
        # Get length of file.
        with open(os.path.join(folder_1, filename)) as f1, open(
                f'{folder_2}/{"Linear" if sample_library == "Lin" else "Cyclic"}/{target}/'
                f'peptide_counts_{sample_name}.csv', 'r') as f2:
            total_seq = len(f1.readlines())/4
            unique_seq = len(f2.readlines())-1
            diversities[target][sample_library][sample[-1]] = unique_seq/total_seq*100
    return diversities


def generate_figures():
    # Plotting of recoveries.
    # plot_recovery_split(recovery_dict_from_excel('TRAP rounds (excel)/TRAP rounds.xlsx', 'Recovery RF1'),
    #                     'RF1 recovery', 'RF1 recovery S1', ymin1=0, ymax1=0.021, step1=0.002,
    #                     ymin2=0.1, ymax2=0.7001, step2=0.1)
    # plot_recovery_split(recovery_dict_from_excel('TRAP rounds (excel)/TRAP rounds.xlsx', 'Recovery RF2'),
    #                     'RF2 recovery', 'RF2 recovery S1', ymin1=0, ymax1=0.021, step1=0.002,
    #                     ymin2=0.1, ymax2=0.7001, step2=0.1)
    # plot_recovery(recovery_dict_from_excel('TRAP rounds (excel)/TRAP rounds 2.xlsx', 'Recovery RF1'),
    #               'RF1 recovery', 'RF1 recovery S2', ymin=0, ymax=0.22001)
    # plot_recovery(recovery_dict_from_excel('TRAP rounds (excel)/TRAP rounds 2.xlsx', 'Recovery RF2'),
    #               'RF2 recovery', 'RF2 recovery S2', ymin=0, ymax=0.22001)

    # Plotting of qPCR diversities.
    docs = [
        'qPCR/qPCR of rPCRs of RF1 lin round 1-5.xlsx',
        'qPCR/qPCR of rPCRs of RF1 cyc round 1-5.xlsx',
        'qPCR/qPCR of rPCRs of RF2 lin round 1-4.xlsx',
        'qPCR/qPCR of rPCRs of RF2 cyc round 1-4.xlsx'
    ]
    titles = [
        'RF1 linear',
        'RF1 cyclic',
        'RF2 linear',
        'RF2 cyclic'
    ]
    cycles = 35
    # plot_diversity_horizontal('Diversity plot', docs, titles, cycles)
    plot_diversity('Epoxy bead clone assay', 'qPCR/Epoxy qpcr.xlsx', 'Epoxy bead clone assay', cycles)
    print('Generated figures.')


if __name__ == '__main__':
    generate_figures()
