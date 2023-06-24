# Plotting functions.
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rc
import numpy as np


def plot_diversity_horizontal(plot_title, docs, titles, cycles):
    # Set the plot colors to VU colors
    plot_colors = {
        'Round 0': '#0077B3',
        'Round 1': '#4FAF48',
        'Round 2': '#E8692D',
        'Round 3': '#8E4DA4',
        'Round 4': '#F2BA2F',
        'Round 5': '#D4CAC8',
        'Round 6': '#575756',
        'Round 7': '#003F6C'
    }
    fig, axes = plt.subplots(1, 4, sharey=True, figsize=(8, 2))
    # Title.
    fig.suptitle('Diversity per round', fontsize=16)

    n_labels = 0
    for ax, doc, title in zip(axes.reshape(-1), docs, titles):
        empty_df = False
        i = -1
        while not empty_df:
            i += 1
            df = pd.read_excel(doc, skiprows=2 + (cycles + 2) * i, nrows=cycles + 1)
            empty_df = df.empty
            if not df.empty:
                ax.plot(df.iloc[:, 1], df.iloc[:, 2], color=plot_colors[f'Round {i}'], label=f'Round {i}')
            ax.set_title(title)
            ax.xaxis.set_ticks(np.arange(5, cycles + 0.0001, 10))
            ax.margins(x=0)
            # ax.set_xlim(0, cycles)
        else:
            if i > n_labels:
                n_labels = i
    # X and Y labels
    fig.text(0.5, 0.04, 'Cycles', ha='center', fontsize=12)
    fig.text(0.02, 0.5, 'Fluorescence (R)', va='center', rotation='vertical', fontsize=12)
    fig.legend([f'Round {n}' for n in range(n_labels)], frameon=False, loc='center right')
    fig.tight_layout()
    fig.subplots_adjust(left=0.10, bottom=0.24, right=0.85)
    plt.savefig(plot_title + '.png', dpi=1200, transparent=True)
    fig.savefig(plot_title + '.pdf', format='pdf')
    plt.show()


def plot_recovery(recoveries, plot_title, ymin, ymax):
    # Initialise the plot.
    fig, axes = plt.subplots(
        1, (len(recoveries)), sharey='row',
        # Make first subplot half the size of the other subplots.
        gridspec_kw={'width_ratios': [1] + [2 for _ in range(len(recoveries) - 1)]}
    )
    # Settings of plot.
    xmin, xmax = 0, 4
    bar_width = 0.8
    padding = bar_width / 2
    # Bar hex colors.
    bar_colors = {
        'Linear +': '#0077B3',
        'Cyclic +': '#E8692D',
        'Linear -': '#DFF2FD',
        'Cyclic -': '#FCD3B6'
    }
    # Title.
    fig.suptitle(plot_title, fontsize=16, y=0.95)
    for i, round_n in enumerate(recoveries):
        # Set the x- and y-axis range.
        axes[i].set_ylim([ymin, ymax])
        # Plot data from dictionary.
        # axes[i].bar(x=list(recoveries[round_n].keys()), height=list(recoveries[round_n].values()),
        #             width=bar_width, align='edge')
        for label, height in recoveries[round_n].items():
            axes[i].bar(x=label, height=height, width=bar_width, align='edge',
                        color=bar_colors[label], label=label, edgecolor='black', linewidth=1.0)
        # Remove x ticks.
        axes[i].set_xticks([])
        # Set round number as the x-axis label.
        axes[i].set_xlabel(round_n)
        # Remove left border and set different x-axis range if it is not the first subplot.
        if i != 0:
            axes[i].spines['left'].set_visible(False)
            axes[i].get_yaxis().set_visible(False)
            axes[i].set_xlim([xmin - padding, xmax + padding])
        # Set different x-axis range and y-axis label if it is the first subplot.
        else:
            axes[i].set_xlim([xmin / 2 - padding / 2, xmax / 2 + padding / 2])
            axes[i].yaxis.set_ticks(np.arange(ymin, ymax, 0.02))
            axes[i].set_ylabel('Recovery (%)')
        # Remove the right border if it is not the last subplot.
        if i != (len(recoveries) - 1):
            axes[i].spines['right'].set_visible(False)
        else:
            axes[i].set_xlim([xmin - padding, xmax])
    # Remove space between subplots.
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.legend(frameon=False)  # , title='Library')
    # Show the plot.
    plt.savefig(plot_title + '.png', dpi=1200, transparent=True)
    fig.savefig(plot_title + '.pdf', format='pdf')
    plt.show()


def plot_recovery_split(recoveries, plot_title, ymin1, ymax1, step1, ymin2, ymax2, step2):
    # Initialise the plot.
    fig, axes = plt.subplots(
        2, (len(recoveries)), sharey='row',
        gridspec_kw={
            # Make first subplot half the width of the other subplots.
            'width_ratios': [1] + [2 for _ in range(len(recoveries) - 1)],
            # Make top part of the plot half the height of the lower part.
            'height_ratios': [1, 2]
        }
    )
    # Settings of plot.
    xmin, xmax = 0, 4
    bar_width = 0.8
    padding = bar_width / 2
    # Bar hex colors.
    bar_colors = {
        'Linear +': '#0077B3',
        'Cyclic +': '#E8692D',
        'Linear -': '#DFF2FD',
        'Cyclic -': '#FCD3B6'
    }
    # Font for figure.
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text')
    # Title.
    fig.suptitle(plot_title, fontsize=16, y=0.95)
    for i, round_n in enumerate(recoveries):
        upper_ax = axes[0, i]
        lower_ax = axes[1, i]
        # Set the x- and y-axis range.
        lower_ax.set_ylim([ymin1, ymax1])
        upper_ax.set_ylim([ymin2, ymax2])
        # Plot data from dictionary.
        # axes[i].bar(x=list(recoveries[round_n].keys()), height=list(recoveries[round_n].values()),
        #             width=bar_width, align='edge')
        for label, height in recoveries[round_n].items():
            lower_ax.bar(x=label, height=height, width=bar_width, align='edge',
                         color=bar_colors[label], label=label, edgecolor='black', linewidth=1.0)
            upper_ax.bar(x=label, height=height, width=bar_width, align='edge',
                         color=bar_colors[label], label=label, edgecolor='black', linewidth=1.0)
        # Remove x ticks.
        lower_ax.set_xticks([])
        upper_ax.set_xticks([])
        # Remove lines in between plots.
        lower_ax.spines['top'].set_visible(False)
        upper_ax.spines['bottom'].set_visible(False)
        # Set round number as the x-axis label.
        lower_ax.set_xlabel(round_n)
        # Todo clean up these if statements, maybe to match?
        #  https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html
        # Remove left border and set different x-axis range if it is not the first subplot.
        if i != 0:
            lower_ax.spines['left'].set_visible(False)
            upper_ax.spines['left'].set_visible(False)
            lower_ax.get_yaxis().set_visible(False)
            upper_ax.get_yaxis().set_visible(False)
            lower_ax.set_xlim([xmin - padding, xmax + padding])
            upper_ax.set_xlim([xmin - padding, xmax + padding])
        # Set different x-axis range and y-axis label if it is the first subplot.
        else:
            lower_ax.set_xlim([xmin / 2 - padding / 2, xmax / 2 + padding / 2])
            upper_ax.set_xlim([xmin / 2 - padding / 2, xmax / 2 + padding / 2])
            lower_ax.yaxis.set_ticks(np.arange(ymin1, ymax1, step1))
            upper_ax.yaxis.set_ticks(np.arange(ymin2, ymax2, step2))
            # lower_ax.set_ylabel('Recovery (%)')
            # Break indicators.
            d = .5
            kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10,
                          linestyle="none", color='k', mec='k', mew=1, clip_on=False)
            upper_ax.plot(0, 0, transform=upper_ax.transAxes, **kwargs)
            lower_ax.plot(0, 1, transform=lower_ax.transAxes, **kwargs)
        # Remove the right border if it is not the last subplot.
        if i != (len(recoveries) - 1):
            lower_ax.spines['right'].set_visible(False)
            upper_ax.spines['right'].set_visible(False)
        else:
            lower_ax.set_xlim([xmin - padding, xmax])
            upper_ax.set_xlim([xmin - padding, xmax])
            upper_ax.legend(frameon=False)  # , title='Library')
            # Break indicators.
            d = .5
            kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10,
                          linestyle="none", color='k', mec='k', mew=1, clip_on=False)
            upper_ax.plot(1, 0, transform=upper_ax.transAxes, **kwargs)
            lower_ax.plot(1, 1, transform=lower_ax.transAxes, **kwargs)
    # Remove space between subplots.
    plt.subplots_adjust(wspace=0, hspace=0.05, left=0.14)
    fig.text(0.04, 0.5, 'Recovery (%)', va='center', rotation='vertical', fontsize=12)
    # Set y-scale to logarithmic and show legend.
    # plt.yscale('log')
    # plt.legend(frameon=False)  # , title='Library')
    # Show the plot.
    plt.savefig(plot_title + '.png', dpi=1200, transparent=True)
    fig.savefig(plot_title + '.eps', format='eps')
    plt.show()
