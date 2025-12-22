#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_capacity_results(csv_file):
    # Read CSV
    df = pd.read_csv(csv_file)

    # Convert capacity to TB
    df['capacity_tb'] = df['capacity'] / 1e12
    df['group_size'] = df['m'] + df['k']

    # Get unique group sizes
    group_sizes = sorted(df['group_size'].unique())

    # Create figure with subplots: 2 rows (Availability, Performability) x N columns (group sizes)
    n_groups = len(group_sizes)
    fig, axes = plt.subplots(2, n_groups, figsize=(5*n_groups, 8))
    fig.suptitle('Availability & Performability by Capacity and EC Configuration', fontsize=14)

    if n_groups == 1:
        axes = axes.reshape(2, 1)

    # Color map for k values
    k_values = sorted(df['k'].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(k_values)))
    color_map = {k: colors[i] for i, k in enumerate(k_values)}
    markers = ['o', 's', '^', 'd', 'v', '<', '>', 'p', '*', 'h']

    for col, group_size in enumerate(group_sizes):
        group_df = df[df['group_size'] == group_size]

        # Plot Availability
        ax1 = axes[0, col]
        for k in sorted(group_df['k'].unique()):
            subset = group_df[group_df['k'] == k].sort_values('capacity_tb')
            m_val = subset['m'].iloc[0]
            marker = markers[k % len(markers)]
            ax1.plot(subset['capacity_tb'], subset['avail_nines'],
                     marker=marker, label=f'm={m_val}, k={k}', color=color_map[k], linewidth=2)
        ax1.axhline(y=6, color='red', linestyle='--', alpha=0.5, label='6 nines')
        ax1.set_xlabel('Capacity (TB)')
        ax1.set_ylabel('Availability (nines)')
        ax1.set_title(f'Group Size = {group_size}')
        ax1.legend(loc='best', fontsize=8)
        ax1.grid(True, alpha=0.3)
        ax1.set_xscale('log')

        # Plot Performability
        ax2 = axes[1, col]
        for k in sorted(group_df['k'].unique()):
            subset = group_df[group_df['k'] == k].sort_values('capacity_tb')
            m_val = subset['m'].iloc[0]
            marker = markers[k % len(markers)]
            ax2.plot(subset['capacity_tb'], subset['credit_avail_nines'],
                     marker=marker, label=f'm={m_val}, k={k}', color=color_map[k], linewidth=2)
        ax2.axhline(y=6, color='red', linestyle='--', alpha=0.5, label='6 nines')
        ax2.set_xlabel('Capacity (TB)')
        ax2.set_ylabel('Performability (nines)')
        ax2.set_title(f'Group Size = {group_size}')
        ax2.legend(loc='best', fontsize=8)
        ax2.grid(True, alpha=0.3)
        ax2.set_xscale('log')

    plt.tight_layout()

    # Save figure
    output_file = csv_file.replace('.csv', '_plot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")

    # Plot 2: Minimum k to achieve 6 nines
    plot_min_k_for_six_nines(df, csv_file)

def plot_min_k_for_six_nines(df, csv_file):
    """Find minimum k to achieve 6 nines for each capacity and group size"""

    capacities = sorted(df['capacity_tb'].unique())
    group_sizes = sorted(df['group_size'].unique())

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Minimum k to Achieve 6 Nines', fontsize=14)

    for idx, (metric, title) in enumerate([('avail_nines', 'Availability'), ('credit_avail_nines', 'Performability')]):
        ax = axes[idx]

        for group_size in group_sizes:
            group_df = df[df['group_size'] == group_size]
            min_k_list = []
            cap_list = []

            for cap in capacities:
                cap_df = group_df[group_df['capacity_tb'] == cap]
                # Find minimum k that achieves >= 6 nines
                passing = cap_df[cap_df[metric] >= 6.0]
                if len(passing) > 0:
                    min_k = passing['k'].min()
                else:
                    min_k = np.nan  # Not achievable
                min_k_list.append(min_k)
                cap_list.append(cap)

            ax.plot(cap_list, min_k_list, marker='o', label=f'Group={group_size}', linewidth=2)

        ax.set_xlabel('Capacity (TB)')
        ax.set_ylabel('Minimum k for 6 nines')
        ax.set_title(f'{title}: Min k to Achieve 6 Nines')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xscale('log')
        ax.set_ylim(0, max(df['k']) + 1)

    plt.tight_layout()

    output_file = csv_file.replace('.csv', '_min_k_plot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Min-k plot saved to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_capacity.py <csv_file>")
        sys.exit(1)

    plot_capacity_results(sys.argv[1])
