#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_inter_results(csv_file):
    # Read CSV
    df = pd.read_csv(csv_file)

    # Filter for inter (network_k > 0) data
    df_inter = df[df['network_k'] > 0].copy()

    # Convert capacity to TB
    df_inter['capacity_tb'] = df_inter['capacity'] / 1e12

    # Group size = network_m + network_k
    df_inter['group_size'] = df_inter['network_m'] + df_inter['network_k']

    # Get unique values
    group_sizes = sorted(df_inter['group_size'].unique())
    capacities = sorted(df_inter['capacity_tb'].unique())
    network_ks = sorted(df_inter['network_k'].unique())

    print(f"Group sizes: {group_sizes}")
    print(f"Capacities (TB): {capacities}")
    print(f"Network k values: {network_ks}")

    # Color map for network_k values
    colors = plt.cm.tab10(np.linspace(0, 1, len(network_ks)))
    color_map = {k: colors[i] for i, k in enumerate(network_ks)}
    markers = ['o', 's', '^', 'd', 'v', '<', '>', 'p']

    # Create figure: 3 rows (Availability, Performability, Rebuild Time) x N columns (group sizes)
    n_groups = len(group_sizes)
    fig, axes = plt.subplots(3, n_groups, figsize=(5*n_groups, 12))
    fig.suptitle('Inter-node EC: Availability, Performability & Rebuild Time by Capacity', fontsize=14, y=0.98)

    if n_groups == 1:
        axes = axes.reshape(3, 1)

    for col, group_size in enumerate(group_sizes):
        group_df = df_inter[df_inter['group_size'] == group_size]

        # Plot Availability
        ax1 = axes[0, col]
        for i, nk in enumerate(sorted(group_df['network_k'].unique())):
            subset = group_df[group_df['network_k'] == nk].sort_values('capacity_tb')
            nm_val = subset['network_m'].iloc[0] if len(subset) > 0 else 0
            marker = markers[i % len(markers)]
            ax1.plot(subset['capacity_tb'], subset['avail_nines'],
                     marker=marker, label=f'network_k={nk}', color=color_map[nk], linewidth=2, markersize=6)
        ax1.axhline(y=6, color='red', linestyle='--', alpha=0.7, label='6 nines')
        ax1.set_xlabel('Capacity (TB)')
        ax1.set_ylabel('Availability (nines)')
        ax1.set_title(f'Group Size = {group_size}')
        ax1.legend(loc='best', fontsize=8)
        ax1.grid(True, alpha=0.3)
        ax1.set_xscale('log')

        # Plot Performability (credit2_avail_nines)
        ax2 = axes[1, col]
        for i, nk in enumerate(sorted(group_df['network_k'].unique())):
            subset = group_df[group_df['network_k'] == nk].sort_values('capacity_tb')
            marker = markers[i % len(markers)]
            ax2.plot(subset['capacity_tb'], subset['credit2_avail_nines'],
                     marker=marker, label=f'network_k={nk}', color=color_map[nk], linewidth=2, markersize=6)
        ax2.axhline(y=6, color='red', linestyle='--', alpha=0.7, label='6 nines')
        ax2.set_xlabel('Capacity (TB)')
        ax2.set_ylabel('Performability (nines)')
        ax2.set_title(f'Group Size = {group_size}')
        ax2.legend(loc='best', fontsize=8)
        ax2.grid(True, alpha=0.3)
        ax2.set_xscale('log')

        # Plot Rebuild Time
        ax3 = axes[2, col]
        for i, nk in enumerate(sorted(group_df['network_k'].unique())):
            subset = group_df[group_df['network_k'] == nk].sort_values('capacity_tb')
            marker = markers[i % len(markers)]
            ax3.plot(subset['capacity_tb'], subset['avg_rebuilding_time'],
                     marker=marker, label=f'network_k={nk}', color=color_map[nk], linewidth=2, markersize=6)
        ax3.set_xlabel('Capacity (TB)')
        ax3.set_ylabel('Avg Rebuild Time (hours)')
        ax3.set_title(f'Group Size = {group_size}')
        ax3.legend(loc='best', fontsize=8)
        ax3.grid(True, alpha=0.3)
        ax3.set_xscale('log')

    plt.tight_layout()

    # Save figure
    output_file = csv_file.replace('.csv', '_inter_plot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.close()

    # Plot 2: Minimum network_k to achieve 6 nines performability
    plot_min_k_for_six_nines(df_inter, csv_file)

def plot_min_k_for_six_nines(df, csv_file):
    """Find minimum network_k to achieve 6 nines for each capacity and group size"""

    capacities = sorted(df['capacity_tb'].unique())
    group_sizes = sorted(df['group_size'].unique())

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Minimum network_k to Achieve 6 Nines (Inter-node EC)', fontsize=14)

    colors = plt.cm.Set1(np.linspace(0, 1, len(group_sizes)))

    for idx, (metric, title) in enumerate([('avail_nines', 'Availability'), ('credit2_avail_nines', 'Performability')]):
        ax = axes[idx]

        for gi, group_size in enumerate(group_sizes):
            group_df = df[df['group_size'] == group_size]
            min_k_list = []
            cap_list = []

            for cap in capacities:
                cap_df = group_df[group_df['capacity_tb'] == cap]
                # Find minimum network_k that achieves >= 6 nines
                passing = cap_df[cap_df[metric] >= 6.0]
                if len(passing) > 0:
                    min_k = passing['network_k'].min()
                else:
                    min_k = np.nan  # Not achievable
                min_k_list.append(min_k)
                cap_list.append(cap)

            ax.plot(cap_list, min_k_list, marker='o', label=f'Group={group_size}',
                    linewidth=2, color=colors[gi], markersize=8)

        ax.set_xlabel('Capacity (TB)')
        ax.set_ylabel('Minimum network_k for 6 nines')
        ax.set_title(f'{title}: Min network_k to Achieve 6 Nines')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xscale('log')
        if df['network_k'].max() > 0:
            ax.set_ylim(0, df['network_k'].max() + 1)

    plt.tight_layout()

    output_file = csv_file.replace('.csv', '_min_k_plot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Min-k plot saved to: {output_file}")
    plt.close()

    # Print table of minimum k values
    print("\n=== Minimum network_k to achieve 6 nines ===")
    print("\nPerformability (credit2_avail_nines):")
    print(f"{'Capacity (TB)':<15}", end="")
    for gs in group_sizes:
        print(f"{'Group='+str(gs):<12}", end="")
    print()

    for cap in capacities:
        print(f"{cap:<15.0f}", end="")
        for gs in group_sizes:
            subset = df[(df['group_size'] == gs) & (df['capacity_tb'] == cap)]
            passing = subset[subset['credit2_avail_nines'] >= 6.0]
            if len(passing) > 0:
                min_k = passing['network_k'].min()
                print(f"{min_k:<12}", end="")
            else:
                print(f"{'N/A':<12}", end="")
        print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_inter_results.py <csv_file>")
        sys.exit(1)

    plot_inter_results(sys.argv[1])
