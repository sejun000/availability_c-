#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_max_target_ratio(csv_file):
    # Read CSV
    df = pd.read_csv(csv_file)

    # Filter for inter (network_k > 0) data
    df = df[df['network_k'] > 0].copy()

    # Convert capacity to TB
    df['capacity_tb'] = df['capacity'] / 1e12

    # Group size = network_m + network_k
    df['group_size'] = df['network_m'] + df['network_k']

    # Get unique values
    group_sizes = sorted(df['group_size'].unique())
    capacities = sorted(df['capacity_tb'].unique())
    network_ks = sorted(df['network_k'].unique())

    print(f"Group sizes: {group_sizes}")
    print(f"Capacities (TB): {capacities}")
    print(f"Network k values: {network_ks}")
    print(f"Total rows: {len(df)}")

    # For each (network_k, capacity, group_size), find max target_perf_ratio achieving 6 nines
    results = []
    for gs in group_sizes:
        for cap in capacities:
            for nk in network_ks:
                subset = df[(df['group_size'] == gs) &
                           (df['capacity_tb'] == cap) &
                           (df['network_k'] == nk)]

                if len(subset) == 0:
                    continue

                # Find max target_perf_ratio where credit_avail_nines >= 6
                passing = subset[subset['credit_avail_nines'] >= 6.0]
                if len(passing) > 0:
                    max_target = passing['target_perf_ratio'].max()
                else:
                    max_target = 0.0  # Can't achieve 6 nines at any target

                results.append({
                    'group_size': gs,
                    'capacity_tb': cap,
                    'network_k': nk,
                    'max_target_ratio': max_target
                })

    result_df = pd.DataFrame(results)
    print("\n=== Max target_perf_ratio for 6 nines performability ===")
    print(result_df.to_string(index=False))

    # Create plots
    n_groups = len(group_sizes)
    fig, axes = plt.subplots(1, n_groups, figsize=(5*n_groups, 5))
    fig.suptitle('Max target_perf_ratio to Achieve 6 Nines Performability', fontsize=14)

    if n_groups == 1:
        axes = [axes]

    colors = plt.cm.tab10(np.linspace(0, 1, len(network_ks)))
    color_map = {k: colors[i] for i, k in enumerate(network_ks)}
    markers = ['o', 's', '^', 'd']

    for col, gs in enumerate(group_sizes):
        ax = axes[col]
        gs_df = result_df[result_df['group_size'] == gs]

        for i, nk in enumerate(network_ks):
            nk_df = gs_df[gs_df['network_k'] == nk].sort_values('capacity_tb')
            if len(nk_df) > 0:
                ax.plot(nk_df['capacity_tb'], nk_df['max_target_ratio'],
                       marker=markers[i % len(markers)], label=f'network_k={nk}',
                       color=color_map[nk], linewidth=2, markersize=8)

        ax.set_xlabel('Capacity (TB)')
        ax.set_ylabel('Max target_perf_ratio for 6 nines')
        ax.set_title(f'Group Size = {gs}')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xscale('log')
        ax.set_ylim(0.35, 1.0)
        ax.axhline(y=0.9, color='gray', linestyle='--', alpha=0.5, label='0.9')

    plt.tight_layout()

    output_file = csv_file.replace('.csv', '_max_target_plot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    plt.close()

    # Also create a combined view by capacity
    fig2, axes2 = plt.subplots(1, len(capacities), figsize=(5*len(capacities), 5))
    fig2.suptitle('Max target_perf_ratio for 6 Nines by Capacity (Higher = Better)', fontsize=14)

    if len(capacities) == 1:
        axes2 = [axes2]

    for col, cap in enumerate(capacities):
        ax = axes2[col]
        cap_df = result_df[result_df['capacity_tb'] == cap]

        x = np.arange(len(network_ks))
        width = 0.2

        for i, gs in enumerate(group_sizes):
            gs_cap_df = cap_df[cap_df['group_size'] == gs]
            values = []
            for nk in network_ks:
                v = gs_cap_df[gs_cap_df['network_k'] == nk]['max_target_ratio'].values
                values.append(v[0] if len(v) > 0 else 0)
            ax.bar(x + i*width, values, width, label=f'Group={gs}')

        ax.set_xlabel('network_k')
        ax.set_ylabel('Max target_perf_ratio')
        ax.set_title(f'Capacity = {cap:.0f} TB')
        ax.set_xticks(x + width * (len(group_sizes)-1) / 2)
        ax.set_xticklabels(network_ks)
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_ylim(0, 1.0)

    plt.tight_layout()

    output_file2 = csv_file.replace('.csv', '_max_target_bar_plot.png')
    plt.savefig(output_file2, dpi=150, bbox_inches='tight')
    print(f"Bar plot saved to: {output_file2}")
    plt.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_max_target_ratio.py <sweep_csv_file>")
        sys.exit(1)

    plot_max_target_ratio(sys.argv[1])
