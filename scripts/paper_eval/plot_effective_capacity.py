#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_effective_capacity(csv_file):
    df = pd.read_csv(csv_file)
    df = df[df['network_k'] > 0].copy()

    # Calculate effective capacity ratio: network_m / (network_m + network_k)
    df['eff_cap_ratio'] = df['network_m'] / (df['network_m'] + df['network_k'])
    df['capacity_tb'] = df['capacity'] / 1e12

    capacities = sorted(df['capacity_tb'].unique())

    # For each unique (eff_cap_ratio, capacity), find max target achieving 6 nines
    results = []
    for _, row in df.drop_duplicates(subset=['eff_cap_ratio', 'capacity']).iterrows():
        ecr = row['eff_cap_ratio']
        cap = row['capacity']

        subset = df[(df['eff_cap_ratio'] == ecr) & (df['capacity'] == cap)]

        passing = subset[subset['credit2_avail_nines'] >= 6.0]
        max_target = passing['target_perf_ratio'].max() if len(passing) > 0 else 0.0

        results.append({
            'eff_cap_ratio': ecr,
            'capacity_tb': cap / 1e12,
            'network_m': row['network_m'],
            'network_k': row['network_k'],
            'max_target_ratio': max_target
        })

    result_df = pd.DataFrame(results)

    # Plot: X = eff_cap_ratio, Y = max target ratio, color = capacity
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, len(capacities)))
    markers = ['o', 's', '^', 'd', 'v', 'p']

    for i, cap in enumerate(capacities):
        cap_df = result_df[result_df['capacity_tb'] == cap].sort_values('eff_cap_ratio')
        ax.plot(cap_df['eff_cap_ratio'], cap_df['max_target_ratio'],
               marker=markers[i % len(markers)], label=f'{cap:.0f} TB',
               color=colors[i], linewidth=2, markersize=8)

    ax.set_xlabel('Effective Capacity Ratio (m / (m+k))', fontsize=12)
    ax.set_ylabel('Max target_perf_ratio for 6 Nines', fontsize=12)
    ax.set_title('Max Performance Target Achieving 6 Nines Performability', fontsize=14)
    ax.legend(title='Capacity', loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.4, 1.0)
    ax.set_ylim(0.4, 1.0)
    ax.axhline(y=0.9, color='red', linestyle='--', alpha=0.5)

    plt.tight_layout()
    output_file = csv_file.replace('.csv', '_eff_ratio_plot.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.close()

    # Print summary table
    print("\n=== Max target_perf_ratio for 6 nines ===")
    print(f"{'Ratio':<8} {'m':<4} {'k':<4}", end="")
    for cap in capacities:
        print(f"{cap:>8.0f}TB", end="")
    print()

    for ecr in sorted(result_df['eff_cap_ratio'].unique(), reverse=True):
        row_df = result_df[result_df['eff_cap_ratio'] == ecr].iloc[0]
        print(f"{ecr:<8.3f} {int(row_df['network_m']):<4} {int(row_df['network_k']):<4}", end="")
        for cap in capacities:
            val = result_df[(result_df['eff_cap_ratio'] == ecr) &
                           (result_df['capacity_tb'] == cap)]['max_target_ratio'].values
            if len(val) > 0 and val[0] > 0:
                print(f"{val[0]:>8.2f}", end="")
            else:
                print(f"{'N/A':>8}", end="")
        print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_effective_capacity.py <sweep_csv_file>")
        sys.exit(1)
    plot_effective_capacity(sys.argv[1])
