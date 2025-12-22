#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_rebuild_ratio_results(csv_file):
    # Read CSV
    df = pd.read_csv(csv_file)

    # ===== Intra (Local parity) plot =====
    df_local = df[df['network_k'] == 0].copy()

    if len(df_local) > 0:
        df_local['group_size'] = df_local['m'] + df_local['k']
        df_local['avg_rebuild_bw_mbps'] = df_local['avg_rebuild_bw'] / (1024 * 1024)

        k_values = sorted(df_local['k'].unique())
        colors = plt.cm.tab10(np.linspace(0, 1, len(k_values)))
        color_map = {k: colors[i] for i, k in enumerate(k_values)}

        fig, ax = plt.subplots(1, 1, figsize=(10, 7))

        for k in k_values:
            df_k = df_local[df_local['k'] == k].sort_values('avg_rebuild_bw_mbps')
            ax.plot(df_k['avg_rebuild_bw_mbps'], df_k['effective_avail_nines'],
                   marker='o', color=color_map[k], label=f'k={k}',
                   alpha=0.7, markersize=8, linewidth=2)

        ax.set_xlabel('Avg Rebuild BW (MB/s)')
        ax.set_ylabel('Normalized Max Bandwidth (nines)')
        ax.set_title('Intra (Local Parity): Rebuild BW vs Normalized Max BW')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, loc='best')

        plt.tight_layout()
        output_file = csv_file.replace('.csv', '_intra_bw_vs_avail.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Intra plot saved to: {output_file}")
        plt.close()
    else:
        print("No intra (local parity) data found")

    # ===== Inter (Network parity) plot =====
    df_network = df[df['network_k'] > 0].copy()

    if len(df_network) > 0:
        df_network['group_size'] = df_network['network_m'] + df_network['network_k']
        df_network['avg_rebuild_bw_mbps'] = df_network['avg_rebuild_bw'] / (1024 * 1024)

        nk_values = sorted(df_network['network_k'].unique())
        colors_net = plt.cm.tab10(np.linspace(0, 1, len(nk_values)))
        color_map_net = {nk: colors_net[i] for i, nk in enumerate(nk_values)}

        fig2, ax2 = plt.subplots(1, 1, figsize=(10, 7))

        for nk in nk_values:
            df_nk = df_network[df_network['network_k'] == nk].sort_values('avg_rebuild_bw_mbps')
            ax2.plot(df_nk['avg_rebuild_bw_mbps'], df_nk['effective_avail_nines'],
                    marker='o', color=color_map_net[nk], label=f'network_k={nk}',
                    alpha=0.7, markersize=8, linewidth=2)

        ax2.set_xlabel('Avg Rebuild BW (MB/s)')
        ax2.set_ylabel('Normalized Max Bandwidth (nines)')
        ax2.set_title('Inter (Network Parity): Rebuild BW vs Normalized Max BW')
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=9, loc='best')

        plt.tight_layout()
        output_file2 = csv_file.replace('.csv', '_inter_bw_vs_avail.png')
        plt.savefig(output_file2, dpi=150, bbox_inches='tight')
        print(f"Inter plot saved to: {output_file2}")
        plt.close()
    else:
        print("No inter (network parity) data found")

    # ===== Combined plot (side by side) =====
    if len(df_local) > 0 and len(df_network) > 0:
        fig3, (ax3, ax4) = plt.subplots(1, 2, figsize=(16, 6))

        # Intra
        for k in k_values:
            df_k = df_local[df_local['k'] == k].sort_values('avg_rebuild_bw_mbps')
            ax3.plot(df_k['avg_rebuild_bw_mbps'], df_k['effective_avail_nines'],
                    marker='o', color=color_map[k], label=f'k={k}',
                    alpha=0.7, markersize=8, linewidth=2)
        ax3.set_xlabel('Avg Rebuild BW (MB/s)')
        ax3.set_ylabel('Normalized Max Bandwidth (nines)')
        ax3.set_title('Intra (Local Parity)')
        ax3.grid(True, alpha=0.3)
        ax3.legend(fontsize=9, loc='best')

        # Inter
        for nk in nk_values:
            df_nk = df_network[df_network['network_k'] == nk].sort_values('avg_rebuild_bw_mbps')
            ax4.plot(df_nk['avg_rebuild_bw_mbps'], df_nk['effective_avail_nines'],
                    marker='o', color=color_map_net[nk], label=f'network_k={nk}',
                    alpha=0.7, markersize=8, linewidth=2)
        ax4.set_xlabel('Avg Rebuild BW (MB/s)')
        ax4.set_ylabel('Normalized Max Bandwidth (nines)')
        ax4.set_title('Inter (Network Parity)')
        ax4.grid(True, alpha=0.3)
        ax4.legend(fontsize=9, loc='best')

        plt.suptitle('Rebuild BW vs Normalized Max Bandwidth', fontsize=14)
        plt.tight_layout()
        output_file3 = csv_file.replace('.csv', '_combined_bw_vs_avail.png')
        plt.savefig(output_file3, dpi=150, bbox_inches='tight')
        print(f"Combined plot saved to: {output_file3}")
        plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_rebuild_ratio.py <csv_file>")
        sys.exit(1)

    plot_rebuild_ratio_results(sys.argv[1])
