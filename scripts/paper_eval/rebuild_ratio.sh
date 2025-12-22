#!/bin/bash
script_dir="$(realpath "$(dirname "${BASH_SOURCE[0]}")/../..")"

seed=42
timestamp=$(date '+%Y%m%d_%H%M%S')
output_file="results/result_rebuild_ratio_${timestamp}.csv"
sweep_csv="results/result_rebuild_ratio_sweep_${timestamp}.csv"

# Fixed parameters
capacity=64000000000000  # 64TB
dwpd=0.2
total_ssds=48
tier_files=("2tier.json")
simulator="$script_dir/build/availability_sim"

# Sweep parameters
rebuild_ratios=(0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
data=(24)  # m + k + l for local parity
parities=(1 2 3 4 5 6)  # --k for local parity
l=0
network_data=(24)  # network_m + network_k
network_parities=(1 2 3 4 5 6)  # network_k

# Local parity (m, k) sweep with varying rebuild_bw_ratio
for t in "${tier_files[@]}"; do
    for rb in "${rebuild_ratios[@]}"; do
        for d in "${data[@]}"; do
            for p in "${parities[@]}"; do
                m=$(($d - $p - $l))
                if [ $m -lt 1 ]; then
                    continue
                fi
                cmd="$simulator --output_file $output_file --m $m --k $p --l $l --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed --rebuild_bw_ratio $rb credit_availability_sweep_csv $sweep_csv"
                echo "========================================"
                echo "Running local: rebuild_ratio=$rb, m=$m, k=$p"
                echo "========================================"
                $cmd
                if [ $? -ne 0 ]; then
                    echo "ERROR: Command failed with exit code $?"
                fi
            done
        done
    done
done

# Network parity (network_m, network_k) sweep with varying rebuild_bw_ratio
for t in "${tier_files[@]}"; do
    for rb in "${rebuild_ratios[@]}"; do
        for nd in "${network_data[@]}"; do
            # First run baseline with network_k=0
            network_m=$nd
            cmd="$simulator --no_result --output_file $output_file --m 1 --k 0 --l 0 --network_m $network_m --network_k 0 --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed --rebuild_bw_ratio $rb"
            echo "========================================"
            echo "Running baseline: rebuild_ratio=$rb, network_m=$network_m"
            echo "========================================"
            $cmd
            if [ $? -ne 0 ]; then
                echo "ERROR: Command failed with exit code $?"
            fi

            # Then run with network_k variations
            for np in "${network_parities[@]}"; do
                network_m=$(($nd - $np))
                if [ $network_m -lt 1 ]; then
                    continue
                fi
                cmd="$simulator --output_file $output_file --m 1 --k 0 --l 0 --network_m $network_m --network_k $np --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed --rebuild_bw_ratio $rb --credit_availability_sweep_csv $sweep_csv"
                echo "========================================"
                echo "Running: rebuild_ratio=$rb, network_m=$network_m, network_k=$np"
                echo "========================================"
                $cmd
                if [ $? -ne 0 ]; then
                    echo "ERROR: Command failed with exit code $?"
                fi
            done
        done
    done
done
