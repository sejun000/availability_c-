#!/bin/bash
script_dir="$(realpath "$(dirname "${BASH_SOURCE[0]}")/../..")"

seed=42
timestamp=$(date '+%Y%m%d_%H%M%S')
output_file="results/result_capacity_${timestamp}.csv"
sweep_csv="results/result_capacity_sweep_${timestamp}.csv"
data=(8 16 24 48) # m + k + l
parities=(1 2 3 4 5 6) # --k
l=0 # --l
capacities=(4000000000000 16000000000000 64000000000000 256000000000000)
dwpd=0.2
total_ssds=48
tier_files=("2tier.json")
simulator="$script_dir/build/availability_sim"

# Local parity (m, k) sweep
for t in "${tier_files[@]}"; do
    for d in "${data[@]}"; do
        for p in "${parities[@]}"; do
            for c in "${capacities[@]}"; do
                m=$(($d - $p - $l))
                cmd="$simulator --output_file $output_file --m $m --k $p --l $l --capacity $c --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed"
                echo "========================================"
                echo "Running: $cmd"
                echo "========================================"
                $cmd
                if [ $? -ne 0 ]; then
                    echo "ERROR: Command failed with exit code $?"
                    echo "Failed command: $cmd"
                fi
            done
        done
    done
done

# Network parity (network_m, network_k) sweep
# m=1, k=0 for local, vary network_m and network_k
network_data=(8 16 24 48) # network_m + network_k
network_parities=(1 2 3 4 5 6) # network_k

for t in "${tier_files[@]}"; do
    for c in "${capacities[@]}"; do
        for nd in "${network_data[@]}"; do
            # First run baseline with network_k=0 (--no_result for staged simulation)
            network_m=$nd
            cmd="$simulator --no_result --output_file $output_file --m 1 --k 0 --l 0 --network_m $network_m --network_k 0 --capacity $c --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed"
            echo "========================================"
            echo "Running baseline (staged): $cmd"
            echo "========================================"
            $cmd
            if [ $? -ne 0 ]; then
                echo "ERROR: Command failed with exit code $?"
                echo "Failed command: $cmd"
            fi

            # Then run with network_k variations
            for np in "${network_parities[@]}"; do
                network_m=$(($nd - $np))
                if [ $network_m -lt 1 ]; then
                    continue
                fi
                cmd="$simulator --output_file $output_file --m 1 --k 0 --l 0 --network_m $network_m --network_k $np --capacity $c --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed --credit_availability_sweep_csv $sweep_csv"
                echo "========================================"
                echo "Running: $cmd"
                echo "========================================"
                $cmd
                if [ $? -ne 0 ]; then
                    echo "ERROR: Command failed with exit code $?"
                    echo "Failed command: $cmd"
                fi
            done
        done
    done
done
