#!/bin/bash
script_dir="$(realpath "$(dirname "${BASH_SOURCE[0]}")/../..")"


seed=42
output_file="results/result_clustered_declustered_parity_$(date '+%Y%m%d_%H%M%S').csv"
data=48 # m + k + l
parities=(1 2 3 4 5) # --k
ls=(0 4 8 12 16 20 24 28 32) # --l
capacity=64000000000000
dwpd=0.2
total_ssds=48
tier_files=("2tier.json")
simulator="$script_dir/build/availability_sim"

for t in "${tier_files[@]}"; do
    for p in "${parities[@]}"; do
        for l in "${ls[@]}"; do
            m=$(($data - $p - $l))
            $simulator --output_file $output_file --m $m --k $p --l $l --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed 
        done
    done
done
