#!/bin/bash
script_dir="$(realpath "$(dirname "${BASH_SOURCE[0]}")/../..")"


seed=42
output_file="results/result_cache_$(date '+%Y%m%d_%H%M%S').csv"
p=2 # --k
rgroups=(1)
cached_ssds=(0 2 4 6 8 10 12)
cache_hit_ratios=(0 0.4151 0.438 0.4609 0.5289 0.592 0.6147)
capacity=64000000000000
dwpds=(0.02 0.066 0.2 0.66 2)
total_ssds=48
tier_files=("2tier.json")
simulator="$script_dir/build/availability_sim"

for t in "${tier_files[@]}"; do
    for rgroup in "${rgroups[@]}"; do
        for dwpd in "${dwpds[@]}"; do
            cs_index=0
            for cs in "${cached_ssds[@]}"; do
                cached_hit_ratio=${cache_hit_ratios[$cs_index]}
                cs_index=$(($cs_index + 1))

                cold_ssds=$(($total_ssds - $cs))
                m=$(($cold_ssds / $rgroup - $p))
                if [ "$m" -lt 1 ]; then
                    continue
                fi
                $simulator --output_file $output_file --m $m --k $p --intra_replicas 2 --cached_ssds $cs --cached_write_ratio $cached_hit_ratio --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed 
            done
            # all tlc configuration
            m=$(($total_ssds / $rgroup - $p))
            $simulator --output_file $output_file --m $m --k $p --intra_replicas 2 --cached_write_ratio 0 --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --seed $seed 
        done
    done
done
