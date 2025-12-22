#!/bin/bash
script_dir="$(realpath "$(dirname "${BASH_SOURCE[0]}")/../..")"


seed=42
output_file="results/result_stripe_$(date '+%Y%m%d_%H%M%S').csv"
data=48 # m + k + l
parities=(1 2 3 4 5 6) # --k
l=0 # --l
stripes=(8 16 24 48 64)
capacity=64000000000000
dwpd=0.2
total_ssds=48
tier_files=("2tier.json")
simulator="$script_dir/build/availability_sim"

for t in "${tier_files[@]}"; do
    for s in "${stripes[@]}"; do
        m=$(($s - "0"))
        $simulator --no_result --output_file $output_file --m 1 --k 0 --l 0 --network_m $m --network_k 0 --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc $single_port --seed $seed 
        for p in "${parities[@]}"; do
            m=$(($s - $p))
            $simulator --output_file $output_file --m 1 --k 0 --l 0 --network_m $m --network_k $p --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc $single_port --seed $seed 
        done
    done
done


parities=(1 2 3 4 5) # --k
l=0 # --l
stripes=(8 16 24 48)

for t in "${tier_files[@]}"; do
    for s in "${stripes[@]}"; do
        for p in "${parities[@]}"; do
            m=$(($s - $p - $l))
            $simulator --output_file $output_file --m $m --k $p --l $l --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc --seed $seed 
        done
    done
done
