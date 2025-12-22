#!/bin/bash
script_dir="$(realpath "$(dirname "${BASH_SOURCE[0]}")/../..")"


seed=42
output_file="results/result_dual_single_$(date '+%Y%m%d_%H%M%S').csv"
data=48 # m + k + l
parities=(1 2 3 4 5 6) # --k
l=0 # --l
capacity=64000000000000
dwpd=0.2
total_ssds=48
tier_files=("2tier.json" "2tier-singleport.json" "3tier.json" "3tier-singleport.json" "2tier-singleport-dualio.json" "3tier-singleport-dualio.json")
simulator="$script_dir/build/availability_sim"

for t in "${tier_files[@]}"; do
    single_port=""
    if [[ "$t" == *"singleport"* ]]; then
        echo "Running for single port tier: $t"
        single_port="--single_port_ssd"
    else
        echo "Running for dual port tier: $t"
    fi
    for p in "${parities[@]}"; do
        m=$(($data - $p - $l))
        $simulator --output_file $output_file --m $m --k $p --l $l --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc $single_port --seed $seed 
    done
done

for t in "${tier_files[@]}"; do
    single_port=""
    if [[ "$t" == *"singleport"* ]]; then
        echo "Running for single port tier: $t"
        single_port="--single_port_ssd"
    else
        echo "Running for dual port tier: $t"
    fi
    m=$(($data - "0"))
    $simulator --no_result --output_file $output_file --m 1 --k 0 --l 0 --network_m $m --network_k 0 --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc $single_port --nprocs 20 --seed $seed 
    for p in "${parities[@]}"; do
        m=$(($data - $p))
        $simulator --output_file $output_file --m 1 --k 0 --l 0 --network_m $m --network_k $p --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc $single_port --nprocs 20 --seed $seed 
    done
done
