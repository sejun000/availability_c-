#!/bin/bash
script_dir="$(realpath "$(dirname "${BASH_SOURCE[0]}")/../..")"


seed=42
output_file="results/result_tiering_dual_single_$(date '+%Y%m%d_%H%M%S').csv"
data=48 # m + k + l
parities=(1 2 3 4 5 6) # --k
l=0 # --l
capacity=64000000000000
dwpd=0.2
total_ssds=48
tier_files=("2tier.json" "2tier-singleport.json" "3tier.json" "3tier-singleport.json")
io_module_mttrs=(0.5 0.5 0.5 0.5)
simulator="$script_dir/build/availability_sim"

i=0
for t in "${tier_files[@]}"; do
    io_module_mttr=${io_module_mttrs[$i]}
    i=$(($i + 1))
    single_port=""
    if [[ "$t" == *"singleport"* ]]; then
        echo "Running for single port tier: $t"
        single_port="--single_port_ssd"
    else
        echo "Running for dual port tier: $t"
    fi
    for p in "${parities[@]}"; do
        m=$(($data - $p - $l))
        $simulator --output_file $output_file --m $m --k $p --l $l --capacity $capacity --config_file $script_dir/json/$t --total_ssds $total_ssds --dwpd $dwpd --qlc $single_port --io_module_mttr $io_module_mttr --seed $seed 
    done
done
