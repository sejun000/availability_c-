#!/bin/bash
# Phase 1: RQ1 (IO Redundancy) + RQ3 (EC Comparison)
# 예상 소요 시간: ~5.1시간 (306 실험 × 1분)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/build"
CONFIG_DIR="$PROJECT_DIR/configs"
RESULT_DIR="$PROJECT_DIR/results/phase1"

mkdir -p "$RESULT_DIR"

# 공통 옵션
BASE_OPTS="--sim_years 5 --n 32 --total_disks 384 --disk_failure_trace $PROJECT_DIR/ssd_failure_times.txt --seed 42"
PERF_RATIOS="0.50 0.60 0.70 0.80 0.90 0.99"

# 결과 파일 초기화
RESULT_FILE="$RESULT_DIR/phase1_results.csv"
echo "rq,config,ec_type,m,k,local_m,local_k,local_n,storage_eff,target_perf,rebuild_ratio,availability,perf_availability,durability,avg_rebuild_time,max_read_bw,max_write_bw,avg_read_bw,avg_write_bw" > "$RESULT_FILE"

# 실험 카운터
total_experiments=306
current=0

run_experiment() {
    local rq=$1
    local config=$2
    local ec_type=$3
    local m=$4
    local k=$5
    local local_m=$6
    local local_k=$7
    local local_n=$8
    local storage_eff=$9
    local perf=${10}
    local rebuild=${11}
    local simulations=${12}
    local extra_opts="${13}"

    current=$((current + 1))
    echo "[$current/$total_experiments] RQ=$rq, EC=$ec_type, m=$m, k=$k, perf=$perf"
    echo "CMD: $BUILD_DIR/availability_sim $CONFIG_DIR/$config $BASE_OPTS --num_sim $simulations --rebuild_bw_ratio $rebuild --degraded_ratio $rebuild --target_perf $perf $extra_opts"

    local output=$("$BUILD_DIR/availability_sim" "$CONFIG_DIR/$config" \
        $BASE_OPTS \
        --num_sim $simulations \
        --rebuild_bw_ratio $rebuild \
        --degraded_ratio $rebuild \
        --target_perf $perf \
        $extra_opts 2>&1)

    # 결과 파싱
    local avail=$(echo "$output" | grep '^  Availability:' | head -1 | awk '{print $2}')
    local perf_avail=$(echo "$output" | grep 'Perf Availability MIN' | head -1 | awk '{print $NF}')
    local durability=$(echo "$output" | grep '^  Durability:' | head -1 | awk '{print $2}')
    local rebuild_time=$(echo "$output" | grep 'Avg rebuild time:' | awk '{print $4}')
    local max_read_bw=$(echo "$output" | grep 'Baseline max READ BW:' | awk '{print $5}')
    local max_write_bw=$(echo "$output" | grep 'Baseline max WRITE BW:' | awk '{print $5}')
    local avg_read_bw=$(echo "$output" | grep 'Avg READ BW (overall):' | awk '{print $5}')
    local avg_write_bw=$(echo "$output" | grep 'Avg WRITE BW (overall):' | awk '{print $5}')

    echo "$rq,$config,$ec_type,$m,$k,$local_m,$local_k,$local_n,$storage_eff,$perf,$rebuild,$avail,$perf_avail,$durability,$rebuild_time,$max_read_bw,$max_write_bw,$avg_read_bw,$avg_write_bw" >> "$RESULT_FILE"
}

echo "=========================================="
echo "Phase 1: Starting experiments"
echo "=========================================="

# ===========================================
# RQ1: IO Module Redundancy (36 experiments)
# k = 1, 2, 3 × perf × config (rebuild=0.10 고정)
# ===========================================
echo ""
echo "=== RQ1: IO Module Redundancy ==="
echo ""

for k in 1 2 3; do
    m=$((24 - k))
    storage_eff=$(echo "scale=3; $m / ($m + $k)" | bc)
    for perf in $PERF_RATIOS; do
        # 1 IO module (1000 simulations)
        run_experiment "RQ1" "16enclosure_distributed_1io.json" "standard" $m $k "" "" "" "$storage_eff" $perf 0.10 1000 "--m $m --k $k"

        # 2 IO modules (10000 simulations)
        run_experiment "RQ1" "16enclosure_distributed_2io.json" "standard" $m $k "" "" "" "$storage_eff" $perf 0.10 10000 "--m $m --k $k"
    done
done

# ===========================================
# RQ3: Erasure Coding Comparison (288 experiments)
# ===========================================
echo ""
echo "=== RQ3: Erasure Coding Comparison ==="
echo ""

# --- Standard EC (k = 1, 2, 3 × m values) ---
# k=1: m=23,16 / k=2: m=22,16 / k=3: m=21,16
echo "--- Standard EC (6 configs × 6 perf = 36 experiments) ---"

# k=1
for m in 23 16; do
    k=1
    storage_eff=$(echo "scale=3; $m / ($m + $k)" | bc)
    for perf in $PERF_RATIOS; do
        run_experiment "RQ3" "16enclosure_distributed_2io.json" "standard" $m $k "" "" "" "$storage_eff" $perf 0.10 10000 "--m $m --k $k"
    done
done

# k=2
for m in 22 16; do
    k=2
    storage_eff=$(echo "scale=3; $m / ($m + $k)" | bc)
    for perf in $PERF_RATIOS; do
        run_experiment "RQ3" "16enclosure_distributed_2io.json" "standard" $m $k "" "" "" "$storage_eff" $perf 0.10 10000 "--m $m --k $k"
    done
done

# k=3
for m in 21 16; do
    k=3
    storage_eff=$(echo "scale=3; $m / ($m + $k)" | bc)
    for perf in $PERF_RATIOS; do
        run_experiment "RQ3" "16enclosure_distributed_2io.json" "standard" $m $k "" "" "" "$storage_eff" $perf 0.10 10000 "--m $m --k $k"
    done
done

# --- LRC (4 configs from test.md table) ---
echo "--- LRC (4 configs × 6 perf = 24 experiments) ---"

# m=24, k=2, local_m=8, local_k=2
run_lrc() {
    local m=$1
    local k=$2
    local local_m=$3
    local local_k=$4

    num_groups=$(echo "$m / $local_m" | bc)
    total_parity=$(echo "$k + $num_groups * $local_k" | bc)
    storage_eff=$(echo "scale=3; $m / ($m + $total_parity)" | bc)

    for perf in $PERF_RATIOS; do
        run_experiment "RQ3" "16enclosure_distributed_2io.json" "lrc" $m $k $local_m $local_k "" "$storage_eff" $perf 0.10 10000 "--ec_type lrc --m $m --k $k --local_m $local_m --local_k $local_k"
    done
}

run_lrc 24 2 8 2
run_lrc 27 2 9 1
run_lrc 26 2 13 2
run_lrc 28 2 14 1

# --- Multi-EC (outer_m = 12, 20 × outer_k = 1, 2 × local_k = 1, 2) ---
# local_n = 24 (고정), outer_m=12: local_m=10, outer_m=20: local_m=20
echo "--- Multi-EC (8 configs × 6 perf = 48 experiments) ---"

local_n=24  # 고정 (enclosure당 디스크 개수)

for outer_m in 12 20; do
    # Set local_m based on outer_m
    if [ $outer_m -eq 12 ]; then
        local_m=10
    else
        local_m=20
    fi

    for outer_k in 1 2; do
        for local_k in 1 2; do
            # Storage efficiency = outer_m / (outer_m + outer_k) * local_m / (local_m + local_k)
            outer_eff=$(echo "scale=6; $outer_m / ($outer_m + $outer_k)" | bc)
            local_eff=$(echo "scale=6; $local_m / ($local_m + $local_k)" | bc)
            storage_eff=$(echo "scale=3; $outer_eff * $local_eff" | bc)

            for perf in $PERF_RATIOS; do
                run_experiment "RQ3" "16enclosure_distributed_2io_multi_ec.json" "multi_ec" $outer_m $outer_k $local_m $local_k $local_n "$storage_eff" $perf 0.10 10000 "--ec_type multi_ec --m $outer_m --k $outer_k --local_n $local_n --local_m $local_m --local_k $local_k --n 16"
            done
        done
    done
done

echo ""
echo "=========================================="
echo "Phase 1 Complete!"
echo "Results saved to: $RESULT_FILE"
echo "Total experiments: $current"
echo "=========================================="
