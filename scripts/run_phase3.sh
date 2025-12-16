#!/bin/bash
# Phase 3: RQ5 (SSD vs HDD)
# 예상 소요 시간: ~1시간 (60 실험 × 1분)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/build"
CONFIG_DIR="$PROJECT_DIR/configs"
RESULT_DIR="$PROJECT_DIR/results/phase3"

mkdir -p "$RESULT_DIR"

# 공통 옵션
BASE_OPTS="--sim_years 5 --num_sim 10000 --n 32 --total_disks 384 --disk_failure_trace $PROJECT_DIR/ssd_failure_times.txt --seed 42"
PERF_RATIOS="0.50 0.60 0.70 0.80 0.90 0.99"
K_VALUES="2 3 4 5 6"

# 결과 파일 초기화
RESULT_FILE="$RESULT_DIR/phase3_results.csv"
echo "rq,disk_type,m,k,storage_eff,target_perf,availability,perf_availability,durability,avg_rebuild_time" > "$RESULT_FILE"

# 실험 카운터
total_experiments=60
current=0

run_experiment() {
    local disk_type=$1
    local config=$2
    local m=$3
    local k=$4
    local storage_eff=$5
    local perf=$6

    current=$((current + 1))
    echo "[$current/$total_experiments] disk=$disk_type, EC($m,$k), perf=$perf"

    local output=$("$BUILD_DIR/availability_sim" "$CONFIG_DIR/$config" \
        $BASE_OPTS \
        --rebuild_bw_ratio 0.10 \
        --degraded_ratio 0.10 \
        --target_perf $perf \
        --m $m --k $k 2>&1)

    # 결과 파싱
    local avail=$(echo "$output" | grep '^  Availability:' | head -1 | awk '{print $2}')
    local perf_avail=$(echo "$output" | grep 'Perf Availability MIN' | head -1 | awk '{print $NF}')
    local durability=$(echo "$output" | grep '^  Durability:' | head -1 | awk '{print $2}')
    local rebuild_time=$(echo "$output" | grep 'Avg rebuild time:' | awk '{print $4}')

    echo "RQ5,$disk_type,$m,$k,$storage_eff,$perf,$avail,$perf_avail,$durability,$rebuild_time" >> "$RESULT_FILE"
}

echo "=========================================="
echo "Phase 3: Starting experiments"
echo "=========================================="

# ===========================================
# RQ5: SSD vs HDD (60 experiments)
# ===========================================
echo ""
echo "=== RQ5: SSD vs HDD ==="
echo ""

# SSD with varying k
echo "--- SSD ---"
for k in $K_VALUES; do
    m=$((24 - k))
    storage_eff=$(echo "scale=3; $m / 24" | bc)
    for perf in $PERF_RATIOS; do
        run_experiment "ssd" "16enclosure_distributed_2io.json" $m $k $storage_eff $perf
    done
done

# HDD with varying k
echo "--- HDD ---"
for k in $K_VALUES; do
    m=$((24 - k))
    storage_eff=$(echo "scale=3; $m / 24" | bc)
    for perf in $PERF_RATIOS; do
        run_experiment "hdd" "16enc_hdd_2io.json" $m $k $storage_eff $perf
    done
done

echo ""
echo "=========================================="
echo "Phase 3 Complete!"
echo "Results saved to: $RESULT_FILE"
echo "=========================================="
