#!/bin/bash
# Phase 2: RQ2 (Hardware Topology) + RQ4 (Encoding Entity)
# 예상 소요 시간: ~48분 (48 실험 × 1분)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/build"
CONFIG_DIR="$PROJECT_DIR/configs"
RESULT_DIR="$PROJECT_DIR/results/phase2"

mkdir -p "$RESULT_DIR"

# 공통 옵션
BASE_OPTS="--sim_years 5 --num_sim 10000 --n 32 --total_disks 384 --disk_failure_trace $PROJECT_DIR/ssd_failure_times.txt --seed 42"
PERF_RATIOS="0.50 0.60 0.70 0.80 0.90 0.99"

# 결과 파일 초기화
RESULT_FILE="$RESULT_DIR/phase2_results.csv"
echo "rq,config,architecture,ec_type,m,k,local_m,local_k,storage_eff,target_perf,availability,perf_availability,durability,avg_rebuild_time,max_read_bw,max_write_bw,avg_read_bw,avg_write_bw" > "$RESULT_FILE"

# 실험 카운터
# RQ2: 4 arch × 2 EC × 6 perf = 48
# RQ4: 4 config × 6 perf = 24
total_experiments=72
current=0

run_experiment() {
    local rq=$1
    local config=$2
    local arch=$3
    local ec_type=$4
    local m=$5
    local k=$6
    local local_m=$7
    local local_k=$8
    local storage_eff=$9
    local perf=${10}
    local extra_opts="${11}"

    current=$((current + 1))
    echo "[$current/$total_experiments] RQ=$rq, arch=$arch, ec=$ec_type, m=$m, k=$k, perf=$perf"
    echo "CMD: $BUILD_DIR/availability_sim $CONFIG_DIR/$config $BASE_OPTS --rebuild_bw_ratio 0.10 --degraded_ratio 0.10 --target_perf $perf $extra_opts"

    local output=$("$BUILD_DIR/availability_sim" "$CONFIG_DIR/$config" \
        $BASE_OPTS \
        --rebuild_bw_ratio 0.10 \
        --degraded_ratio 0.10 \
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

    echo "$rq,$config,$arch,$ec_type,$m,$k,$local_m,$local_k,$storage_eff,$perf,$avail,$perf_avail,$durability,$rebuild_time,$max_read_bw,$max_write_bw,$avg_read_bw,$avg_write_bw" >> "$RESULT_FILE"
}

echo "=========================================="
echo "Phase 2: Starting experiments"
echo "=========================================="

# ===========================================
# RQ2: Hardware Topology (48 experiments)
# 4 architectures × 2 EC types × 6 perf
# EC types: 30+2 Standard, (14+1)*2+2 LRC
# ===========================================
echo ""
echo "=== RQ2: Hardware Topology ==="
echo ""

# EC Configurations
# Standard EC: m=30, k=2 (storage eff = 30/32 = 0.9375)
# LRC: m=28, k=2, local_m=14, local_k=1 (storage eff = 28/32 = 0.875)

declare -A ARCH_CONFIGS
ARCH_CONFIGS["2tier"]="16enclosure_distributed_2io.json"
ARCH_CONFIGS["3tier"]="16enc_3tier_2io.json"
ARCH_CONFIGS["vastdata"]="16enc_vastdata_2io.json"
ARCH_CONFIGS["purestorage"]="16enc_purestorage_2io.json"

for arch in "2tier" "3tier" "vastdata" "purestorage"; do
    config="${ARCH_CONFIGS[$arch]}"
    echo "--- Architecture: $arch ---"

    # Standard EC (30+2)
    for perf in $PERF_RATIOS; do
        run_experiment "RQ2" "$config" "$arch" "standard" 30 2 "" "" "0.9375" $perf "--m 30 --k 2"
    done

    # LRC ((14+1)*2+2): m=28, k=2, local_m=14, local_k=1
    for perf in $PERF_RATIOS; do
        run_experiment "RQ2" "$config" "$arch" "lrc" 28 2 14 1 "0.875" $perf "--ec_type lrc --m 28 --k 2 --local_m 14 --local_k 1"
    done
done

# ===========================================
# RQ4: Encoding Entity (24 experiments)
# Using Standard EC (30+2) for comparison
# ===========================================
echo ""
echo "=== RQ4: Encoding Entity ==="
echo ""

# Controller encoding (baseline) - 16 IO modules
for perf in $PERF_RATIOS; do
    run_experiment "RQ4" "16enc_encoding_controller.json" "2tier" "standard" 30 2 "" "" "0.9375" $perf "--m 30 --k 2"
done

# IO module encoding - 8 IO modules
for perf in $PERF_RATIOS; do
    run_experiment "RQ4" "8enc_encoding_iomodule.json" "2tier" "standard" 30 2 "" "" "0.9375" $perf "--m 30 --k 2"
done

# IO module encoding - 16 IO modules
for perf in $PERF_RATIOS; do
    run_experiment "RQ4" "16enc_encoding_iomodule.json" "2tier" "standard" 30 2 "" "" "0.9375" $perf "--m 30 --k 2"
done

# IO module encoding - 32 IO modules
for perf in $PERF_RATIOS; do
    run_experiment "RQ4" "32enc_encoding_iomodule.json" "2tier" "standard" 30 2 "" "" "0.9375" $perf "--m 30 --k 2"
done

echo ""
echo "=========================================="
echo "Phase 2 Complete!"
echo "Results saved to: $RESULT_FILE"
echo "Total experiments: $current"
echo "=========================================="
