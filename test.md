# Storage System Availability & Performance Evaluation

## 실험 공통 설정
- **Disks**: 384 (16 enclosures × 24 SSDs)
- **EC Group Size (n)**: 32 disks per group
- **Simulations**: 10,000 (`--num_sim 10000`, 단 single IO module은 1,000)
- **Duration**: 5 years (`--sim_years 5`)
- **AFR**: 3.27% (`--disk_failure_trace ssd_failure_times.txt`)
- **target_perf sweep**: 0.50, 0.60, 0.70, 0.80, 0.90, 0.99 (6개, `--target_perf`)
- **rebuild_bw_ratio / degraded_ratio**: 0.10 고정 (모든 RQ)

### 시뮬레이터 Argument 참조
```bash
./build/availability_sim <config.json> \
    --total_disks 384 \
    --n 32 \
    --sim_years 5 \
    --num_sim 10000 \
    --disk_failure_trace ssd_failure_times.txt \
    --target_perf 0.50 \
    --rebuild_bw_ratio 0.10 \
    --degraded_ratio 0.10 \
    --seed 42 \
    --m 20 --k 4 \
    --ec_type standard|lrc|multi_ec \
    --local_m 6 --local_k 1 --local_n 4
```

**Note**: `--seed` 옵션으로 재현 가능한 시뮬레이션 수행. 같은 seed를 사용하면 availability/durability는 동일하고, target_perf에 따라 perf_availability만 변경됨.

## 공통 측정 방법
모든 RQ에서 다음을 측정:
1. **Availability**: 데이터 접근 가능 여부 (data loss 없음)
2. **Performance Availability**: target_perf 이상의 성능을 유지하는 시간 비율 (Perf Availability MIN)
3. **Durability**: 데이터 손실 없이 유지되는 비율
4. **Avg Rebuild Time**: 평균 rebuild 시간 (hours)

### 시뮬레이터 출력 예시
```
Availability: 1.0
Durability: 1.0
Perf Availability MIN (>=50%): 0.9999916666666667
Avg rebuild time: 66.5711 hours
```

### 공통 그래프 패턴
- **Graph Type A**: X축 = target_perf (0.5~0.99), Y축 = Performance Availability (nines)
- **Graph Type B**: X축 = 독립변수, Y축 = "6 nines perf availability를 달성하는 최대 target_perf"

---

## RQ1: IO Module Redundancy의 가치는 어느 정도인가?

### 가설
Dual IO module 구성은 single IO module 대비 performance availability를 크게 향상시키며,
특히 높은 target_perf_ratio에서 그 차이가 극대화될 것이다.

### 실험 설계
| Config | IO Modules | Switch | 비교 포인트 |
|--------|-----------|--------|------------|
| 16enc_1io | 1 | 2 | Baseline (SPOF 존재) |
| 16enc_2io | 2 | 2 | IO module redundancy 효과 |

**Sweep 변수:**
- k (parity): 1, 2, 3 (m = 24-k)
- target_perf: 0.50, 0.60, 0.70, 0.80, 0.90, 0.99
- rebuild_bw_ratio / degraded_ratio: 0.10 고정

**실험 횟수**: 3 k × 6 perf × 2 config = **36개**
- 1 IO module: `--num_sim 1000` (failure 빈도가 높아 적은 simulation으로도 충분)
- 2 IO modules: `--num_sim 10000`

### 그래프 설계

#### Graph 1.1: Performance Availability vs Target Perf Ratio
- X축: Target Performance Ratio (0.5 ~ 0.99)
- Y축: Performance Availability (nines)
- Lines: 1 IO module vs 2 IO modules (rebuild_ratio=0.1 고정)
- **목적**: IO redundancy가 높은 target_perf에서 얼마나 큰 차이를 만드는가

#### Graph 1.2: 6-Nines Achievable Performance vs Rebuild Ratio
- X축: rebuild_bw_ratio (0.02 ~ 0.20)
- Y축: 6 nines availability를 달성하는 최대 target_perf_ratio
- Lines: 1 IO module vs 2 IO modules
- **목적**: rebuild ratio가 증가할수록 IO redundancy의 가치가 어떻게 변하는가

### 예상 Finding
> **F1**: Dual IO module은 target_perf_ratio ≥ 0.8에서 availability를 2+ nines 향상시킨다.
> **F2**: Single IO module은 target_perf_ratio > 0.7에서 6 nines를 달성하기 어렵다.

---

## RQ2: Hardware Topology가 Availability에 미치는 영향은?

### 가설
계층이 깊어질수록 failure domain이 세분화되어 availability가 향상되지만,
component 수 증가로 인한 failure rate 상승과 trade-off가 존재한다.

### 실험 설계
| Architecture | Topology | Config | 특징 |
|-------------|----------|--------|------|
| 2-tier | switch → io_module → SSD | 16enclosure_distributed_2io.json | Simple, fewer components |
| 3-tier | switch → controller → io_module → SSD | 16enc_3tier_2io.json | Controller adds isolation |
| VastData-like | switch → controller → leaf_switch → SSD | 16enc_vastdata_2io.json | Distributed architecture |
| PureStorage-like | switch → blade → io_module → SSD | 16enc_purestorage_2io.json | Blade-based consolidation |

**EC Configuration:**
- **Standard EC (30+2)**: m=30, k=2, storage_eff=93.75%
- **LRC ((14+1)×2+2)**: m=28, k=2, local_m=14, local_k=1, storage_eff=87.5%

**Sweep 변수:**
- target_perf: 0.50, 0.60, 0.70, 0.80, 0.90, 0.99
- rebuild_bw_ratio / degraded_ratio: 0.10 고정

**실험 횟수**: 4 arch × 2 EC × 6 perf = **48개**

### 그래프 설계

#### Graph 2.1: Performance Availability vs Target Perf Ratio
- X축: Target Performance Ratio (0.5 ~ 0.99)
- Y축: Performance Availability (nines)
- Lines: 각 architecture별 선
- **목적**: 어떤 architecture가 가장 높은 availability를 달성하는가

#### Graph 2.2: 6-Nines Achievable Performance by Architecture
- X축: Architecture (categorical)
- Y축: 6 nines availability를 달성하는 최대 target_perf_ratio
- Bars: 각 architecture별 막대
- **목적**: architecture 선택이 achievable performance에 미치는 영향

### 예상 Finding
> **F3**: 3-tier는 controller 장애 시 더 작은 failure domain을 가져 2-tier 대비 availability가 향상된다.
> **F4**: VastData 구조는 높은 target_perf_ratio에서 가장 높은 availability를 제공한다.

---

## RQ3: Erasure Coding 방식이 Performance-Durability Trade-off에 미치는 영향은?

### 가설
LRC와 Multi-EC는 standard EC 대비 rebuild traffic을 줄여 performance availability를 향상시키지만,
storage overhead가 증가한다. Storage efficiency와 achievable performance 사이의 Pareto frontier가 존재한다.

### 실험 설계

#### Standard EC (6 configs × 6 perf = 36 experiments)
| k | m values | Storage Efficiency |
|---|----------|-------------------|
| 1 | 23, 16 | 95.8%, 94.1% |
| 2 | 22, 16 | 91.7%, 88.9% |
| 3 | 21, 16 | 87.5%, 84.2% |

#### LRC (4 configs × 6 perf = 24 experiments)
| m | k | local_m | local_k | Storage Efficiency |
|---|---|---------|---------|-------------------|
| 24 | 2 | 8 | 2 | 75.0% |
| 27 | 2 | 9 | 1 | 87.1% |
| 26 | 2 | 13 | 2 | 81.3% |
| 28 | 2 | 14 | 1 | 87.5% |

#### Multi-EC (8 configs × 6 perf = 48 experiments)
| outer_m | outer_k | local_m | local_k | local_n | Storage Efficiency |
|---------|---------|---------|---------|---------|-------------------|
| 12 | 1, 2 | 10 | 1, 2 | 24 | 다양 |
| 20 | 1, 2 | 20 | 1, 2 | 24 | 다양 |

**총 EC config: 6 + 4 + 8 = 18개**
**Sweep 변수:** target_perf × 6, rebuild_bw_ratio = 0.10 고정
**실험 횟수**: 36 + 24 + 48 = **108개**

### 그래프 설계

#### Graph 3.1: Storage Efficiency vs 6-Nines Achievable Performance (핵심 그래프)
- X축: Storage Efficiency (50% ~ 100%)
- Y축: 6 nines perf availability를 달성하는 최대 target_perf_ratio
- Points: Standard EC (파란), LRC (초록), Multi-EC (빨강)
- Pareto frontier line 표시
- **목적**: 같은 storage efficiency에서 어떤 EC 방식이 더 높은 performance를 달성하는가

#### Graph 3.2: Performance Availability Curves by EC Type
- X축: Target Performance Ratio (0.5 ~ 0.99)
- Y축: Performance Availability (nines)
- Lines: 대표적인 config들 (비슷한 storage efficiency에서 비교)
- **목적**: 비슷한 storage efficiency에서 curve 형태가 어떻게 다른가

#### Graph 3.3: Durability vs Storage Efficiency
- X축: Storage Efficiency (50% ~ 100%)
- Y축: Durability (nines)
- Points: 각 EC type별
- **목적**: Durability 관점에서의 trade-off 확인

#### Graph 3.4: Rebuild Traffic vs Local Group Size (Multi-EC)
- X축: Local group size (local_n)
- Y축: Average rebuild traffic
- Lines: 각 outer config별
- **목적**: local group size가 rebuild traffic에 미치는 영향

### 예상 Finding
> **F5**: 동일 storage efficiency에서 Multi-EC는 standard EC 대비 6-nines achievable perf가 5~15%p 높다.
> **F6**: Multi-EC의 local_n이 클수록 rebuild traffic이 증가하지만 durability도 향상된다.
> **F7**: LRC는 local_m이 클수록 local rebuild 비율이 높아져 performance가 향상된다.
> **F8**: Pareto frontier에서 Multi-EC가 70~90% efficiency 구간을 지배한다.

---

## RQ4: Encoding Entity (Controller vs IO Module)가 Write Performance에 미치는 영향은?

### 가설
IO module에서 encoding 시 cross-traffic이 발생하여 write bandwidth가 감소한다.
IO module 수가 많을수록 cross-traffic overhead가 줄어든다.

### 실험 설계
| Config | Encoding Entity | IO Modules | 예상 Write BW 감소 |
|--------|----------------|------------|------------------|
| 16enc_encoding_controller.json | Controller | 16 | 0% (baseline) |
| 8enc_encoding_iomodule.json | IO Module | 8 | ~25% |
| 16enc_encoding_iomodule.json | IO Module | 16 | ~12% |
| 32enc_encoding_iomodule.json | IO Module | 32 | ~6% |

**EC Configuration:** Standard EC (30+2): m=30, k=2

**Sweep 변수:**
- target_perf: 0.50, 0.60, 0.70, 0.80, 0.90, 0.99
- rebuild_bw_ratio / degraded_ratio: 0.10 고정

**실험 횟수**: 6 perf × 4 config = **24개**

### 그래프 설계

#### Graph 4.1: Write Performance Availability Comparison
- X축: Target Performance Ratio (0.5 ~ 0.99)
- Y축: Performance Availability (nines) - Write 기준
- Lines: Controller vs IO module (8, 16, 32 modules)
- **목적**: encoding entity가 write performance availability에 미치는 영향

#### Graph 4.2: Cross-traffic Overhead vs IO Module Count
- X축: IO Module 개수 (8, 16, 32)
- Y축: Write bandwidth 감소율 (%)
- Line: 실측값, Dashed: 이론값 (2/N)
- **목적**: cross-traffic overhead가 이론과 얼마나 일치하는가

### 예상 Finding
> **F9**: IO module encoding은 N개 모듈에서 약 2/N의 cross-traffic overhead를 발생시킨다.
> **F10**: 16개 이상의 IO module에서는 cross-traffic overhead가 12% 미만으로 감소한다.

---

## RQ5: SSD vs HDD의 Rebuild Time이 Availability에 미치는 영향은?

### 가설
HDD의 낮은 I/O bandwidth로 인해 rebuild time이 길어지고,
이 기간 동안 추가 장애 발생 확률이 높아져 durability가 크게 감소한다.

### 실험 설계
| Disk Type | Config | Read BW | Write BW | Capacity |
|-----------|--------|---------|----------|----------|
| SSD | 16enclosure_distributed_2io.json | 7 GB/s | 4 GB/s | 64 TB |
| HDD | 16enc_hdd_2io.json | 200 MB/s | 200 MB/s | 16 TB |

**Sweep 변수:**
- target_perf: 0.50, 0.60, 0.70, 0.80, 0.90, 0.99
- k (parity): 2, 3, 4, 5, 6 (m = 24-k)
- rebuild_bw_ratio / degraded_ratio: 0.10 고정

**실험 횟수**: 6 perf × 2 disk × 5 k = **60개**

### 그래프 설계

#### Graph 5.1: Performance Availability - SSD vs HDD
- X축: Target Performance Ratio (0.5 ~ 0.99)
- Y축: Performance Availability (nines)
- Lines: SSD (k=4), HDD (k=4), HDD (k=6)
- **목적**: 같은 parity에서 SSD와 HDD의 availability 차이

#### Graph 5.2: Required Parity for 6-Nines Durability
- X축: Disk Type (SSD, HDD)
- Y축: 6 nines durability를 달성하는 최소 k값
- Bars: 각 disk type별 필요 k값
- **목적**: HDD에서 SSD 수준의 durability를 얻으려면 얼마나 더 많은 parity가 필요한가

### 예상 Finding
> **F11**: HDD 기반 시스템은 SSD 대비 rebuild time이 약 35배 길어 durability가 3+ nines 감소한다.
> **F12**: HDD 환경에서 SSD k=4와 유사한 durability를 달성하려면 k=6 이상이 필요하다.

---

## 실험 요약 및 소요 시간

| RQ | 실험 횟수 | 1분/실험 기준 |
|----|----------|--------------|
| RQ1 | 36 | 36분 |
| RQ2 | 48 | 48분 |
| RQ3 | 108 | 1.8시간 |
| RQ4 | 24 | 24분 |
| RQ5 | 60 | 1시간 |
| **총합** | **276** | **~4.6시간** |

### Phase별 소요 시간
| Phase | RQ | 실험 횟수 | 소요 시간 |
|-------|-----|----------|----------|
| **Phase 1** | RQ1 + RQ3 | 144 | **~2.4시간** |
| Phase 2 | RQ2 + RQ4 | 72 | ~1.2시간 |
| Phase 3 | RQ5 | 60 | ~1시간 |

---

## 실험 우선순위

### Phase 1: Core Findings (필수)
1. **RQ1** - IO redundancy 효과 (기본 질문)
2. **RQ3** - EC 방식 비교 (핵심 contribution, Graph 3.1이 main figure)

### Phase 2: Architecture Comparison
3. **RQ2** - Hardware topology 비교
4. **RQ4** - Encoding entity 영향

### Phase 3: Extended Analysis
5. **RQ5** - SSD vs HDD

---

## 실험 실행 스크립트

### Phase 1: RQ1 + RQ3
```bash
./scripts/run_phase1.sh
```

### Phase 2: RQ2 + RQ4
```bash
./scripts/run_phase2.sh
```

### Phase 3: RQ5
```bash
./scripts/run_phase3.sh
```

---

## 결과 저장 위치

```
results/
├── phase1/phase1_results.csv   # RQ1 + RQ3
├── phase2/phase2_results.csv   # RQ2 + RQ4
└── phase3/phase3_results.csv   # RQ5
```
