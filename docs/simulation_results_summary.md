# Availability Simulator - 작업 요약 및 결과

## 날짜: 2025-12-11

---

## 1. 구현된 주요 기능

### 1.1 커스텀 EC 그룹 매핑
- JSON에서 `ec_groups` 필드로 직접 디스크-그룹 매핑 지정 가능
- 파일: `include/erasure_coding.hpp`, `src/erasure_coding.cpp`

```json
"ec_groups": [
    [0, 48, 96, 144, 192, 240, 288, 336, ...],  // 그룹 0
    [3, 51, 99, 147, 195, 243, 291, 339, ...]   // 그룹 1
]
```

### 1.2 분산 배치 (Distributed Placement)
- 8개 enclosure에서 각 3개씩 디스크를 가져와 24개 디스크 그룹 구성
- enclosure 장애 시 각 EC 그룹에 최대 3개 디스크만 영향

### 1.3 프로파일링 옵션
- CMakeLists.txt에 `ENABLE_PROFILING` 옵션 추가
- 사용법: `cmake .. -DENABLE_PROFILING=ON`

### 1.4 성능 가용성 (Performance Availability)
- `--target_perf` 옵션으로 성능 기반 가용성 측정
- EC 그룹 중 target_perf 비율 이상 정상인 시간 비율 계산

---

## 2. 테스트 구성

### 2.1 하드웨어 구성
- 8개 NVMe Enclosure
- 각 Enclosure에 48개 SSD = 총 384 SSDs
- 4개 IO Module, 2개 Switch
- SSD ARR: 3.27% (연간 교체율)

### 2.2 테스트된 EC 구성

| EC 방식 | m | k | local_m | local_k | n | 효율 |
|---------|---|---|---------|---------|---|------|
| Standard EC | 22 | 2 | - | - | 24 | 91.7% |
| Standard EC | 20 | 4 | - | - | 24 | 83.3% |
| LRC | 20 | 2 | 10 | 1 | 24 | 83.3% |

---

## 3. 시뮬레이션 결과

### 3.1 Sequential vs Distributed 매핑 비교 (EC 22+2)

| 매핑 방식 | Availability (nines) | Durability (nines) | 실행 시간 |
|----------|---------------------|-------------------|----------|
| Sequential | 5.77 | 10.13 | 60초 |
| Distributed | 5.76 | 10.14 | 35초 |

**결론**: 가용성은 유사, 분산 매핑이 2배 빠름 (상태 공간 감소)

### 3.2 LRC vs Standard EC 비교 (n=24, 분산 배치)

| EC 방식 | Availability (nines) | Durability (nines) | 효율 |
|---------|---------------------|-------------------|------|
| LRC (m=20, k=2, local 10+1) | 6.02 | 13.0 | 83% |
| Standard EC (m=20, k=4) | 9.30 | 13.5 | 83% |

**결론**: 같은 효율에서 Standard EC가 availability 3+ nines 더 높음

### 3.3 성능 가용성 비교 (Performance Availability)

#### target_perf = 0.50 (50% 이상 성능 유지 시간)

| EC 방식 | Availability | Perf Avail (>=50%) | Avg Rebuild Time |
|---------|-------------|-------------------|------------------|
| LRC (20+2, local 10+1) | 6.02 nines | **100%** (13.0 nines) | 29.5시간 |

#### target_perf = 0.90 (90% 이상 성능 유지 시간)

| EC 방식 | Availability | Perf Avail (>=90%) |
|---------|-------------|-------------------|
| LRC (20+2, local 10+1) | 6.03 nines | **96.5%** (1.46 nines) |
| EC 20+4 | 9.44 nines | 89.1% (0.96 nines) |

#### target_perf = 0.95 (95% 이상 성능 유지 시간)

| EC 방식 | Availability | Perf Avail (>=95%) |
|---------|-------------|-------------------|
| LRC (20+2, local 10+1) | 6.03 nines | **74.1%** (0.59 nines) |
| EC 20+4 | 13.0 nines | 56.0% (0.36 nines) |

**핵심 인사이트**:
1. EC 20+4는 데이터 가용성(Availability)이 훨씬 높음 (9-13 nines)
2. LRC는 성능 가용성(Perf Avail)이 훨씬 높음 (target=0.9에서 96.5% vs 89.1%)
3. LRC의 로컬 rebuild가 성능 영향을 줄이는 데 효과적
4. target_perf=0.5에서 LRC는 항상 50% 이상 성능 유지 (data loss 0건)

---

## 4. 설정 파일 위치

```
/home/sejun000/availability_c-/configs/
├── 2tier.json                      # 기본 1 enclosure 설정
├── 8enclosure_distributed.json     # 8 enclosure, EC 22+2, 분산 배치
├── 8enclosure_lrc_distributed.json # 8 enclosure, LRC, 분산 배치
├── 8enclosure_ec20_4_distributed.json # 8 enclosure, EC 20+4, 분산 배치
└── ... (기타 설정 파일)
```

---

## 5. 명령어 예시

### 기본 실행
```bash
./availability_sim ../../configs/8enclosure_distributed.json \
    --total_disks 384 --m 22 --k 2 --nprocs 40 --num_sim 40000
```

### LRC 실행
```bash
./availability_sim ../../configs/8enclosure_lrc_distributed.json \
    --total_disks 384 --m 20 --k 2 --n 24 \
    --ec_type lrc --local_m 10 --local_k 1 \
    --nprocs 40 --num_sim 40000
```

### 성능 가용성 측정
```bash
./availability_sim ../../configs/8enclosure_ec20_4_distributed.json \
    --total_disks 384 --m 20 --k 4 \
    --target_perf 0.9 \
    --nprocs 40 --num_sim 40000
```

### 프로파일링 빌드
```bash
rm -rf CMakeCache.txt CMakeFiles
cmake .. -DENABLE_PROFILING=ON
make -j4
./availability_sim config.json --nprocs 1 --num_sim 5000
gprof ./availability_sim gmon.out > profile.txt
```

---

## 6. 코드 변경 요약

### 6.1 erasure_coding.hpp
- `ErasureCodingScheme` 클래스에 커스텀 매핑 지원 추가
- `initialize_with_custom_groups()` 메서드 추가
- `disk_to_group_`, `group_to_disks_` 멤버 추가

### 6.2 erasure_coding.cpp
- `get_disk_group()`, `get_stripe_disks()`, `get_total_groups()` 커스텀 매핑 지원
- `get_group_disks()` 새 메서드 추가

### 6.3 simulation_core.cpp
- JSON에서 `ec_groups` 파싱 및 `initialize_with_custom_groups()` 호출
- 디버그 메시지 추가 ("Using custom EC group mapping")

### 6.4 CMakeLists.txt
- `ENABLE_PROFILING` 옵션 추가 (-pg 플래그)

---

## 7. 향후 작업 제안

1. **total_disks를 JSON에 포함**: CLI 대신 JSON 설정 파일에서 읽기
2. **자동 분산 배치**: `ec_group_assignment: "distributed"` 옵션으로 자동 생성
3. **Enclosure 장애 시뮬레이션**: 전체 enclosure 장애 시나리오 테스트
4. **스레드 간 memoization 테이블 공유**: 성능 최적화 (현재 각 스레드 독립)
5. **Accelerated Testing (Importance Sampling)**: 고장률 가속 시뮬레이션

---

## 8. Accelerated Testing 기법 (검토 중)

### 8.1 개요
- Data loss가 거의 발생하지 않아 (13 nines durability) rare event 추정이 어려움
- AFR을 인위적으로 높여서 시뮬레이션 후, likelihood ratio로 보정하는 방법

### 8.2 Importance Sampling
```
1. AFR 10x로 시뮬레이션 (빠르게 수렴)
2. 각 failure event에 weight = (original_AFR / accelerated_AFR)^n 적용
3. 원래 AFR에서의 durability 추정
```

### 8.3 유효성 조건
| 조건 | 설명 |
|------|------|
| ✅ 유효 | 고장 메커니즘이 동일할 때 (같은 failure mode) |
| ✅ 유효 | 가속 인자(Acceleration Factor)가 검증되었을 때 |
| ❌ 무효 | 새로운 failure mode가 나타날 정도로 극단적 가속 |

### 8.4 안전한 가속 범위
- 현재 AFR 3.27% 기준
- **10x (AFR 33%)**: 아마 안전 - 단일 디스크 고장 후 rebuild 중 추가 고장이 주요 failure mode
- **30x+ (AFR 100%)**: 위험 - 동시다발 고장, rebuild 병목 등 새로운 failure mode 발생 가능

### 8.5 새로운 Failure Mode 예시 (과도한 가속 시)
1. **Rebuild 병목**: 동시에 너무 많은 디스크 rebuild → 네트워크/CPU 포화
2. **Correlated timing**: 고장 시점이 몰림 → 동시 고장 과대평가
3. **Spare 고갈**: 교체용 디스크 부족 (현재 무한 spare 가정)

---

## 9. 결론

- **데이터 안정성 우선**: Standard EC (20+4) 선택 - 9+ nines availability
- **성능 우선**: LRC (20+2, local 10+1) 선택 - 96.5% 성능 가용성
- **분산 배치**: enclosure 장애 시 영향 최소화, 시뮬레이션 속도 향상
