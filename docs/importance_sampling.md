# Importance Sampling for Rare Event Simulation

## 문제: Rare Event 시뮬레이션의 어려움

Storage 시스템의 durability 측정 시, data loss는 매우 드물게 발생합니다:
- AFR 3%인 디스크 48개, EC(22,2) 구성
- Data loss 확률: 연간 ~10^-9 수준
- 10만 번 시뮬레이션으로도 data loss 0건 → durability 측정 불가

## 해결책: Importance Sampling

**핵심 아이디어**: failure rate를 인위적으로 높여서 시뮬레이션하고, 나중에 수학적으로 보정

### 기본 원리

원래 구하고 싶은 기대값:
```
E[X] = ∫ x·f(x) dx    (f = 원래 분포)
```

f에서 샘플링하기 어려우니, 다른 분포 g를 사용:
```
E[X] = ∫ x·f(x) dx
     = ∫ x·(f(x)/g(x))·g(x) dx
     = E_g[X · w(x)]

여기서 w(x) = f(x)/g(x) = likelihood ratio = importance weight
```

Monte Carlo 추정:
```
E[X] ≈ Σ(x_i · w_i) / Σ(w_i)
```

## Disk Failure에 적용

### Exponential Distribution

디스크 failure는 exponential distribution을 따름:
```
f(t) = λ · e^(-λt)    where λ = 1/MTTF = AFR/8760
```

### 가속 (Acceleration)

k배 가속 시:
```
원래: f(t) = λ · e^(-λt)
가속: g(t) = kλ · e^(-kλt)
```

Likelihood ratio:
```
w(t) = f(t) / g(t) = (1/k) · e^(λt(k-1))
```

### 단순화: 왜 (1/k)^n 만 사용해도 되는가?

1. **Event-based 관점**

   시간 T 동안 n번 failure가 발생할 확률:
   ```
   원래: P_f(n) = (λT)^n · e^(-λT) / n!      (Poisson)
   가속: P_g(n) = (kλT)^n · e^(-kλT) / n!

   Ratio = P_f(n) / P_g(n) = (1/k)^n · e^((k-1)λT)
   ```

2. **상수 상쇄**

   `e^((k-1)λT)` 항은:
   - 시뮬레이션 기간 T에만 의존
   - 모든 샘플에 동일하게 적용

   비율 계산 시:
   ```
   결과 = Σ(x_i · w_i · C) / Σ(w_i · C) = Σ(x_i · w_i) / Σ(w_i)
   ```
   → 상수 C가 상쇄됨!

3. **결론**
   ```
   w = (1/k)^n

   where:
     k = acceleration_factor
     n = disk failure 횟수
   ```

## 구현

### 1. Failure Event 생성 시

```cpp
// 가속된 MTTF 사용
double effective_mttf = disk_mttf / acceleration_factor;
std::exponential_distribution<double> dist(1.0 / effective_mttf);
double failure_time = current_time + dist(rng);

// Weight 업데이트
importance_weight *= (1.0 / acceleration_factor);
```

### 2. 결과 집계

```cpp
// 각 시뮬레이션 결과
struct SimulationResult {
    double data_loss_ratio;
    double importance_weight;  // = (1/k)^n
};

// Weighted average
double weighted_data_loss = 0;
double sum_weights = 0;

for (auto& result : results) {
    weighted_data_loss += result.data_loss_ratio * result.importance_weight;
    sum_weights += result.importance_weight;
}

double avg_data_loss = weighted_data_loss / sum_weights;
double durability = 1.0 - avg_data_loss;
```

### 3. Effective Sample Size (ESS)

Importance sampling의 품질 지표:
```cpp
double sum_w = Σ(w_i);
double sum_w2 = Σ(w_i^2);

ESS = sum_w^2 / sum_w2
ESS_ratio = ESS / N   // N = 총 샘플 수
```

- ESS_ratio > 0.5: 좋음
- ESS_ratio > 0.1: 적당함
- ESS_ratio < 0.1: 가속이 너무 과함, variance 높음

## 가속 계수 선택 가이드

| Acceleration | AFR (원래 3%) | 권장 여부 | 비고 |
|--------------|---------------|-----------|------|
| 1x | 3% | - | 가속 없음 (baseline) |
| 5x | 15% | 권장 | 안전한 가속 |
| 10x | 30% | 권장 | 일반적으로 좋은 선택 |
| 20x | 60% | 주의 | ESS 확인 필요 |
| 30x+ | 90%+ | 위험 | Failure mode 변경 가능 |

### 주의사항

1. **Failure mode 변경**: 너무 높은 가속은 비현실적인 failure 패턴 생성
   - 예: rebuild 중 3-4개 추가 failure (실제로는 거의 불가능)

2. **Non-disk failure**: io_module, switch 등은 가속하지 않음
   - 이미 disk보다 failure rate가 높음
   - 가속하면 disk failure 기반 data loss를 왜곡

3. **Empirical CDF**: 실제 failure trace 사용 시 가속 적용 어려움
   - 이 경우 가속 비활성화 권장

## 검증 방법

1. **가속 없이 baseline 측정** (가능한 경우)
2. **여러 가속 계수로 테스트**
   ```
   k=1:  durability = 0.999999xxx (수렴 느림)
   k=5:  durability = 0.999999yyy (빠르게 수렴)
   k=10: durability = 0.999999zzz (더 빠름)
   ```
3. **결과가 일관되면 OK**
4. **ESS_ratio > 0.1 확인**

## 참고 문헌

- Heidelberger, P. (1995). "Fast Simulation of Rare Events in Queueing and Reliability Models"
- Juneja, S., & Shahabuddin, P. (2006). "Rare-Event Simulation Techniques"
- Rubino, G., & Tuffin, B. (2009). "Rare Event Simulation using Monte Carlo Methods"
