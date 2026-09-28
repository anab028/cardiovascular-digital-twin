# Benchmark results

Synthetic, single-seed proof of concept. Lower MAE is better. Confidence intervals bootstrap patients, not time points.

| Scenario | Model | MAE (mmHg), 95% CI | Trajectory coverage | Band width (mmHg) |
|---|---|---|---|---|
| standard | persistence | 11.71 [9.65, 14.03] | 87% | 62.77 |
| standard | physiology | 0.82 [0.71, 0.95] | 97% | 5.48 |
| standard | data_driven | 6.72 [5.26, 8.54] | 90% | 33.99 |
| standard | hybrid | 0.47 [0.40, 0.54] | 100% | 3.54 |
| sparse | persistence | 9.48 [7.60, 11.78] | 93% | 62.77 |
| sparse | physiology | 0.86 [0.74, 0.98] | 90% | 5.48 |
| sparse | data_driven | 6.33 [4.49, 8.48] | 90% | 33.99 |
| sparse | hybrid | 0.58 [0.47, 0.69] | 87% | 3.54 |
| noisy | persistence | 11.06 [8.39, 13.91] | 90% | 62.77 |
| noisy | physiology | 0.94 [0.79, 1.11] | 77% | 5.48 |
| noisy | data_driven | 7.76 [5.90, 9.92] | 80% | 33.99 |
| noisy | hybrid | 0.70 [0.56, 0.86] | 73% | 3.54 |
| intervention | persistence | 46.48 [42.11, 51.41] | 0% | 62.77 |
| intervention | physiology | 1.29 [1.12, 1.46] | 57% | 5.48 |
| intervention | data_driven | 39.11 [35.42, 42.62] | 0% | 33.99 |
| intervention | hybrid | 0.72 [0.60, 0.86] | 63% | 3.54 |
| shift | persistence | 9.67 [7.94, 11.64] | 97% | 62.77 |
| shift | physiology | 2.39 [2.14, 2.69] | 20% | 5.48 |
| shift | data_driven | 5.04 [4.10, 5.96] | 97% | 33.99 |
| shift | hybrid | 1.67 [1.45, 1.93] | 17% | 3.54 |

Coverage targets the complete discrete forecast trajectory at nominal 90% under exchangeability. Shifted scenarios have no coverage guarantee. Bands are calibrated against noiseless synthetic truth, not future noisy measurements.

Matched two-element generator, known HR/SV: mean relative parameter errors {'R': 0.0022913717398945896, 'C': 0.028250529724102285}.

The three-element generator creates intentional model mismatch. Its fitted two-element parameters are effective parameters, not necessarily the generating physiological values.
