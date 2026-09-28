# Methods and interpretation

## Dynamics and numerical method

Pressure is in mmHg, time in seconds, volume in mL. Resistance R and proximal impedance Z have units mmHg·s/mL; compliance C is mL/mmHg. Venous pressure is fixed to 5 mmHg. Initial capacitor pressure is 85 mmHg in generated patients.

For constant inflow on one time step, the exact capacitor update is:

`Pc[n+1] = a Pc[n] + (1-a)(Pv + R Q[n])`, where `a = exp(-dt/(RC))`.

This is exact for zero-order-held inflow, not for the continuous sinusoidal signal. The reference test checks a truly constant flow against its analytic solution. Arterial pressure is `Pc + Z Q`. At the intervention boundary, the proximal term responds immediately; the capacitor responds on the next step.

Synthetic parameters are sampled independently: R 0.8–1.5, C 0.9–2.2, SV 55–90 mL, HR 55–100 bpm. These are illustrative ranges, not a clinical population model. Standard Z is 0.005–0.025; shifted Z is 0.035–0.075. Standard future inflow multipliers are 0.8–1.2; the intervention stress test uses 1.5. Noise is independent Gaussian with SD 2 mmHg, or 5 mmHg in the noisy scenario.

## Personalization and information boundaries

Only the first 8 seconds of noisy, incomplete observations enter fitting. R, C, and initial pressure are fit with bounded least squares, treating HR and SV as known. The fit assumes Z = 0 even though the forecasting generator uses Z > 0. Thus fitted R/C are effective quantities under model mismatch.

A separate parameter experiment uses a matched Z = 0 generator and known flow scale to quantify recoverability. The pressure-only experiment analytically illustrates non-identifiability when SV is unknown; it does not fit a Bayesian posterior.

Missing history is interpolated only inside the observed history. Twenty-four evenly spaced history samples form ML features, along with HR, SV, planned inflow multiplier, forecast horizon, sine/cosine phase. The hybrid additionally receives fitted R/C and the physics prediction. Targets are future simulated pressure for the data model and future pressure minus physics prediction for the hybrid. Extra Trees has 80 trees and minimum leaf size six; parameters are fixed before test evaluation.

Each patient contributes one intervention trajectory. There is no repeated-patient leakage. The generated split manifest records every identifier. Test scenarios use independent virtual patients; comparisons between models within a scenario use the same patients, while comparisons between scenarios are not paired.

## Uncertainty

For each calibration patient, calculate maximum absolute error across the entire forecast grid. The conformal radius is the ceil((n+1)(1-alpha))-th order statistic, or infinity if that rank exceeds n. Predictions receive constant symmetric bands. This targets simultaneous coverage on the discrete trajectory, under exchangeability of calibration and test patients, and conditional on the fitted model. It is not continuous-time coverage and not a posterior distribution over parameters.

Calibration uses noise-free simulator truth as the outcome. This is possible in simulation; it is not a deployable clinical calibration procedure. No noisy future measurement coverage is claimed. Standard calibration is deliberately reused under shifts to reveal failures. No distribution-free guarantee is asserted for shifted populations.

MAE is first averaged over time within each patient, then across patients. The 95% intervals resample patients 1,000 times; they do not quantify uncertainty across model-training seeds. Small cohort coverage is an empirical proportion and will vary.

## Reproducibility

The CLI seed controls patient generation, noise, random missingness, the tree models, and the bootstrap. The output stores package versions and cohort sizes. Re-running with a different seed is a new experiment. For publication-level conclusions, run multiple seeds and larger independent cohorts, preserve every run, and predefine model-selection rules.
