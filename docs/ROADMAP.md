# Research roadmap

The implemented benchmark is a starting point. These extensions are not completed:

1. Independent generator: integrate CVSim/RCVSIM with a documented export adapter, license review, and fixed input/output unit checks. Repeat tests across simulator families.
2. Observation study: compare BP+HR against BP+HR+flow, varying observation windows and sampling rates. Use profile likelihood or posterior sampling to quantify parameter ambiguity.
3. Dynamic correction: compare residual trees with neural state-space or neural differential-equation corrections, while testing whether they preserve stability and physiological constraints.
4. Better study design: multiple training seeds, larger patient cohorts, paired stress tests, realistic correlated missingness, and sensitivity to parameter distributions.
5. Clinical bridge: seek appropriately licensed physiological waveform data with measured flow or credible reference measurements. Validate factual forecasts before considering treatment-response claims.
6. Treatment modeling: define an explicit intervention mechanism and evaluate identifiability and causal assumptions. Inflow multipliers in this prototype are not medication effects.

A useful research contribution would be a reproducible account of which measurements resolve ambiguity and which apparent improvements disappear under generator mismatch. Do not infer novelty from the project title.
