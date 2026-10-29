# FitzHugh_Nagumo validation data

The obs_data files hold the synthetic benchmark of Ramsay et al. (2007), as shipped with circulatory_autogen (`resources/FitzHugh_Nagumo_obs_data.json`; the operands are renamed to this repo's test vessel `mod`).

- **Benchmark:** J. O. Ramsay, G. Hooker, D. Campbell and J. Cao (2007), *Parameter estimation for differential equations: a generalized smoothing approach*, J. R. Statist. Soc. B 69(5), 741–796.
- **True values** (a, b, c) = (0.2, 0.2, 3.0), with V(0) = −1 and R(0) = 1.
- **Noise:** std 0.1 on V and 0.05 on R, sampled every 0.2 over 0–60.
- **Calibration** fits 0–40 from the poor initial guess (0.8, 0.9, 2.0), which sits in a local minimum of the least-squares surface. It then predicts 40–60.
- **Limitation:** this is a calibration and identifiability benchmark on synthetic data. It is not a comparison with physical measurements.
- **Noise:** the shipped series look noise-free (the declared std is 0.1 / 0.05, but the values are exact model output), so recovery is almost exact.

The items below have no species or preparation: the data are synthetic model output.

## mod/V

- Kind: derived (synthetic: FitzHugh-Nagumo output at (a, b, c) = (0.2, 0.2, 3.0))
- Role: target
- Source: ramsay2007parameter, the benchmark as shipped with circulatory_autogen (`resources/FitzHugh_Nagumo_obs_data.json`)
- Obtained: copied from the fixture; 201 samples every 0.2 over 0–40
- Sigma: 0.1, the benchmark's declared noise std (the values themselves look noise-free, above).

## mod/R

- Kind: derived (synthetic: FitzHugh-Nagumo output at (a, b, c) = (0.2, 0.2, 3.0))
- Role: target
- Source: ramsay2007parameter, the benchmark as shipped with circulatory_autogen (`resources/FitzHugh_Nagumo_obs_data.json`)
- Obtained: copied from the fixture; 201 samples every 0.2 over 0–40
- Sigma: 0.05, the benchmark's declared noise std.

## V_validation

- Kind: derived (synthetic, as mod/V)
- Role: held out
- Source: ramsay2007parameter, as mod/V
- Obtained: 301 samples every 0.2 over 0–60; the calibration fits 0–40, so 40–60 is the prediction
- Sigma: 0.1, as mod/V.

## R_validation

- Kind: derived (synthetic, as mod/R)
- Role: held out
- Source: ramsay2007parameter, as mod/R
- Obtained: 301 samples every 0.2 over 0–60; the calibration fits 0–40, so 40–60 is the prediction
- Sigma: 0.05, as mod/R.

## V_max_validation

- Kind: derived (the maximum of V_validation)
- Role: held out
- Source: V_validation (cam_testing.calibrate.feature_items)
- Sigma: V_validation's std at that sample.

## V_min_validation

- Kind: derived (the minimum of V_validation)
- Role: held out
- Source: V_validation (cam_testing.calibrate.feature_items)
- Sigma: V_validation's std at that sample.

## V_mean_validation

- Kind: derived (the mean of V_validation)
- Role: held out
- Source: V_validation (cam_testing.calibrate.feature_items)
- Sigma: propagated, sqrt(Σ std²) / n = 0.1 / √301.

## R_max_validation

- Kind: derived (the maximum of R_validation)
- Role: held out
- Source: R_validation (cam_testing.calibrate.feature_items)
- Sigma: R_validation's std at that sample.

## R_min_validation

- Kind: derived (the minimum of R_validation)
- Role: held out
- Source: R_validation (cam_testing.calibrate.feature_items)
- Sigma: R_validation's std at that sample.

## R_mean_validation

- Kind: derived (the mean of R_validation)
- Role: held out
- Source: R_validation (cam_testing.calibrate.feature_items)
- Sigma: propagated, sqrt(Σ std²) / n.
