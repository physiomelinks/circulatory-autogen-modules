# FitzHugh_Nagumo validation data

The obs_data files hold the synthetic benchmark of Ramsay et al. (2007), as shipped with circulatory_autogen (`resources/FitzHugh_Nagumo_obs_data.json`; the operands are renamed to this repo's test vessel `mod`).

- **Benchmark:** J. O. Ramsay, G. Hooker, D. Campbell and J. Cao (2007), *Parameter estimation for differential equations: a generalized smoothing approach*, J. R. Statist. Soc. B 69(5), 741–796.
- **True values** (a, b, c) = (0.2, 0.2, 3.0), with V(0) = −1 and R(0) = 1.
- **Noise:** std 0.1 on V and 0.05 on R, sampled every 0.2 over 0–60.
- **Calibration** fits 0–40 from the poor initial guess (0.8, 0.9, 2.0), which sits in a local minimum of the least-squares surface. It then predicts 40–60.
- **Limitation:** this is a calibration and identifiability benchmark on synthetic data. It is not a comparison with physical measurements.
- **Noise:** the shipped series look noise-free (the declared std is 0.1 / 0.05, but the values are exact model output), so recovery is almost exact.
