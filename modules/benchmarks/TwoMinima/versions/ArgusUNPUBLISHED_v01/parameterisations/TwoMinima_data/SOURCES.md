# TwoMinima calibration check

`TwoMinima__nn_obs_data.json` is circulatory_autogen's two-minima fixture (`resources/TwoMinima_obs_data.json`, with operands renamed to this repo's test vessel `mod`).

- **Data:** one constant item, the steady-state average of x = 4 ± 0.5.
- **Why two minima:** x settles to a², so a = +2 and a = −2 fit equally well.
- **What the test does:** calibrates with a local optimiser (L-BFGS-B, `sp_minimize`) from a = +1 and from a = −1 and requires |a| = 2 within 5% from each side, with the fitted steady state within 2 std of 4.

This checks calibration and identifiability, not physical validation. The validated range recorded for `a` is the envelope [−2, 2], but only the two endpoints were actually validated.
