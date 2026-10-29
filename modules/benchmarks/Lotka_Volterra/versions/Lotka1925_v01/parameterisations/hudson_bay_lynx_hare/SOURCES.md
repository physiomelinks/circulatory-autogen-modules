# Lotka_Volterra validation data

`hudson_bay_lynx_hare.csv` holds the lynx and snowshoe hare pelts (thousands) traded by the Hudson's Bay Company from 1900 to 1920.

- **File taken from:** the Stan example models: https://github.com/stan-dev/example-models/blob/master/knitr/lotka-volterra/hudson-bay-lynx-hare.csv. Its header says it came from P. Howard, *Modeling Basics* (Texas A&M), http://www.math.tamu.edu/~phoward/m442/modbasics.pdf.
- **Original records:** C. G. Hewitt (1921), *The Conservation of the Wild Life of Canada*, Scribner.
- **Published fit (used by `validation_test_baseline`):** B. Carpenter (2018), *Predator-Prey Population Dynamics: the Lotka-Volterra model in Stan*, https://mc-stan.org/learn-stan/case-studies/lotka-volterra-predator-prey.html.
  - Posterior means, with the 10% and 90% quantiles:

    | Parameter | Mean | 10% | 90% |
    |---|---|---|---|
    | α | 0.545 | 0.465 | 0.630 |
    | β | 0.028 | 0.022 | 0.033 |
    | γ | 0.803 | 0.692 | 0.926 |
    | δ | 0.024 | 0.020 | 0.029 |
    | initial hare | 33.956 | 30.415 | 37.630 |
    | initial lynx | 5.933 | 5.273 | 6.614 |

  - The measurement error is lognormal with σ ≈ 0.25.
  - In the model, hare is the prey `x` and lynx is the predator `y`. Time is in years from 1900.

**Sigma of every item:** 0.25 × the value (checked in the obs_data). This matches the σ ≈ 0.25 of
Carpenter's lognormal error, but nothing here records that it was chosen for that reason
(unconfirmed). The scalar features' sigmas follow cam_testing.calibrate.feature_items.

## mod/x

- Kind: measured (pelt counts, thousands, a proxy for the hare population)
- Role: target
- Source: hewitt1921conservation, via carpenter2018predator (the Stan example's `hudson-bay-lynx-hare.csv`)
- Species/preparation: snowshoe hare; Hudson's Bay Company pelt records, 1900–1915
- Obtained: table values, the Hare column of `hudson_bay_lynx_hare.csv` (16 years)
- Sigma: 0.25 × value (above).

## mod/y

- Kind: measured (pelt counts, thousands, a proxy for the lynx population)
- Role: target
- Source: hewitt1921conservation, via carpenter2018predator (the Stan example's `hudson-bay-lynx-hare.csv`)
- Species/preparation: Canada lynx; Hudson's Bay Company pelt records, 1900–1915
- Obtained: table values, the Lynx column of `hudson_bay_lynx_hare.csv` (16 years)
- Sigma: 0.25 × value (above).

## x_validation

- Kind: measured (as mod/x)
- Role: held out
- Source: hewitt1921conservation, via carpenter2018predator
- Species/preparation: snowshoe hare; pelt records, 1900–1920
- Obtained: the Hare column, all 21 years; calibration fits 1900–1915, so 1915–1920 is the prediction
- Sigma: 0.25 × value (above).

## y_validation

- Kind: measured (as mod/y)
- Role: held out
- Source: hewitt1921conservation, via carpenter2018predator
- Species/preparation: Canada lynx; pelt records, 1900–1920
- Obtained: the Lynx column, all 21 years; calibration fits 1900–1915, so 1915–1920 is the prediction
- Sigma: 0.25 × value (above).

## x_max_validation

- Kind: derived (the maximum of x_validation)
- Role: held out
- Source: x_validation (cam_testing.calibrate.feature_items)
- Sigma: x_validation's std at that year.

## x_mean_validation

- Kind: derived (the mean of x_validation)
- Role: held out
- Source: x_validation (cam_testing.calibrate.feature_items)
- Sigma: propagated, sqrt(Σ std²) / n.

## y_max_validation

- Kind: derived (the maximum of y_validation)
- Role: held out
- Source: y_validation (cam_testing.calibrate.feature_items)
- Sigma: y_validation's std at that year.

## y_min_validation

- Kind: derived (the minimum of y_validation)
- Role: held out
- Source: y_validation (cam_testing.calibrate.feature_items)
- Sigma: y_validation's std at that year.

## y_mean_validation

- Kind: derived (the mean of y_validation)
- Role: held out
- Source: y_validation (cam_testing.calibrate.feature_items)
- Sigma: propagated, sqrt(Σ std²) / n.
