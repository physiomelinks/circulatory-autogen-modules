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
