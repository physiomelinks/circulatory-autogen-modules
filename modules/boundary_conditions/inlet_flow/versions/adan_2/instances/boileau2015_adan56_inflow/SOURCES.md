# open_loop validation data

Both files are published benchmark inflow waveforms from

- E. Boileau, P. Nithiarasu, P. J. Blanco, L. O. Müller, F. E. Fossan, L. R. Hellevik, W. P. Donders,
  W. Huberts, M. Willemet and J. Alastruey (2015), *A benchmark study of numerical schemes for one-dimensional
  arterial blood flow modelling*, Int. J. Numer. Meth. Biomed. Engng. 31(10), e02732, doi:10.1002/cnm.2732.

They were taken (unchanged, re-written as CSV with columns `t` [s], `Q` [m3/s]) from the openBF repository,
which redistributes the benchmark set-ups built from the paper's supplementary material:

| file | openBF source | benchmark case | period |
|---|---|---|---|
| `boileau2015_adan56_inflow.csv` | `models/boileau2015/adan56/adan56_inlet.dat` | ADAN56 aortic-root inflow, 22 points | 1.0 s |
| `boileau2015_ibif_inflow.csv` | `models/boileau2015/ibif/ibif_inlet.dat` | iliac (aortic) bifurcation inflow, 100 points | 1.1 s |

URL: https://github.com/INSIGNEO/openBF (fetched 2026-09-26).

Checks done while preparing the review:

- The `adan_flow_interp_BC_type` knot table is identical to the ADAN56 inlet data (zeros in the module where
  the data have 1e-10 m3/s).
- The `adan_flow_fourier_BC_type` series is a 21-harmonic fit of the same data (max deviation 3.9e-5 m3/s,
  about 7 % of the 5.7e-4 m3/s peak; mean 112.4 vs 112.9 mL/s).
- The `aorticbif_flow_BC_type` series reproduces the iliac-bifurcation data to round-off when T = 1.1 s
  (mean 7.9853e-6 m3/s); with T = 1 s it is time-compressed.

Limitation: these are prescribed input waveforms, so the baseline tests check that the module reproduces
the published inflow, not that a model prediction matches a measurement.
