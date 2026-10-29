# Validation data provenance: microvasculature_network

No data files; the one proposed baseline uses scalar targets written inline in the spec.

- **capillary / pp_micro, scalar target** (status: proposed). Source: Williams SA et al. (1988), *Dynamic measurement of human capillary blood pressure*, Clin Sci 74:507–512, doi:10.1042/cs0740507. The values come from the abstract.
  - Measured mean nailfold capillary pressures: 37.7 ± 3.7 mmHg at the arteriolar limb, 19.4 ± 1.0 mmHg at the apex and 14.6 ± 0.5 mmHg at the venular limb. The ± is as reported and is probably the SEM.
  - The test holds the pp capillary between the limb pressures: u_in = 37.7 mmHg = 5026 Pa and u_out = 14.6 mmHg = 1946.5 Pa. It compares the steady mid-segment pressure u with the apex value (z ≤ 2).
  - The model gives 26.1 mmHg (the mean of the two limb pressures), so the target fails. See the review questions.
