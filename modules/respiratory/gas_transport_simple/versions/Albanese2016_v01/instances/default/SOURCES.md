# gas_exchange validation data

All values are proposed and need the owner's confirmation.

## `pulmonary_GE_normal_blood_gases.csv` (pulmonary_GE baseline, status `proposed`)

This file holds constant target values from t = 150 to 180 s, after the harness has reached a quasi-steady state:
- P_O2 = 100 mmHg (13332 Pa)
- P_CO2 = 40 mmHg (5333 Pa)

**Source:** A. Albanese (2014), PhD thesis, Columbia University, doi:10.7916/D8JQ0Z52, Table 2-12, "Normal Value" column. That table gives:
- arterial P_O2 100 mmHg and P_CO2 40 mmHg
- alveolar P_O2 104 mmHg and P_CO2 40 mmHg
- mixed venous P_O2 40 mmHg and P_CO2 46 mmHg
- mixed venous O2 content 15 ml/dl and CO2 content 53 ml/dl

The thesis attributes the partial pressures to its reference [74]. In its bibliography, [74] is a coronary CT paper, so this is probably a numbering slip. The values are the standard textbook normals (West, *Respiratory Physiology: The Essentials*, 8th ed., 2008, is the thesis's reference [70]).

**How it is used:**
- The harness feeds pulmonary_GE with mixed-venous blood held at 15 ml/dl O2 and 53 ml/dl CO2, taken from the same table.
- The model's capillary P_O2_p and P_CO2_p are compared with the arterial targets using a relative RMS error (nrmse on constant data), with a threshold of 0.1.
- The capillary value comes before shunt mixing. With a 1.7% shunt, arterial P_O2 would be a few mmHg lower.
- The result depends on the simple_lung_bg parameters in the harness, which are library values (see the lung review).

## `severinghaus1979_odc.csv` (candidate for gas_transport_simple / pulmonary_GE, status `pending`)

This file holds points on the standard human O2 dissociation curve. They are computed from eq. (1) of J. W. Severinghaus (1979), J Appl Physiol 46(3):599-602, doi:10.1152/jappl.1979.46.3.599:

    S = 1 / (23400 / (P^3 + 150 P) + 1)        (P in mmHg)

The abstract states that this equation fits the standard curve (37 C, pH 7.40, P_CO2 40 mmHg) to within ±0.0055 saturation. The tabulated standard curve itself was not accessible, so the points come from the equation, not from the table.

**Offline check:** Spencer/C_sat at P_CO2 = 40 mmHg, using the library constants, is within 0.03 saturation of these points. The largest deviation, 0.029, is at P_O2 = 20-30 mmHg; above 50 mmHg the deviation is below 0.007.

**Not set up yet:** the baseline test compares time series, and this is an x-y curve (framework need).
