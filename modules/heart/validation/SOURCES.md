# heart validation candidates

No validation data files have been added yet: every component's validation status is `pending`. The candidates below were identified for the owner to confirm.

## Normal adult values (baseline / calibrate for heart_simple, _wCont, _Ca_input, _LVprop, _nonstiff, _new_valve)

- **Maceira et al. (2006), J Cardiovasc Magn Reson 8(3):417-426, doi:10.1080/10976640600572889** (`Maceira2006LV`).
  - Table 8 (SSFP CMR, n = 120: 60 men, 60 women; mean +- SD):
    - LV EDV: 142 +- 21 ml (men 156 +- 21, women 128 +- 21)
    - LV ESV: 47 +- 10 ml
    - SV: 95 +- 14 ml
    - EF: 67 +- 4.6 %
  - Table 1: brachial cuff pressure 125 +- 7 / 73 +- 5 mmHg for men aged 20-29.
- **Maceira et al. (2006), Eur Heart J 27(23):2879-2888, doi:10.1093/eurheartj/ehl336** (`Maceira2006RV`). The values were read from the reproduction in Kawel-Boehm et al. (2020), Table 9, because the primary article was paywalled:
  - RV EDV: 163 +- 27 ml (men), 127 +- 24 ml (women)
  - RV ESV: 57 +- 17 ml (men), 44 +- 15 ml (women)
  - RVEF: 66 +- 7 %
- **Kawel-Boehm et al. (2020), J Cardiovasc Magn Reson 22:87, doi:10.1186/s12968-020-00683-3** (`KawelBoehm2020`).
  - Table 2: LV cardiac output 5.6 +- 1.1 l/min (men), 4.5 +- 0.9 l/min (women).

**How they would be used:** the comparison is of scalar features of the periodic state (max/min of q_lv and q_rv, max/min of aortic_root/u, mean aortic flow), taken after the 20 s pre-time.

**Why they are not set up:** `validation_test_baseline` compares time series from t = 0, and the calibrate evaluation runs with pre_time 0. Both would include the start-up transient, so this needs a framework change first.

**Preliminary comparison** (nominal parameters, test harness): see the component reasons in `heart_tests.yaml`. For example, heart_simple gives LVEDV 203 ml, ESV 81 ml, EF 60 % and CO 7.2 l/min.

## heart_devel / heart_new_valve

The candidate is the owner's BioBeat subject data used by circulatory_autogen `new_valve_p_est`:
- r_ao, r_lvot and A_aov
- T = 0.725 s and T_vc = 0.331 s
- gamma, from the population identification of 07/09/22

These data are not public. The owner should confirm whether they can be used.

## heart_ASD

No source has been chosen yet. A moderate secundum ASD data set (Qp/Qs, RV/LV volumes) is needed.
