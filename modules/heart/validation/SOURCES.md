# heart validation candidates

No data files are needed: the proposed baseline uses scalar `targets` in `heart_tests.yaml` (status `proposed`). They are set on heart vp_wCont, vp, vp_Ca, nonstiff and new_valve, and evaluated on the logged 5 s after a 20 s pre_time. Maceira 2006 Tables 1 and 8 were checked directly in the paper PDF (scmr.org). The Kawel-Boehm 2020 and Maceira RV numbers below come from a literature search and were not re-verified, so they are not used as targets.

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

**Targets used:** EDV, ESV, SV and EF from Table 8 (pass at 2 SD). Aortic systolic 115-135 and diastolic 68-79 mmHg: these ranges span the Table 1 age-decile means for men and women, but Table 1 gives brachial cuff pressures, so they are only a proxy for central aortic pressure.

**Result at nominal parameters** (heart_simple): SV 122 ml, EF 60 % and aortic 131/77 mmHg pass. EDV 203 ml and ESV 81 ml fail. heart_new_valve additionally misses SV (125 ml) and diastolic pressure (67 mmHg).

## heart_devel / heart_new_valve

The candidate is the owner's BioBeat subject data used by circulatory_autogen `new_valve_p_est`:
- r_ao, r_lvot and A_aov
- T = 0.725 s and T_vc = 0.331 s
- gamma, from the population identification of 07/09/22

These data are not public. The owner should confirm whether they can be used.

## heart_ASD

No source has been chosen yet. A moderate secundum ASD data set (Qp/Qs, RV/LV volumes) is needed.
