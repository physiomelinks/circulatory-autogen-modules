# chamber/vv review: Liang 2009 data in progress (2026-10-12)

Not yet placed in the lv parameterisation. Liang et al. 2009 (MBEC 47:743), p. 751, rendered at 300 dpi:
- fig4a_raw.png: Fig. 4a; digitised.json "4a": the normal (thick black) LV pressure every 0.01 s over one
  beat (peak 108.1 mmHg at 0.21 s, minimum 1.6 at 0.45 s, ~7 mmHg at 0 and 1 s).
- fig4e_raw.png: Fig. 4e, the normal LV P-V loop: ESV 41.8 ml, EDV 117.6 ml (its vertical edges).
- digitise_liang.py: the script (axis calibration from tick marks).

Decisions (user, 2026-10-12): parameterisations lv/rv/la/ra (default = lv) from Liang Table 2 (p. 748) and
Korakianitis & Shi 2006 J Biomech Table 1 (V0: lv 5, rv 10, la 4, ra 4 ml), mmHg converted with 133.322368 Pa;
q_init lv = EDV from Fig. 4e, rv/la/ra = their V0; the omitted viscoelastic S dV/dt term of Liang Eq. 9 as a
bold departure note; validation of lv in a minimal LV loop against Fig. 4a/4e. Then point the 68 system-model
chamber records (17 each of lv/rv/la/ra) at their chamber's parameterisation.
