# lv validation and calibration data

All of it is from liang2009multi (Liang et al. 2009, MBEC 47:743), Fig. 4 (p. 751), the normal case,
digitised (source_figures/ holds the panels). These curves are the output of Liang's own model (a
55-artery 1-D tree coupled to a lumped heart), not measurements: the items are model-to-model
comparisons, with no species or preparation. The calibration set-up is in lv_obs_data.json's comment.

## aortic_pressure

- Kind: derived (Liang et al.'s simulated ascending-aortic pressure, normal case)
- Role: target
- Source: liang2009multi
- Figure: Fig. 4b (p. 751), thick black curve; source_figures/liang2009_fig4b.png
- Obtained: digitised every 0.01 s into `liang2009_fig4b_aortic_pressure.csv` (axes from the tick marks)
- Sigma: 2 mmHg, the digitising uncertainty.

## lv_end_diastolic_volume

- Kind: derived (Liang et al.'s simulated P-V loop, normal case)
- Role: target
- Source: liang2009multi
- Figure: Fig. 4e (p. 751), right edge of the normal loop; source_figures/liang2009_fig4e.png
- Obtained: digitised, 117.6 ml
- Sigma: 2 ml, the digitising uncertainty.

## lv_end_systolic_volume

- Kind: derived (Liang et al.'s simulated P-V loop, normal case)
- Role: target
- Source: liang2009multi
- Figure: Fig. 4e (p. 751), left edge of the normal loop; source_figures/liang2009_fig4e.png
- Obtained: digitised, 41.8 ml
- Sigma: 2 ml, the digitising uncertainty.

## lv_pressure

- Kind: derived (Liang et al.'s simulated LV pressure, normal case)
- Role: held out
- Source: liang2009multi
- Figure: Fig. 4a (p. 751), thick black curve; source_figures/liang2009_fig4a.png
- Obtained: digitised every 0.01 s into `liang2009_fig4a_lv_pressure.csv`
- Sigma: 2 mmHg, the digitising uncertainty.

## LV peak systolic pressure

- Kind: derived (from the Fig. 4a curve)
- Role: validation
- Source: liang2009multi
- Figure: Fig. 4a: 108.14 mmHg
- Obtained: digitised
- Sigma: the target is the range 97.33–118.95 mmHg (±10 % of 108.14; the spec records the range, not why it is 10 %: unconfirmed).

## LV minimum (early diastolic) pressure

- Kind: derived (from the Fig. 4a curve)
- Role: validation
- Source: liang2009multi
- Figure: Fig. 4a: 1.56 mmHg
- Obtained: digitised
- Sigma: the target is the range 0–5 mmHg (the spec records the range, not its rationale: unconfirmed).

## LV end-diastolic pressure

- Kind: derived (from the Fig. 4a curve)
- Role: validation
- Source: liang2009multi
- Figure: Fig. 4a: 7.29 mmHg at the start of contraction
- Obtained: digitised
- Sigma: the range 5.29–9.29 mmHg, ±2 mmHg, the digitising uncertainty as in the calibration.
