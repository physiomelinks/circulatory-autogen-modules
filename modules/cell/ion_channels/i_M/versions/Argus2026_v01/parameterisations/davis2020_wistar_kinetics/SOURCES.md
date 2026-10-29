# davis2020_wistar_kinetics calibration data

Both items come from the representative Wistar trace of davis2020downregulation Fig. 3C
(source_figures/davis2020_fig3c.png), digitised and normalised by tools/i_M_davis2020_kinetics.py. Each
relaxation is mapped onto the gate w with the module's fixed w_inf, so the data constrain tau_w only.

## w_deactivation_at_-55mV

- Kind: measured (one cell's raw current, normalised onto w)
- Role: target
- Source: davis2020downregulation; the protocol from davis2021transcriptomic p. 172
- Figure: Fig. 3C, the representative Wistar trace (deactivation at −55 mV)
- Species/preparation: Wistar rat stellate neuron, n = 1; perforated patch, room temperature; raw current (not XE-991-subtracted); hold −25 mV, 1 s step to −55 mV
- Obtained: digitised column by column, normalised to the relaxation's own end points, mapped onto w (w_inf(−25) = 0.7311, w_inf(−55) = 0.1192)
- Primary read: single-exponential fit of the digitised current, tau 74 ± 4 ms, amplitude 180 pA
- Sigma: the residual SD about that fit (0.052 of the relaxation) in w units; the first sample, inside the voltage jump, has std 1 (not fitted).

## w_reactivation_at_-25mV

- Kind: measured (one cell's raw current, normalised onto w)
- Role: target
- Source: davis2020downregulation; the protocol from davis2021transcriptomic p. 172
- Figure: Fig. 3C, the representative Wistar trace (reactivation on return to −25 mV)
- Species/preparation: Wistar rat stellate neuron, n = 1; perforated patch, room temperature; raw current (not XE-991-subtracted); return from −55 mV to −25 mV
- Obtained: digitised column by column, normalised to the relaxation's own end points, mapped onto w (w_inf(−25) = 0.7311, w_inf(−55) = 0.1192)
- Primary read: single-exponential fit of the digitised current, tau 159 ± 20 ms, amplitude 203 pA
- Sigma: the residual SD about that fit (0.089 of the relaxation) in w units; the first 3 samples, inside the voltage jump, have std 1 (not fitted).
