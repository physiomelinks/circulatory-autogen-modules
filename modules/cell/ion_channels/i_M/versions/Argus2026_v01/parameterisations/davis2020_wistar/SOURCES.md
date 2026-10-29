# davis2020_wistar calibration data

The protocol and the model's measure are in davis2020_wistar_obs_data.json's comment; the source figure
and the thesis pages defining the measure are in source_figures/.

## XE991_sensitive_I_M_deactivation

- Kind: measured (the current density); the conversion to a current uses an assumed Cm = 30 pF (not reported by Davis; the SN model's soma capacitance)
- Role: target
- Source: davis2020downregulation; the measure's definition from davis2021transcriptomic pp. 172–173
- Figure: Fig. 3C, the Wistar column (source_figures/davis2020_fig3c.png); thesis Fig. 4.10C is the same data
- Species/preparation: Wistar rat stellate neurons; perforated patch, room temperature; 10 µM XE-991-sensitive current; hold −25 mV, 1 s step to −55 mV
- Obtained: table value from the figure legend (−4.032 pA/pF, median, n = 11) × 30 pF = −0.12096 nA
- Primary read: −4.032 pA/pF (median, n = 11, legend)
- Sigma: the Fig. 3C error bar, digitised: 0.34 pA/pF × 30 pF = 0.0102 nA. The 11 points give SEM 0.34 pA/pF, which is what the bar shows, although the legend says median ± IQR.

Not corrected for the liquid junction potential (24.3 mV), as in the source.
