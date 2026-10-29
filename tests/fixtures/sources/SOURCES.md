# enterocyte, whole cell

The data the `whole_cell` step of the calibration workflow fits (from the enterocyte repo, converted to
the SOURCES.md shape of modules/README.md, "Data sources"). One heading per `data_item_name` in
`enterocyte_obs_data.json`; value and sigma are in that file.

## u_e_rest

- Kind: measured
- Role: target
- Source: Rose RC, Schultz SG (1971). Studies on the electrical potential profile across rabbit ileum. J Gen Physiol 57:639–663. doi:10.1085/jgp.57.6.639
- Species/preparation: rabbit ileum, mucosal solution without sugar
- Primary read: apical ψ_mc −36 mV (abstract); basolateral about −40.5 mV as tabulated by Afshar et al. 2019 (secondary)
- Sigma: 3 mV stands for one membrane across the two readings.

The serosal composition has not been checked against the paper's methods.

## c_Na_rest — 15 ± 4 mM

- Kind: assumed (not measured in an enterocyte)
- Role: target
- Source: ion-selective microelectrode reviews: cytosolic Na 10–20 mM across mammalian epithelial and cardiac cells
- Sigma: ± 4 mM covers that range.

The weakest datum in the step, and the one to replace if an enterocyte measurement can be found.

## u_e_depolarisation_with_sugar

- Kind: measured
- Role: target
- Source: rose1971studies p. 648
- Species/preparation: rabbit ileum; 20 mM D-glucose added to the mucosal side
- Primary read: +6.0 ± 0.5 mV, n = 39 (pooled with L-alanine); glucose alone, digitised from Fig. 8, +5.8 ± 0.7, n = 15
- Sigma: widened to 1.5 mV for the pooling, the chamber concentration, and the cells' −32 mV rest against this model's −38 mV.

## na_time_to_90pc_after_sugar

- Kind: measured
- Role: target
- Source: Grasset, Gunter-Smith & Schultz 1983; Hudson & Schultz 1984; Lapointe, Hudson & Schultz 1986
- Species/preparation: Necturus small intestine

Used for the time course only: Necturus Na:sugar coupling is about 1:1 against this model's 2:1.

## c_Na_rise_with_sugar   (proposed)

- Kind: measured
- Role: target
- Source: the same Necturus papers: cytosolic Na rises from about 12 to 20 mM with luminal sugar.

## v_NKA_ratio_with_sugar   (proposed)

- Kind: measured
- Role: target
- Source: the same Necturus papers: Na transport rises three- to fourfold with luminal sugar.

## c_Na_peak_sugar

- Kind: measured
- Role: held out
- Source: the same Necturus papers: cytosolic Na peaks near 20 mM with luminal sugar.
- Species/preparation: Necturus small intestine

A prediction to compare against, not a target.

## u_e_with_sugar — no data

- Kind: derived: term of u_e_depolarisation_with_sugar

The membrane potential with sugar, which that item subtracts u_e_rest from. Weight 0, sigma 1e15.

## c_Na_sugar

- Kind: derived: term of c_Na_rise_with_sugar

## v_NKA_rest

- Kind: derived: term of v_NKA_ratio_with_sugar

No one measured the absolute resting pump rate in Necturus — the papers give the fold change.

## v_NKA_sugar
