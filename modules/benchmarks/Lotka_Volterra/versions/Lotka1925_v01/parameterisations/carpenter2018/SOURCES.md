# carpenter2018 baseline data

The baseline (a fit check, not an independent validation: Carpenter fitted these parameters to the same
data) compares the model at Carpenter's posterior means with `../hudson_bay_lynx_hare/hudson_bay_lynx_hare.csv`.
Where the file came from: `../hudson_bay_lynx_hare/SOURCES.md`.

## Hare

- Kind: measured (pelt counts, thousands, a proxy for the hare population)
- Role: validation
- Source: hewitt1921conservation, via carpenter2018predator (the Stan example's `hudson-bay-lynx-hare.csv`)
- Species/preparation: snowshoe hare; Hudson's Bay Company pelt records, 1900–1920
- Obtained: table values, the Hare column (model output `x`)

## Lynx

- Kind: measured (pelt counts, thousands, a proxy for the lynx population)
- Role: validation
- Source: hewitt1921conservation, via carpenter2018predator (the Stan example's `hudson-bay-lynx-hare.csv`)
- Species/preparation: Canada lynx; Hudson's Bay Company pelt records, 1900–1920
- Obtained: table values, the Lynx column (model output `y`)
