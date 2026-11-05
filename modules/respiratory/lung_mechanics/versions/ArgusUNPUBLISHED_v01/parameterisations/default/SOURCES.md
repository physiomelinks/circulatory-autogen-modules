# lung validation data

No data files are set up yet; both lung validations are `pending`.

Candidate for simple_lung_bg: the quiet-breathing targets Albanese used to tune the model. They come from A. Albanese (2014), PhD thesis, Columbia University, doi:10.7916/D8JQ0Z52, section 2.4.3 and Table 2-5, which cite J. B. West, *Respiratory Physiology: The Essentials*, 8th ed., 2008:
- FRC 2.4 L
- tidal volume about 500 ml at 12 breaths/min
- end-expiratory pleural pressure -5 cmH2O

With the library values the model gives an end-expiratory lung volume of 14.0 L and VT of 0.88 L. With Albanese's C_A, C_b and q_A_us it gives 2.48 L and 0.67 L.
