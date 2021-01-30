---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: R
  language: R
  name: ir
---

# 2021-01-21 Phosphate Sensor Calibration

- Calibration of phosphate sensor for the samples run from 01/17-01/19
- Also received a new tube of sensor from the same lot (2198260E), so included that in the testing
    - For this tube, I aliquoted into 6 tubes with 33 uL, then had 45 uL leftover that I used for the calibration
- Wanted to get a sense for variations between calibrations with the same lot and different dilutions, so I did not mix these two tubes, but instead staggered them so that I had titrations of phosphate concentration read out with both. Then I could build curves with each set separately and both sets together and compare.
- One mistake I made was not freeze-thawing the two sensor aliquotes the same number of times. I.e. both tubes were delivered frozen, on dry ice. But the "sensor 1" was thawed, aliquoted, and re-frozen with LN2 on 01/17, while the "sensor 2" was thawed for the first time on 01/21 and aliquoted, but then immediately diluted for use in the calibration. I probably should have done one additional freeze-thaw cycle for the "sensor 2" aliquot, since that matches what is done to the aliquots that are used in the GAP assays.  
- Also, it should be noted that the dilutions are slightly different, both from each other and from what is used in the GAP assays. For the GAP assays, I thaw a 33 uL aliquot and dilute 1:4 with 99 uL assay buffer, to 125 uM. For the calibration I diluted instead to 40 uM ( a 40:500 or 1:12.5 dilution), so that I would have a larger volume to pipette (50 uL of 40 uM sensor instead of 19.36 uL of 125 uM sensor).
    - For "sensor 1" (aliquoted 11/17, from lot 21928260E), this dilution was 33 uL + 379.5 uL buffer (19.5 uL with P20 + 360 uL with P1000). 33 uL diluted to 412.5 with 397.5 uL is a 12.5x dilution
    - For "sensor 2" (aliquoted 11/21, from lot 21928260E), this dilution was 45 uL + 517.5 uL buffer (17.5 uL with P20 + 500 uL with P1000). 45 uL diluted to 562.5 with 517.5 uL is a 12.5x dilution
- The buffer was the same assay buffer from GAP assays: 40 mM HEPES pH 7.5, 100 mM NaCl, 4 mM MgCl2, 1 mM DTT
- The phosphate samples were made using 100 mM dibasic potassium phosphate (K2HPO4), prepared previously by Tina.
- The dilutions were from 15 uM down to 0 uM, and the expected linear regime according to Tina is 60% of the sensor concentration (so up to 12 uM of phosphate for 20 uM of sensor)
- The gain was measured at 57, as has been down for all GAP assays in this project, and the signal measured using excitation at 425nm and emission at 457nm
- plate set-up at `20210121_sensor_calib.xslx`.
    - One thing to note is that I did not have enough of sensor 1 for sample 14, which was supposed to have 14 uM phosphate. I left this one last intentionally in case I would run out. I did have enough of sensor 1 for this sample, but before I realized this I had removed the buffer and phosphate from well B4. So I remade the sample in B6, and that change is reflected in the excel table


```{code-cell} r
:tags: [remove-stderr]

library(tidyverse)
library(readxl)
read_excel('./20210121_sensor_calib.xlsx')
```
