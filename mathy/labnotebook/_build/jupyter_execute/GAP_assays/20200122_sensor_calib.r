# 2020-01-22 Phosphate Sensor Calibration

library(tidyverse)

d <- read_delim('../../data/GAP_assay/20210121_Phosphate_Calibration.txt',
                delim='\t', col_types=cols(), skip=46,
                locale = readr::locale(encoding = "windows-1252"))

head(d)

library(tidyverse)
library(readxl)

read_excel('../GAP_assays//20210121_sensor_calib.xlsx')

install.packages('readxl')

