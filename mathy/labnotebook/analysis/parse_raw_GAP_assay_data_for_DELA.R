#### import libraries
library(tidyverse)
library(viridis)
library(lubridate)
library(ggrepel)
library(glue)


# first, compile index of all the sample conditions from the various days
import_index_table <- function(date) {
       filename <- glue('../GAP_assays/{date}/index_table_temp.txt')
       read_delim(filename, delim='\t', col_types=cols()) %>%
              mutate('column' = as.numeric(substr(well, 2, nchar(well))),
                     'GAP_conc' = ifelse(column %% 2 == 0, 0, 1),
                     'sensor_conc' = 20,
                     'date' = as.character(date)) %>% 
              select(well, 'sample' = prot, 'conc' = fin.Ran.conc, GAP_conc, sensor_conc, date)
}

index <- 
    data.frame('data_file' = list.files(path = '../../', pattern = '*_GAP_assay.txt', recursive = TRUE)) %>% 
    mutate(date = str_extract(data_file, '202[0-9]+')) %>% 
    left_join(x = map_dfr(.$date, import_index_table), y=., by = 'date') %>%
    mutate(data_file = paste0('../../', data_file))

# prepare calibration information, from Tina and from `20210121_sensor_calib.ipynb`
calibration <- tibble('sensor_lot' = c('1853001A', '2198260E'),  # took the values for Lot 1853001A from Tina's parse_raw_ R script
                      'sensor_conc' = c(20, 20), 
                      'cal_slope' = c(0.0003713, 0.0004629023),
                      'cal_intercept' = c(0.1462813, 0.2113335605))

# prepare i/o
dir_for_DELA <- 'parsed_data_for_DELA'
dir.create(file.path(dir_for_DELA), showWarnings = FALSE)
parsed_dataset_file <- paste0(dir_for_DELA, '/', 'DMS_GAP_data_parsed.txt')
intrinsic_hydrolysis_outfile <- paste0(dir_for_DELA, '/', 'DMS_intrinsic_hydrolysis_data.txt')
calculated_product_conc_file <- paste0(dir_for_DELA, '/', 'calculated_product_conc.txt')

## READ and FORMAT the raw ASCII data from Synergy H1
####### a function to read in all the biotek files, gather the data and reformat the time columns
read_and_gather <- function(file) {
  print(file)
  raw_in <-  read_lines(file)
  first_data_row <- grep(raw_in, pattern = 'Time\tT')
  last_data_row <- grep(raw_in, pattern = 'Results')
  if (length(last_data_row) > 0) {
    n_rows_to_read <- last_data_row - first_data_row - 2
    data_in <- read_delim(file, delim='\t', col_names = T, skip = (first_data_row - 1),
                          col_types=cols(), n_max = n_rows_to_read,
                          locale = locale(encoding = 'windows-1250'))
  } else {
    data_in <-  read_delim(file, delim='\t', col_names = T, skip = (first_data_row - 1),
                           col_types=cols(), locale = locale(encoding = 'UTF-8'))
  }
  data_gathered <- data_in %>% 
    select(., -2) %>% 
    gather(., key = well, value = fluorescence, -Time) %>% 
    mutate('Time' = as.integer(as.numeric(hms(Time))))  #### hms gives a warning when parsing 00:00:00
  data <- index %>% 
    filter(., data_file == file) %>% 
    select(-data_file) %>% 
    inner_join(., data_gathered, by = "well") %>% 
    mutate('row' = str_sub(well, 1, 1), "column" = str_sub(well, 2)) %>% 
    mutate('condition' = str_c(date, sample, well, conc, GAP_conc, sensor_conc, sep = '-')) %>% 
    filter(! is.na(fluorescence))
  no_GAP_blank_start <- data %>% filter(GAP_conc == 0) %>% 
    select(sample, date, conc, sensor_conc, Time, fluorescence) %>% 
    filter(Time < 100) %>% 
    group_by(sample, date, conc, sensor_conc) %>% 
    summarize('mean_starting_no_GAP_fluorescence' = mean(fluorescence)) %>% 
    ungroup()
  data <- data %>% filter(conc != 0) %>% 
    inner_join(., no_GAP_blank_start, by = c('sample', 'date', 'conc', 'sensor_conc')) %>% 
    mutate('corrected_fluorescence' = fluorescence - mean_starting_no_GAP_fluorescence) %>%
    mutate('sensor_lot' = case_when(date == 20201117 ~ '1853001A', TRUE ~ '2198260E')) %>%  
    inner_join(., calibration, by = c('sensor_conc', 'sensor_lot')) %>% 
    group_by(sample, date, conc, sensor_conc) %>% 
    mutate('product_conc' = as.integer(corrected_fluorescence)*cal_slope + cal_intercept) %>% 
    filter(! is.na(fluorescence)) %>% 
    ungroup()
  return(data)
}

# combine all the files into one tibble

files <- index %>% 
  pull(data_file) %>% unique()
### read in the data files, join them with the info from the index file and make them tidy
( dataset <- do.call("bind_rows", lapply(files, FUN = read_and_gather)) )
write_delim(dataset, parsed_dataset_file, delim='\t')


# plot individual files for DELA fitting
dataset %>% 
    select(condition, Time, fluorescence) %>% 
    mutate('group' = condition) %>% 
    group_by(group) %>% 
    nest() %>% 
    pwalk(~write_tsv(x = .y, file = paste0(dir_for_DELA, '/', .x, '.txt'))) %>% 
    unnest(data) %>% ungroup()


noGAP_blanks <- dataset %>% filter(GAP_conc == 0)

dataset <- dataset %>% filter(GAP_conc != 0)

final_product_conc <- dataset %>% 
  group_by(sample, well, date, conc, condition, row) %>% 
  summarize("final_product_conc" = max(product_conc, na.rm = T)) %>% 
  write_delim(calculated_product_conc_file, delim='\t')

# final_product_conc %>% filter(date == '20201117')   ##### filter by date of experiment to get data for today's fit

### intrinsic hydrolysis
intrinsic_hydrolysis <- noGAP_blanks %>% 
  separate(col = sample, into = c('PE', 'mutant'), remove = F) %>% 
  select(sample, mutant, conc, condition, sensor_conc, Time, product_conc,fluorescence, corrected_fluorescence, row, column, date) %>% 
  inner_join(., final_product_conc, by = c('sample', 'row', 'date', 'conc')) %>% 
  select('condition' = condition.x, everything(), -condition.y)
fit_intrinsic_hydrolysis <- function(data) {
  lm.fit <- lm(product_conc ~ Time, data = data)
  data$slope <- lm.fit$coefficients[2]
  data$intercept <- lm.fit$coefficients[1]
  data$rel_rate <- data$slope / data$final_product_conc # uM/s/uM [s-1]
  return(data)
}
intrinsic_hydrolysis <- intrinsic_hydrolysis %>% 
  group_by(condition) %>% 
  do(fit_intrinsic_hydrolysis(.)) %>% 
  ungroup()

# intrinsic_hydrolysis %>% 
#   filter(mutant %in% c('T34E')) %>% 
#   select(mutant, intercept, slope, rel_rate, final_product_conc, conc, date, condition) %>% 
#   unique() %>% 
#   ggplot(aes(x = rel_rate, y = final_product_conc, color = as.character(condition), shape = mutant)) + geom_point()

write_delim(intrinsic_hydrolysis, intrinsic_hydrolysis_outfile, delim='\t')
