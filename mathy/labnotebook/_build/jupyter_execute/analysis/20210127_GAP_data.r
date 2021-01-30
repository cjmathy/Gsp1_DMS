library(tidyverse)
library(readxl)
library(lubridate)

# set theme for plotting

# Set font sizes
SMALL_SIZE = 10
MEDIUM_SIZE = 11
BIG_SIZE = 12

# SMALL_SIZE = 6
# MEDIUM_SIZE = 7
# BIG_SIZE = 8

theme_custom <- theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = BIG_SIZE),
    axis.title = element_text(size = BIG_SIZE),
    axis.text = element_text(size = BIG_SIZE),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
#     legend.position = 'right',
#     legend.spacing.y = unit(0.01, 'cm'),
#     legend.box = 'horizontal',
#     legend.box.just = 'left',
#     legend.text = element_text(size = SMALL_SIZE),
#     legend.title = element_text(size = SMALL_SIZE),
#     legend.margin = margin(t = 0, unit='cm'),
    axis.line = element_line(size = 0.1),
    strip.text.x = element_text(size = SMALL_SIZE)
  )

# options(repr.plot.width=20, repr.plot.height=5)

10.69340436/0.2728471828

df

domnegs = c('F28V','H32D','H32I','G35D','H50N','P51A',
            'P51E','F54A','F54K','N156A','F159L','F163L')
tolerant = c('F28Y','H32L','H32E','G35A','H50R','P51G',
             'F54W','N156W','F159W','F163Y')


# read in GAP data
df <- read_delim('GAP_computed_kcat_Km/GAP_table_raw.txt', delim='\t', col_types=cols())
WT <- filter(df, mutant=='WT')

df %>% 
    mutate('label' = ifelse(mutant %in% domnegs, 'DN', 'TOL')) %>% 
    ggplot(aes(mutant, mean_kcat_Km, fill=label)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes(ymin = mean_kcat_Km - se, ymax =  mean_kcat_Km + se), width = 0.25) +
    scale_fill_manual(values = c('salmon','black')) +
    theme_custom +
    coord_cartesian(ylim = c(0, 100))

df_fit6 <- read_csv('../../../Data/6gen_fitness_current.csv', col_types=cols()) %>% 
    select(mutant, score)


df %>% 
    mutate('label' = ifelse(mutant %in% domnegs, 'DN', 'TOL')) %>% 
    left_join(df_fit6, by='mutant') %>% 
    ggplot(aes(score, mean_kcat_Km, color=label)) +
    geom_point(size=5) +
    geom_text(aes(label=mutant), hjust=0.75, vjust=1.5) +
    scale_color_manual(values = c('salmon','black')) +
    theme_custom

