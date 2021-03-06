library(tidyverse)
SubstrateDilute <- function(stock, work, volume, Mw) {
  molar.stock <- 1000*(stock/Mw)
  dilutionX <- molar.stock/work
  protein.V <- round(volume/dilutionX,1)
  buffer.V <- round(volume - protein.V,1)
  return(data.frame(protein.V, buffer.V, "dilution" = round(dilutionX, 2)))
}
outdir <- "~/gdrive/gsp1_dms/mathy/labnotebook/GAP_assays/20210119/"   ### the output you need will be in a reaction_mix_temp.txt file
GAPs <- c("SpRNA1")
reaction.volume <- 100 # ul
Ran.premix.volume <- 110
### Set up  experiment parameters  ### you edit this box of parameters for each experiment
#################################################
#################################################
plate <- 1
row <- "G"
date <- "20210119"
Rans <- c('F159L', 'F159W', 'F163L','F159L', 'F159W', 'F163L')
GAP.stock <- 10  # nM
GAP.working <- 1 # nM
Rans.stock.conc <- c(164, 44.8, 233, 210, 42.3, 267) # uM
#################################################
#################################################
sensor <- c('sensor')
# sensor.stock.conc <- 500
sensor.stock.conc <- 125
final.Ran.conc <- list()
final.Ran.conc[[1]] = c(5, 5)
final.Ran.conc[[2]] = c(5, 5)
final.Ran.conc[[3]] = c(5, 5)
final.Ran.conc[[4]] = c(5, 5)
final.Ran.conc[[5]] = c(5, 5)
final.Ran.conc[[6]] = c(5, 5)

GAP.V <- round(reaction.volume/(GAP.stock/GAP.working), 2)
Ran.premix.to.add <- reaction.volume - GAP.V

reaction.mix.table <- data.frame()
index.table <- data.frame()
well_count <- 0
for (i in 1:length(Rans)) {
  prot <- Rans[i]
  for (j in seq_along(final.Ran.conc[[i]])) {
    well_count <- well_count + 1
    fin.Ran.conc <- final.Ran.conc[[i]][j]
    if (fin.Ran.conc < 3) {
      sensor.working <- 10 # uM
    } else if (fin.Ran.conc < 10) {
      sensor.working <- 20
    } else if (fin.Ran.conc < 30) {
      sensor.working <- 50
    } else {
      sensor.working <- 70
    }
    Ran.premix.conc <- round(fin.Ran.conc * (1 + ( (reaction.volume - Ran.premix.to.add) / reaction.volume )), 2)
    stock.conc <- Rans.stock.conc[i]
    Ran.volume <- round(Ran.premix.volume/(stock.conc/Ran.premix.conc), 2)
    
    sensor.premix.conc <- round(sensor.working * (1 + ( (reaction.volume - Ran.premix.to.add) / reaction.volume )), 2)
    #sensor.V = reaction.volume/(sensor.stock.conc/sensor.working)
    sensor.V = round( Ran.premix.volume / (sensor.stock.conc / sensor.premix.conc), 2)
    
    
    
    Ran.premix.buffer.volume <- round(Ran.premix.volume - Ran.volume - sensor.V) 
    reaction.mix.table <- rbind(reaction.mix.table, data.frame(plate, "well" = paste0(row, well_count), prot, fin.Ran.conc, Ran.volume, 
                                                               sensor.V, Ran.premix.buffer.volume, 
                                                               GAP.stock, GAP.working, 
                                                               Ran.premix.to.add, GAP.V, date))
    index.table <- rbind(index.table, data.frame(plate, prot, "well" = paste0(row, well_count), fin.Ran.conc, GAP.working, date, "cutoff_time" = 3000))
  }
}
### (j+((i-1)*6))
reaction.mix.table
write.table(reaction.mix.table, file = paste0(outdir, "reaction_mix_temp.txt"), quote = F, row.names = F, sep = "\t")
(Ran.volume.sum <- reaction.mix.table %>%
    group_by(prot) %>% 
    summarise( "total_Ran_V" = sum(Ran.volume)) )
(sensor.volume.sum <- reaction.mix.table %>%
    group_by(prot) %>% 
    summarise( "total_diluted_125mM_sensor_V" = 1.05 * sum(sensor.V)) )   ### prepare 5 % extra diluted sensor
index.table
write.table(index.table, file = paste0(outdir, "index_table_temp.txt"), quote = F, row.names = F, sep = "\t")


