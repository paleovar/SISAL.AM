### packages to load
library(tidyverse)

## packages that have to be installed but should not be loaded 
# plyr

#####################################
####### example for plotting ########
#####################################

final.plot(eID = c(1,2,3,4), fp_input = '~/SISAL_Data/sisalv2_20200130/', fp_output = '~/SISAL_Data/plots/')

### eID: vector of entity_id's for which plots are to be executed
### fp_input: file path to SISAL files
### fp_output: file path to folder where plots are to be saved




#########################################
###### functions for plotting ###########
#########################################

#' function to load sisal files
#' @param fiel_path file path to folder where SISAL csv files are saved 
#' @param prefix prefix to SISAL csv file names; no prefix vor SISALv1b and SISALv2 files required
#' 
#' @return SISAL database object
read.SISAL.files <- function(file_path, 
                             prefix = '' 
){
  
  # read in files
  composite_link_entity <- read.csv(file.path(file_path, paste('/',prefix, 'composite_link_entity.csv', sep = '')), header = T,stringsAsFactors = F)
  d13C <- read.csv(file.path(file_path,paste('/',prefix, 'd13c.csv',sep='')),header = T, stringsAsFactors = F)
  d13C <- rename(d13C, iso_std_d13C = iso_std)
  d18O <- read.csv(file.path(file_path,paste('/',prefix, 'd18o.csv', sep ='')),header = T, stringsAsFactors = F)
  d18O <- rename(d18O, iso_std_d18O = iso_std)
  dating_lamina <- read.csv(file.path(file_path,paste('/',prefix, 'dating_lamina.csv', sep = '')), header = T, stringsAsFactors = F) %>% 
    mutate_at(vars(depth_lam, lam_thickness, lam_age_uncert_pos, lam_age_uncert_neg), as.numeric)
  dating <- read.csv(file.path(file_path,paste('/',prefix, 'dating.csv',sep = '')), header = T, stringsAsFactors = F) %>%
    mutate_at(vars(depth_dating, dating_thickness,min_weight, max_weight, uncorr_age, uncorr_age_uncert_pos, uncorr_age_uncert_neg, starts_with('X'), corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
  entity_link_reference <- read.csv(file.path(file_path,paste('/',prefix, 'entity_link_reference.csv', sep = '')), header =T, stringsAsFactors = F)
  entity <- read.csv(file.path(file_path,paste('/',prefix, 'entity.csv', sep = '')), header = T, stringsAsFactors = F) %>% mutate_at(vars(cover_thickness, distance_entrance), as.numeric)
  gap <- read.csv(file.path(file_path,paste('/',prefix, 'gap.csv', sep = '')), header = T, stringsAsFactors = F)
  hiatus <- read.csv(file.path(file_path,paste('/',prefix, 'hiatus.csv', sep ='')), header = T, stringsAsFactors = F)
  notes <- read.csv(file.path(file_path,paste('/',prefix, 'notes.csv', sep = '')), header = T, stringsAsFactors = F)
  original_chronology <- read.csv(file.path(file_path,paste('/',prefix, 'original_chronology.csv', sep = '')), header = T, stringsAsFactors = F) %>% mutate_at(vars(interp_age_uncert_pos, interp_age_uncert_neg), as.numeric)
  reference <- read.csv(file.path(file_path,paste('/',prefix, 'reference.csv', sep = '')), header = T, stringsAsFactors = F)
  sample <- read.csv(file.path(file_path,paste('/',prefix, 'sample.csv', sep = '')), header = T, stringsAsFactors = F) %>% mutate_at(vars(sample_thickness, depth_sample), as.numeric)
  sisal_chronology <- read.csv(file.path(file_path,paste('/',prefix, 'sisal_chronology.csv', sep = '')), header = T, stringsAsFactors = F) %>% mutate_at(vars(everything()), as.numeric)
  site <- read.csv(file.path(file_path,paste('/',prefix, 'site.csv', sep = '')), header = T, stringsAsFactors = F) %>% mutate_at(vars(elevation), as.numeric)
  
  # correct for depth_ref 'from base'
  entity_from_base <- entity %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
  sample_from_base <- sample %>% filter(entity_id %in% entity_from_base$entity_id) %>% select(entity_id,depth_sample) %>% group_by(entity_id) %>% summarise(max = max(depth_sample))
  
  dating_from_base <- full_join(dating, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_dating, NA_real_)) %>% 
    mutate(depth_dating = if_else(!is.na(depth_conv), depth_conv, depth_dating)) %>%
    select(-depth_conv) %>% arrange(., depth_dating, .by_group = T)
  
  sampling_from_base <- full_join(sample, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>% 
    mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
    select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)
  
  # create sisal object
  d <- list()
  d$composit_link_entity <- composite_link_entity
  d$d13C <- d13C
  d$d18O <- d18O
  d$dating_lamina <- dating_lamina
  d$dating <- dating_from_base
  d$entity_link_reference <- entity_link_reference
  d$entity <- entity
  d$gap <- gap
  d$hiatus <- hiatus
  d$notes <- notes
  d$original_chronology <- original_chronology
  d$reference <- reference
  d$sample <- sampling_from_base
  d$sisal_chronology <- sisal_chronology
  d$site <- site
  
  return(d)
}


#' function to plot age models
#' 
#' @param chrono RData object of one SISAL entity
#' @param .entity_id entity_id of the SISAL entity
#' @param dating. SISAL csv dating file
#' @param entity_name entity name of the SISAL entity
#' 
#' @return x_lim for further plotting and plots all available age models in one graph
plot_am <- function(chrono,.entity_id, dating., entity_name){
  not.age <- dating. %>% filter(entity_id == .entity_id) %>% filter(date_used == 'no' | date_used == 'unknown') %>%
    mutate(corr_age_uncert = 0.5*(corr_age_uncert_pos + corr_age_uncert_neg)) %>%
    mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
                          if_else(date_type == "TIMS", 8,
                                  if_else(date_type == "ICP-MS U/Th Other", 4,
                                          if_else(date_type == "Alpha U/Th", 5,
                                                  if_else(date_type ==  "U/Th unspecified",7, 
                                                          if_else(date_type == 'Event; start of laminations' |
                                                                    date_type == 'Event; end of laminations', 3, 
                                                                  if_else(date_type == 'C14', 12, NA_real_)))))))) 
  
  age <- dating. %>% filter(entity_id == .entity_id) %>% filter( date_used == 'yes') %>%
    mutate(corr_age_uncert = 0.5*(corr_age_uncert_pos + corr_age_uncert_neg)) %>% 
    mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
                          if_else(date_type == "TIMS", 8,
                                  if_else(date_type == "ICP-MS U/Th Other", 4,
                                          if_else(date_type == "Alpha U/Th", 5,
                                                  if_else(date_type ==  "U/Th unspecified",7, NA_real_))))))
  

  color <- RColorBrewer::brewer.pal(12,"Paired")[c(1:5,8:10)]
  
  n <- rep(0,8)
  if(!plyr::empty(data.frame(chrono$origAM$interp_age))){n[1]<- 1}
  if(!plyr::empty(data.frame(chrono$copRa$copRa_age))){n[2]<- 2}
  if(!plyr::empty(data.frame(chrono$StalAge$StalAge_age))){n[3]<- 3}
  if(!plyr::empty(data.frame(chrono$linInterp$lin_interp_age))){n[4]<- 4}
  if(!plyr::empty(data.frame(chrono$linReg$lin_reg_age))){n[5]<- 5}
  if(!plyr::empty(data.frame(chrono$Bacon$Bacon_age))){n[6]<- 6}
  if(!plyr::empty(data.frame(chrono$Bchron$Bchron_age))){n[7]<- 7}
  if(!plyr::empty(data.frame(chrono$OxCal$OxCal_age))){n[8]<- 8}
  
  n <- n[which(n !=0)]
  
  matplot(x = cbind(chrono$origAM$interp_age,chrono$copRa$copRa_age, chrono$StalAge$StalAge_age, chrono$linInterp$lin_interp_age, chrono$linReg$lin_reg_age, chrono$Bacon$Bacon_age, 
                    chrono$Bchron$Bchron_age, chrono$OxCal$OxCal_age), 
          y= chrono$origAM$depth_sample, 
          col = color[n], lty = 1, type = 'l', lwd = 1.5, 
          xlim = range(not.age$corr_age, not.age$corr_age +not.age$corr_age_uncert, not.age$corr_age -not.age$corr_age_uncert, 
                       age$corr_age, age$corr_age +age$corr_age_uncert, age$corr_age -age$corr_age_uncert, 
                       chrono$origAM$interp_age, 
                       chrono$copRa$copRa_age, 
                       chrono$StalAage$StalAge_age,
                       chrono$linInterp$lin_interp_age,
                       chrono$linReg$lin_reg_age,
                       chrono$Bacon$Bacon_age,
                       chrono$Bchron$Bchron_age,
                       chrono$OxCal$OxCal_age,
                       na.rm=TRUE),
          ylim = c(max(range(chrono$origAM$depth_sample, not.age$depth_dating, age$depth_dating, na.rm = T), na.rm = T),min(range(chrono$origAM$depth_sample, not.age$depth_dating, age$depth_dating, na.rm = T), na.rm = T)), xlab = '', ylab = '')
  mtext(side = 1, line = 2, text = 'Median age [yr BP]')
  mtext(side = 2, line = 2, text = 'Depth from top [mm]')
  
  if(!plyr::empty(data.frame(not.age))){
    points(x = not.age$corr_age, y=not.age$depth_dating, lty = 2, col = 'red', pch = not.age$.pch)
    arrows(not.age$corr_age-not.age$corr_age_uncert, not.age$depth_dating, not.age$corr_age+not.age$corr_age_uncert, not.age$depth_dating, length=0.05, angle=90, code=3, col = 'red')
  }
  
  points(x = age$corr_age, y=age$depth_dating, lty = 2, col = 'black', pch = age$.pch)
  arrows(age$corr_age-age$corr_age_uncert, age$depth_dating, age$corr_age+age$corr_age_uncert, age$depth_dating, length=0.05, angle=90, code=3, col = 'black')
  
  if (any(age$date_type=="Event; actively forming")){
    abline(h=age$depth_dating[which(age$date_type=="Event; actively forming")],lty=3,col="orange")
  }
  if (!plyr::empty(data.frame(chrono$hiatus))) {
    abline(h =  chrono$hiatus$depth_sample, col = 'grey', lty = 3)
  }
  
  datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
  pch.datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
  
  if(plyr::empty(data.frame(not.age))){
    legend("topright",legend = c(paste(datetype, '- used')),
           pch = c(pch.datetype),cex=1,
           col=c(rep('black', length(datetype))),ncol=1,title="Date Type",bty="n")
  }
  
  
  
  if(!plyr::empty(data.frame(not.age))){
    not.datetype <- (not.age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
    pch.not.datetype <- (not.age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
    legend("topright",legend = c(paste(datetype, '- used'), paste(not.datetype, '- not used/unkown')),
           pch = c(pch.datetype, pch.not.datetype),cex=1,
           col=c(rep('black', length(datetype)), rep('red', length(not.datetype))),ncol=1,title="Date Type",bty="n")
  }
  
  mtext(line = 0, text = 'Age-depth model')
  
  return(range(not.age$corr_age, not.age$corr_age +not.age$corr_age_uncert, not.age$corr_age -not.age$corr_age_uncert, 
               age$corr_age, age$corr_age +age$corr_age_uncert, age$corr_age -age$corr_age_uncert, 
               chrono$origAM$interp_age, 
               chrono$copRa$copRa_age, 
               chrono$StalAge$StalAge_age,
               chrono$linInterp$lin_interp_age,
               chrono$linReg$lin_reg_age,
               chrono$Bacon$Bacon_age,
               chrono$Bchron$Bchron_age,
               chrono$OxCal$OxCal_age,
               na.rm=TRUE))
}

#' function to plot age uncertainties as interquartile ranges
#' 
#' @param chrono RData object of one SISAL entity
#' @param .entity_id entity_id of the SISAL entity
#' @param dating. SISAL csv dating file
#' @param entity_name entity name of the SISAL entity
#' 
#' @return plots all available age uncertainties in one graph
plot_iqr <- function(chrono,.entity_id, dating., entity_name){
  age <- dating. %>% filter(entity_id == .entity_id) %>% filter( date_used == 'yes') %>%
    mutate(corr_age_uncert = 0.5*(corr_age_uncert_pos + corr_age_uncert_neg)) 
  
  color <- RColorBrewer::brewer.pal(12,"Paired")[c(1:5,8:10)]
  
  n <- rep(0,7)
  if(!plyr::empty(data.frame(chrono$copRa$copRa_age))){n[1]<- 1}
  if(!plyr::empty(data.frame(chrono$StalAge$StalAge_age))){n[2]<- 2}
  if(!plyr::empty(data.frame(chrono$linInterp$lin_interp_age))){n[3]<- 3}
  if(!plyr::empty(data.frame(chrono$linReg$lin_reg_age))){n[4]<- 4}
  if(!plyr::empty(data.frame(chrono$Bacon$Bacon_age))){n[5]<- 5}
  if(!plyr::empty(data.frame(chrono$Bchron$Bchron_age))){n[6]<- 6}
  if(!plyr::empty(data.frame(chrono$OxCal$OxCal_age))){n[7]<- 7}
  
  n <- n[which(n !=0)]
  
  if(.entity_id == 1){
    matplot(x = cbind(chrono$copRa$copRa_age_uncert_pos + chrono$copRa$copRa_age_uncert_neg, chrono$StalAge$StalAge_age_uncert_pos + chrono$StalAge$StalAge_age_uncert_neg, 
                      chrono$linInterp$lin_interp_age_uncert_pos + chrono$linInterp$lin_interp_age_uncert_neg, chrono$linReg$lin_reg_age_uncert_pos + chrono$linReg$lin_reg_age_uncert_neg,
                      chrono$Bacon$Bacon_age_uncert_pos + chrono$Bacon$Bacon_age_uncert_neg, chrono$Bchron$Bchron_age_uncert_pos+chrono$Bchron$Bchron_age_uncert_neg,
                      chrono$OxCal$OxCal_age_uncert_pos + chrono$OxCal$OxCal_age_uncert_neg),
            y= chrono$origAM$depth_sample, col = color[n], lty = 1, type = 'l', lwd = 1, 
            xlim = c(10,max(chrono$copRa$copRa_age_uncert_pos + chrono$copRa$copRa_age_uncert_neg, chrono$StalAge$StalAge_age_uncert_pos + chrono$StalAge$StalAge_age_uncert_neg, 
                         chrono$linInterp$lin_interp_age_uncert_pos + chrono$linInterp$lin_interp_age_uncert_neg, chrono$linReg$lin_reg_age_uncert_pos + chrono$linReg$lin_reg_age_uncert_neg,
                         chrono$Bacon$Bacon_age_uncert_pos + chrono$Bacon$Bacon_age_uncert_neg, chrono$Bchron$Bchron_age_uncert_pos+chrono$Bchron$Bchron_age_uncert_neg,
                         chrono$OxCal$OxCal_age_uncert_pos + chrono$OxCal$OxCal_age_uncert_neg, na.rm=TRUE)),
            ylim = c(max(range(chrono$origAM$depth_sample, age$depth_dating, na.rm = T), na.rm = T),min(range(chrono$origAM$depth_sample, age$depth_dating, na.rm = T), na.rm = T)),
            xlab = '', ylab = '', log = 'x')
  } else {
    matplot(x = cbind(chrono$copRa$copRa_age_uncert_pos + chrono$copRa$copRa_age_uncert_neg, chrono$StalAge$StalAge_age_uncert_pos + chrono$StalAge$StalAge_age_uncert_neg, 
                      chrono$linInterp$lin_interp_age_uncert_pos + chrono$linInterp$lin_interp_age_uncert_neg, chrono$linReg$lin_reg_age_uncert_pos + chrono$linReg$lin_reg_age_uncert_neg,
                      chrono$Bacon$Bacon_age_uncert_pos + chrono$Bacon$Bacon_age_uncert_neg, chrono$Bchron$Bchron_age_uncert_pos+chrono$Bchron$Bchron_age_uncert_neg,
                      chrono$OxCal$OxCal_age_uncert_pos + chrono$OxCal$OxCal_age_uncert_neg),
            y= chrono$origAM$depth_sample, col = color[n], lty = 1, type = 'l', lwd = 1, 
            xlim = range(chrono$copRa$copRa_age_uncert_pos + chrono$copRa$copRa_age_uncert_neg, chrono$StalAge$StalAge_age_uncert_pos + chrono$StalAge$StalAge_age_uncert_neg, 
                         chrono$linInterp$lin_interp_age_uncert_pos + chrono$linInterp$lin_interp_age_uncert_neg, chrono$linReg$lin_reg_age_uncert_pos + chrono$linReg$lin_reg_age_uncert_neg,
                         chrono$Bacon$Bacon_age_uncert_pos + chrono$Bacon$Bacon_age_uncert_neg, chrono$Bchron$Bchron_age_uncert_pos+chrono$Bchron$Bchron_age_uncert_neg,
                         chrono$OxCal$OxCal_age_uncert_pos + chrono$OxCal$OxCal_age_uncert_neg, na.rm=TRUE),
            ylim = c(max(range(chrono$origAM$depth_sample, age$depth_dating, na.rm = T), na.rm = T),min(range(chrono$origAM$depth_sample, age$depth_dating, na.rm = T), na.rm = T)),
            xlab = '', ylab = '', log = 'x')
  }
  mtext(side = 1, line = 2, text = 'IQR [yr]')
  abline(h = age$depth_dating, lty = 3, col = 'black')
  
  if (any(age$date_type=="Event; actively forming")){
    abline(h=age$depth_dating[which(age$date_type=="Event; actively forming")],lty=3,col="orange")
  }
  if (!plyr::empty(data.frame(chrono$hiatus))) {
    abline(h =  chrono$hiatus$depth_sample, col = 'grey', lty = 3)
  }
  
  
  mtext(line = 0, text = 'Age IQR')
  
}

#' function to plot existing isotope measurements
#'
#' @param chrono RData object of one SISAL entity
#' @param .entity_id entity_id of the SISAL entity
#' @param entity_name entity name of the SISAL entity
#' @param x.lim x-axis limits 
#' 
#' @return plots all available isotope measurements in one graph
plot_isotopes <- function(chrono, .entity_id, entity_name, x.lim){
  iso<- chrono$proxy
  if(!all(is.na(iso$interp_age))){
    # print('1')
    if(all(is.na(iso$d18O_measurement))){
      matplot(x = iso$interp_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
              xlim = x.lim,
              ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
      mtext(side = 1, line = 2, text = 'Original AM age [yr BP]')
      axis(side = 4)
      mtext(side = 4, line = 3, text = expression(paste(delta^{13},'C [\u2030]')))
      legend('topleft', legend = c(expression(paste(delta^{13},'C'))), col = c('goldenrod2'),
             lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      
    } else {
      matplot(x = iso$interp_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
              xlim = x.lim,
              ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
      mtext(side = 2, line = 2, text = expression(paste(delta^{18},'O [\u2030]')))
      mtext(side = 1, line = 2, text = 'Original AM age [yr BP]')
      if(!all(is.na(iso$d13C_measurement))){
        par(new = T)
        plot(x = iso$interp_age, y = iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5,
             xlim = x.lim, ylim = range(iso$d13C_measurement, na.rm = T), xlab = NA, ylab = NA, axes = F)
        axis(side = 4)
        mtext(side = 4, line = 3, text = expression(paste(delta^{13},'C [\u2030]')))
        legend('topleft', legend = c(expression(paste(delta^{18},'O')),expression(paste(delta^{13},'C'))), col = c('dodgerblue','goldenrod2'),
               lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      } else {
        legend('topleft', legend = c(expression(paste(delta^{18},'O'))), col = c('dodgerblue'),
               lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      }
    }
    
  }else{
    if(all(is.na(iso$d18O_measurement))){
      if(!plyr::empty(data.frame(chrono$copRa))){
          matplot(x = chrono$copRa$copRa_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'copRa median age [yr BP]')
      } else if(!all(is.na(chrono$linInterp$lin_interp_age))) {
        matplot(x = chrono$linInterp$lin_interp_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Linear interp. median age [yr BP]')
      } else if(!all(is.na(chrono$Bchron$Bchron_age))) {
        matplot(x = chrono$Bchron$Bchron_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Bchron median age [yr BP]')
      } else if(!all(is.na(chrono$Bacon$Bacon_age))) {
        matplot(x = chrono$Bacon$Bacon_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Bacon median age [yr BP]')
      }
      axis(side = 4)
      mtext(side = 4, line = 3, text = expression(paste(delta^{13},'C [\u2030]')))
      legend('topleft', legend = c(expression(paste(delta^{13},'C'))), col = c('goldenrod2'),
             lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
    } else {
      if(!all(is.na(chrono$copRa$copRa_age))){
          matplot(x = chrono$copRa$copRa_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                  xlim = x.lim,
                  ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'copRa median age [yr BP]')
        } else if(!all(is.na(chrono$linInterp$lin_interp_age))) {
        matplot(x = chrono$linInterp$lin_interp_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Linear interp. median age [yr BP]')
      } else if(!all(is.na(chrono$linReg$lin_reg_age))) {
        matplot(x = chrono$linReg$lin_reg_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Linear reg. median age [yr BP]')
      } else if(!all(is.na(chrono$Bchron$Bchron_age))) {
        matplot(x = chrono$Bchron$Bchron_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Bchron median age [yr BP]')
      } else if(!all(is.na(chrono$Bacon$Bacon_age))) {
        matplot(x = chrono$Bacon$Bacon_age, y= iso$d18O_measurement, col = 'dodgerblue', lty = 1, type = 'l', lwd = 1.5, 
                xlim = x.lim,
                ylim = range(iso$d18O_measurement, na.rm = T), xlab = '', ylab = '')
        mtext(side = 1, line = 2, text = 'Bacon median age [yr BP]')
      }
      
      mtext(side = 2, line = 2, text = expression(paste(delta^{18},'O [\u2030]')))
      if(!all(is.na(iso$d13C_measurement))){
        par(new = T)
        if(!all(is.na(chrono$copRa$copRa_age))){
            plot(x = chrono$copRa$copRa_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
                 xlim = x.lim,
                 ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
            mtext(side = 1, line = 2, text = 'copRa median age [yr BP]')
          } else if(!all(is.na(sc$linear_age))) {
          plot(x = chrono$linInterp$lin_interp_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
               xlim = x.lim,
               ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'Linear interp. median age [yr BP]')
        } else if(!all(is.na(chrono$linReg$lin_reg_age))) {
          plot(x = chrono$linReg$lin_reg_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
               xlim = x.lim,
               ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'Linear reg. median age [yr BP]')
        } else if(!all(is.na(chrono$Bchron$Bchron_age))) {
          plot(x = chrono$Bchron$Bchron_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
               xlim = x.lim,
               ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'Bchron median age [yr BP]')
        } else if(!all(is.na(chrono$Bacon$Bacon_age))) {
          plot(x = chrono$Bacon$Bacon_age, y= iso$d13C_measurement, col = 'goldenrod2', lty = 1, type = 'l', lwd = 1.5, 
               xlim = x.lim,
               ylim = range(iso$d13C_measurement, na.rm = T), xlab = '', ylab = '')
          mtext(side = 1, line = 2, text = 'Bacon median age [yr BP]')
        }
        axis(side = 4)
        mtext(side = 4, line = 3, text = expression(paste(delta^{13},'C [\u2030]')))
        legend('topleft', legend = c(expression(paste(delta^{18},'O')),expression(paste(delta^{13},'C'))), col = c('dodgerblue','goldenrod2'),
               lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      } else {
        legend('topleft', legend = c(expression(paste(delta^{18},'O'))), col = c('dodgerblue'),
               lty = 1, title = 'Information',bty = 'n', lwd = 1.5)
      }
    }
    
  }
  
  
  mtext(line = 0, text = 'Isotopes')
  
}

#' function to plot dating information for detrital control
#'
#' @param chrono RData object of one SISAL entity
#' @param .entity_id entity_id of the SISAL entity
#' @param dating. SISAL csv dating file
#' 
#' @return plots the available U/Th dating information in one graph
plot_dating <- function(chrono,.entity_id, dating.){
  age <- dating. %>% filter(entity_id == .entity_id) %>% filter( date_used == 'yes' & date_type != 'Event; hiatus' & date_type != 'Event; actively forming') %>%
    mutate_at(vars(corr_age, corr_age_uncert_pos, corr_age_uncert_neg,uncorr_age, uncorr_age_uncert_pos, uncorr_age_uncert_neg), as.numeric) %>% 
    mutate(.pch = if_else(date_type == "MC-ICP-MS U/Th", 2,
                          if_else(date_type == "TIMS", 8,
                                  if_else(date_type == "ICP-MS U/Th Other", 4,
                                          if_else(date_type == "Alpha U/Th", 5,
                                                  if_else(date_type ==  "U/Th unspecified",7, NA_real_))))))
  
  if(!all(is.na(age$uncorr_age))){
    if(!all(is.na(age$uncorr_age_uncert_pos))){
      plot(x = age$uncorr_age, y= age$corr_age, col = 'black', pch = age$.pch, type = 'p', 
           xlim = range(age$uncorr_age + age$uncorr_age_uncert_pos, age$uncorr_age - age$uncorr_age_uncert_neg,
                        na.rm=TRUE),
           ylim = range(age$corr_age + age$corr_age_uncert_pos,age$corr_age - age$corr_age_uncert_neg,
                        na.rm = T), xlab = '', ylab = '')
      mtext(side = 1, line = 2, text = 'Uncorrected age [yr BP]')
      mtext(side = 2, line = 2, text = 'Corrected age [yr BP]')
      arrows(age$uncorr_age,age$corr_age-age$corr_age_uncert_neg, age$uncorr_age,age$corr_age+age$corr_age_uncert_pos, length=0.05, angle=90, code=3, col = 'black')
      arrows(age$uncorr_age + age$uncorr_age_uncert_pos,age$corr_age, age$uncorr_age - age$uncorr_age_uncert_neg,age$corr_age, length=0.05, angle=90, code=3, col = 'black')
      
      datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
      pch.datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
      
      abline(a = 0, b = 1, lty = 2)
      
      legend("topleft",legend = c(paste(datetype,'age'), 'Angle bissectrice'),
             pch = c(pch.datetype, NA),cex=1, lty = c(rep(NA,length(datetype)), 2),
             col=c(rep('black', length(datetype))),ncol=1,title="Information",bty="n")
    } else {
      plot(x = age$uncorr_age, y= age$corr_age, col = 'black', pch = age$.pch, type = 'p', 
           xlim = range(age$uncorr_age,
                        na.rm=TRUE),
           ylim = range(age$corr_age + age$corr_age_uncert_pos,age$corr_age - age$corr_age_uncert_neg,
                        na.rm = T), xlab = '', ylab = '')
      mtext(side = 1, line = 2, text = 'Uncorrected age [yr BP]')
      mtext(side = 2, line = 2, text = 'Corrected age [yr BP]')
      arrows(age$uncorr_age,age$corr_age-age$corr_age_uncert_neg, age$uncorr_age,age$corr_age+age$corr_age_uncert_pos, length=0.05, angle=90, code=3, col = 'black')
      
      datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(date_type))$date_type
      pch.datetype <- (age %>% filter(!(date_type %in% c("Event; hiatus","Event; actively forming"))) %>% distinct(.pch))$.pch
      
      abline(a = 0, b = 1, lty = 2)
      
      legend("topleft",legend = c(paste(datetype,'age'), 'Angle bissectrice'),
             pch = c(pch.datetype, NA),cex=1, lty = c(rep(NA,length(datetype)), 2),
             col=c(rep('black', length(datetype))),ncol=1,title="Information",bty="n")
    }
    
    mtext(line = 0, text = 'Evaluation of the detrital contribution')
  } else {
    plot.new()
    legend("center", legend = 'No uncorrected age', col = 'black', bty = 'n')
  }
}


#' function to produce overview plot 
#' 
#' @param eID vector of entity_id's for which plots are to be executed, i.e. c(1,2,3,4)
#' @param fp_input file path to sisal files, i.e. '~/SISAL_Data/sisalv2_20200130/'
#' @param fp_output file path to folder where plots are to be saved, i.e. '~/SISAL_Data/plots/'
#' 
#' @return creates the final overview plot of all available dating and age model information
final.plot<- function(eID,
                      fp_input, 
                      fp_output){ 
  
  eId <- eID
  sisal <- read.SISAL.files(fp_input,'')
  sisal_new <- left_join(sisal$sample %>% select(entity_id, sample_id, depth_sample), sisal$original_chronology, by = 'sample_id') %>% 
    left_join(., sisal$sisal_chronology, by = 'sample_id') %>% ungroup() %>% left_join(., sisal$d18O, by = 'sample_id') %>% left_join(., sisal$d13C, by = 'sample_id') %>% left_join(.,sisal$hiatus, by = 'sample_id')
  
  for(i in eId) {
    
    print(i)
    
    entity_name <- (sisal$entity %>% filter(entity_id == i))$entity_name
    graphics.off()
    
    ## build RData Object
    df_fil <- list()
  
    df_fil$origAM <- sisal_new %>% filter(entity_id == i) %>% select(interp_age, depth_sample, age_model_type)
    df_fil$linReg <- sisal_new %>% filter(entity_id == i) %>% select(lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg)
    df_fil$linInterp <- sisal_new %>% filter(entity_id == i) %>% select(lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg)
    df_fil$copRa <- sisal_new %>% filter(entity_id == i) %>% select(copRa_age, copRa_age_uncert_pos, copRa_age_uncert_neg)
    df_fil$StalAge <- sisal_new %>% filter(entity_id == i) %>% select(StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg)
    df_fil$Bacon <- sisal_new %>% filter(entity_id == i) %>% select(Bacon_age, Bacon_age_uncert_pos, Bacon_age_uncert_neg)
    df_fil$Bchron <- sisal_new %>% filter(entity_id == i) %>% select(Bchron_age, Bchron_age_uncert_pos, Bchron_age_uncert_neg)
    df_fil$OxCal <- sisal_new %>% filter(entity_id == i) %>% select(OxCal_age, OxCal_age_uncert_pos, OxCal_age_uncert_neg)
    df_fil$proxy <- sisal_new %>% filter(entity_id == i) %>% select(sample_id, depth_sample, interp_age, d13C_measurement, d18O_measurement) 
    df_fil$hiatus <- sisal_new %>% filter(entity_id == i & hiatus == 'H') %>% select(depth_sample)
    
    
    cairo_pdf(paste(fp_output,'/',i,'-',entity_name,'.pdf',sep = ''),12, 9)
    layout(matrix(c(1,1,2,3,3,4,5,5,5), nrow=3, ncol=3, byrow = T))
    par(mar = c(4,3,3,4), oma = c(1,3,2,1))
    x.lim <- plot_am(df_fil, i, sisal$dating, entity_name)
    if(any(c(!is_empty(data.frame(df_fil$Bacon)),!is_empty(data.frame(df_fil$Bchron)),!is_empty(data.frame(df_fil$copRa)),!is_empty(data.frame(df_fil$StalAge)),!is_empty(data.frame(df_fil$linInterp)),
             !is_empty(data.frame(df_fil$linReg)),!is_empty(data.frame(df_fil$OxCal))))){
      plot_iqr(df_fil, i, sisal$dating, entity_name)
      mtext(paste(i,'-', entity_name, sep = ''), side = 3, line = -2, outer = TRUE)
    } else {
      plot.new()
      mtext(paste(i,'-', entity_name, sep = ''), side = 3, line = -2, outer = TRUE)
      plot.new()
      plot.new()
    }
    plot_isotopes(df_fil, i, entity_name, x.lim)
    plot_dating(df_fil, i, sisal$dating)
    plot.new()
    
    orig_am <- df_fil$origAM %>% filter(age_model_type != 'NA') %>% distinct(age_model_type)
    
    legend('center',legend=c(paste('original AM:', orig_am),'copRa median age','StalAge mean age','lin. interp. median age','lin. reg. median age',
                             'Bacon median age','Bchron median age', 'OxCal median age', "Event; actively forming", "Event; hiatus"),
           lty=c(1,1,1,1,1,1,1,1,3,3),cex=1.5,
           col = c(RColorBrewer::brewer.pal(12,"Paired")[c(1:5,8:10)],'orange','grey'),
           title=paste("Age model Information"),bty = "n",  ncol = 3, lwd = 4)
    
    dev.off()
    
    
  }
  
}


