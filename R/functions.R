#' Import libraries
#'
#' @import tidyverse
#' @import rbacon
#' @import Bchron
#' @import clam
#' @import tidyr
#' @import tibble
#' @import Hmisc
#' @import plotrix
#' @import rlist
#' @import ggplot2
#' @importFrom Rdpack reprompt
#' @export filter_SISAL
#' @export AM_SISAL
#' @export merge_SISAL_chrono
#' @export eval
#' @export plot_sisal_eval
#'
NULL


#' Filter SISAL db for records of class I to IV.
#'
#' class I: Records that run without any additional manipulation: at most tractable reversals, no hiatuses, U-series dates \\
#' class II: Records that at most consist of tractable reversals, only include U-series dates but have at least one hiatus \\
#' class III: Any record that is not U-series dated \\
#' class IV: Records that cannot be evaluated, since they have too few or poor inconsistent dates (i.e. non-tractable reversals) \\
#'
#' @param m Number of class to filter for; m in[1,2,3,4]
#' @return A vector of all records belonging to the chosen class.
filter_SISAL <- function(m, dating. = SISAL.AM::dating, sample. = SISAL.AM::sample, entity. = SISAL.AM::entity, hiatus. = SISAL.AM::hiatus){
  dating_tb <- dating. %>% left_join(.,entity., by = "entity_id") %>% filter(entity_status == 'current') %>%
    mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
  sample_tb <- sample. %>% left_join(., hiatus., by = "sample_id") %>% left_join(., entity., by = "entity_id") %>% filter(entity_status == 'current') %>%
    mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample), as.numeric)

  hiatus_tb <- sample_tb %>% filter(hiatus == 'H') %>% distinct(entity_id)
  reversals <- dating_tb %>% filter(date_used == 'yes') %>%
    select(entity_id, depth_dating, corr_age, corr_age_uncert_pos, corr_age_uncert_neg) %>%
    group_by(entity_id) %>%
    arrange(depth_dating, .by_group = T) %>%
    mutate(diff = lead(corr_age)-corr_age) %>%
    mutate(reversal = if_else(!is.na(diff) & diff < 0, TRUE, FALSE)) %>%
    group_by(entity_id) %>%
    arrange(depth_dating, .by_group = TRUE) %>%
    mutate(tractable = if_else(reversal & (abs(corr_age-lead(corr_age)) < (corr_age_uncert_pos + lead(corr_age_uncert_neg))), TRUE, FALSE))

  reversal_tractable <- reversals %>% select(entity_id, tractable) %>% group_by(entity_id) %>% filter(tractable) %>% distinct(entity_id)

  reversal_nontractable <- reversals %>% filter(reversal) %>%filter(!(entity_id %in% reversal_tractable$entity_id)) %>% distinct(entity_id)
  nr_dates <- dating_tb %>% filter(date_used == 'yes') %>% select(entity_id, corr_age) %>% group_by(entity_id) %>%count() %>% filter(n>=3)
  no_depth_sample <- sample_tb %>% group_by(entity_id) %>% summarise(depths = if_else(all(is.na(depth_sample)), FALSE, TRUE)) %>% filter(!depths)

  UTh_dates <- dating_tb %>% filter(date_used == 'yes') %>% filter(date_type == "MC-ICP-MS U/Th" | date_type == "TIMS" | date_type ==  "ICP-MS U/Th Other" | date_type == "Alpha U/Th" | date_type ==  "U/Th unspecified")  %>% distinct(entity_id)
  not_UTh_dates <- dating_tb %>% filter(date_used == 'yes') %>% filter(date_type == 'Event; start of laminations' | date_type == 'Event; end of laminations' | date_type == 'C14' | date_type =='Multiple methods' | date_type =='other') %>%
    distinct(entity_id) %>% filter(!(entity_id %in% UTh_dates$entity_id))


  if(m == 1){
    run <- dating_tb %>% distinct(entity_id) %>%
      filter(entity_id %in% nr_dates$entity_id) %>%
      filter(!(entity_id %in% no_depth_sample$entity_id)) %>%
      filter(!(entity_id %in% hiatus_tb$entity_id)) %>%
      filter(!(entity_id %in% not_UTh_dates$entity_id)) %>%
      filter(!(entity_id %in% reversal_nontractable$entity_id)) %>%
      arrange(., entity_id)
  }

  if(m == 2){
    run <- dating_tb %>% distinct(entity_id) %>%
      filter(entity_id %in% nr_dates$entity_id) %>%
      filter(!(entity_id %in% no_depth_sample$entity_id)) %>%
      filter(entity_id %in% hiatus_tb$entity_id) %>%
      filter(!(entity_id %in% not_UTh_dates$entity_id)) %>%
      filter(!(entity_id %in% reversal_nontractable$entity_id)) %>%
      arrange(., entity_id)
  }

  if(m == 3){
    run <- not_UTh_dates %>% arrange(., entity_id)
  }

  if(m == 4){
    run <- dating_tb %>% distinct(entity_id) %>%
      filter(!(entity_id %in% nr_dates$entity_id) | entity_id %in% reversal_nontractable$entity_id)
  }

  return(run)
}

#' Write age model input data
#'
#' @param entid entity_id
#' @param bacon Logical. A parameter specifiying if Bacon AM is to be executed.
#' @param bchron Logical. A parameter specifiying if Bchron AM is to be executed.
#' @param stalage Logical. A parameter specifiying if StalAge AM is to be executed. #'
#' @param linInterp Logical. A parameter specifiying if lin. Interp. AM is to be executed.
#' @param linReg Logical. A parameter specifiying if lin. Reg. AM is to be executed.
#' @param working_directory File path, where the age model files are to be saved.
#' @param site. SISAL site table
#' @param entity. SISAL entity table
#' @param entity_link_reference. SISAL entity_link_reference table
#' @param reference. SISAL reference table
#' @param notes. SISAL notes table
#' @param sample. SISAL sample table
#' @param hiatus. SISAL hiatus table
#' @param gap. SISAL gap table
#' @param original_chronology. SISAL original_chronology table
#' @param sisal_chronology. SISAL sisal_chronology table
#' @param d13C. SISAL d13C table
#' @param d18O. SISAL d18O table
#' @return file_name. The combined entity_id and entity_name. Writes the files modifed and prepared to run the AM.
#' @example write_files(237, site_tb, dating_tb, sample_tb, bacon = T, '~/Documents/runAM)
write_files <- function(entid, bacon = F, bchron = F, stalage = F, linInterp = F, linReg = F, dating. = SISAL.AM::dating, working_directory, site. = SISAL.AM::site, entity. = SISAL.AM::entity, entity_link_reference. = SISAL.AM::entity_link_reference, reference. = SISAL.AM::reference, notes. = SISAL.AM::notes,
                        sample. = SISAL.AM::sample, hiatus. = SISAL.AM::hiatus, gap. = SISAL.AM::gap, original_chronology. = SISAL.AM::original_chronology, sisal_chronology. = SISAL.AM::sisal_chronology, d13C. = SISAL.AM::d13C, d18O. = SISAL.AM::d18O){

  site_tb <- left_join(site., entity., by = 'site_id') %>% left_join(., entity_link_reference., by = 'entity_id') %>%
    left_join(., reference., by = 'ref_id') %>% left_join(., notes., by = 'site_id') %>% mutate_at(vars(site_id, entity_id), as.numeric)
  dating_tb <- dating. %>%
    mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric) %>%
    group_by(entity_id) %>%mutate(laminar_dated = if_else((entity_id %in% dating_lamina$entity_id), 'yes', 'no')) %>% ungroup()
  sample_tb <- plyr::join_all(list(sample.,hiatus., gap., original_chronology., sisal_chronology., d13C., d18O.), by = 'sample_id', type = 'left', match = 'all') %>%
    mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, COPRA_age,
                   COPRA_age_uncert_pos, COPRA_age_uncert_neg, linear_age, linear_age_uncert_pos, linear_age_uncert_neg, d13C_measurement,
                   d13C_precision, d18O_measurement, d18O_precision), as.numeric)

  entity_from_base <- site_tb %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
  sample_from_base <- sample_tb %>% filter(entity_id %in% entity_from_base$entity_id) %>% select(entity_id,depth_sample) %>% group_by(entity_id) %>% summarise(max = max(depth_sample))

  dating_from_base <- full_join(dating_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>%
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_dating, NA_real_)) %>%
    mutate(depth_dating = if_else(!is.na(depth_conv), depth_conv, depth_dating)) %>%
    select(-depth_conv) %>% arrange(., depth_dating, .by_group = T)

  sampling_from_base <- full_join(sample_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>%
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>%
    mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
    select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)

  dating_tb_new <- dating_from_base %>% group_by(entity_id) %>%arrange(depth_dating) %>%
    mutate(depth_dating_new = depth_dating/10,
      corr_age_uncert = (corr_age_uncert_pos + corr_age_uncert_neg)/4,
      thickness_new = if_else(is.na(dating_thickness), 0.5, dating_thickness/20)) %>% ungroup() %>%
    mutate(calib_curve_new = case_when(
      calib_used == 'NULL' ~ 'normal',
      calib_used == 'unknown' ~ 'unknown',
      calib_used == "not calibrated" ~ 'ask again',
      calib_used == "INTCAL13 NH" ~ 'intcal13')) %>%
    mutate(cc = case_when(
      date_type == 'C14' ~ 1,
      date_type != 'C14' ~ 0 ))

  sample_tb_new <- sampling_from_base %>% group_by(entity_id) %>% mutate(depth_sample_new = depth_sample/10) %>% ungroup()

  original_chrono <- sample_tb_new %>% filter(entity_id == entid) %>%
    arrange(., depth_sample) %>%
    mutate(interp_age = if_else(is.na(interp_age) & !is.na(COPRA_age), COPRA_age, interp_age),
           interp_age_uncert_pos =  if_else(is.na(interp_age_uncert_pos) & !is.na(COPRA_age_uncert_pos), COPRA_age_uncert_pos, interp_age_uncert_pos),
           interp_age_uncert_neg =  if_else(is.na(interp_age_uncert_neg) & !is.na(COPRA_age_uncert_neg), COPRA_age_uncert_neg, interp_age_uncert_neg),
           age_model_type = if_else(is.na(age_model_type), rep("COPRA", length(interp_age)), age_model_type)) %>%
    select(sample_id, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, age_model_type)

  entity_name <- unique((site_tb %>% filter(entity_id == entid))$entity_name)
  file_name <- paste(entid,"-",entity_name, sep = '')
  setwd(working_directory)
  dir.create(file.path(working_directory,file_name))
  setwd(file.path(working_directory, file_name))
  #write.csv(unique((site_tb %>% filter(entity_id == entid))$notes), "notes.txt", row.names = F, col.names = F)
  write.csv(site_tb %>% filter(entity_id == entid) %>% select(site_id, site_name, latitude, longitude, entity_id, entity_name, depth_ref, speleothem_type, contact, data_DOI_URL, ref_id, citation, publication_DOI, notes), 'info.csv', row.names = F, col.names = T)
  write.csv(sample_tb_new %>% filter(entity_id == entid & hiatus =='H') %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'hiatus.csv', row.names = F, col.names = T)
  write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(depth_dating, corr_age, corr_age_uncert_pos, corr_age_uncert_neg, date_type) %>% arrange(., depth_dating), 'used_dates.csv', row.names = F, col.names = T)
  write.csv(dating_tb_new %>% filter(entity_id == entid  & date_type != 'Event; hiatus') %>% filter(date_used == 'no' | date_used == 'unknown') %>% select(depth_dating, corr_age) %>% arrange(., depth_dating), 'not_used_dates.csv', row.names = F, col.names = T)
  write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample, interp_age, d13C_measurement, d18O_measurement) %>% arrange(., depth_sample), "proxy_data.csv", row.names = F, col.names = T)
  write.csv(original_chrono, "original_chronology.csv", row.names = F, col.names = T)

  if (bacon) {
    setwd(file.path(working_directory, file_name))
    dir.create('Bacon_runs')
    setwd(file.path(getwd(),'/Bacon_runs'))
    dir.create(file_name)
    setwd(file.path(working_directory, file_name,'Bacon_runs', file_name))

    hiatus_bacon <- sample_tb_new %>% filter(entity_id == entid & hiatus =='H') %>% arrange(depth_sample) %>% mutate(depth_sample_bacon = depth_sample/10)  %>% select(sample_id, depth_sample_bacon)
    id <- sample_tb_new %>% filter(entity_id == entid & is.na(hiatus) & is.na(gap)) %>% select(sample_id, depth_sample_new) %>% arrange(., depth_sample_new)

    write.csv(id$sample_id, 'sample_id.csv', row.names = F, col.names = F)
    write.table(id$depth_sample_new, paste(entid,'-',entity_name,'_depths.txt', sep = ''), row.names = F, col.names = F)
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating_new, cc) %>% arrange(., depth_dating_new), paste(entid,'-',entity_name,'.csv', sep = ''), row.names = F, col.names = T)
    write.csv(hiatus_bacon,'hiatus_bacon.csv', row.names = F, col.names = T)
  }

  if (bchron){
    setwd(file.path(working_directory, file_name))
    dir.create('Bchron')
    setwd(file.path(getwd(),'/Bchron'))

    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>%
                select(dating_id, corr_age, corr_age_uncert, depth_dating_new, thickness_new, calib_curve_new) %>% arrange(., depth_dating_new), 'ages.csv', row.names = F, col.names = T)
    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample_new) %>% arrange(., depth_sample_new), 'depths.csv', row.names = F, col.names = T)
  }

  if (stalage) {
    setwd(file.path(working_directory, file_name))
    dir.create('StalAge')
    setwd(file.path(getwd(),'/StalAge'))

    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating) %>% arrange(., depth_dating), 'ages.csv', row.names = F, col.names = T)
    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
  }

  if (linInterp) {
    setwd(file.path(working_directory, file_name))
    dir.create('linInterp')
    setwd(file.path(getwd(),'/linInterp'))

    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>%
                arrange(., depth_dating),'ages.csv', row.names = F, col.names = T)
  }

  if (linReg) {
    setwd(file.path(working_directory, file_name))
    dir.create('linReg')
    setwd(file.path(getwd(),'/linReg'))

    write.csv(sample_tb_new %>% filter(entity_id == entid & is.na(hiatus) & is.na(gap)) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'depths.csv', row.names = F, col.names = T)
    write.csv(sample_tb_new %>% filter(entity_id == entid) %>% select(sample_id, depth_sample) %>% arrange(., depth_sample), 'id.csv', row.names = F, col.names = T)
    write.csv(dating_tb_new %>% filter(entity_id == entid & date_used == 'yes' & date_type != 'Event; hiatus') %>% select(dating_id, corr_age, corr_age_uncert, depth_dating, date_type) %>%
                arrange(., depth_dating),'ages.csv', row.names = F, col.names = T)
  }

  return(file_name)
}

#' Add artificial hiatus date to the input dating files for StalAge, Bchron and lin. Interp.
#'
#' @param data Input dating information.
#' @param hiatus_tb Table containing the sample_id and depth of each hiatus.
#' @param stalage Logical. TRUE if the artificial hiatus ages is to be added to the Stalage input file.
#' @param bchron Logical. TRUE if the artificial hiatus ages is to be added to the Bchron input file.
#' @param linInterp Logical. TRUE if the artificial hiatus ages is to be added to the lin. Interp. input file.
#' @return The modified date file.
#' @example add_hiatus(dating_tb, hiatus, bchron=T)
add_hiatus <- function(data, hiatus_tb, stalage = F, bchron =F, linInterp=F) {

  age <- unlist(data[,2])
  depth_dating <- unlist(data[,4])
  #uncert <- unlist(data[,3])
  x_out <- unlist(hiatus_tb$depth_sample)

  e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')

  if(linInterp){
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus_tb$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus_tb$depth_sample)),depth_dating = e$x, date_type = rep('Hiatus',length(hiatus_tb$depth_sample)))
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)))
  }

  if(stalage){
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus_tb$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus_tb$depth_sample)),depth_dating = e$x)
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)))
  }

  if(bchron){
    x_out <- x_out/10
    e <- approx(x = depth_dating, y = age, xout = x_out, method = 'linear')
    h <- data.frame(dating_id = hiatus_tb$sample_id,corr_age = e$y, corr_age_uncert = rep(NA, length(hiatus_tb$depth_sample)),depth_dating_new = e$x, thickness_new = rep(NA, length(hiatus_tb$depth_sample)), calib_curve_new = 'normal')
    new <- rbind(data, h)
    new <- new %>% arrange(., corr_age) %>% mutate(corr_age_uncert = if_else(is.na(corr_age_uncert), ((lead(corr_age_uncert)+lead(corr_age)-corr_age)+(corr_age-(lag(corr_age)-lag(corr_age_uncert))))/2, as.double(corr_age_uncert)),
                                                   thickness_new = if_else(is.na(thickness_new), 0.01, as.double(thickness_new)),
                                                   calib_curve_new = if_else(is.na(calib_curve_new), 'normal', calib_curve_new))

  }

  return(new)
}

#' Determine median and quantiles for Bacon MC ensemble.
#'
#' @param depth_eval Sample depths for which interpolated age is wanted.
#' @param hiatus_tb Table containing the hiatus depths.
#' @param bacon_mcmc Bacon MC ensemble
#' @param q1 quantile 1
#' @param q2 quantile 2
#' @return Table containing the sample depths, median ages and uncertainties.
#'
#' @example get_bacon_median_quantile(depth_tb, hiatus_depth, bacon_mcmc, 0.05, 0.95)
get_bacon_median_quantile <- function(depth_eval, hiatus_tb, bacon_mcmc,q1 = 0.05, q2=0.95) {
  #bacon_mcmc <- sapply(depth_eval, Bacon.Age.d)
  bacon_age <- apply(bacon_mcmc,2,median)
  bacon_quantile <- apply(bacon_mcmc, 2, function(x){quantile(x, probs = c(q1,q2), na.rm = T)})

  data <- cbind(depth_eval, bacon_age, bacon_age_uncert_pos = bacon_quantile[2,]-bacon_age,
                bacon_age_uncert_neg = bacon_age - bacon_quantile[1,])
  h <- data.frame(depth_eval = hiatus_tb$depth_sample_bacon, bacon_age = replicate(dim(hiatus_tb)[1],NA),
                  bacon_age_uncert_pos = replicate(dim(hiatus_tb)[1],NA),
                  bacon_age_uncert_neg = replicate(dim(hiatus_tb)[1],NA))
  data <- rbind(data, h)
  data <- data[order(data[,1]),]

  return(data)
}

#' Linear interpolation for single MC simulation.
#'
#' @param data Table containing single age ensemble and dating depths
#' @param depth_eval Sample depths.
#' @param hiatus_tb Table conatining hiatus depths and sample_id's
#' @return Interpolated ages for sample depths.
#' @example get_lin_interp(cbind(depth_dating, age_ensemble[,1]), depth_sample, hiatus)
get_lin_interp <- function(data, depth_eval, hiatus_tb) {

  age <- unlist(data[,2])
  depth_dating <- unlist(data[,1])
  x_out <- unlist(depth_eval)

  e <- approxExtrap(x = depth_dating, y = age, xout = x_out)
  linInterp <- as_tibble(data.frame(depth_eval = unlist(e$x), lin_interp_age = unlist(e$y)))

  if (!plyr::empty(data.frame(hiatus_tb))) {
    linInterp <- linInterp %>% rowwise() %>% mutate(lin_interp_age = if_else(depth_eval %in% hiatus_tb, NA_real_, lin_interp_age))
  }

  return(linInterp)
}

#' Linear regression slopes and intecepts for single MC simulation
#'
#' @param data Table containing single age ensemble and dating depths
#' @param hiatus_tb Table conatining hiatus depths and sample_id's
#' @return Table containing slopes and interceptions for each section.
linear_regression <- function(data, hiatus_tb){ # data = c("depth","age")

  # initialize
  j <- length(hiatus_tb)
  m <- list(list())

  # calculate slope and intercept for each section between hiati

  idx <- data[,1] < hiatus_tb[1]
  m[[1]] <- lm(data[idx,2]~data[idx,1])
  if(is.null(dim(data[idx,]))) {m[[1]] <- NA}

  if (j>=2) {
    for (i in seq(from = 2, to = j)){
      idx <- (data[,1] < hiatus_tb[i]) & (data[,1] > hiatus_tb[i-1])
      m[[i]] <- lm(data[idx,2]~data[idx,1])
      if(is.null(dim(data[idx,]))) {m[[i]] <- NA}
    }

    idx <- data[,1] > hiatus_tb[length(hiatus_tb)]
    m[[j+1]] <- lm(data[idx,2]~data[idx,1])
    if(is.null(dim(data[idx,]))) {m[[j+1]] <- NA}

  } else {
    idx <- data[,1] > hiatus_tb[length(hiatus_tb)]
    m[[j+1]] <- lm(data[idx,2]~data[idx,1])
    if(is.null(dim(data[idx,]))) {m[[j+1]] <- NA}
  }

  # retrun data
  return(m)
}

#' Interpolate ages using lin. reg. for single MC simulation; with hiatus.
#'
#' @param m Table containing slopes and interceptions for each section.
#' @param depth_eval Sample depths.
#' @param hiatus_tb Table containing hiatus depths and sample_id's.
#' @return Lin. regression fitted ages for sample depths.
lin_reg_ages <- function(m, depth_eval, hiatus_tb) {

  # initialize
  d <- length(unlist(depth_eval))
  j <- length(m)
  lin_reg_age <- replicate(d, 0)

  for (i in seq(1,j)) {
    if(!is.na(m[[i]])) {
      for (k in c(1,2)) {
        if (is.na(m[[i]]$coefficients[[k]])) {
          m[[i]]$coefficients[[k]] = 0
        }
      }
    }
  }

  # add column with ages to depth table
  depth <- cbind(depth_eval, lin_reg_age)

  idx <- depth[,1] < hiatus_tb[1]
  if(is.na(m[[1]])){
    depth[idx,2] <- NA
  } else {
    depth[idx,2]<- m[[1]]$coefficients[[1]] + depth[idx,1]*m[[1]]$coefficients[[2]]
  }


  if (j>2) {
    for (i in seq(2,length(hiatus_tb))){
      idx <- (depth[,1] < hiatus_tb[i]) & (depth[,1] > hiatus_tb[i-1])
      #depth[idx,2]<- m[[i]]$coefficients[[1]] + depth[idx,1]*m[[i]]$coefficients[[2]]
      if(is.na(m[[i]])){
        depth[idx,2] <- NA
      } else {
        depth[idx,2]<- m[[i]]$coefficients[[1]] + depth[idx,1]*m[[i]]$coefficients[[2]]
      }
    }

    idx <- depth[,1] > hiatus_tb[length(hiatus_tb)]
    #depth[idx,2]<- m[[length(hiatus)+1]]$coefficients[[1]] + depth[idx,1]*m[[length(hiatus)+1]]$coefficients[[2]]
    if(is.na(m[[1]])){
      depth[idx,2] <- NA
    } else {
      depth[idx,2]<- m[[length(hiatus_tb)+1]]$coefficients[[1]] + depth[idx,1]*m[[length(hiatus_tb)+1]]$coefficients[[2]]
    }
  } else {
    idx <- depth[,1] > hiatus_tb[length(hiatus_tb)]
    #depth[idx,2]<- m[[length(hiatus)+1]]$coefficients[[1]] + depth[idx,1]*m[[length(hiatus)+1]]$coefficients[[2]]
    if(is.na(m[[1]])){
      depth[idx,2] <- NA
    } else {
      depth[idx,2]<- m[[length(hiatus_tb)+1]]$coefficients[[1]] + depth[idx,1]*m[[length(hiatus_tb)+1]]$coefficients[[2]]
    }
  }
  data <- data.frame(depth_eval = depth[,1], lin_reg_age = depth[,2])
  #data <- data[-which(depth_eval == hiatus),]
  h <- data.frame(depth_eval = hiatus_tb, lin_reg_age = replicate(length(hiatus_tb),NA))
  data <- rbind(data, h)
  data <- data[order(data[,1]),]

  return(data)

}

#' Interpolate ages using lin. reg. for single MC simulation; no hiatus.
#'
#' @param data Table containing dates and dating depths.
#' @param depth_eval Sample depths.
#' @return Lin. regression fitted ages for sample depths.
lin_reg_no_hiatus <- function(data, depth_eval) {
  m_lR<- lm(data[,2]~data[,1])

  age_lR <- replicate(length(unlist(depth_eval)), 0)
  depth <- cbind(depth_eval, age_lR)
  depth[,2]<- m_lR[1]$coefficients[[1]] + depth[,1]*m_lR[1]$coefficients[[2]]
  return(depth)

}

#' Determine median and quantiles for MC ensemble; not Bacon.
#'
#' @param upd Table of Mc ensemble.
#' @param q1 quantile 1
#' @param q2 quantile 2
#' @return Table containing the sample depths, median ages and uncertainties.
#'
#' @example get_median_quantile(mcmc, 0.05, 0.95)
get_median_quantiles <- function(upd,q1,q2){
  age_median <- apply(upd, 1, median)
  age_sd <- apply(upd, 1, function(x){quantile(x, probs = c(q1,q2), na.rm = T)})

  return(cbind(age_median, t(age_sd)))
}

#' Generate lin. reg. ensemble.
#'
#' @param N Number of interations.
#' @param hiatus_tb Table containing hiatus depths and hiatus sample_ids.
#' @param depth_dating Dating depths.
#' @param age_ensemble Mc ensemble of dating table.
#' @param depth_sample Sample depths.
#' @return Lin. reg. ensemble of N interations for sample depths.
mc_linReg <- function(N, hiatus_tb, depth_dating, age_ensemble, depth_sample) {
  #sample <- rep(NA, N)
  for (i in 1:N){
    if (i == 1) {
      if (plyr::empty(data.frame(hiatus_tb))){
        age_mc <- lin_reg_no_hiatus(cbind(depth_dating, age_ensemble[i,]), depth_sample)
      }else{
        m <- linear_regression(cbind(depth_dating, age_ensemble[i,]), hiatus_tb) # hiatus = hiatus[[1]]
        age_mc <- linear_regression_ages(m, depth_sample, hiatus_tb)
      }
      sample_ensemble <- age_mc
    } else {
      if (plyr::empty(data.frame(hiatus_tb))){
        age_mc <- lin_reg_no_hiatus(cbind(depth_dating, age_ensemble[i,]), depth_sample)
      } else {
        m <- linear_regression(cbind(depth_dating, age_ensemble[i,]), hiatus_tb)
        age_mc <- linear_regression_ages(m, depth_sample, hiatus_tb)
      }
      sample_ensemble <- cbind(sample_ensemble,age_mc[,2])
    }
  }
  return(sample_ensemble)
}

#' Generate lin. interp. ensemble.
#'
#' @param N Number of interations.
#' @param hiatus_tb Table containing hiatus depths and hiatus sample_ids.
#' @param depth_dating Dating depths.
#' @param age_ensemble Mc ensemble of dating table.
#' @param depth_sample Sample depths.
#' @return Lin. interp. ensemble of N interations for sample depths.
mc_linInt <- function(N, hiatus_tb, depth_dating, age_ensemble, depth_sample){
  for(j in 1:N){
    if(j == 1){
      age_mc <- get_lin_interp(cbind(depth_dating, age_ensemble[j,]), depth_sample,hiatus_tb)
      sample_ensemble <- age_mc
      #print(dim(age))
    } else {
      age_mc <- get_lin_interp(cbind(depth_dating, age_ensemble[j,]), depth_sample,hiatus_tb)
      #print(dim(age))
      sample_ensemble <- cbind(sample_ensemble,age_mc[,2])
    }
  }
  return(sample_ensemble)

}

#' Generate age ensembles depending on AM method.
#'
#' @param linReg Logical. TRUE if MC ensemble is to be generated for lin. reg.
#' @param linInterp Logival. TRUE if MC ensemble is to be generated for lin. interp.
#' @param age Dates.
#' @param age_error Uncertainties to input dates.
#' @param N Number of interations.
#'
#' @return Lin. reg. ensemble of N interations for sample depths.
mc_ensemble <- function(linReg = F, linInterp = F, age, age_error, N, working_directory, file_name) { # N number of MC simulations
  number <- 0
  d <- 0
  age_ensemble_final <- NA

  if(linReg) {
    age_ensemble <- apply(cbind(age,age_error), 1, function(x) rnorm(n = N, mean = x[1], sd = x[2])) # calculate N deviates for each age
    age_ensemble_final <- age_ensemble
    return(age_ensemble_final)
  }


  if(linInterp){
    while(d < N){
      #print('hello')
      k<- N-d
      #print(k)

      number <- number+1
      age_ensemble <- apply(cbind(age,age_error), 1, function(x) rnorm(n = k, mean = x[1], sd = x[2])) # calculate N deviates for each age

      if(k==1){
        age_ensemble_diff <- diff(age_ensemble)
        run <- T

        if (any(age_ensemble_diff < 0 & !is.na(age_ensemble_diff))){
          run <- F
          #print('k=1')
        }

      } else {
        age_ensemble_diff <- apply(age_ensemble, 1, diff) # each column contains the derivatives for one MC run -> N columns
        del <- NULL #age_ensemble_copy <- age_ensemble
        run <- T
        #if (k ==2) {print('k=2')}
        for (i in seq(1,k)) {
          #print(age_ensemble_diff[,i])
          if(any(age_ensemble_diff[,i] <0)){
            if (is.null(dim(age_ensemble)[1])) {
              #age_ensemble <- age_ensemble[-i,]
              run <- F
            } else {
              del <- c(del,i)
              #print(age_ensemble)
            }

          }

        }
        if(!is.null(del)){age_ensemble <- age_ensemble[-del,]}
        #if(k==1 && any(diff(age_ensemble)<0)){age_ensemble<-NULL}

      }

      if (run) {
        if (k == N){
          age_ensemble_final <- age_ensemble
        } else {
          age_ensemble_final <- rbind(age_ensemble_final, age_ensemble)
        }
      }

      run <- T
      d <- dim(age_ensemble_final)[1]
      #diff_final <- apply(age_ensemble_final, 1, diff)
      #print(any(diff_final<0))
      if(is.null(d)){d <- 1}
      if(number>2000 && d > 100){return(age_ensemble_final)
        break
      } else if(number > 2000 && d < 100) {
        setwd(file.path(working_directory, file_name, '/linInterp'))
        write.csv(c(number,d),'mc_fail.csv', row.names = F)
        stop('ERROR: too many iterations, check data!')}
      #print(c('Dim lI:', dim(age_ensemble_linInt_final)[1]))
      #print(c('Dim lR:', dim(age_ensemble_linReg_final)[1]))
    }
    return(age_ensemble_final) # ensemble
  }

}


#' Braodly scan dating input for reversals/outliers and increase errors at tagged depths.
#'
#'#' This is a modified version of scan_fine() to work with the input data of the lin. interp.
#'
#' @param dating_tb Table containing dates, errors and depths.
#' @return Modified dating file.
#'
#' @references Scholz, D. and Hoffmann, D. L., Quaternary Geochronology 6, 369-382 (2011)
scan_lin_interp<-function(dating_tb) {

  library(Hmisc)
  #library(rpanel)
  library(plotrix)

  dating_id <- dating_tb$dating_id
  age <- dating_tb$corr_age
  depth <- dating_tb$depth_dating
  error <- dating_tb$corr_age_uncert
  date_type <- dating_tb$date_type

  hilf<-0

  test<-array(NA, c(length(depth), length(depth)))	#Test-Array f?r Screening


  age<-age[order(depth)]		#Sortiert die Vektoren aufsteigend nach der Tiefe
  error<-error[order(depth)]
  depth<-depth[order(depth)]


  for (i in 1:length(depth)) {	#Schleife zum Testen auf Inversionen. Wenn Inversion, schreibe 1 in Matrix.

    j<-i+1

    while (j<=length(depth)) {

      if (age[j]+error[j]<age[i]-error[i]) test[i, j]<-1 else test[i, j]<-0	#wenn 2-sigma Fehler nicht ?berlappen, schreibe 1 ins Feld

      j<-j+1

    }

  }

  #print(test)

  max<-1		#Variable f?r das Maximum der Inversionen
  max_pos<-0	#Variable f?r den Punkt des Maximums


  repeat {

    for (i in 1:length(depth)) {

      if (sum(test[i,], na.rm=TRUE)>max) {	#Gehe Zeilen durch ... ## NA z?hlen nicht da na.rm = T

        max<-sum(test[i,], na.rm=TRUE)
        max_pos<-i

      }

      if (sum(test[,i], na.rm=TRUE)>max) {	#Gehe Spalten durch ...

        max<-sum(test[,i], na.rm=TRUE)
        max_pos<-i

      }

    }

    if (max>1) {

      if (max(test[max_pos,], na.rm=TRUE)==1) { ## was soll es sonst sein????? es muss mindestens 2 einser haben, da max>1

        hilf<-age[max_pos]

        for (i in (max_pos+1):length(depth)) {

          if (is.na(test[max_pos, i])) next ### wann kann das auftreten ???? wir gehen hier reihen durch

          if (test[max_pos, i]==1 && age[i]<hilf) hilf<-age[i]	#; print(hilf)

          error[max_pos]<-age[max_pos]-hilf

        }

      }

      if (max(test[, max_pos], na.rm=TRUE)==1) {

        hilf<-age[max_pos]

        for (i in 1:(max_pos-1)) {

          if (is.na(test[i, max_pos])) next

          if (test[i, max_pos]==1 && age[i]>hilf) hilf<-age[i]	#; print(hilf)

          error[max_pos]<-hilf-age[max_pos]

        }

      }

      test[max_pos,]<-NA
      test[,max_pos]<-NA


    }

    if (max==1) break

    max<-1

  }


  test[is.na(test)]<-0			#Ersetze NA's in Matrix durch Nullen.


  depth<-depth[!is.na(age)]		#L?sche rausgeworfene Punkte
  error<-error[!is.na(age)]
  age<-age[!is.na(age)]


  Daten<-data.frame(dating_id = dating_id, corr_age = age, corr_age_uncert = error, depth_dating=depth, date_type = date_type)	#Speichere Punkte in Dataframe

  return(Daten)

}

#' Fine scan dating input for reversals/outliers and increase errors at tagged depths.
#'
#' This is a modified version of scan_fine() to work with the input data of the lin. interp.
#'
#' @param dating_tb Table containing dates, errors and depths.
#' @return Modified dating file.
#'
#' @references Scholz, D. and Hoffmann, D. L., Quaternary Geochronology 6, 369-382 (2011)
scan_fine_lin_interp<-function(dating_tb){

  #attach(Daten)
  dating_id <- dating_tb$dating_id
  age <- dating_tb$corr_age
  depth <- dating_tb$depth_dating
  error <- dating_tb$corr_age_uncert
  date_type <- dating_tb$date_type


  Anteil<-0							#Vektoren f?r das Testen der Fits
  counter<-0
  part<-0

  part[1]<-1							#Definition des Beteiligungs-Vektors
  part[length(depth)]<-1
  part[2]<-2
  part[length(depth)-1]<-2

  for (i in 3:(length(depth)-2)) part[i]<-3

  #print(part)


  test<-0								#Pruefvariable f?r ?bergeordneten Test
  pruef<-0							#Pruefvariable f?r andere Tests


  while (test==0) {						#Schleife f?r Iteration bis Fehler alle passen

    test<-1								#Setzen der allgemeinen Pr?fvariable auf 1

    for (i in 1:length(depth)) Anteil[i]<-0				#Nullsetzen von Anteil
    for (i in 1:length(depth)) counter[i]<-0			#Nullsetzen von counter


    for (i in 1:(length(depth)-2)) {

      #print(i)


      fit<-lm(age[i:(i+2)] ~ depth [i:(i+2)], weights=1/error[i:(i+2)])	#Fit ?ber 3-Punkt-Intervall
      attr(fit$coefficients, "names")<-NULL

      #curve(fit$coefficients[2]*x+fit$coefficients[1], from=depth[i], to=depth[i+2], add=TRUE, col=i)	#Plotten des Fits

      age_fit<-0							#Bestimmung des Alters des Fits an den jeweiligen Punkten
      age_fit[i:(i+2)]<-fit$fitted.values


      pruef<-0							#Null-Setzen der Pr?fvariable

      for (k in i:(i+2)) {						#Test ob linearer Fit m?glich

        if (age_fit[k]>age[k]+error[k] || age_fit[k]<age[k]-error[k]) pruef<-pruef+1	#wenn Punkt au?erhalb der Fehlergrenzen, z?hle Pruefvariable um 1 hoch

      }

      #print (pruef)

      if (pruef>0) {

        counter[i:(i+2)]<-counter[i:(i+2)]+1				#wenn der Fit nicht geht, setze counter um 1 hoch

      } else {

        if (fit$coefficients[2]<0) {					#wenn die Steigung des Fits negativ ist

          if (slope(depth[i:(i+2)], age[i:(i+2)], error[i:(i+2)])<0.2) counter[i:(i+2)]<-counter[i:(i+2)]+1	#wenn weniger als 30% der Fits eine positive Steigung haben, setze counter um 1 hoch

        }

      }

    }

    #print(counter)

    Anteil<-counter/part

    #print(Anteil)

    for (i in 1:length(depth)) {				#gehe Anteil nach 1ern durch

      if (Anteil[i]==1) {					#wenn Anteil=1, dann vergr??ere Fehler um 10% und setze allgemeine Pr?fvariable auf 1

        error[i]<-error[i]*1.1
        test<-0

      }

    }

    #errbar(depth, age, age+error, age-error, add=TRUE)		#Plotten der Daten

  }

  #x11()
  #errbar(depth, age, age+error, age-error)			#Plotten der Daten
  #title(main="Age data screened for minor outliers")

  Daten<-data.frame(dating_id = dating_id, corr_age = age, corr_age_uncert = error, depth_dating=depth, date_type = date_type)

  #detach(Daten)

  return(Daten)

}
