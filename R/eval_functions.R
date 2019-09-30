#' Build SISAL chronology.
#'
#' Builds the SISAL chronology for successfully executed age models. Supply runFile to method and chronology is build for the included entity ids.
#'
#' @param runFile Tabel containing entity ids and working_directory where age model chrnologies are saved.
#' @param entity SISAL v1b entity csv file
#' @param sample SISAL v1b sample csv file
#' @param dating SISAL v1b dating csv file
#'
#' @return SISAL chronolgy containg NA for the missing  age models
#' @example
#' m <- merge_sisal_chrono(runFile)
#' r <- as.data.frame(m[1])
#' s <- as.data.frame(m[2])
merge_SISAL_chrono <- function(runfile, .dating = SISAL.AM::dating, .sample = SISAL.AM::sample, .entity = SISAL.AM::entity){

  entity_from_base <-  .entity %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
  sample_from_base <- .sample %>% filter(entity_id %in% entity_from_base$entity_id) %>% select(entity_id,depth_sample) %>% group_by(entity_id) %>% summarise(max = max(depth_sample))
  sampling_from_base <- full_join(.sample, sample_from_base, by = 'entity_id') %>%
    mutate_at(vars(max),as.numeric) %>%
    group_by(entity_id) %>%
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>%
    mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
    select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)
  dating_from_base <- full_join(.dating, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>%
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_dating, NA_real_)) %>%
    mutate(depth_dating = if_else(!is.na(depth_conv), depth_conv, depth_dating)) %>%
    select(-depth_conv) %>% arrange(., depth_dating, .by_group = T)

  se <- left_join(sampling_from_base, .entity, by = 'entity_id')


  SISAL_chronology <- data.frame(entity_id = se$entity_id, sample_id = se$sample_id, depth_sample = se$depth_sample,
                                     lin_reg_age = rep(NA_real_, length(sample$sample_id)),
                                     lin_reg_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                     lin_reg_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                     lin_interp_age = rep(NA_real_, length(sample$sample_id)),
                                     lin_interp_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                     lin_interp_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                     StalAge_age = rep(NA_real_, length(sample$sample_id)),
                                     StalAge_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                     StalAge_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                     bacon_age = rep(NA_real_, length(sample$sample_id)),
                                     bacon_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                     bacon_age_uncert_neg = rep(NA_real_, length(sample$sample_id)),
                                     bchron_age = rep(NA_real_, length(sample$sample_id)),
                                     bchron_age_uncert_pos = rep(NA_real_, length(sample$sample_id)),
                                     bchron_age_uncert_neg = rep(NA_real_, length(sample$sample_id)))

  for(i in runFile$entity_id) {
    print(i)
    y <- runFile %>% filter(entity_id == i)
    entity_name <- (entity %>% filter(entity_id == i))$entity_name
    file_name <- paste(i,"-",entity_name, sep = '')

    if(y$linReg){
      setwd(file.path(y$working_directory, file_name, '/linReg'))
      if(file.exists('linReg_chronology.csv')){
        lR_chrono <- read.csv('linReg_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
        names(lR_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(lR_chrono, SISAL_chronology_new, by = 'sample_id') %>%
          mutate(lin_reg_age = if_else(sample_id %in% lR_chrono$sample_id, age, lin_reg_age),
                 lin_reg_age_uncert_pos = if_else(sample_id %in% lR_chrono$sample_id, uncert_pos, lin_reg_age_uncert_pos),
                 lin_reg_age_uncert_neg = if_else(sample_id %in% lR_chrono$sample_id, uncert_neg, lin_reg_age_uncert_neg)) %>%
          select(-age, -uncert_pos, -uncert_neg)
      } else {y$linReg <- F}
    }

    if(y$linInterp){
      setwd(file.path(y$working_directory, file_name, '/linInterp'))
      if(file.exists('linInt_chronology.csv')){
        lI_chrono <- read.csv('linInt_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
        names(lI_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(lI_chrono, SISAL_chronology_new, by = 'sample_id')  %>%
          mutate(lin_interp_age = if_else(sample_id %in% lI_chrono$sample_id, age, lin_interp_age),
                 lin_interp_age_uncert_pos = if_else(sample_id %in% lI_chrono$sample_id, uncert_pos, lin_interp_age_uncert_pos),
                 lin_interp_age_uncert_neg = if_else(sample_id %in% lI_chrono$sample_id, uncert_neg, lin_interp_age_uncert_neg)) %>%
          select(-age, -uncert_pos, -uncert_neg)
      } else { y$linInterp <- F}
    }

    if(y$Bchron){
      setwd(file.path(y$working_directory, file_name, '/Bchron'))
      if(file.exists('bchron_chronology.csv')){
        bchron_chrono <- read.csv('bchron_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
        names(bchron_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(bchron_chrono, SISAL_chronology_new, by = 'sample_id') %>%
          mutate(bchron_age = if_else(sample_id %in% bchron_chrono$sample_id, age, bchron_age),
                 bchron_age_uncert_pos = if_else(sample_id %in% bchron_chrono$sample_id, uncert_pos, bchron_age_uncert_pos),
                 bchron_age_uncert_neg = if_else(sample_id %in% bchron_chrono$sample_id, uncert_neg, bchron_age_uncert_neg)) %>%
          select(-age, -uncert_pos, -uncert_neg)
      } else {y$Bchron <- F}
    }

    if(y$StalAge){
      setwd(file.path(y$working_directory, file_name, '/StalAge'))
      if(file.exists('StalAge_chronology.csv')){
        stalage_chrono <- read.csv('StalAge_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
        names(stalage_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(stalage_chrono, SISAL_chronology_new, by = 'sample_id')  %>%
          mutate(StalAge_age = if_else(sample_id %in% stalage_chrono$sample_id, age, StalAge_age),
                 StalAge_age_uncert_pos = if_else(sample_id %in% stalage_chrono$sample_id, uncert_pos, StalAge_age_uncert_pos),
                 StalAge_age_uncert_neg = if_else(sample_id %in% stalage_chrono$sample_id, uncert_neg, StalAge_age_uncert_neg)) %>%
          select(-age, -uncert_pos, -uncert_neg)
      } else {y$StalAge <- F}
    }

    if(y$Bacon){
      setwd(file.path(y$working_directory, file_name, '/Bacon_runs'))
      if(file.exists('bacon_chronology.csv')){
        bacon_chrono <- read.csv('bacon_chronology.csv', header = T, colClasses = c('numeric', 'numeric', 'numeric', 'numeric')) #
        names(bacon_chrono) <- c('sample_id', 'age', 'uncert_pos', 'uncert_neg')
        SISAL_chronology_new <- full_join(bacon_chrono, SISAL_chronology_new, by = 'sample_id') %>%
          mutate(bacon_age = if_else(sample_id %in% bacon_chrono$sample_id, age, bacon_age),
                 bacon_age_uncert_pos = if_else(sample_id %in% bacon_chrono$sample_id, uncert_pos, bacon_age_uncert_pos),
                 bacon_age_uncert_neg = if_else(sample_id %in% bacon_chrono$sample_id, uncert_neg, bacon_age_uncert_neg)) %>%
          select(-age, -uncert_pos, -uncert_neg)
      } else {y$Bacon <- F}
    }

  }
  return(list(runFile, SISAL_chronology))
}

#' Evaluate SISAL chronology.
#'
#' Evaluates SISAL chronology following 4 checks:
#' 1: necessary: No reversals (AM fails if this check fails)
#' 2: necessary: Flexibility - AM has to folllow clear growth rate changes (AM fails if this checl fails)
#' 3: sufficient: Uncertainty increases between dates and at hiatuses
#' 4: sufficient: Increasing absolute uncertainty with age
#'
#' @param sisal_chrono SISAL chrnology that is meant to be evaluted
#' @param .dating SISAL v1b dating csv file
#' @param runFile runFile containing the entity ids of the executed age models
#'
#' @return 5 files: the evauated reversal check, the flexibility check, the 1st uncertainty check, the 2nd uncertainty check and a file containing the full evaluation of all 4 checks. Each passed check gives 1 point.
#' The more points a AM recieved the better the model.
#' @example evaluation <- eval(s,r)
eval <- function(sisal_chrono, runFile,.dating = SISAL.AM::dating){

  s_new_fil <- sisal_chrono %>% filter(entity_id %in% runFile$entity_id)

  dt <- .dating %>% filter(entity_id %in% s_new_fil$entity_id & date_type != 'Event; hiatus' & date_used == 'yes') %>% mutate_at(vars(depth_dating), as.numeric)

  d <- data.frame(entity_id = dt$entity_id,
                  sample_id = dt$dating_id,
                  depth_sample = dt$depth_dating,
                  lin_reg_age = dt$corr_age,
                  lin_reg_age_uncert_pos = dt$corr_age_uncert_pos,
                  lin_reg_age_uncert_neg = dt$corr_age_uncert_neg,
                  lin_interp_age = dt$corr_age,
                  lin_interp_age_uncert_pos = dt$corr_age_uncert_pos,
                  lin_interp_age_uncert_neg = dt$corr_age_uncert_neg,
                  bacon_age = dt$corr_age,
                  bacon_age_uncert_pos = dt$corr_age_uncert_pos,
                  bacon_age_uncert_neg = dt$corr_age_uncert_neg,
                  bchron_age = dt$corr_age,
                  bchron_age_uncert_pos = dt$corr_age_uncert_pos,
                  bchron_age_uncert_neg = dt$corr_age_uncert_neg,
                  StalAge_age = dt$corr_age,
                  StalAge_age_uncert_pos = dt$corr_age_uncert_pos,
                  StalAge_age_uncert_neg = dt$corr_age_uncert_neg) %>%
    mutate_at(vars(everything()), as.character) %>%
    mutate_at(vars(everything()), as.numeric)

  dt_count <- dt %>% group_by(entity_id) %>% count()

  d_fil <- d
  dt_new <- dt %>% as_tibble() %>% mutate(depth_new = (lead(depth_dating)-depth_dating)/2+depth_dating)
  d_new <- data.frame(sample_id = dt_new$dating_id, depth_new = dt_new$depth_new)

  SISAL_eval_rev <- s_new_fil %>% ungroup() %>% group_by(entity_id) %>%
    mutate(., diff_lr = lead(lin_reg_age) - lin_reg_age,
           diff_li = lead(lin_interp_age) - lin_interp_age,
           diff_bacon = lead(bacon_age) - bacon_age,
           diff_bchron = lead(bchron_age) - bchron_age,
           diff_stalage = lead(StalAge_age) - StalAge_age) %>%
    mutate(rev_lr = if_else( diff_lr< 0, FALSE, TRUE),
           rev_li = if_else( diff_li< 0, FALSE, TRUE),
           rev_bacon = if_else( diff_bacon< 0, FALSE, TRUE),
           rev_bchron = if_else( diff_bchron< 0, FALSE, TRUE),
           rev_StalAge = if_else( diff_stalage< 0, FALSE, TRUE)) %>%
    group_by(entity_id) %>%
    summarise(lR_reversal = if_else(all(is.na(rev_lr)),-1, if_else(any(!(rev_lr), na.rm = T), 0, 1)),
              lI_reversal = if_else(all(is.na(rev_li)),-1, if_else(any(!(rev_li), na.rm = T), 0, 1)),
              bacon_reversal = if_else(all(is.na(rev_bacon)),-1, if_else(any(!(rev_bacon), na.rm = T), 0, 1)),
              bchron_reversal = if_else(all(is.na(rev_bchron)),-1, if_else(any(!(rev_bchron), na.rm = T), 0, 1)),
              stalage_reversal = if_else(all(is.na(rev_StalAge)),-1, if_else(any(!(rev_StalAge), na.rm = T), 0, 1)))


  SISAL_eval_fit <- bind_rows(s_new_fil, d_fil) %>% ungroup() %>% group_by(entity_id) %>% dplyr::arrange(., depth_sample, .by_group =T) %>%
    mutate(rev_lr = if_else(abs(lin_reg_age-lead(lin_reg_age)) < (2*lin_reg_age_uncert_pos + 2*lead(lin_reg_age_uncert_neg)), TRUE, FALSE),
           rev_li = if_else(abs(lin_interp_age-lead(lin_interp_age)) < (2*lin_interp_age_uncert_pos + 2*lead(lin_interp_age_uncert_neg)), TRUE, FALSE),
           rev_bacon = if_else(abs(bacon_age-lead(bacon_age)) < (2*bacon_age_uncert_pos + 2*lead(bacon_age_uncert_neg)), TRUE, FALSE),
           rev_bchron = if_else(abs(bchron_age-lead(bchron_age)) < (2*bchron_age_uncert_pos + 2*lead(bchron_age_uncert_neg)), TRUE, FALSE),
           rev_StalAge = if_else(abs(StalAge_age-lead(StalAge_age)) < (2*StalAge_age_uncert_pos + 2*lead(StalAge_age_uncert_neg)), TRUE, FALSE)) %>%
    filter(sample_id %in% d_fil$sample_id) %>%
    summarise(lR_count = sum(!rev_lr, na.rm = T),
              lI_count = sum(!rev_li, na.rm = T),
              bacon_count = sum(!rev_bacon, na.rm = T),
              bchron_count = sum(!rev_bchron, na.rm = T),
              stalage_count = sum(!rev_StalAge, na.rm = T)) %>%
    left_join(., dt_count, by = "entity_id") %>% group_by(entity_id) %>% summarise(lR_fitness = if_else(lR_count <= 0.7*n, 1, 0),
                                                                                   lI_fitness = if_else(lI_count <= 0.7*n, 1, 0),
                                                                                   bacon_fitness = if_else(bacon_count <= 0.7*n, 1, 0),
                                                                                   bchron_fitness = if_else(bchron_count <= 0.7*n, 1, 0),
                                                                                   stalage_fitness = if_else(stalage_count <= 0.7*n, 1, 0))



  SISAL_IQR_check <- bind_rows(s_new_fil, d_fil) %>% ungroup() %>% group_by(entity_id) %>% dplyr::arrange(., depth_sample, .by_group =T) %>% left_join(., d_new, by = "sample_id") %>%
    mutate(depth_new_d = sapply(depth_new, function(x) depth_sample[order(abs(x - depth_sample))][1])) %>%
    mutate(depth_new_d = if_else(!is.na(depth_new), depth_new_d,  NA_real_)) %>% nest() %>%
    mutate(filter = purrr::map(.$data, function(x) filter(x, !is.na(depth_new_d))$depth_new_d)) %>%
    mutate(proxies = purrr::map2(.$data, .$filter, function(x, y) filter(x, depth_sample %in% y & is.na(depth_new) & is.na(depth_new_d)) %>% mutate(
      lR_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 ,
      lI_iqr =  (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2,
      bacon_iqr =  (bacon_age_uncert_pos+bacon_age_uncert_neg)/2,
      bchron_iqr =  (bchron_age_uncert_pos+bchron_age_uncert_neg)/2,
      StalAge_iqr =  (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2) %>% select(sample_id, depth_sample, lR_iqr, lI_iqr, bacon_iqr, bchron_iqr, StalAge_iqr))) %>%
    mutate(filter2 = purrr::map(.$data, function(x) filter(x, !is.na(depth_new))$sample_id)) %>%
    mutate(dating =  purrr::map2(.$data, .$filter2, function(x, y) filter(x, sample_id %in% y) %>% mutate(
      date_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2) %>% select(sample_id, depth_new_d, date_iqr))) %>%
    mutate(final = purrr::map2(.$proxies, .$dating, function(x,y) left_join(y %>% rename(depth_join = depth_new_d),x %>% rename(depth_join = depth_sample), by= 'depth_join'))) %>%
    unnest(final) %>%
    mutate(lR = if_else((date_iqr + lead(date_iqr))/2<lR_iqr, 1,0),
           lI = if_else((date_iqr + lead(date_iqr))/2<lI_iqr , 1,0),
           bacon = if_else((date_iqr + lead(date_iqr))/2< bacon_iqr , 1,0),
           bchron = if_else((date_iqr + lead(date_iqr))/2<bchron_iqr , 1,0),
           stalage = if_else((date_iqr + lead(date_iqr))/2<StalAge_iqr, 1,0)) %>% group_by(entity_id) %>%
    summarise(general = n(),
              lR_iqr = if_else(sum(lR, na.rm = T)>= 0.6*general, 1,0),
              lI_iqr = if_else(sum(lI, na.rm = T)>= 0.6*general, 1, 0),
              bacon_iqr= if_else(sum(bacon, na.rm = T)>= 0.6*general,1,0),
              bchron_iqr = if_else(sum(bchron, na.rm = T)>= 0.6*general,1,0),
              StalAge_iqr = if_else(sum(stalage, na.rm = T)>= 0.6*general,1,0)) %>%
    select(entity_id, lR_iqr, lI_iqr, copRa_iqr, bacon_iqr, bchron_iqr, StalAge_iqr)




  SISAL_IQR <- s_new_fil %>% ungroup() %>% group_by(entity_id) %>% arrange(., depth_sample, .by_group =T) %>%
    mutate(lR_iqr = (lin_reg_age_uncert_pos+lin_reg_age_uncert_neg)/2 ,
           lI_iqr =  (lin_interp_age_uncert_pos+lin_interp_age_uncert_neg)/2,
           bacon_iqr =  (bacon_age_uncert_pos+bacon_age_uncert_neg)/2,
           bchron_iqr =  (bchron_age_uncert_pos+bchron_age_uncert_neg)/2,
           StalAge_iqr =  (StalAge_age_uncert_pos+StalAge_age_uncert_neg)/2)

  SISAL_young <- SISAL_IQR %>% group_by(entity_id) %>% mutate(general = n()/3) %>%
    slice(.,1:general) %>% group_by(entity_id) %>% summarise(lR_young = mean(lR_iqr, na.rm = T),
                                                             lI_young = mean(lI_iqr, na.rm = T),
                                                             bacon_young = mean(bacon_iqr, na.rm = T),
                                                             bchron_young = mean(bchron_iqr, na.rm = T),
                                                             stalage_young = mean(StalAge_iqr, na.rm = T))

  SISAL_old <- SISAL_IQR %>% group_by(entity_id) %>% mutate(general = n()/3) %>%
    slice(.,general*2:n()) %>% group_by(entity_id) %>% summarise(lR_old = mean(lR_iqr, na.rm = T),
                                                                 lI_old = mean(lI_iqr, na.rm = T),
                                                                 bacon_old = mean(bacon_iqr, na.rm = T),
                                                                 bchron_old = mean(bchron_iqr, na.rm = T),
                                                                 stalage_old = mean(StalAge_iqr, na.rm = T))

  SISAL_IQR_check2 <- left_join(SISAL_young, SISAL_old, by = "entity_id") %>% mutate(lR_iqr = if_else(lR_young < lR_old , 1, 0),
                                                                                     lI_iqr = if_else(lI_young < lI_old, 1,0),
                                                                                     bacon_iqr = if_else(bacon_young < bacon_old,1,0),
                                                                                     bchron_iqr = if_else(bchron_young < bchron_old,1,0),
                                                                                     StalAge_iqr = if_else(stalage_young < stalage_old, 1,0)) %>%
    mutate(lR_iqr_2 = if_else(is.na(lR_young) & is.na(lR_old), -1, lR_iqr),
           lI_iqr_2 = if_else(is.na(lI_young) & is.na(lI_old), -1, lI_iqr),
           bacon_iqr_2 = if_else(is.na(bacon_young) & is.na(bacon_old), -1, bacon_iqr),
           bchron_iqr_2 = if_else(is.na(bchron_young) & is.na(bchron_old),-1,bchron_iqr),
           StalAge_iqr_2 = if_else(is.na(stalage_young) & is.na(stalage_old),-1,StalAge_iqr)) %>%
    select(entity_id, lR_iqr_2, lI_iqr_2, bacon_iqr_2, bchron_iqr_2, StalAge_iqr_2)



  SISAL_eval <- left_join(SISAL_eval_rev, SISAL_eval_fit, by = "entity_id") %>% left_join(., SISAL_IQR_check, by = "entity_id") %>%
    left_join(., SISAL_IQR_check2, by = "entity_id") %>%
    mutate(lR = if_else(lR_reversal == -1, -1, if_else(lR_reversal == 0, 0, 1 + lR_fitness + lR_iqr +lR_iqr_2)),
           lI = if_else(lI_reversal == -1,-1,if_else(lI_reversal == 0, 0, 1 + lI_fitness + lI_iqr + lI_iqr_2)),
           bacon = if_else(bacon_reversal == -1, -1, if_else(bacon_reversal == 0, 0, 1 + bacon_fitness + bacon_iqr + bacon_iqr_2)),
           bchron = if_else(bchron_reversal == -1, -1, if_else(bchron_reversal == 0, 0, 1+ bchron_fitness + bchron_iqr + bchron_iqr_2)),
           stalage = if_else(stalage_reversal == -1,-1, if_else(stalage_reversal == 0, 0, 1 + stalage_fitness + StalAge_iqr + StalAge_iqr_2))) %>%
    mutate(success = if_else(lR>0 | lI >0 | bacon >0 | bchron>0 | stalage >0, 'success', '-1')) %>%
    select(entity_id, lR, lI,copRa, bacon, bchron, stalage, success)

  return(list(SISAL_eval_rev, SISAL_eval_fit, SISAL_IQR_check, SISAL_IQR_check2, SISAL_eval))
}


#' Heatmap plot of the evaluated SISAL chronology.
#'
#' Plots the evaluated chrnology with ggplot.
#'
#' @param sisal_eval Evaluation for each entity and age model
#' @return Plot.
#'
#' @example
#' m <- merge_sisal_chrono(runFile)
#' r <- as.data.frame(m[1])
#' s <- as.data.frame(m[2])
#' evaluation <- eval(s,r)
#' plot <- plot_sisal_eval(as.data.frame(evaluation[5]))
#' plot
plot_sisal_eval <- function(sisal_eval){
  names(sisal_eval) <- c('entity_id', 'lin. Reg.', 'lin. Interp.','Bacon', 'Bchron', 'StalAge')
  data <- as.data.frame(gather(data = sisal_eval, key = Class, value = Fit, c(-1)))
  eval.heatmap <- ggplot(data = data, mapping = aes(x = Class, y = factor(entity_id), fill = as.character(Fit))) +
    geom_tile()  +
    facet_grid( ~ Class,switch = 'x',scales = 'free_x', space = 'free_x') +
    theme(plot.title = element_text(hjust = 0.5),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          strip.placement = 'outside',
          strip.background = element_rect(fill = 'white',color = 'black')) +
    xlab(label = 'Age Model') +
    ylab(label = 'Entity ID') +
    ggtitle(label = 'Filter for Reversals') +
    scale_fill_manual(guide = guide_legend(title = 'Fit yes-no'),values = c('0'='indianred', '1'='skyblue', '-1' = 'black'))

  return(eval.heatmap)
}
