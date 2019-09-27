#' Age model function for linear regression.
#'
#' @param working_directory File path, where input files are saved.
#' @param file_name name of the entity_id.
#' @param N Number of iterations.
#' @return Saves MC ensembel, lin_reg_chronology and AM plot.
runLinReg <- function(working_directory, file_name, N = 2000){

  print('---------------- Read in data -------------')
  setwd(file.path(working_directory, file_name, '/linReg'))
  dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
  depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))
  id <- read.csv('id.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))

  setwd(file.path(working_directory, file_name))
  hiatus_tb <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)

  sample <- data.frame(sample_id = id$sample_id, depth_eval = id$depth_sample)

  print('------- MC Simulations----------')
  mc_runs <- mc_ensemble(linReg = T, age = dating_tb$corr_age, age_error = dating_tb$corr_age_uncert, N = 2000,  working_directory = working_directory, file_name = file_name)

  print('--------lin Reg ------------')
  N <- dim(mc_runs)[1]
  lr <- mc_linReg(N, hiatus_tb$depth_sample, dating_tb$depth_dating, mc_runs, depth_sample$depth_sample)

  lr <- merge(sample, lr, by = 'depth_eval', all.x = T, all.y = T)
  lr <- select(lr, sample_id, depth_eval, everything())

  print('---------------save data Lin Reg --------------')
  setwd(file.path(working_directory, file_name, '/linReg'))
  write.table(as.matrix(mc_runs),'dating_mc_linReg_ensemble.txt', col.names = F, row.names = F)
  write.table(as.matrix(lr),'mc_linReg_ensemble.txt', col.names = F, row.names = F)

  print('--------------median and quantiles--------------')
  stats <- get_median_quantiles(lr[,3:N+2], q1 = 0.05, q2 = 0.95)
  age_median <- stats[,1]
  age_sd_low <- stats[,2]
  age_sd_high <- stats[,3]
  lin_reg <- cbind(lr$sample_id,age_median, age_sd_high-age_median, age_median - age_sd_low)
  colnames(lin_reg) <- c('sample_id', 'lin_reg_age', 'lin_reg_age_uncert_pos', 'lin_reg_age_uncert_neg')
  write.csv(lin_reg, 'linReg_chronology.csv',  row.names = F, sep =',')

  pdf('final_age_model.pdf', 6,4)
  matplot(x = age_median, y= lr$depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1,
          ylim = c(max(lr$depth_eval),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=age_sd_high,y=lr$depth_eval, lty = 2, col = 'red')
  lines(x=age_sd_low,y = lr$depth_eval, lty = 2, col = 'red')
  points(x = dating_tb$corr_age, y=dating_tb$depth_dating, lty = 2, col = 'orange', pch = 4)
  arrows(dating_tb$corr_age-dating_tb$corr_age_uncert, dating_tb$depth_dating, dating_tb$corr_age+dating_tb$corr_age_uncert, dating_tb$depth_dating, length=0.05, angle=90, code=3, col = 'orange')
  if (!plyr::empty(data.frame(hiatus_tb))) {
    abline(h = hiatus_tb$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()
}

#' Age model function for linear interpolation.
#'
#' @param working_directory File path, where input files are saved.
#' @param file_name name of the entity_id.
#' @param N Number of iterations.
#' @return Saves MC ensembel, lin_interp_chronology and AM plot.
runLinInterp <- function(working_directory, file_name, N = 2000){

  print('---------------- Read in data -------------')
  setwd(file.path(working_directory, file_name, '/linInterp'))
  dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
  depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))

  setwd(file.path(working_directory, file_name))
  hiatus_tb <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)

  sample <- data.frame(sample_id = depth_sample$sample_id, depth_eval = depth_sample$depth_sample)

  if (!plyr::empty(data.frame(hiatus_tb))) {
    dating_tb <- add_hiatus(dating_tb, hiatus_tb, linInterp =T)
    write.csv(dating_tb, 'hiatus_dates_interp.csv', row.names = F, sep = ',')
  }

  dating_tb <- scan_lin_interp(dating_tb)
  dating_tb <- scan_fine_lin_interp(dating_tb)

  setwd(file.path(working_directory, file_name, '/linInterp'))
  write.csv(dating_tb, 'new_dating_tb.csv', row.names = F)

  print('------- MC Simulations----------')
  mc_runs <- mc_ensemble(linInterp = T, age = dating_tb$corr_age, age_error = dating_tb$corr_age_uncert, N = 2000, working_directory = working_directory, file_name = file_name)

  print('------------- lin Interp -----------')
  N <- dim(mc_runs)[1]
  li <- mc_linInt(N, hiatus_tb$depth_sample, dating_tb$depth_dating, mc_runs, depth_sample$depth_sample)
  li <- merge(sample, li, by = 'depth_eval', all.x = T, all.y = T)
  li <- select(li, sample_id, depth_eval, everything())

  print('---------------save data Lin Interp --------------')
  setwd(file.path(working_directory, file_name, '/linInterp'))
  write.table(as.matrix(mc_runs[[2]]),'dating_mc_linInt_ensemble.txt', col.names = F, row.names = F)
  write.table(as.matrix(li),'mc_linInt_ensemble.txt', col.names = F, row.names = F)

  print('--------------median and quantiles--------------')
  stats <- get_median_quantiles(li[,3:N+2], q1 = 0.05, q2 = 0.95)
  age_median <- stats[,1]
  age_sd_low <- stats[,2]
  age_sd_high <- stats[,3]
  lin_interp <- cbind(li$sample_id, age_median, age_sd_high-age_median, age_median - age_sd_low)
  colnames(lin_interp) <- c('sample_id', 'lin_interp_age', 'lin_interp_age_uncert_pos', 'lin_interp_age_uncert_neg')
  write.csv(lin_interp, 'linInt_chronology.csv', row.names = F, sep = ',')

  pdf('final_age_model.pdf', 6,4)
  matplot(x = age_median, y= li$depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1,
          ylim = c(max(li$depth_eval),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=age_sd_high,y=li$depth_eval, lty = 2, col = 'red')
  lines(x=age_sd_low,y = li$depth_eval, lty = 2, col = 'red')
  points(x = dating_tb$corr_age, y=dating_tb$depth_dating, lty = 2, col = 'orange', pch = 4)
  arrows(dating_tb$corr_age-dating_tb$corr_age_uncert, dating_tb$depth_dating, dating_tb$corr_age+dating_tb$corr_age_uncert, dating_tb$depth_dating, length=0.05, angle=90, code=3, col = 'orange')
  if (!plyr::empty(data.frame(hiatus_tb))) {
    abline(h = hiatus_tb$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()
}

#' Age model function for Bacon.
#'
#' @param working_directory File path, where input files are saved.
#' @param file_name name of the entity_id.
#' @param N Number of iterations.
#' @return Saves MC ensembel, bacon_chronology and AM plot.
runBacon <- function(working_directory, file_name, postbomb = 0, cc = 0 ){
  setwd(file.path(working_directory, file_name, '/Bacon_runs/', file_name))
  depth_eval <- matrix(read.table(paste(file_name,'_depths.txt', sep = ''), col.names = ''))[[1]]
  sample_id <-  read.csv('sample_id.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric'))
  hiatus_tb <- read.csv('hiatus_bacon.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  core <- read.csv(paste(file_name,'.csv', sep = ''), header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))


  accMean<- sapply(c(1, 2, 5), function(x) x * 10^(-1:2))
  ballpacc <- lm(core[,2] * 1.1 ~ core[, 4])$coefficients[2]
  ballpacc <- abs(accMean - ballpacc)
  ballpacc <- ballpacc[ballpacc > 0]
  accMean <- accMean[order(ballpacc)[1]]

  k <- seq(floor(min(depth_eval, na.rm = T)), ceiling(max(depth_eval, na.rm = T)), by = 5)
  if (k < 10) {
    thickness <- pretty(5 * (k/10), 10)
    thickness <- min(thickness[thickness > 0])
  } else if (k > 20) {
    thickness <- max(pretty(5 * (k/20)))
  }

  j <- 2000
  tho <- c()


  setwd(file.path(working_directory, file_name))
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  print('#------------ run bacon ---------------#')
  if (dim(hiatus_tb)[1] == 0) {
    tryCatch({Bacon(core = file_name, depths.file = TRUE, thick = thickness, acc.mean = accMean, postbomb = postbomb, cc = cc, suggest = F, ask = F, ssize = j, th0 = tho)},
             error = function(e){write.table(x=paste("ERROR in Bacon:",conditionMessage(e)),file = 'bacon_error.txt')})
  } else {
    tryCatch({Bacon( core = file_name, depths.file = TRUE, thick = thickness, acc.mean = accMean ,postbomb = postbomb, hiatus.depths = hiatus_tb$depth_sample_bacon, cc =cc, suggest = F, ask = F, ssize = j, th0 = tho)},
             error = function(e){write.table(x=paste("ERROR in Bacon:",conditionMessage(e)),file = 'bacon_error.txt')})
  }

  print('--------------save data -----------------')
  bacon_mcmc <- sapply(depth_eval, Bacon.Age.d)
  bacon_age <- get_bacon_median_quantile(depth_eval, hiatus_tb, bacon_mcmc)
  bacon_mcmc <- rbind(depth_eval, bacon_mcmc)
  bacon_mcmc <- t(bacon_mcmc)
  bacon_mcmc <- cbind(sample_id, bacon_mcmc)

  h <- cbind(hiatus_tb, matrix(NA,nrow = dim(hiatus_tb)[1], ncol = dim(bacon_mcmc)[2]-2))
  names(h) <- names(bacon_mcmc)

  bacon_mcmc <- rbind(bacon_mcmc, h)
  bacon_mcmc <- bacon_mcmc[order(bacon_mcmc[,2]),]

  sample_id <- bacon_mcmc[,1]

  setwd(file.path(working_directory, file_name,'/Bacon_runs'))
  write.table(bacon_mcmc, 'mc_bacon_ensemble.txt', col.names = F, row.names = F)
  write.csv(cbind(sample_id, bacon_age[,2:4]),'bacon_chronology.csv' ,col.names = T, row.names = F)

  pdf('final_age_model.pdf', 6,4)
  matplot(x = bacon_age[,2], y= bacon_age[,1]*10, col = 'black', lty = 1, type = 'l', lwd = 1,
          ylim = c(max(depth_eval*10),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=bacon_age[,3]+bacon_age[,2],y=bacon_age[,1]*10, lty = 2, col = 'red')
  lines(x=bacon_age[,2]-bacon_age[,4],y = bacon_age[,1]*10, lty = 2, col = 'red')
  points(x = core$corr_age, y=core$depth_dating_new*10, lty = 2, col = 'orange', pch = 4)
  arrows(core$corr_age-core$corr_age_uncert, core$depth_dating_new*10, core$corr_age+core$corr_age_uncert, core$depth_dating_new*10, length=0.05, angle=90, code=3, col = 'orange')
  if (!plyr::empty(data.frame(hiatus_tb))) {
    abline(h = hiatus_tb$depth_sample_bacon*10, col = 'grey', lty = 2)
  }
  dev.off()

}

#' Age model function for Bchron.
#'
#' @param working_directory File path, where input files are saved.
#' @param file_name name of the entity_id.
#' @param N Number of iterations.
#' @return Saves MC ensembel, bchron_chronology and AM plot.
runBchron <- function(working_directory,file_name){ # sample ids missing
  setwd(file.path(working_directory, file_name, '/Bchron'))
  dating_tb <- read.csv('ages.csv', header = T, stringsAsFactors = F)
  depth_sample <- read.csv('depths.csv', header = T, stringsAsFactors = F, colClasses = c('numeric', 'numeric'))

  depth_eval <- depth_sample$depth_sample # in cm !!!!!!

  setwd(file.path(working_directory, file_name))
  hiatus_tb <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)

  #dating_tb <- dating_tb %>% add_row(.,dating_id =1,corr_age =3000,corr_age_uncert= 60,depth_dating_new= 0,thickness_new =0.1,calib_curve_new='normal', .before = 1)
  if (!plyr::empty(data.frame(hiatus_tb))) {
    dating_tb <- add_hiatus(dating_tb, hiatus_tb, bchron = T)
    write.csv(dating_tb, 'hiatus_dates_bchron.csv', row.names = F, sep = ',')
  }

  run <- Bchronology(ages=dating_tb$corr_age,
                     ageSds=dating_tb$corr_age_uncert,
                     positions = dating_tb$depth_dating_new,
                     positionThicknesses = dating_tb$thickness_new,
                     calCurves = dating_tb$calib_curve_new,
                     ids = dating_tb$dating_id,
                     predictPositions = depth_eval,
                     jitterPositions = T)

  mcmc <- run$thetaPredict
  bchron_age <- apply(mcmc,2,median)
  bchron_quantile <- apply(mcmc, 2, function(x){quantile(x, probs = c(0.05,0.95), na.rm = T)})

  d <- cbind(depth_eval, bchron_age,
             bchron_age_uncert_pos = bchron_quantile[2,]-bchron_age,
             bchron_age_uncert_neg = bchron_age - bchron_quantile[1,])
  d <- data.frame(d)

  if (!plyr::empty(data.frame(hiatus_tb))) {
    hiatus_new <- hiatus_tb %>% mutate(depth_sample = depth_sample/10)
    d <- d %>% mutate(bchron_age = if_else(depth_eval %in% hiatus_new$depth_sample, NA_real_, bchron_age),
                      bchron_age_uncert_pos = if_else(depth_eval %in% hiatus_new$depth_sample, NA_real_, bchron_age_uncert_pos),
                      bchron_age_uncert_neg = if_else(depth_eval %in% hiatus_new$depth_sample, NA_real_, bchron_age_uncert_neg))
  }

  setwd(file.path(working_directory, file_name, '/Bchron'))
  write.table(mcmc, 'bchron_ensemble.txt', col.names = F, row.names = F)
  b_chrono <- data.frame(sample_id = depth_sample$sample_id,
                         bchron_age = d[,2],
                         bchron_age_uncert_pos = d[,3],
                         bchron_age_uncert_neg = d[,4])
  write.csv(b_chrono,'bchron_chronology.csv' ,col.names = T, row.names = F)

  x_min <- min(b_chrono$bchron_age, na.rm = T)
  x_max <- max(b_chrono$bchron_age, na.rm = T)

  pdf('final_age_model.pdf', 6,4)
  matplot(x = b_chrono$bchron_age, y= d$depth_eval*10, col = 'black', lty = 1, type = 'l', lwd = 1,
          ylim = c(max(d$depth_eval*10),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]', xlim = c(3000, x_max))
  lines(x=b_chrono$bchron_age - b_chrono$bchron_age_uncert_neg,y=d$depth_eval*10, lty = 2, col = 'red')
  lines(x=b_chrono$bchron_age_uncert_pos + b_chrono$bchron_age,y = d$depth_eval*10, lty = 2, col = 'red')
  points(x = dating_tb$corr_age, y=dating_tb$depth_dating_new*10, lty = 2, col = 'orange', pch = 4)
  arrows(dating_tb$corr_age-dating_tb$corr_age_uncert, dating_tb$depth_dating_new*10, dating_tb$corr_age+dating_tb$corr_age_uncert, dating_tb$depth_dating_new*10, length=0.05, angle=90, code=3, col = 'orange')
  if (!plyr::empty(data.frame(hiatus_tb))) {
    abline(h = hiatus_tb$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()
}

#' Age model function for StalAge
#'
#' @param working_directory File path, where input files are saved.
#' @param file_name name of the entity_id.
#' @param N Number of iterations.
#' @return Saves MC ensembel, StalAge_chronology and AM plot.
runStalAge <- function(working_directory, file_name){

  #Einlesen der Daten
  print('#-----------READ IN DATA-------------#')
  setwd(file.path(working_directory, file_name))
  hiatus_tb <- read.csv('hiatus.csv', header = T, stringsAsFactors = F,  colClasses = c('numeric', 'numeric'))
  unknown_age <- read.csv('not_used_dates.csv', header = T)

  setwd(file.path(working_directory, file_name, '/StalAge'))
  Daten_orig <- read.csv("ages.csv", header = T)
  Daten_orig <- Daten_orig %>% mutate(corr_age_uncert = 2*corr_age_uncert) %>% select(dating_id, corr_age, corr_age_uncert, depth_dating)

  if (!plyr::empty(data.frame(hiatus_tb))) {
    Daten_orig <- add_hiatus(Daten_orig, hiatus_tb, stalage=T)
    write.csv(Daten_orig, 'hiatus_dates_StalAge.csv', row.names = F, sep = ',')
  }

  names(Daten_orig) <- c('dating_id', 'age','error', 'depth')

  depths <- read.csv("depths.csv", header = T)
  names(depths) <- c('sample_id', 'dft')



  #StalAge Modellierung
  print('#----------SCAN DATA FOR REVERSALS ----------------#')
  Daten <- scan(Daten_orig$depth, Daten_orig$age, Daten_orig$error)
  Daten <- scan_fine(Daten)

  print('#----------FIT DATA ----------------#')
  fit <- age_model(Daten, depths$dft, Daten_orig)
  names(fit) <- c('depth_sample', 'StalAge_age','StalAge_age_uncert_pos', 'StalAge_age_uncert_neg')

  print('#----------WRITE DATA TO CSV FILES----------------#')
  sample_id <- depths$sample_id
  depth_eval <- depths$dft
  StalAge_age <- fit$StalAge_age
  StalAge_age_uncert_pos <- fit$StalAge_age_uncert_pos - StalAge_age
  StalAge_age_uncert_neg <- StalAge_age - fit$StalAge_age_uncert_neg

  d <- cbind(sample_id, depth_eval, StalAge_age, StalAge_age_uncert_pos, StalAge_age_uncert_neg)
  d <- data.frame(d)

  if(!plyr::empty(data.frame(hiatus_tb))){
    d <- d %>% rowwise() %>% mutate(StalAge_age = if_else(depth_eval %in% hiatus_tb$depth_sample, NA_real_, StalAge_age),
                                    StalAge_age_uncert_pos = if_else(depth_eval %in% hiatus_tb$depth_sample, NA_real_, StalAge_age_uncert_pos),
                                    StalAge_age_uncert_neg = if_else(depth_eval %in% hiatus_tb$depth_sample, NA_real_, StalAge_age_uncert_neg))
  }


  #chronology <- cbind(sample_id,StalAge_age, StalAge_age_uncert_pos - StalAge_age, StalAge_age - StalAge_age_uncert_neg)
  #colnames(chronology) <- c('sample_id','StalAge_age', 'StalAge_age_uncert_pos', 'StalAge_age_uncert_neg')

  setwd(file.path(working_directory, file_name, '/StalAge'))
  write.csv(d, file="StalAge_results.csv",
            col.names = T, row.names = F)
  write.table(d[,c(1,3,4,5)], file="StalAge_chronology.csv", col.names = T,sep = ',' ,row.names = F)
  write.csv(Daten, file = 'StalAge_date_used.csv', col.names = T, row.names = F)

  pdf('final_age_model.pdf', 6,4)
  matplot(x = StalAge_age, y= depth_eval, col = 'black', lty = 1, type = 'l', lwd = 1,
          ylim = c(max(depth_eval),0), xlab = 'Age [yrs BP]', ylab = 'Depth from top [mm]')
  lines(x=StalAge_age_uncert_pos+StalAge_age,y=depth_eval, lty = 2, col = 'red')
  lines(x=StalAge_age - StalAge_age_uncert_neg,y = depth_eval, lty = 2, col = 'red')
  points(x = Daten_orig$age, y=Daten_orig$depth, lty = 2, col = 'orange', pch = 4)
  arrows(Daten_orig$age-Daten_orig$error, Daten_orig$depth, Daten_orig$age+Daten_orig$error, Daten_orig$depth, length=0.05, angle=90, code=3, col = 'orange')
  if (!plyr::empty(data.frame(hiatus_tb))) {
    abline(h = hiatus_tb$depth_sample, col = 'grey', lty = 2)
  }
  dev.off()

}

#' Run age models for chosen entities.
#'
#' @param prefix prefix to SISAL files
#' @param run File specifying which entitis and age models for each entity to run, where to save the files and where the SISAL files are saved.
#' @return Runs the age models.
#'
#' @example
#' setwd('~')
#' m <- filter_SISAL(1) %>% filter(entity_id < 10)
#' AM_SISAL('',m)
#'
AM_SISAL <- function(prefix, run, file_path = file.path(getwd(),"SISAL_Age_Model"), b=T, bc=T, stalage=T, lI=T, lR=T, default = T) {

  if(file_path == file.path(getwd(),"SISAL_Age_Model")){
    dir.create(file.path(getwd(),"SISAL_Age_Model"))
  }

    runFile <- tibble(entity_id = run$entity_id, Bacon = rep(b, length(run$entity_id)), Bchron = rep(bc, length(run$entity_id)),
                      StalAge = rep(stalage, length(run$entity_id)), linInterp = rep(lI, length(run$entity_id)), linReg = rep(lR, length(run$entity_id)),
                      working_directory = rep(file_path, length(run$entity_id)))

  j <- 1

  sapply(1: dim(runFile)[1],
         function(j, x){
           y <- x[j,]
           run_SISAL_chrono(runFile, entid = y$entity_id, working_directory = y$working_directory, bacon = y$Bacon, bchron = y$Bchron, stalage = y$StalAge, linInterp = y$linInterp, linReg = y$linReg, j. = j)
           graphics.off()
           gc()
         },
         x = runFile
  )
}

#' Run AM for single entity_id.
#'
#' @param entid entity_id
#' @param data List of site_tb, dating_tb and ample_tb.
#' @param working_directory File path to where input files are saved.
#' @inheritParams write_files
#' @return Executes the specified age models, marks in a sepearte file the successfully and failed age models, writes an error file and plots the summary.
run_SISAL_chrono <- function(runFile, entid, working_directory, bacon, bchron, stalage, linInterp, linReg,j.) {
  file_name <- write_files(entid, bacon, bchron, stalage, linInterp, linReg, working_directory)

  err <- NULL
  tryCatch({
    i = 7
    runLinReg(working_directory, file_name)
    print('LinReg done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in linReg:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in linReg:",conditionMessage(e), "\n"))
    runFile[j.,i] <- F
    linReg <<- F})

  tryCatch({
    i = 5
    runLinInterp(working_directory, file_name)
    print('LinInterp done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in linInterp:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in linInterp:",conditionMessage(e), "\n"))
    runFile[j.,i] <- F
    linInterp <<- F})

  tryCatch({
    i = 2
    runBacon(working_directory, file_name)
    print('Bacon done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in Bacon:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in Bacon:",conditionMessage(e), "\n"))
    runFile[j.,i] <- F
    bacon <<- F})

  tryCatch({
    i = 3
    runBchron(working_directory, file_name)
    print('Bchron done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in Bchron:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in Bchron new:",conditionMessage(e), "\n"))
    runFile[j.,i] <- F
    bchron <<- F})

  tryCatch({
    i = 4
    runStalAge(working_directory, file_name)
    print('StalAge done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in StalAge:",conditionMessage(e), "\n"))
    #error = function(e){print(paste("ERROR in StalAge:",conditionMessage(e), "\n"))
    runFile[j.,i] <- F
    stalage <<- F})

  setwd(file.path(working_directory, file_name))
  print('wd changed')

  tryCatch({
    plotAMC(working_directory, file_name, b = bacon, bc = bchron, stalage = stalage, linreg = linReg, lininterp = linInterp)
    print('plot done')
    graphics.off()},
    error = function(e){err <<- append(err, paste("ERROR in plot:",conditionMessage(e), "\n"))})
  #error = function(e){print(paste("ERROR in plot:",conditionMessage(e), "\n"))})

  write.table(err, 'errors.txt', row.names = F, col.names = F)

  setwd(file.path(working_directory))
  write.csv(runFile, 'runFile.csv', row.names = F, col.names = T)

  print(paste(file_name, 'done'))
}

