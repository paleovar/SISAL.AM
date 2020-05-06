#' Plot summary plot.
#'
#' @param working_directory File path from where data is loaded and plot is to be saved.
#' @param file_name Name of the folder for which the plot is to be generated.
#' @param b Logical. TRUE if Bacon is to be plotted.
#' @param bc Logical. TRUE if Bchron is to be plotted.
#' @param stalage Logical. TRUE if Stalage is to be plotted.
#' @param linreg Logical. TRUE if lin. reg. is to be plotted.
#' @param lininterp Logical. TRUE if lin. interp. is to be plotted.
#' @return Plots all successfully age models and confidence intervals. Marks modified and artificial dates.
plotAMC <- function(working_directory, file_name, b, bc, stalage , linreg, lininterp){
  
  setwd(file.path(working_directory,file_name))
  proxy <- read.csv('proxy_data.csv', header = T)
  hiatus_tb <- read.csv('hiatus.csv', header = T)
  original_chrono <- read.csv('original_chronology.csv', header = T)
  unknown_age <- read.csv('not_used_dates.csv', header = T)
  used_age <- read.csv('used_dates.csv', header = T)
  #used_age <- read.csv('hiatus_dates.csv', header = T)
  
  depth_sample <- proxy$depth_sample
  original_age <- original_chrono$interp_age
  
  if(all(is.na(depth_sample))){stop('No depth_sample!!!')}
  
  date_type <- unique(used_age$date_type)
  age_model <- unique(original_chrono$age_model_type)
  
  
  if(stalage){ # in mm
    setwd(file.path(working_directory, file_name, '/StalAge'))
    StalAge <- read.csv('StalAge_chronology.csv', header = T)
    StalAge_age <- StalAge$StalAge_age
    StalAge_iqr <- StalAge$StalAge_age_uncert_pos + StalAge$StalAge_age_uncert_neg
    col_StalAge <- 'forestgreen'
  }
  #print('stalage')
  
  if(b){ # in cm
    setwd(file.path(working_directory, file_name, '/Bacon_runs'))
    bacon <- read.csv('bacon_chronology.csv', header = T)
    bacon_age <- bacon$bacon_age
    bacon_iqr <- bacon$bacon_age_uncert_pos + bacon$bacon_age_uncert_neg
    col_bacon <- 'palevioletred4'
    #col_bacon <- 'yellow2'
  }
  #print('bacon')
  
  if(bc){ # in cm
    setwd(file.path(working_directory, file_name,'/Bchron'))
    bchron <- read.csv('bchron_chronology.csv', header = T)
    bchron_age <- bchron$bchron_age
    bchron_iqr <- bchron$bchron_age_uncert_pos + bchron$bchron_age_uncert_neg
    col_bchron <- 'sienna1'
  }
  #print('bchron')
  
  if(linreg){ # in mm
    setwd(file.path(working_directory, file_name, '/linReg'))
    linReg <- read.csv('linReg_chronology.csv', header = T)
    lin_reg_age <- linReg$lin_reg_age
    lin_reg_iqr <- linReg$lin_reg_age_uncert_pos + linReg$lin_reg_age_uncert_neg
    
    col_linReg <- 'springgreen2'
    
  }
  
  if(lininterp){
    setwd(file.path(working_directory, file_name, '/linInterp'))
    linInt <- read.csv('linInt_chronology.csv', header = T)
    lin_interp_age <- linInt$lin_interp_age
    lin_interp_iqr <- linInt$lin_interp_age_uncert_pos + linInt$lin_interp_age_uncert_neg
    
    #depth_sample_new <- sort(c(depth_sample, hiatus))
    
    #col_linReg <- 'springgreen2'
    col_linInt <- 'skyblue'
  }
  #print('mc')
  
  age <- cbind(depth_sample, original_age)
  iqr <- cbind(depth_sample)
  
  color <- 'black'
  color_iqr <- NULL
  
  if(linreg){
    age <- cbind(age, lin_reg_age)
    iqr <- cbind(iqr, lin_reg_iqr)
    
    color <- c(color, col_linReg)
    color_iqr <- c(color_iqr, col_linReg)
  }
  
  if(lininterp){
    age <- cbind(age, lin_interp_age)
    iqr <- cbind(iqr, lin_interp_iqr)
    
    color <- c(color, col_linInt)
    color_iqr <- c(color_iqr, col_linInt)
  }
  
  if(stalage){
    
    age <- cbind(age, StalAge_age)
    iqr <- cbind(iqr, StalAge_iqr)
    
    color <- c(color, col_StalAge)
    color_iqr <- c(color_iqr, col_StalAge)
    
  }
  
  if(b){
    age <- cbind(age, bacon_age)
    iqr <- cbind(iqr, bacon_iqr)
    
    color <- c(color, col_bacon)
    color_iqr <- c(color_iqr, col_bacon)
  }
  
  
  if(bc){
    age <- cbind(age, bchron_age)
    iqr <- cbind(iqr, bchron_iqr)
    
    color <- c(color, col_bchron)
    color_iqr <- c(color_iqr, col_bchron)
  }
  
  
  leg <- colnames(age)[2:length(colnames(age))]
  leg_iqr <- colnames(iqr)[2:length(colnames(iqr))]
  
  
  x_lim <- max(age[,2:ncol(age)], na.rm = T)/1000
  if(!plyr::empty(data.frame(unknown_age)) && min(unknown_age$corr_age)< min(age[,2:ncol(age)], na.rm = T)){x_lim_bottom <- min(unknown_age$corr_age)/1000} else {x_lim_bottom <- min(age[,2:ncol(age)], na.rm = T)/1000}
  y_lim <- max(age[,1])+3
  if(x_lim/10 > 1) {
    steps_x <- round(x_lim/10)
  } else {steps_x <- round(x_lim/10, digits = 2)}
  #steps_x <- round(x_lim/10)
  steps_y <- round(y_lim/5, digits = 0)
  
  setwd(file.path(working_directory, file_name))
  pdf(paste('am_iqr_',file_name,'.pdf', sep = ''),9.5, 12.75)
  #pdf(paste('iqr_',fname,'.pdf', sep = ''),9.5, 11)
  layout(matrix(c(1,2,3), nrow=3, ncol=1, byrow = TRUE), widths = lcm(20), heights = c(rep(lcm(11),2),5))
  #layout(matrix(c(1,2), nrow=2, ncol=1, byrow = TRUE), widths = lcm(22), heights = c(rep(lcm(16),1),5))
  par(mar = c(5,5,3,1), oma = c(1,2,2,1))
  plot_AM(age,cbind(unknown_age$corr_age, unknown_age$depth_dating),
          cbind(used_age$corr_age, used_age$depth_dating, used_age$corr_age_uncert_pos, used_age$corr_age_uncert_neg),
          hiatus_tb$depth_sample, date_type, x_lim, x_lim_bottom, y_lim, steps_x,
          steps_y, cex_legend = 0.2, color, file_name)
  
  plot_IQR(iqr, age, used_age, leg_iqr, cex_legend = 0.4, color_iqr, file_name)
  
  plot.new()
  legend("top",legend = c(paste(date_type,"age - used"),"Age - not used", paste('Used Age model:',age_model)), pch = c(4,25,1),col = c("orange", "black", 'black') ,bty='n', ncol = 1, title = "AM date types", cex = 1.5, pt.cex = 4, title.adj = 0.5)
  #legend("top",legend = c("Age - used"), pch = c(4),col = c("orange") ,bty='n', title = "AM date types", cex = 1.4, pt.cex = 4, title.adj = 0.5)
  #legend('bottom', legend= c(leg, 'Hiatus'),
  #       lwd = 3,col = c(color, 'grey'), lty = c(rep(1,length(age)-1),2),bty ='n', title = "Additional AM information", cex = 2, ncol=2)
  legend('bottom', legend= c(leg, 'Hiatus'),
         lwd = 3,col = c(color,'grey'), lty = c(1,1,1,1,2),bty ='n', title = "Additional AM information", cex = 1.5, ncol=2)
  
  dev.off()
  
}

#' Plots successfully executed AM.
#'
#' @param age Table containing dates, uncertainies and dating depths of used ages.
#' @param unknown Table containing dates, uncertainies and dating depths of not used ages.
#' @param hiatuses Table containing hiatus depths.
#' @param x_lim Upper limit of x axis.
#' @param x_lim_bottom Lower limit of x axis.
#' @param y_lim Upper limit of y axis.
#' @param color Color for each AM
#' @param steps_x Ticks steps of x axis
#' @param steps_y Ticks steps of y axis.
#' @param file_name Entity_id and entity_name.
plot_AM <- function(age, unknown, used, hiatuses, x_lim = 50,x_lim_bottom, y_lim  = 600, steps_x = 10,
                    steps_y = 100,cex_legend ,color ,file_name) {
  j <- ncol(age)
  
  age[,seq(2,j)] <- age[,seq(2,j)]/1000
  unknown[,1] <- unknown[,1]/1000 # to ka
  used[,c(1,3,4)] <- used[,c(1,3,4)]/1000 # to ka
  
  matplot(age[, seq(2,j)], age[,1], type = 'l', lty = 1, col = color, xlab = '', ylab = '', lwd = 2, axes = FALSE, ylim = c(y_lim,0), xlim = c(x_lim_bottom, x_lim))
  #points(used[,2], used[,1], pch = 4, col = "orange", cex = 2, lwd = 2)
  points(used[,1], used[,2], pch = 4, col = "black", cex = 2, lwd = 2)
  arrows(used[,1]-used[,4], used[,2], used[,1]+used[,3], used[,2], length=0.05, angle=90, code=3, col = 'black')
  
  if(length(hiatuses) != 0) {
    abline(h = hiatuses, col = 'grey', lty = 2, lwd = 1.5)
  }
  
  if (length(unknown) != 0) {
    points(unknown[,1],unknown[,2], pch = 25, col ="red", cex = 2)
  }
  
  axis(side = 1,  lwd = 2, cex.axis = 2.2, at = seq(x_lim_bottom, x_lim, by=steps_x), lwd.ticks = 2,labels = seq(x_lim_bottom,x_lim, by = steps_x), padj = 0.5)
  axis(side = 2,  lwd = 2, lwd.ticks = 2, cex.axis = 2.2, labels = seq(0,y_lim, by = steps_y), at= seq(0,y_lim, by = steps_y))#, padj = -0.5)
  
  mtext(side = 3, text = paste('Age-depth model:', file_name), line = 0, cex = 2)
  mtext(side = 2, text = 'Depth from top [mm]', line = 3, cex  =1.7)
}

#' Plots IQRs of successfully executed AM.
#'
#' @param iqr Table containing IQRs of successfully executed AM.
#' @param depths Sample depths.
#' @param color_iqr Color for each AM
#' @param used Table containing dating depths.
#' @param file_name Entity_id and entity_name.
plot_IQR <- function(iqr,depths,used, cex_legend, color_iqr, file_name) {
  j <- ncol(iqr)
  
  iqr <- iqr/1000
  #age[,seq(3,j+1)] <- age[,seq(3,j+1)]/1000
  
  y_lim <- max(age[,1]) + 10
  x_lim <- max(apply(iqr[,seq(2,j)], 2, function(x){max(x, na.rm = T)})) +1
  y_steps <- round(y_lim/5)
  x_steps <- round(x_lim/10)
  
  matplot(iqr[,seq(2,j)], age[,1], type = 'l', lty = 1, col = color, xlab = '', ylab = '', lwd = 2, axes = FALSE,ylim = c(y_lim,0) )
  #legend('topright', legend= legend, col = color, lwd = 2,lty = 1, bty ='n', title = "Information")
  abline(h = used$depth_dating, col = 'grey', lty = 2, lwd = 1.5)
  #points(x = used/1000, y = rep(0, length(used)), pch = 4, col = "orange", cex = 2, lwd = 2)
  
  axis(side = 1,  at = seq(0, x_lim, by = x_steps), labels = seq(0,x_lim, by = x_steps), lwd = 2, lwd.ticks = 2, cex.axis = 2.2, padj = 0.5)
  axis(side = 2,  lwd = 2, lwd.ticks = 2, cex.axis = 2.2, labels = seq(0, round(y_lim),by = y_steps), at= seq(0,round(y_lim), by = y_steps))#, padj = -0.5)
  #axis(side = 2,  lwd = 2, lwd.ticks = 2, cex.axis = 2, labels = seq(0,1.5, by = 0.25), at= seq(0,1.5, by = 0.25))#, padj = -0.5)
  
  mtext(side = 3, text = paste('Age uncertainty quantification:', file_name), line = 1, cex = 2)
  mtext(side = 2, text = 'Depth from top [mm]', line = 4, cex = 1.7)
  mtext(side = 1, text = 'Interquartile range [kyrs]', line = 3, cex  =1.7)
}