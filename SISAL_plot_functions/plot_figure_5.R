################################
####### plot Figure 5 ##########
################################

plot_figure_5 <- function(eID = 342,
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
    layout(matrix(c(1,1,2,3,3,4), nrow=2, ncol=3, byrow = T))
    par(mar = c(4,3,3,4), oma = c(1,3,2,1))
    x.lim <- plot_am(chrono, i, sisal$dating, entity_name, fig5 = T)
    if(any(c(!is_empty(data.frame(df_fil$Bacon)),!is_empty(data.frame(df_fil$Bchron)),!is_empty(data.frame(df_fil$copRa)),!is_empty(data.frame(df_fil$StalAge)),!is_empty(data.frame(df_fil$linInterp)),
             !is_empty(data.frame(df_fil$linReg)),!is_empty(data.frame(df_fil$OxCal))))){
      plot_iqr(df_fil, i, sisal$dating, entity_name, fig5 = T)
    } else {
      plot.new()
      plot.new()
    }
    plot_isotopes_figure_5(df_fil, x.lim)
    #plot_dating(df_fil, i, sisal$dating)
    plot.new()

    orig_am <- df_fil$origAM %>% filter(age_model_type != 'NA') %>% distinct(age_model_type)

    legend('center',legend=c(paste('original AM:', orig_am),'lin. interp. median age','lin. reg. median age', 'Bchron median age',
                             'Bacon median age', 'OxCal median age','copRa median age','StalAge mean age'),
           lty=c(1,1,1,1,1,1,1,1,3,3),cex=1,
           col = c('black','royalblue3','skyblue','forestgreen','yellowgreen','lightcoral', 'hotpink4','grey'),
           bty = "n",  ncol = 1, lwd = 4)

    dev.off()


  }

}

plot_isotopes_figure_5 <- function(chrono, x.lim){
  iso<- chrono$proxy

  color = c('black','hotpink4', 'grey','royalblue3', 'skyblue', 'yellowgreen', 'forestgreen', 'lightcoral')

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
          y = iso$d18O_measurement,
          col = color,
            #c('black',color[8],'grey',color[c(2,1,3,4,5)]),
          lty = 1, type = 'l', lwd = 1.5,
          xlim = x.lim,
          ylim = c(max(range(iso$d18O_measurement, na.rm = T), na.rm = T),min(range(iso$d18O_measurement, na.rm = T), na.rm = T)), xlab = '', ylab = '')
  mtext(side = 1, line = 2, text = 'Age [yr BP]')
  mtext(side = 2, line = 2, text = expression(paste(delta^{18},'O [\u2030]')))

  mtext(line = 0, text = expression(bold(Isotopes)), adj = 0)

}
