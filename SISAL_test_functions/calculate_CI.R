###############################################
######## Calculate CI sisal v2 30.01.2020 ##########
###############################################

###########################
######### example eID 342 #########
###########################

### read in ensembels and calculate 95% CI and median ages ###
ens_LR <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-linReg.RData'))
m_LR <- apply(ens_LR[,-c(1,2)],1 , function(x){median(unlist(x))}) # medium ages
q_LR <- apply(ens_LR[,-c(1,2)],1 , function(x){quantile(unlist(x), probs = c(0.025, 0.975))}) # 2.5% and 97.5% quatiles
LR_chrono <- data.frame(sample_id = ens_LR[,1], depth_sample = ens_LR[,2],
                        lin_reg_age = m_LR, lin_reg_age_uncert_pos = q_LR[2,]-m_LR, lin_reg_age_uncert_neg = m_LR-q_LR[1,]) # create data frame

ens_LI <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-linInterp.RData'))
m_LI <- apply(ens_LI[,-c(1,2)],1 , function(x){median(unlist(x))})
q_LI <- apply(ens_LI[,-c(1,2)],1 , function(x){quantile(unlist(x), probs = c(0.025, 0.975))})
LI_chrono <- data.frame(sample_id = ens_LI[,1], depth_sample = ens_LI[,2],
                        lin_interp_age = m_LI, lin_interp_age_uncert_pos = q_LI[2,]-m_LI, lin_interp_age_uncert_neg = m_LI-q_LI[1,])

ens_copRa <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-copRa.RData'))
m_copRa <- apply(ens_copRa[,-c(1,2)],1 , function(x){median(unlist(x))})
q_copRa <- apply(ens_copRa[,-c(1,2)],1 , function(x){quantile(unlist(x), probs = c(0.025, 0.975))})
copRa_chrono <- data.frame(sample_id = ens_copRa[,1], depth_sample = ens_copRa[,2],
                        copRa_age = m_copRa, copRa_age_uncert_pos = q_copRa[2,]-m_copRa, copRa_age_uncert_neg = m_copRa-q_copRa[1,])

ens_bacon <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-Bacon.RData'))
m_bacon <- apply(ens_bacon[,-c(1,2)],1 , function(x){median(unlist(x))})
q_bacon <- apply(ens_bacon[,-c(1,2)],1 , function(x){quantile(unlist(x), probs = c(0.025, 0.975))})
bacon_chrono <- data.frame(sample_id = ens_bacon[,1], depth_sample = ens_bacon[,2],
                        bacon_age = m_bacon, bacon_age_uncert_pos = q_bacon[2,]-m_bacon, bacon_age_uncert_neg = m_bacon-q_bacon[1,])

ens_bchron <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-Bchron.RData'))
m_bchron <- apply(ens_bchron[,-c(1,2)],1 , function(x){median(unlist(x))})
q_bchron <- apply(ens_bchron[,-c(1,2)],1 , function(x){quantile(unlist(x), probs = c(0.025, 0.975))})
bchron_chrono <- data.frame(sample_id = ens_bchron[,1], depth_sample = ens_bchron[,2],
                        bchron_age = m_bchron, bchron_age_uncert_pos = q_bchron[2,]-m_bchron, bchron_age_uncert_neg = m_bchron-q_bchron[1,])
