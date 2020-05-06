###############################################
######## Test CI sisal v2 30.01.2020 ##########
###############################################


### read in chronology ###
sc <- read.csv('/stacydata/data/ariana/Data/SISAL_chronologies/sisal_chrono_4.csv', stringsAsFactors = F) %>% mutate_at(vars(everything()), as.numeric)


###########################
##### test for eID 342 ####
###########################

### filter sisal chronology for eID 342 ###
sc_342 <- sc %>% filter(entity_id == 342) 
print(sc_342[1,])

### read in ensembels and calculate CI and median ages ###
ens_LR <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-linReg.RData'))
m_LR <- median(unlist(ens_LR[1,-c(1,2)]))
q_LR <- quantile(unlist(ens_LR[1,-c(1,2)]), probs = c(0.05, 0.95))
print(c(m_LR, m_LR-q_LR[1],q_LR[2]-m_LR))

ens_LI <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-linInterp.RData'))
m_LI <- median(unlist(ens_LI[1,-c(1,2)]))
q_LI <- quantile(unlist(ens_LI[1,-c(1,2)]), probs = c(0.05, 0.95))
print(c(m_LI, m_LI-q_LI[1],q_LI[2]-m_LI))

ens_copRa <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-copRa.RData'))
m_copRa <- median(unlist(ens_copRa[1,-c(1,2)]))
q_copRa <- quantile(unlist(ens_copRa[1,-c(1,2)]), probs = c(0.05, 0.95))
print(c(m_copRa, m_copRa-q_copRa[1],q_copRa[2]-m_copRa))

ens_bacon <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-Bacon.RData'))
m_bacon <- median(unlist(ens_bacon[1,-c(1,2)]))
q_bacon <- quantile(unlist(ens_bacon[1,-c(1,2)]), probs = c(0.05, 0.95))
print(c(m_bacon, m_bacon-q_bacon[1],q_bacon[2]-m_bacon))

ens_bchron <- get(load('/stacywork/ariana/SISALv2_ensembles/342-KNI-51-H-Bchron.RData'))
m_bchron <- median(unlist(ens_bchron[1,-c(1,2)]))
q_bchron <- quantile(unlist(ens_bchron[1,-c(1,2)]), probs = c(0.05, 0.95))
print(c(m_bchron, m_bchron-q_bchron[1],q_bchron[2]-m_bchron))



