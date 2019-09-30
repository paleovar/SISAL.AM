####  Sample script

wd <- '...' #choose a directory to work in

eId <- filter_SISAL(1) %>%
  filter(entity_id < 20) # run a selection of records of class 1

AM_SISAL('',eID, file_path = wd) # run age models and create folder to save them in; creates folder in current working directory

run <- read.csv(file = file.path(wd, "runFile.csv"),header = T)
m <- merge_SISAL_chrono(run)

r <- as.data.frame(m[1]) # checked runFile; if error is thrown after the chronology was built, the age model is marked as failed in the previous runFile
s <- as.data.frame(m[2]) # built chronolgy

# save files
write.csv(s, file = file.path(wd,"chronology.csv"), row.names = F)
write.csv(r, file = file.path(wd, "runFile_final.csv"), row.names = F)

evaluation <- eval(s,r) # evaluates sisal chronolgy for entity ids included in r
plot <- plot_sisal_eval(as.data.frame(evaluation[5])) # plot with ggplot the evaluated chronology

plot + ggsave(filename = file.path(wd,"evaluation.pdf")) # saves the plot to the working directory
