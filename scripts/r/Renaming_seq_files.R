rm(list=ls())
library(stringr)
args= commandArgs()


##NB this was done after renaming the files with the python script
setwd(args[1])
setwd('~/renamed_NON_COLLAPSED copy/')
dir.create('renamed')
files <- list.files(include.dirs = F)


newname <- c()
for(i in 1:length(files)){
  if(length(grep('GC', files[i]))>0){
    if(length(grep('BLANK', files[i]))>0){
      nom <- gsub('Galaxy12__Filter_sequences_by_length_on_data_11_.GC.EC.7330.Dave.2017.4.','', files[i])
      nom <- gsub('Galaxy12__Filter_sequences_by_length_on_data_11_.GC.EC.7330.Dave.2017.2.','', files[i])
      nom <- gsub('Galaxy12__Filter_sequences_by_length_on_data_11_.GC.EC.7330.Dave.2017.3.','', files[i])
      
      newname <- c(newname, paste('renamed/galaxy_', nom, sep = ''))
    }else{
      newname <- c(newname, paste('renamed/galaxy_',str_split(files[i], pattern = '\\.')[[1]][8],'.fasta', sep = ''))
    }
  }else{
    newname <- c(newname, paste('renamed/galaxy_', str_split(files[i], pattern = '\\.')[[1]][2], '.fasta', sep =''))
    # print(files[i])
    # print(newname)
  }
  if(i>1){
    if(newname[i]==newname[i-1]){
      print(newname[i-1])
      print(newname[i])
    }
  }
  
  file.copy(from=files[i], to = newname[i])
  #print(files[i])
  #print(newname)
}
