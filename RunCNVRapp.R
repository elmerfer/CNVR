library(BiocParallel)
library(ExomeDepth)
library(plyr)
# library(CNVR)

.DB.FULL.NAME <- "/media/respaldo4t/FLENI/CNVs/ONCO/CNVDB/CNVDB.CLINICAL_MGEN_"
.LIBRARY <- "TWIST" #"SSV6" #
.DB.FULL.NAME <- paste0(.DB.FULL.NAME,.LIBRARY,".RDS")


db <- CNVR::LoadDB(.DB.FULL.NAME)

dir.path <- rstudioapi::selectDirectory(path = "/media/respaldo4t/FLENI/CNVs/ONCO/PACIENTES_T" )

if(!is.null(dir.path)){
   bam.file <- list.files(dir.path,full.names = T)
   bam.file<- bam.file[stringr::str_detect(basename(bam.file),".BAM|.bam")]
   bam.file<- bam.file[!stringr::str_detect(basename(bam.file),".bai|BAM.|trim")]
   if(!file.exists(bam.file)){
      stop(paste0("ERROR: el archivo ", bam.file," NO encontrado"))
   }
   if(!c(unlist(stringr::str_split(basename(bam.file),"_"))[1] %in% colnames(db$DB))){
      cat(paste0("\nAdding the subject ", dirname(dir.path)," into the database\n"))
      db <- AddNewReference(db,sbjsBamFiles = bam.file )   
   }else{
      cat(paste0("\nThe subject ", dirname(dir.path)," is already in the database\n"))
   }
   cat(paste0("\nNow searching for CNVs on subject ",dirname(dir.path)))
    sbj <- CNVR::CNVcallFromDB(db, unlist(stringr::str_split(basename(bam.file),"_"))[1])
    CNVReport(cnv = sbj[[1]], dir.path)
}else{
  cat("Not patient selected")
}
