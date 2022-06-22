library(BiocParallel)
library(ExomeDepth)
library(plyr)
library(CNVR)
if(Sys.info()["user"]=="biomolecular"){
   data.path <- "~/DATA/NGS/ONCO/"
}else{
   data.path <- "/media/respaldo4t/FLENI/CNVs/ONCO"
}
.DB.FULL.NAME <- file.path(data.path,"/CNVDB/CNVDB.CLINICAL_MGEN_")
.LIBRARY <- "TWIST" #"SSV6" #
.DB.FULL.NAME <- paste0(.DB.FULL.NAME,.LIBRARY,".RDS")


db <- CNVR::LoadDB(.DB.FULL.NAME)

dir.path <- rstudioapi::selectDirectory(path = file.path(data.path,"PACIENTES/" ))

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
    filexlsx <- CNVR::CNVReport(cnv = sbj[[1]], dir.path)
    cat("\n-------------------\n")
    cat(paste0("\nFile created ....:", ifelse(is.null(filexlsx),"NOT CREATED", basename(filexlsx)),"\n"))
    
}else{
  cat("Not patient selected")
}
