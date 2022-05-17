#' LoadDB
#' @param path a path directory to the DB
#' @export
#' @return an object of class CNVRdb, derived from a list with the following slots
#' It should hold the following slots:
#' DB : (data.matrix) the exon counts data matrix of the reference subjects (ExS)
#' exons : (data.frame) annotation information chromosome, start, end , strand , ...
#' path : (character) full path where the library is saved
#' GenomeDBversion : the genome DB version (see \code{\link[GenomeDB]{GetGenome}}) (all the samples should be processed by the same GenomeDBversion)
#' Aligner : the used aligner (all the samples should be processed by the same Aligner)
LoadDB <- function(path){
  if(file.exists(path)==TRUE){
    db <- readRDS(path)
    db$path <- path
    # db <- list(DB=db, exons = db[,1:4],path=path, GenomeDBversion = NULL, Aligner = NULL)
  }else{
    cat("db not found")#la base de datos no existe
    db <- list(DB=NULL,exons=NULL,path=path, GenomeDBversion = NULL, Aligner = NULL)
    class(db) <- c("CNVRdb",class(db))
  }
  if(!.Validate(db)) stop("not valid DB (LoadDB")
  return(db)
}

#' .Validate
#' it validate the existence of the data base as well as its correct class
#' this is an internal function, not for users
#' @param db an CNVRdb object
#' @return 
#' it returns FALSE if all the elements are NULL, TRUE otherwise or stop if the db is not CNVRdb
.Validate <- function(db){
  if(all(unlist(lapply(db, is.null)))){
    return(FALSE)
  }
  if(!any(stringr::str_detect(class(db),"CNVRdb"))){
    stop("it is not a CNVRdb class")
  }
  return(TRUE)
}

#' SaveDB
#' It saves the reference data base of exon counts for CNV detection
#' It should hold the following slots:
#' DB : (data.matrix) the exon counts data matrix of the reference subjects (ExS)
#' exons : (data.frame) annotation information chromosome, start, end , strand , ...
#' path : (character) full path where the library is saved
#' GenomeDBversion : the genome DB version (see \code{\link[GenomeDB]{GetGenome}}) (all the samples should be processed by the same GenomeDBversion)
#' Aligner : the used aligner (all the samples should be processed by the same Aligner)
#' @param db an Obeject CNVRdb
#' 
#' @export 
SaveDB <- function(db){
  if(!.Validate(db)) stop("bad db format (SaveDB")
  if(any(stringr::str_detect(class(db),"CNVRdb"))==FALSE){
    stop("it is not a CNVRdb class object")
  }
  if(!is.null(db$path)){
    saveRDS(db,db$path)  
  }else{
    db$path <- file.path(getwd(),"CNVRdb.RData")
    saveRDS(db,db$path)
    cat(paste0("DB saved at ",db$path, "pls change to a better location"))
  }
  
}

#' AddNewReference
#' This function will create or increase the Exon counts data base.
#' @param db exon counts database 
#' @param sbjsBamFiles path to subjects
#' @param exonsBed it is required that the first three columns should be chromosome, start and end
#' @param saveTo (character) default NULL. The path to save the database. It is only valid if the database do not exists
#' @export
AddNewReference <- function(db, sbjsBamFiles, exonsBed= NULL, saveTo=NULL){
  if(is.null(db$DB)){## la base de datos no existe
    if(is.null(exonsBed)){
      stop("exons gtf annotation should be provided as a data.table")
    }
    ##lee la información de alineamiento desde el BAM. La almacena en la base de datos para 
    ##tener en cuenta que deben utilizarse el mismo alineador
    ##si la base de datos no existe,m, le asigna la del primer sujeto
    hd <- .ReadBamHeader(sbjsBamFiles[1])
    db$GenomeDBversion <- hd$GenomeDBversion
    db$Aligner <- hd$ProgramName[hd$ProgramName %in% c("bwa","subread","STAR")]
    
    db$exons <- exonsBed
  }else{## la base de datos existe
    exonsBed <- db$exons
    ##check names
    fn <- unlist(lapply(sbjsBamFiles, function(sbj){
      unlist(stringr::str_split(basename(sbj),"_"))[1] 
    }))
    cual <- fn %in% colnames(db$DB)
    if(all(cual==TRUE)){
      stop("All subjects in DB")
    }
    #filter out those files already present in the database
    sbjsBamFiles <- sbjsBamFiles[!cual]
  } 
  ## check for GenomeDB and Aligner of BAM files
  df.info <- plyr::ldply(sbjsBamFiles, function(bfile){
    hd <- .ReadBamHeader(bfile)
    data.frame(GenomeDBversion=hd$GenomeDBversion, Aligner=hd$ProgramName[hd$ProgramName %in% c("bwa","subread","STAR")])
  })
  rownames(df.info) <- basename(sbjsBamFiles)
  
  # if(!all(df.info$GenomeDBversion == db$GenomeDBversion)) stop("All GenomeDBversions should be the same")
  #esto esta comentado porque no son iguales siempre (REVISAR)
  
  ##esto tiene mucha mas sentido, no mezclar alineadores
  if(!all(df.info$Aligner == db$Aligner)) stop("All Aligners should be the same")
  
   counts <- bplapply(sbjsBamFiles, function(sbj){
    # sbj <- sbjsBamFiles[1]
    bh <- .ReadBamHeader(sbj)
    if(all(stringr::str_detect(bh$CHR,"chr"))){
      if(!all(stringr::str_detect(db$exons$chromosome,"chr"))){
        db$exons$chromosome <- paste0("chr",db$exons$chromosome)
      }
    }else{
      if(all(stringr::str_detect(db$exons$chromosome,"chr"))){
        db$exons$chromosome <- stringr::str_remove_all(db$exons$chromosome,"chr")
      }
    }
    
    cts <- ExomeDepth::getBamCounts(bed.frame = db$exons,
                        bam.files = sbj, 
                        include.chr = F, ##assuming that it comes from Aligner (RSubread)
                        referenceFasta = NULL)
    colnames(cts)[ncol(cts)] <- unlist(stringr::str_split(basename(sbj),"_"))[1] 
    attr(cts[,colnames(cts)[ncol(cts)]],"BamHeader") <- bh
     return(cts)
    },BPPARAM = MulticoreParam(workers = length(sbjsBamFiles)))
  ##---- fin aparallel
   names(counts) <- sbjsBamFiles
  if(length(counts)>1){
    counts <- cbind(counts[[1]],lapply(2:length(counts), function(x) {
      a <- counts[[x]][,5, drop=T]
      dim(a) <- c(length(a),1)
      colnames(a) <- colnames(counts[[x]])[5]
      return(a)
    }))
  }else{
    counts <- counts[[1]]
  }
  
  if(is.null(db$DB)){#si la DB no existe, la crea
    db$DB <- counts[,-c(1:4)]
    db$GenomeDBversion <- unlist(df.info[basename(sbjsBamFiles),"GenomeDBversion"])##ojo aqui, porque pueden haber sido procesados con distintos genomas
    ##esto estaba pensado mas que nada para un solo tipo de pipeline. que seria lo esperable de ahroa en mas.
    db$Aligner <- unlist(df.info[basename(sbjsBamFiles),"Aligner"])
    class(db) <- c("CNVRdb",class(db))
    
  }else{##si existe agrega el sujeto
    db$DB <- cbind(db$DB,counts[,-c(1:4)])
    db$GenomeDBversion <- c(db$GenomeDBversion,df.info[basename(sbjsBamFiles),"GenomeDBversion"])
    colnames(db$DB) <- c(colnames(db$DB[-ncol(db$DB)]),colnames(counts)[5])
  }
  
  # colnames(db$DB)[4] <- "name"
  
  if(!is.null(saveTo) & is.null(db$path)){## it will overwrite the last path
    db$path <- saveTo
    SaveDB(db)
  }else{
    if(!is.null(db$path)){
      SaveDB(db)
    }
  }
  # class(db) <- c("CNVRdb",class(db))
  return(db)
}

#' CNVcall  
#' This function execute the CNV call based on ExomeDepth. 
#' It assumes that the Subject is not present on the current Exon Count DB.
#' If the subjects is already in the DB, it will stop.
#' @param db exon counts database
#' @param sbjsBamFiles subject BAM files
#' @param reference the reference type 
#' @export
CNVcall <- function(db, sbjsBamFiles, reference = c("auto","all")){
  .Validate(db)
  fn <- unlist(lapply(sbjsBamFiles, function(sbj){
    unlist(stringr::str_split(basename(sbj),"_"))[1] 
  }))
  cual <- fn %in% colnames(db$DB)
  if(all(cual==TRUE)){
    stop("All subjects in DB")
  }
  
  sbjsBamFiles <- sbjsBamFiles[!cual]
  
  reference <- match.arg(tolower(reference)[1],c("auto","all"))
  
  # exonsBed <- db$DB[,1:4]
  
  
  # my.matrix <- as.matrix( cnvDB[, my.choice$reference.choice, drop = FALSE])
  
  # print(my.choice$reference.choice)
  
  
  
   subjs.calls <-bplapply(sbjsBamFiles, function(sbj){
     # sbj <- sbjsBamFiles
     bh <- .ReadBamHeader(sbj)
     if(db$GenomeDBversion != bh$GenomeDBversion) stop()
     if(db$Aligner != bh$ProgramName) stop()
     if(bh$ProgramName == "subread"){
       db$exons$chromosome <- stringr::str_remove_all(db$exons$chromosome,"chr")
     }
     Genome <- ifelse(stringr::str_detect(bh$GenomeDBversion,"GRCh38"),"GRCh38","hg19")
     
    cts <- getBamCounts(bed.frame = db$exons,
                        bam.files = sbj, 
                        include.chr = F, ##assuming that it comes from Aligner (RSubread)
                        referenceFasta = NULL)
    colnames(cts)[ncol(cts)] <- unlist(stringr::str_split(basename(sbj),"_"))[1] 
    
    if(reference == "auto"){
      my.choice <- select.reference.set (test.counts = cts[,-c(1:4)],
                                         reference.counts = data.matrix(db$DB) ,
                                         bin.length = (db$exons$end - db$exons$start)/1000,
                                         n.bins.reduced = 10000)  
    }else{
      my.choice <- list()
      my.choice$reference.choice <- colnames(ref.set)
    }
    my.reference.set <- as.matrix(db$DB[,my.choice$reference.choice, drop=FALSE])
    my.reference.selected <- apply(X=my.reference.set,MAR=1,FUN=sum)
    all.exons <- new('ExomeDepth',
                     test = cts[,5],
                     reference = my.reference.selected,
                     formula = 'cbind(test, reference) ~ 1')
    
    all.exons <- CallCNVs(x = all.exons,
                          transition.probability = 10^-4,
                          chromosome = db$exons$chromosome,
                          start = db$exons$start,
                          end = db$exons$end,
                          name = db$exons$name)
    
    
    exons.gencode.granges <- GenomicRanges::GRanges(seqnames = db$exons$chromosome,
                                                    IRanges::IRanges(start=db$exons$start,end=db$exons$end),
                                                    names = db$exons$name)
    
    all.exons <- AnnotateExtra(x = all.exons,reference.annotation = exons.gencode.granges,
                               min.overlap = 0.0001,
                               column.name = 'genecode')  
    
    all.exons@annotations$freq <- all.exons@test/(all.exons@reference + all.exons@test)
    all.exons@annotations$ratio <- all.exons@annotations$freq/all.exons@expected                                                           
    all.exons@annotations$LR <- log2(all.exons@annotations$ratio)
    all.exons@annotations$type <- NA
    all.exons@annotations$type[all.exons@annotations$LR > 0] <- "Gain"
    all.exons@annotations$type[all.exons@annotations$LR < 0] <- "Loss"
    all.exons@annotations$type[!c(all.exons@annotations$freq > 0)] <- NA
    all.exons@annotations$type[!c(all.exons@annotations$ratio > 0)] <- NA
    all.exons@annotations$counts.weighted <- (all.exons@test + all.exons@reference) * all.exons@expected
    #cuidado aqui, que siempre sea GRCh38
    
    all.exons@CNV.calls$loc <- unlist(lapply(unique(all.exons@CNV.calls$chromosome),function(x){
      ctr <- CNVR::Centromere[[Genome]]$left[CNVR::Centromere[[Genome]]$chr==x]
      loc <- unlist(lapply(all.exons@CNV.calls$end[all.exons@CNV.calls$chromosome==x], function(xx){
        ifelse(xx<ctr,"p","q")
      }))
      return(loc)
    }))
    all.exons@CNV.calls$size <- all.exons@CNV.calls$end - all.exons@CNV.calls$start
    all.exons@CNV.calls$Genes <- .FormatGenecodeColumn(all.exons@CNV.calls$genecode)
    attr(all.exons@CNV.calls,"BamHeader") <- bh
      return(all.exons)
    },BPPARAM = MulticoreParam(workers = length(sbjsBamFiles)))

  names(subjs.calls) <- basename(sbjsBamFiles)
  return(subjs.calls)
  
}

#' CNVcallFromDB  
#' This function execute the CNV call based on ExomeDepth. The subjects should be already in the Exon DB count database \code{\link{LoadDB}}
#' @usage 
#' CNVcallFromDB(db, dbSubjects, reference = c("auto","all"))
#' It assumes that the Subject is not present on the current Exon Count DB.
#' If the subjects is already in the DB, it will stop.
#' @param db exon counts database
#' @param dbSubjects the subjects to test
#' @param reference the reference type 
#' @details 
#' The subjects should be already in the DB by calling \code{\link{AddNewReference}}. Then the test subject should 
#' be identified by it colnames(db$DB).
#' if moretha one subject is referred, then all the subjects in the DB, except the one tested, 
#' will be used to build the reference according to the "reference" type used.
#' if reference == auto , then the refernce will be chosen by correlation
#' if reference == all , all the remainders subjects will be summed up to build the reference
#' @export
#' @examples 
#' \dontrun{
#' DB <- LoadDB("path/to/db/RDS")
#' DB <- AddNewReference(db=DB, 
#' sbjsBamFiles = "/..../PATIENTS/27040/27040_1.fastq.gz.subread.BAM",save=TRUE)
#' CNVs <- CNVcallFromDB(DB, "27040")
#' }
CNVcallFromDB <- function(db, dbSubjects, reference = c("auto","all")){
  .Validate(db)
  sbjsBamFiles <- dbSubjects
  fn <- unlist(lapply(sbjsBamFiles, function(sbj){
    unlist(stringr::str_split(basename(sbj),"_"))[1] 
  }))
  cual <- fn %in% colnames(db$DB)
  if(any(cual==FALSE)){
    cat("The following subjects are not in DB :\n" , paste0(fn[cual==FALSE],"\n"))
  }
  
  sbjsBamFiles <- sbjsBamFiles[cual]
  
  reference <- match.arg(tolower(reference)[1],c("auto","all"))
  
  # exonsBed <- db$DB[,1:4]
  
  
  # my.matrix <- as.matrix( cnvDB[, my.choice$reference.choice, drop = FALSE])
  
  # print(my.choice$reference.choice)
  
  ## No tengo info sobre el BamHeader cuando el paciente esta en la base de datos, debo agregar esa data alli.
  
  subjs.calls <-bplapply(sbjsBamFiles, function(sbj){
    cat("\nProcessing : ", sbj,"\n")
    # sbj <- sbjsBamFiles
    bh <- attr(db$DB[,sbj],"BamHeader")
    # if(db$GenomeDBversion != bh$GenomeDBversion) stop()
    # if(db$Aligner != bh$ProgramName) stop()
    # if(bh$ProgramName == "subread"){
    #   db$exons$chromosome <- stringr::str_remove_all(db$exons$chromosome,"chr")
    # }
    Genome <- ifelse(stringr::str_detect(db$GenomeDBversion,"GRCh38"),"GRCh38","hg19")[which(colnames(db$DB) %in% sbj)]
    
    # cts <- getBamCounts(bed.frame = db$exons,
    #                     bam.files = db$exons[,], 
    #                     include.chr = F, ##assuming that it comes from Aligner (RSubread)
    #                     referenceFasta = NULL)
    # colnames(cts)[ncol(cts)] <- unlist(stringr::str_split(basename(sbj),"_"))[1] 
    
    if(reference == "auto"){
      my.choice <- select.reference.set (test.counts = data.matrix(db$DB[,sbj]),
                                         reference.counts = data.matrix(db$DB[,colnames(db$DB)!=sbj]) ,
                                         bin.length = (db$exons$end - db$exons$start)/1000,
                                         n.bins.reduced = 10000)  
    }else{
      my.choice <- list()
      my.choice$reference.choice <- colnames(ref.set)
    }
    my.reference.set <- as.matrix(db$DB[,my.choice$reference.choice, drop=FALSE])
    my.reference.selected <- apply(X=my.reference.set,MAR=1,FUN=sum)
    all.exons <- new('ExomeDepth',
                     test = db$DB[,sbj],
                     reference = my.reference.selected,
                     formula = 'cbind(test, reference) ~ 1')
    
    all.exons <- CallCNVs(x = all.exons,
                          transition.probability = 10^-4,
                          chromosome = db$exons$chromosome,
                          start = db$exons$start,
                          end = db$exons$end,
                          name = db$exons$name)
    
    
    exons.gencode.granges <- GenomicRanges::GRanges(seqnames = db$exons$chromosome,
                                                    IRanges::IRanges(start=db$exons$start,end=db$exons$end),
                                                    names = db$exons$name)
    
    all.exons <- AnnotateExtra(x = all.exons,reference.annotation = exons.gencode.granges,
                               min.overlap = 0.0001,
                               column.name = 'genecode')  
    
    all.exons@annotations$freq <- all.exons@test/(all.exons@reference + all.exons@test)
    all.exons@annotations$ratio <- all.exons@annotations$freq/all.exons@expected                                                           
    all.exons@annotations$LR <- log2(all.exons@annotations$ratio)
    all.exons@annotations$type <- NA
    all.exons@annotations$type[all.exons@annotations$LR > 0] <- "Gain"
    all.exons@annotations$type[all.exons@annotations$LR < 0] <- "Loss"
    all.exons@annotations$type[!c(all.exons@annotations$freq > 0)] <- NA
    all.exons@annotations$type[!c(all.exons@annotations$ratio > 0)] <- NA
    all.exons@annotations$counts.weighted <- (all.exons@test + all.exons@reference) * all.exons@expected
    #cuidado aqui, que siempre sea GRCh38
    
    all.exons@CNV.calls$loc <- unlist(lapply(unique(all.exons@CNV.calls$chromosome),function(x){
      ctr <- CNVR::Centromere[[Genome]]$left[CNVR::Centromere[[Genome]]$chr==x]
      loc <- unlist(lapply(all.exons@CNV.calls$end[all.exons@CNV.calls$chromosome==x], function(xx){
        ifelse(xx<ctr,"p","q")
      }))
      return(loc)
    }))
    all.exons@CNV.calls$size <- all.exons@CNV.calls$end - all.exons@CNV.calls$start
    all.exons@CNV.calls$Genes <- CNVR:::.FormatGenecodeColumn(all.exons@CNV.calls$genecode)
    attr(all.exons@CNV.calls,"BamHeader") <- bh
    attr(all.exons@reference,"DBreference") <- my.choice$reference.choice
    # CNVReport(all.exons)
    # print("PASEEEEE")
    return(all.exons)
  },BPPARAM = MulticoreParam(workers = length(sbjsBamFiles)))
  
  names(subjs.calls) <- basename(sbjsBamFiles)
  return(subjs.calls)
  
}

#' CNVReport
#' @usage CNVReport(cnv)
#' @param cnv an ExomeDepth CNV call object \code{\link[ExomeDepth]{CallCNVs}}
#' @param sbjPath the subject path directory
#' @export
#' @return 
#' It creates an Excel file in the subject directory. The file is named "subject_CNVs.xlsx"
#' @examples 
#' \dontrun{
#' DB <- LoadDB("path/to/db/RDS")
#' DB <- AddNewReference(db=DB, 
#' sbjsBamFiles = "/..../PATIENTS/27040/27040_1.fastq.gz.subread.BAM",save=TRUE)
#' CNVs <- CNVcallFromDB(DB, "27040")
#' CNVReport(CNVs[[1]])
#' }
CNVReport <- function(cnv,sbjPath){
  bh <- attr(cnv@CNV.calls, "BamHeader")
  if(is.null(bh)){
    aligment.code <- "UNKNOWN"
    subject<- ""
    warning("header information not found")
  }else{
    aligment.code <- stringr::str_split(bh$Code[stringr::str_detect(bh$Code,"_1.fastq|_R1.fastq|fq")]," ")
    subject<- basename(unlist(aligment.code)[which(stringr::str_detect(unlist(aligment.code),"_1.fastq|_R1.fastq|fq"))[1]])
  }
  
  if(missing(sbjPath)){
    sbjPath <- paste0(subject,"_CNVs.xlsx")
  }else{
    sbjPath <- file.path(sbjPath,paste0(subject,"_CNVs.xlsx"))
  }
  cnorder <-c("chromosome","id","type","Genes","nexons","loc","size","BF","reads.ratio","reads.expected","reads.observed", "start.p","end.p","start","end","genecode")
  openxlsx::write.xlsx(list(CNV=cnv@CNV.calls[,cnorder],Code=aligment.code), sbjPath, overwrite = T)
}

