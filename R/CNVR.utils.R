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
    cat("db not found")
    db <- list(DB=NULL,exons=NULL,path=path, GenomeDBversion = NULL, Aligner = Aligner)
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
#' @param db exon counts database 
#' @param sbjsBamFiles path to wubjects
#' @param exonsBed it is required that the first three columns should be chromosome, start and end
#' @param save TRUE if you want to save the Database
#' @export
AddNewReference <- function(db, sbjsBamFiles, exonsBed= NULL, saveTo=NULL){
  if(is.null(db$DB)){## la base de datos no existe
    if(is.null(exonsBed)){
      stop("exons gtf annotation should be provided as a data.table")
    }
    hd <- .ReadBamHeader(sbjsBamFiles[1])
    db$GenomeDBversion <- hd$GenomeDBversion
    db$Aligner <- hd$Aligner
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
    data.frame(GenomeDBversion=hd$GenomeDBversion, Aligner=hd$ProgramName)
  })
  
  if(!all(df.info$genomeDBversion == db$genomeDBversion)) stop("All GenomeDBversions should be the same")
  if(!all(df.info$Aligner == db$Aligner)) stop("All Aligners should be the same")
  
  counts <- bplapply(sbjsBamFiles, function(sbj){
    bh <- .ReadBamHeader(sbj)
    if(bh$ProgramName == "subread"){
      db$exons$chromosome <- stringr::str_remove_all(db$exons$chromosome,"chr")
    }
  # sbj <- sbjsBamFiles[1]
    cts <- getBamCounts(bed.frame = db$exons,
                        bam.files = sbj, 
                        include.chr = F, ##assuming that it comes from Aligner (RSubread)
                        referenceFasta = NULL)
    colnames(cts)[ncol(cts)] <- unlist(stringr::str_split(basename(sbj),"_"))[1] 
    return(cts)
   },BPPARAM = MulticoreParam(workers = length(sbjsBamFiles)))
  
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
  
  if(is.null(db$DB)){
    db$DB <- counts[,-c(1:4)]
    db$GenomeDBversion <- df.info$GenomeDBversion[1]
    db$Aligner <- df.info$Aligner[1]
    class(db) <- c("CNVRdb",class(db))
    
  }else{
    db$DB <- cbind(db$DB,counts[,-c(1:4)])
  }
  
  # colnames(db$DB)[4] <- "name"
  
  if(!is.null(saveTo)){## it will overwrite the last path
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
    
    all.exons@CNV.calls$loc <- unlist(lapply(unique(all.exons@CNV.calls$chromosome),function(x){
      ctr <- Centromere$GRCh38$left[Centromere$GRCh38$chr==x]
      loc <- unlist(lapply(all.exons@CNV.calls$end[all.exons@CNV.calls$chromosome==x], function(xx){
        ifelse(xx<ctr,"p","q")
      }))
      return(loc)
    }))
    all.exons@CNV.calls$size <- all.exons@CNV.calls$end - all.exons@CNV.calls$start
    all.exons@CNV.calls$Genes <- .FormatGenedodeColumn(all.exons@CNV.calls$genecode)
    attr(all.exons@CNV.calls,"BamHeader") <- bh
      return(all.exons)
    },BPPARAM = MulticoreParam(workers = length(sbjsBamFiles)))

  names(subjs.calls) <- basename(sbjsBamFiles)
  return(subjs.calls)
  
}


