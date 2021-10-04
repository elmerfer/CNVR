#' LoadDB
#' @param path a path direcotry to the DB
#' @export
#' @return an object of class CNVRdb, derived from a list with the following slots
#' DB : a data.table with chromosome star end name subject1 subject2 .....
#' path : path
#' GenomeDBversion 
LoadDB <- function(path){
  if(file.exists(path)==FALSE){
    db <- readRDS(path)
    db$path <- path
  }else{
    db <- list(DB=NULL,exons=NULL,path=path, GenomeDBversion = NULL)
  }
  class(db) <- c("CNVRdb",class(db))
  return(db)
}

.Validate <- function(db){
  if(!any(stringr::str_detect(class(db),"CNVRdb"))){
    stop("it is not a CNVRdb class")
  }
}

#' SaveDB
#' @param db an Obeject CNVRdb
#' @export 
SaveDB <- function(db){
  if(all(stringr::str_detect(class(db),"CNVRdb"))==FALSE){
    stop("it is not a CNVRdb class object")
  }
    
  saveRDS(db,db$path)
}

#' AddNewReference
#' @param db exon counts database 
#' @param sbjsBamFiles path to wubjects
#' @param genomeDBversion genome version see \code{\link[GenomeDB]{GetGenome}}
#' @param exonsBed it is required that the first three columns should be chromosome, start and end
#' @param save TRUE if you whant to save the Database
#' @export
AddNewReference <- function(db, sbjsBamFiles, genomeDBversion = NULL, exonsBed= NULL, save=TRUE){
  if(is.null(db$DB)){## la base de datos no existe
    if(is.null(exonsBed)){
      stop("exons gtf annotation should be provided as a data.table")
    }
    if(is.null(genomeDBversion)){
      stop("genomeDBversion should be provided , pls see GenomeDB library")
    }
    db$GenomeDBversion <- genomeDBversion
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
    sbjsBamFiles <- sbjsBamFiles[!cual]
  } 
  
  counts <- bplapply(sbjsBamFiles, function(sbj){
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
    db$DB <- counts  
  }else{
    db$DB <- cbind(db$DB,counts[,-c(1:4)])
  }
  
  colnames(db$DB)[4] <- "name"
  
  if(save){
    SaveDB(db)
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
  
  exonsBed <- db$DB[,1:4]
  
  
  # my.matrix <- as.matrix( cnvDB[, my.choice$reference.choice, drop = FALSE])
  
  # print(my.choice$reference.choice)
  
  
  
   subjs.calls <-bplapply(sbjsBamFiles, function(sbj){
    # sbj <- sbjsBamFiles
    cts <- getBamCounts(bed.frame = db$exons,
                        bam.files = sbj, 
                        include.chr = F, ##assuming that it comes from Aligner (RSubread)
                        referenceFasta = NULL)
    colnames(cts)[ncol(cts)] <- unlist(stringr::str_split(basename(sbj),"_"))[1] 
    
    if(reference == "auto"){
      my.choice <- select.reference.set (test.counts = cts[,-c(1:4)],
                                         reference.counts = data.matrix(db$DB[,-c(1:4)]) ,
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
     return(all.exons)
   },BPPARAM = MulticoreParam(workers = length(sbjsBamFiles)))

  names(subjs.calls) <- basename(sbjsBamFiles)
  return(subjs.calls)
  
}


