##' .ReadBamHeader
##' @param bamFile bam file full path
##' @return a list with the following slots:
##' CHR: the list of chromosomes and sequence names in the bam file
##' ProgramName : the alignment program (PG:PN see [BAM header](https://samtools.github.io/hts-specs/SAMv1.pdf))
##' ProgramVersion : the version of the aligner
##' Code : executed source line code
##' GenomeDBversion : the genome version used (see GenomeDB)
##' GenomeDBpath : the path to genome version
##' 
.ReadBamHeader <- function(bamFile){
  hd <- Rsamtools::scanBamHeader(bamFile)
  ##program name and version
  if(!any(stringr::str_detect(hd,"subread"))){
    return(.ReadMODApyBamHeader(hd))
  }
  ##aqui si usa Rsubread
  pname <- stringr::str_remove_all(grep(pattern = "PN", unlist(hd),value=TRUE),"PN:")
  pversion <- stringr::str_remove_all(grep(pattern = "VN", unlist(hd),value=TRUE),"VN:")
  bam.chrs <- names(hd[1]$targets)
  CL <- grep(pattern = "CL", unlist(hd),value=TRUE)
  id.genome<-which(stringr::str_detect(unlist(stringr::str_split(CL," ")),"-i"))       
  gpath <- stringr::str_remove_all(unlist(stringr::str_split(CL," "))[id.genome+1],"\"")
  genome <- stringr::str_remove_all(basename(gpath),"\"")
  return(list(CHR=bam.chrs,ProgramName = pname, ProgramVersion=pversion, Code =CL, GenomeDBversion=genome, GenomeDBpath = gpath))
}

.ReadMODApyBamHeader <- function(hdBam){
  hdt <- hdBam[[1]]$text
  ProgramName <- plyr::ldply(hdt[(names(hdt) %in% "@PG")],function(x){
    x<-unlist(x)
    pname <- stringr::str_remove_all(x[which(stringr::str_detect(x,"ID:"))[1]],"ID:")
    version <- stringr::str_remove_all(x[which(stringr::str_detect(x,"VN:"))[1]],"VN:")
    code <- stringr::str_remove_all(x[which(stringr::str_detect(x,"CL:"))[1]],"CL:")
    return(data.frame(ProgramName=pname, Version=version,Code=code))
  })
  
  id.aligner <- which(ProgramName$ProgramName %in% c("bwa","subread","STAR"))
  aligner<- ProgramName$ProgramName[id.aligner]
  version <- ProgramName$Version[id.aligner]
  code <- ProgramName$Code[id.aligner]
  get.bwa.reference <- function(CL){
    id.genome<-which(stringr::str_detect(unlist(stringr::str_split(CL," ")),"fa|fasta"))
    gpath <- stringr::str_remove_all(unlist(stringr::str_split(CL," "))[id.genome],"\"")
    gpath<- gpath[!stringr::str_detect(gpath,".fastq")]
    genome <- stringr::str_remove_all(basename(gpath),"\"")
    return(list(genome=genome, gpath=gpath))
  }
  reference <- switch(aligner,
                      "bwa" = get.bwa.reference(code),
                      "subread")
  
  
  pname <- unlist(ProgramName$ProgramName)
  pversion <- unlist(ProgramName$version)
  bam.chrs <- names(hdBam[[1]]$targets)
  CL <- unlist(ProgramName$Code)
  return(list(CHR=bam.chrs,ProgramName = pname, ProgramVersion=pversion, Code =CL, GenomeDBversion=reference$genome, GenomeDBpath = reference$gpath))
}
##'.FormatGenedodeColumn
##'internal function to format the output of CNVcalls 
##'@param genecode the column of the CNVcall object
.FormatGenecodeColumn <- function(genecode){
  unlist(lapply(genecode, 
                   function(x){
                     genes <- as.character(unique(lapply(unlist(stringr::str_split(x,",")),function(xv){unlist(stringr::str_split(xv,"_"))[1]})))
                     data.frame(genes=paste0(genes,collapse = "|"))
                   }
  ))
}