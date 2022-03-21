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
  pname <- stringr::str_remove_all(grep(pattern = "PN", unlist(hd),value=TRUE),"PN:")
  pversion <- stringr::str_remove_all(grep(pattern = "VN", unlist(hd),value=TRUE),"VN:")
  bam.chrs <- names(hd[1]$targets)
  CL <- grep(pattern = "CL", unlist(hd),value=TRUE)
  id.genome<-which(stringr::str_detect(unlist(stringr::str_split(CL," ")),"-i"))       
  gpath <- stringr::str_remove_all(unlist(stringr::str_split(CL," "))[id.genome+1],"\"")
  genome <- stringr::str_remove_all(basename(gpath),"\"")
  return(list(CHR=bam.chrs,ProgramName = pname, ProgramVersion=pversion, Code =CL, GenomeDBversion=genome, GenomeDBpath = gpath))
}

##'.FormatGenedodeColumn
##'internal function to format the output of CNVcalls 
##'@param genecode the column of the CNVcall object
.FormatGenedodeColumn <- function(genecode){
  unlist(lapply(genecode, 
                   function(x){
                     genes <- as.character(unique(lapply(unlist(stringr::str_split(x,",")),function(xv){unlist(stringr::str_split(xv,"_"))[1]})))
                     data.frame(genes=paste0(genes,collapse = "|"))
                   }
  ))
}