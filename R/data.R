 # centromere.hg19 <- openxlsx::read.xlsx("/home/elmer/Elmer/Rstudio/libraries/CNVR/Centromere.hg19.xlsx")
 # centromere.hg38 <- openxlsx::read.xlsx("/home/elmer/Elmer/Rstudio/libraries/CNVR/Centromere.hg38.xlsx")
 # colnames(centromere.hg19) <-colnames(centromere.hg38) <- c("chr","left","right")
 # Centromere <- list(hg19=centromere.hg19, GRCh38=centromere.hg38)
 # save(Centromere, file="./data/Centromere.RData")

#' Centromere
#' 
#' Information regardin to Centromere position for hg19 and GRCh38
#' 
#' Is a list with two data.frames with three columns "chr", left , right (the last two identifies the position of the base partition. left corresponde to the p. and right to q)
#' @docType data
#' @keywords datasets
#' @name Centromere
#' @usage data(Centromere)
#' @format A data frame with 61 rows (neoantigens) and 10 variables (neoantigens info)
NULL