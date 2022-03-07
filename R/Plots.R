




# anno <- x@annotations
#' GetCNVsAnnotation
#' internal function
#' @param chr (chromosome) any of chr1..chr22
#' @param x a CNV object from ExonDepth
#' @param countThreshold (default = 10). The minimum expected exon count
#' @return an annotation data frame for the specified chromosome to be used for plotting
#' @export
GetCNVsAnnotation <- function(chr, x, countThreshold =10){
  anno <- x@annotations
  selected <- which(anno$chromosome == chr & (x@test + x@reference) * 
                      x@expected > countThreshold)
  
  if (length(selected) == 0) {
    warning("No exon seem to be located in the requested region. It could also be that the read count is too low for these exons? In any case no graph will be plotted.")
    return(NULL)
  }
  anno <- anno[selected, ]
  anno$expected <- x@expected[selected]
  anno$freq <- x@test[selected]/(x@reference[selected] + 
                                   x@test[selected])
  anno$middle <- 0.5 * (anno$start + anno$end)
  anno$ratio <- anno$freq/anno$expected
  anno$test <- x@test[selected]
  anno$reference <- x@reference[selected]
  anno$total.counts <- anno$test + anno$reference
  if (length(x@phi) == 1) 
    anno$phi <- x@phi
  else anno$phi <- x@phi[selected]
  anno <- cbind(anno,plyr::ldply(1:nrow(anno), function(i){
    c(my.min.norm.prop=qbetabinom(p = 0.025, size = anno$total.counts[i], 
                                  phi = anno$phi[i], prob = anno$expected[i]),
      my.max.norm.prop=qbetabinom(p = 0.975, size = anno$total.counts[i], 
                                  phi = anno$phi[i], prob = anno$expected[i]) )
  }))
  
  anno$my.min.norm.prop <- anno$my.min.norm/anno$total.counts
  anno$my.max.norm.prop <- anno$my.max.norm/anno$total.counts
  
  return(list(Exons=anno, CNVs=subset(x@CNV.calls, chromosome == chr)))
}

GetAllExons <- function(cnvCalls, countThreshold =10){
  chrs <- unique(cnvCalls$CNVs@CNV.calls$chromosome)
  chr.exons <- bplapply(chrs, GetCNVsAnnotation, x=cnvCalls$CNVs, BPPARAM = bpparam())
  names(chr.exons) <- chrs
  return(invisible(chr.exons))
}

.PlotCNVchromosome <- function(annot,thLength=500, geneTable){
  ##assume chrome.bed
  # annot <- chr19
  # thLength <- 500
  # annot <- CHRi
  CNVs <- annot$CNVs
  CNVs$width <- (CNVs$end-CNVs$start)
  chr <- stringr::str_replace_all(annot$Exons$chromosome[1],"chr","Chr")
  CNVs <- subset(CNVs,width > thLength)
  
  annot <- annot$Exons
  ordm <- order(annot$middle)
  puntos <- data.frame(x=annot$middle[ordm], y=annot$ratio[ordm])
  puntos$colscale <- puntos$y
  puntos$colscale[puntos$colscale>2] <- 2
  # puntos$Type <- NA
  # puntos$Type[puntos$y > annot$my.max.norm.prop] <- "Gain"
  # puntos$Type[puntos$y < annot$my.min.norm.prop] <- "Loss"
  puntos$ls <- annot$my.max.norm.prop[ordm]/annot$expected[ordm]
  puntos$li <- annot$my.min.norm.prop[ordm]/annot$expected[ordm]
  # lsup <- loess(ls~x,puntos, span = 0.3)
  # linf <- loess(li ~x, puntos, span = 0.3)
  # df.l<- data.frame(x=puntos$x,y=lsup$fitted)
  # df.li<- data.frame(x=puntos$x,y=linf$fitted)
  p <- ggplot(puntos, aes(x,y,colour=colscale)) + geom_point(size=0.5) + scale_color_gradient2(midpoint = 1, low = "red", high = "blue") + 
    geom_smooth(method="loess",span=0.1, se=F,col="black", size=0.5,linetype = "dashed") + 
    geom_smooth(data = puntos,  method="loess",span=0.1, se=F, aes(x=x,y=ls),col="aquamarine4", size=0.5,linetype = "dashed") +
    geom_smooth(data = puntos,  method="loess",span=0.1, se=F, aes(x,y=li),col="aquamarine4", size=0.5,linetype = "dashed")# +
  
  # geom_vline(xintercept = chrom.bed[chr,"mid.point"], linetype = "dashed" ) + xlab("Chr. position") + ylab("log ratio")
  if(nrow(CNVs)>0){
    if(any(CNVs$reads.ratio>1)){
      p <- p + geom_segment(data=subset(CNVs,reads.ratio>1), aes(x = start, y = reads.ratio, xend = end, yend = reads.ratio), col="blue", size=2)
    }
    if(any(CNVs$reads.ratio<1)){
      p <- p + geom_segment(data=subset(CNVs,reads.ratio<1), aes(x = start, y = reads.ratio, xend = end, yend = reads.ratio), col="red", size=2)
    }
  }
  if(missing(geneTable)==FALSE){
    chr <- unique(annot$chromosome)
    genes <- subset(geneTable, chromosome == chr)
    if(nrow(genes)>0){
      exons <- plyr::ldply(genes$Gene,function(x){
        ret <- subset(annot, stringr::str_detect(name, paste0("^",x,"_")))
        ret$gene <- x
        ret
      })
      require(dplyr)
      labels <- exons %>% group_by(gene) %>% summarise(x=min(middle),y=min(ratio))
      p <- p + geom_point(data=exons, mapping=aes(x=middle,y=ratio, shape=gene), colour = "black") + 
        geom_text(data=as.data.frame(labels), aes(x,y, label=gene),angle=45,colour="black")
    }
    
  }
  p <- p + labs(title = paste0("Chr ",unique(annot$chromosome)),x="base position",y="Ratio") + theme(legend.position = "none")
  # print(p)
  
  return(invisible(p))
}