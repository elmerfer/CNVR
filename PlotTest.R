
geneTable <- openxlsx::read.xlsx("./data/GenesTargetCNVs.xlsx")
geneTable$chromosome <- stringr::str_remove_all(geneTable$chromosome,"[a-z]")
View(geneTable)

pls <- bplapply(1:22, function(xc){
  print(xc)
   
  obj <- GetCNVsAnnotation(chr = xc  , x=sbj[[1]])
  # attr(obj$CNVs,"BamHeader")
  pp<-CNVR:::.PlotCNVchromosome(obj, geneTable = geneTable)
} , BPPARAM = bpparam())


library(ggpubr)

#pdf(file="27040.pdf",paper="a4r",height = 14,width = 14)
ggarrange(plotlist = pls, nrow=4, ncol=6)
# dev.off()

obj <- GetCNVsAnnotation(chr = 10  , x=CallsBP[[1]])
pp<-CNVR:::.PlotCNVchromosome(obj, geneTable = geneTable)
pp



foo <- function(genecode){
  unlist(lapply(genecode, function(gline){
    df <- ldply(unlist(stringr::str_split(gline,",")),function(x) unlist(stringr::str_split(x,"_")))
    df <- df %>% group_by(V1) %>% summarise(paste0(V2,collapse=";"))
    df<-as.data.frame(df)
    paste0(unlist(lapply(1:nrow(df), function(x){
      
      paste0(df[x,1],"(",paste0(df[x,2],collapse = ";"),")")
      
    })),collapse = "-")    
  }))

}

CallsBP[[2]]@CNV.calls$GenesInCNVRegion <- foo(CallsBP[[2]]@CNV.calls$genecode)
CHRi$CNVs$nl <- foo(CHRi$CNVs$genecode)
head(CHRi$CNVs$nl)
openxlsx::write.xlsx(CallsBP[[2]]@CNV.calls[,-13], "Test_38148.xlsx", overwrite = T)
print(pls[[1]])
CHRi <- GetCNVsAnnotation(chr=9 , x=CallsBP[[2]])
pl<-CNVR:::.PlotCNVchromosome(CHRi, geneTable = geneTable)
print(.PlotC(annot=CHRi,thLength=500, geneTable=geneTable))
  
  
chr <- unique(CHRi$Exons$chromosome)
genes <- subset(geneTable, chromosome == chr)
exons <- plyr::ldply(genes$Gene,function(x){
  ret <- subset(CHRi$Exons, stringr::str_detect(name, paste0("^",x,"_")))
  ret$gene <- x
  ret
})
p + geom_point(data=exons, mapping=aes(x=middle,y=ratio, shape=gene), colour = "black") + 
  geom_text(data=as.data.frame(labels), aes(x,y, label=gene),angle=45,colour="black")

ldply(exons$gene, function(x){
  dplyr::summarise(.data=exons,x=unique(gene),y=min(ratio),groups_by = gene)
})
library(dplyr)
labels <- exons %>% group_by(gene) %>% summarise(x=min(middle),y=min(ratio))
colnames(labels)
