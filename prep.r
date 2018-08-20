# install the packages
source("http://bioconductor.org/biocLite.R")
biocLite(c("Rsubread", "edgeR", "limma", "DESeq2"))

library(Rsubread)
library(edgeR)
library(limma)
library(DESeq2)

args = commandArgs(trailingOnly=TRUE) 

bamfiles <- list.files(pattern="*.bam")
annfile <- args[1]
fc <- featureCounts(files=bamfiles, annot.ext=annfile, isGTFAnnotationFile=TRUE,
    GTF.featureType="gene",GTF.attrType="gene_id", nthreads=4)
write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],fc$counts,stringsAsFactors=FALSE),
    file="subread.counts.txt",quote=FALSE,sep="\t",row.names=FALSE)

y <- DGEList(counts=fc$counts)
save(y, file="EdgeR.RData")
