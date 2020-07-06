## --------------------------------------
## quick DE analysis, including all genes
## --------------------------------------

myref<-"Pax"

if(myref=="Pax"){
    setwd("/Users/doro/Documents/Style/analysis/Pax/DE/")
}else if(myref=="Pex"){
    setwd("/Users/doro/Documents/Style/analysis/Pex/DE/")
}else{
    stop("reference must be 'Pax' or 'Pex'")
}

library("DESeq2")
library("apeglm") ## for lfcShrink

## featureCounts output
if(myref=="Pax"){
    countsfile<-"/Users/doro/Documents/Style/analysis/Pax/DE/featureCounts_results.txt"
}else{
    countsfile<-"/Users/doro/Documents/Style/analysis/Pex/DE/featureCounts_results.txt"
}


## read raw data
## -------------

cts<-read.table(countsfile,header=TRUE)
rownames(cts)<-cts$Geneid
colnames(cts)<-gsub("^.*bam\\.","",colnames(cts))

anno<-cts[,2:6] ## annotation info
cts<-cts[,-c(1:6)] ## count data

saminfo<-data.frame(condition=gsub("_[0-9]+","",colnames(cts)),
                    species=gsub("(LS|MS|SS).*$","",colnames(cts)),
                    stage=gsub("^.*(LS|MS[0-9]|SS).*$","\\1",colnames(cts)),
                    replicate=gsub("^.*_([0-9]+)$","\\1",colnames(cts)),
                    row.names=colnames(cts)) ## sample info

all(colnames(cts)==rownames(saminfo))


## clean-up annotations
## chr
tmp<-apply(anno[,1,drop=FALSE],1,function(x) length(unique(as.vector(unlist(strsplit(x,";")))))==1)
if(sum(tmp==TRUE)==nrow(anno)){
    chr<-apply(anno[,1,drop=FALSE],1,function(x) unique(as.vector(unlist(strsplit(x,";")))))
    chr<-gsub("Peax|Peex302","",chr)
}else{
    print("Not all chromosome names are equal for meta-feature")
}
## start
start<-apply(anno[,2,drop=FALSE],1,function(x) as.vector(unlist(strsplit(x,";")))[1])
## end
end<-apply(anno[,3,drop=FALSE],1,function(x) tail(as.vector(unlist(strsplit(x,";"))),1L))
## strand
tmp<-apply(anno[,4,drop=FALSE],1,function(x) length(unique(as.vector(unlist(strsplit(x,";")))))==1)
if(sum(tmp==TRUE)==nrow(anno)){
    strand<-apply(anno[,4,drop=FALSE],1,function(x) unique(as.vector(unlist(strsplit(x,";")))))
}else{
    print("Not all strands are equal for meta-feature")
}
## length
length<-anno[,5]

anno<-data.frame(chr=chr,gene.start=as.integer(start),
                 gene.end=as.integer(end),gene.strand=strand,
                 gene.length=as.integer(length),
                 stringsAsFactors=FALSE)
rm(tmp,chr,start,end,strand,length)


## ---------------------
## DE analysis per stage
## ---------------------


## make DESeq data set
## -------------------

## select three out of the AxLS, but AxLS_5, AxLS_6
##sample(paste0("AxLS_",1:4),3) ## "AxLS_3" "AxLS_2" "AxLS_1"
## also exclude ExLS_1, ExLS_3, ExLS_6,

allstages<-c("SS","MS1","MS2","LS")

for(mystage in allstages){


    myset<-c("Ax","Ex")
    mysub<-grepl(paste0("^(",myset[1],"|",myset[2],")$"),saminfo$species) &
        grepl(paste0("^",mystage,"$"),saminfo$stage) &
            !is.element(rownames(saminfo),paste0("ExLS_",c(1,3,6))) &
                !is.element(rownames(saminfo),paste0("AxLS_",c(4,5,6)))
    if(sum(!is.na(unique(saminfo$species[mysub])))!=2){
        stop("require two conditions")
    }

    ## sum(colnames(cts[,mysub])==rownames(saminfo[mysub,]))==sum(mysub)

    dds<-DESeqDataSetFromMatrix(countData=cts[,mysub],
                                colData=saminfo[mysub,],
                                design= ~ species)


    ## add meta data

    mcols(dds)<-DataFrame(mcols(dds),anno)

    ## important to drop factor levels that become irrelevant
    dds$condition<-droplevels(dds$condition)
    dds$species<-droplevels(dds$species)
    dds$stage<-droplevels(dds$stage)
    dds$replicate<-droplevels(dds$replicate)

    ## can be useful to always have specific factor first
    if(myref=="Pax"){
        dds$species<-relevel(dds$species,"Ax")
    }else{
        dds$species<-relevel(dds$species,"Ex")
    }

    ## estimate size factors
    dds<-estimateSizeFactors(dds)

    ## DE analysis
    ## -----------

    myAlpha<-0.05

    dds <- DESeq(dds)
    res <- results(dds, alpha = myAlpha, pAdjustMethod="BH")
    ## (pAdjustMethod="BH" seems to be default)

    ## NOTE: alpha needs to be set to the FDR value used later on
    ##       (i.e., res$padj < alpha);
    ##       this is important for independent filtering (as described
    ##       on pp. 22-23 of http://www.bioconductor.org/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf

    ## DESeq function should be applied to raw counts

    ## Log fold change shrinkage for visualization and ranking
    ## -------------------------------------------------------
    if(myref=="Pax"){
        resLFC <- lfcShrink(dds=dds, res=res, coef="species_Ex_vs_Ax",
                            type="apeglm")
    }else{
        resLFC <- lfcShrink(dds=dds, res=res, coef="species_Ax_vs_Ex",
                            type="apeglm")
    }
    ## other shrinkage estimators available, "apeglm" superior
    ## (cite 10.1093/bioinformatics/bty895 if used)


    ## save results
    ## ------------
    tmp<-data.frame(geneid=rownames(res),anno[,1:4],as.data.frame(res))
    tmp$LFC.log2FoldChange<-resLFC$log2FoldChange
    tmp$LFC.lfcSE<-resLFC$lfcSE
    write.csv(tmp,row.names=FALSE,file=paste0("results_",mystage,".csv"),
              quote=FALSE)

    rm(tmp,dds,res,resLFC)

}
