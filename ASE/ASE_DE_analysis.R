## ---------------------------------------------------------------------
## combined analysis ASE and DE to get first set of good candidate genes
## ---------------------------------------------------------------------

## candidate genes, for example, when the following criteria are met:
## - ASE log2 - fold change >=|1| (med SNP counts DESeq2)
## - significant ASE (geneiase)
## - no significant heterogeneity (MBASED)
## - DE log2 - fold change >=|1| (featureCounts DESeq2)
## - significant DE (DESeq2)
## - relatively highly expressed (adjust for gene length: RPKM)
##      (instead requiring minimum count number)
## - all changes in same direction
## - readCorresp okay
## - ortholog fulfills same criteria

## - when not all criteria are met, genes can nevertheless be interesting
## ---------------------------------------------------------------------

## source("/Users/doro/Dropbox/Style/scripts/ASE_DE_analysis.R")

## =================================================================

setwd("/Users/doro/Documents/Style/analysis/")

mystage<-"MS1" ## do this for each of ("MS1","MS2","LS")
myref<-"Pax" ## genome used for mapping

## =================================================================


library("DESeq2")
library("apeglm") ## for lfcShrink
library("pheatmap") ## for 'pheatmap'
library("MBASED")
library("NOISeq")
library("scales") ## for transparent colors


## =================================================================
## FUNCTIONS
## =================================================================

## convert two-sided p-values to one-sided z-scores
## -> currently returns both a "left" and a "right" value,
##    but could add choice for one as function parameter
twop2onez<-function(pvals,pdirection,replzero=NULL,random=TRUE){

    if(is.vector(pvals)){
        pvals<-as.matrix(pvals,ncol=1)
    }
    if(is.vector(pdirection)){
        pdirection<-as.matrix(pdirection,ncol=1)
    }

    if(sum(dim(pvals)!=dim(pdirection))>0){
        stop("pvals and pdirection need to be of same dimension")
    }

    ## replace zero's in p values
    if(is.null(replzero)){
        replzero<-mean(c(0,min(pvals[!is.na(pvals) & pvals!=0])))
    }
    pvals[!is.na(pvals) & pvals==0]<-replzero

    ## compute one-sided p-values based on pdirection
    effdir<-pdirection
    effdir[!is.na(pdirection) & pdirection<0]<-1
    effdir[!is.na(pdirection) & pdirection>0]<-0
    if(random==TRUE){
        effdir[!is.na(pdirection) & pdirection==0]<-sample(c(0,1),
                                        sum(pdirection==0,na.rm=TRUE),
                                        replace=TRUE)
    }else{
        effdir[!is.na(pdirection) & pdirection==0]<-NA
    }

    p.left<-pvals/2
    p.left[!is.na(effdir) & effdir==1]<-(1-p.left)[!is.na(effdir) & effdir==1]
    p.left[is.na(effdir)]<-NA
    p.right<-pvals/2
    p.right[!is.na(effdir) & effdir==0]<-(1-p.right)[!is.na(effdir) & effdir==0]
    p.right[is.na(effdir)]<-NA

    ## compute z-scores
    z.left<-apply(p.left,2,function(x) qnorm(x))
    z.right<-apply(p.right,2,function(x) qnorm(x))
    rownames(z.left)<-rownames(pvals)
    rownames(z.right)<-rownames(pvals)

    return(list(z.left=z.left,z.right=z.right))
}
##z.nom.med.logfold2<-twop2onez(p.nom,med.logfold2)


## sumz method ("Stouffer’s method", without weights)
## Z_s=(sum{i=1:i=k}(z_i))/sqrt(k); k=nsamples
sumz<-function(zvals){

    if(!is.vector(zvals)){
        stop("zvals needs to be a vector")
    }

    if(sum(!is.na(zvals))==0){
        zs<-NA
    }else{
        zs<-sum(zvals,na.rm=TRUE)/sqrt(sum(!is.na(zvals)))
    }
    return(zs)
}
##sumz(z.left[1,])



## convert one-sided z-scores back to two-sided p-values
onez2twop<-function(zvals){

    if(!is.vector(zvals)){
        stop("zvals needs to be a vector")
    }

    onep<-numeric(length(zvals))
    ## have to change 'lower.tail' when passing 0
    lower<-!is.na(zvals) & zvals<0
    upper<-!is.na(zvals) & zvals>=0
    onep[lower]<-pnorm(zvals[lower],lower.tail=TRUE)
    onep[upper]<-pnorm(zvals[upper],lower.tail=FALSE)
    onep[is.na(zvals)]<-NA
    twop<-onep*2
    return(twop)
}

## =================================================================


samples<-list(MS1=list(
                  sp1=c("AxMS1_1","AxMS1_2","AxMS1_3"),
                  sp2=c("ExMS1_1","ExMS1_2","ExMS1_3"),
                  hyb=c("AxExMS1_1","AxExMS1_2","AxExMS1_3")),
              MS2=list(
                  sp1=c("AxMS2_1","AxMS2_2","AxMS2_3"),
                  sp2=c("ExMS2_1","ExMS2_2","ExMS2_3"),
                  hyb=c("AxExMS2_1","AxExMS2_2","AxExMS2_3")),
              LS=list(
                  sp1=c("AxLS_1","AxLS_2","AxLS_3"),
                  sp2=c("ExLS_2","ExLS_4","ExLS_5"),
                  hyb=c("AxExLS_1","AxExLS_2","AxExLS_3"))
              )

if(myref=="Pax"){
    sp1<-samples[[which(names(samples)==mystage)]]$sp1
    sp2<-samples[[which(names(samples)==mystage)]]$sp2
    contrast<-"P. exserta/P. axillaris"
}else if(myref=="Pex"){
    sp1<-samples[[which(names(samples)==mystage)]]$sp2
    sp2<-samples[[which(names(samples)==mystage)]]$sp1
    contrast<-"P. axillaris/P. exserta"
}else{
    print(paste("unknown reference",myref))
}
hyb<-samples[[which(names(samples)==mystage)]]$hyb




##snps<-read.table(paste0(myref,"/SNPs.txt"),header=TRUE,as.is=TRUE,sep="\t")

## -------------------------------------------------------------------

## ====
## ASE
## ====

myAlpha<-0.05 ## significance level


bias<-numeric(length(hyb))
names(bias)<-hyb

for(i in 1:length(hyb)){
    bias[i]<-read.table(paste0(myref,"/ASE/",hyb[i],"/","bias.txt"))[1,1]
}
## this is average alt/(alt+ref)


## (1) run DESeq2 for each set
## ===========================

## This is only for computing ASE log2 - fold changes, which
## are not output by geneiase nor mbased

## read raw data
## -------------

tmp<-vector("list",length(hyb))

for(i in 1:length(hyb)){
    tmp[[i]]<-read.table(paste0(myref,"/ASE/",hyb[i],"/","input.tab"),
                         header=TRUE,as.is=TRUE)
    colnames(tmp[[i]])[3:4]<-paste0(hyb[i],c(".alt",".ref"))
}
## these should all have the same dimensions
##identical(tmp[[1]]$snp.id,tmp[[2]]$snp.id)

cts<-tmp[[1]][,3:4]
for(i in 2:length(hyb)){
    cts<-cbind(cts,tmp[[i]][,3:4])
}
rownames(cts)<-paste0(tmp[[1]]$gene,".",tmp[[1]]$snp.id)

## sample info
saminfo<-data.frame(sample=rep(hyb,each=2),
                    altref=rep(c("alt","ref"),length(hyb)))
rownames(saminfo)<-colnames(cts)

## make DESeq data set
## -------------------

dds<-DESeqDataSetFromMatrix(countData=cts,
                            colData=saminfo,
                            design= ~ sample + altref)

## can be useful to always have specific factor first
dds$altref<-relevel(dds$altref,"ref") ## consistent with "Ax" (Pax ref)

## ## assign size factors as in Combs & Fraser 2018
## ## <https://doi.org/10.1371/journal.pgen.1007631>
## ## <https://github.com/TheFraserLab/HybridSliceSeq/blob/master/DESeq.R>
## dds<-estimateSizeFactors(dds)
## sf<-sizeFactors(dds)
## alt.sf<-as.numeric(sf[seq(1, length(sf), 2)])
## ref.sf<-as.numeric(sf[seq(2, length(sf), 2)])
## sf[seq(1, length(sf), 2)] <- ref.sf + alt.sf
## sf[seq(2, length(sf), 2)] <- ref.sf + alt.sf
## sizeFactors(dds) <- sf


## assign size factors based on general idea as in Combs & Fraser 2018
## <https://doi.org/10.1371/journal.pgen.1007631>
## <https://github.com/TheFraserLab/HybridSliceSeq/blob/master/DESeq.R>
## >>>>> BUT: ADJUST FOR BIAS <<<<<
dds<-estimateSizeFactors(dds)
sf<-sizeFactors(dds)
alt.sf<-as.numeric(sf[seq(1, length(sf), 2)])
ref.sf<-as.numeric(sf[seq(2, length(sf), 2)])
sf[seq(1, length(sf), 2)] <- (ref.sf + alt.sf) * bias
sf[seq(2, length(sf), 2)] <- (ref.sf + alt.sf) * (1 - bias)
sizeFactors(dds) <- sf



## DE analysis
## -----------

dds <- DESeq(dds)
res <- results(dds, alpha = myAlpha)

## Log fold change shrinkage for visualization and ranking
##resultsNames(dds)
resLFC <- lfcShrink(dds=dds, res=res, coef="altref_alt_vs_ref", type="apeglm")
## other shrinkage estimators available, "apeglm" superior
## (cite 10.1093/bioinformatics/bty895 if used)


## compute median per gene
gene<-as.vector(matrix(unlist(strsplit(rownames(res),".",fixed=TRUE)),
                       ncol=2,byrow=TRUE)[,1])

med.resLFC.log2fold<-ave(resLFC$log2FoldChange,
                         as.factor(gene),FUN=median)
med.resLFC.log2fold<-med.resLFC.log2fold[!duplicated(gene)]
names(med.resLFC.log2fold)<-gene[!duplicated(gene)]

## :plot:
pdf(file=paste0(myref,"/",mystage,"_ASE_DESeq2.pdf"),height=5,width=5)
hist(resLFC$log2FoldChange,breaks=100,main=mystage,
     xlab="Log2_lfc Alt/Ref (ASE)")
abline(v=0,col="red")
dev.off()
## negative values: ref allele -> there is some reference bias
## >>>>> LARGELY GONE WITH TWEAKED SF ESTIMATION
## -----------------------------------------------------------

cat("quantiles med.resLFC.log2fold:\n")
cat(paste(c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1),
          quantile(med.resLFC.log2fold,
                   probs=c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1),
                   na.rm=TRUE),
          sep=":\t",collapse="\n"),"\n")
cat("------\n")



## (2) read geneiase output
## ========================

## This is for getting ASE p-values, combined over all hybrids

gase<-vector("list",length(hyb))

for(i in 1:length(hyb)){
    gase[[i]]<-read.table(paste0(myref,"/ASE/",hyb[i],"/","geneiase.out"),
                          header=TRUE,as.is=TRUE)
    names(gase)[i]<-hyb[i]
}
## these should all have the same dimensions
##identical(gase[[1]]$feat,gase[[2]]$feat)

## !!! adjust order of geneiase outfiles in gase
## ----
for(i in 1:length(hyb)){
    neword<-match(names(med.resLFC.log2fold),gase[[i]]$feat)
    gase[[i]]<-gase[[i]][neword,]
}
##sum(gase[[1]]$feat==names(med.resLFC.log2fold))==length(med.resLFC.log2fold)

## get raw p-values from geneiase
p.nom<-gase[[1]]$p.nom
for(i in 2:length(hyb)){
    p.nom<-cbind(p.nom,gase[[i]]$p.nom)
}
rownames(p.nom)<-gase[[1]]$feat
colnames(p.nom)<-hyb

## !!! replace "0" in p.nom by small value
## (e.g., by half of smallest non-zero value)
## [0 means that p<Nsimulations]
replzero<-mean(c(0,min(p.nom[p.nom!=0])))
p.nom[p.nom==0]<-replzero


## compute one-sided z-scores with med.logfold, computed
## as medium logfold change for each SNP in a gene

## get average logfold change per gene from cts (log2((alt+1)/(ref+1)))
logfold<-log2((cts[,seq(1,ncol(cts),2)] + 1)/(cts[,seq(2,ncol(cts),2)] + 1))
gene<-as.vector(matrix(unlist(strsplit(rownames(logfold),".",fixed=TRUE)),
                       ncol=2,byrow=TRUE)[,1])
med.logfold<-apply(logfold,2,function(x) ave(x,as.factor(gene),FUN=median))
## retain unique genes
med.logfold<-med.logfold[!duplicated(gene),]
rownames(med.logfold)<-gene[!duplicated(gene)]
colnames(med.logfold)<-hyb

## ## remove some of the bias in med.logfold by scaling first
## med.logfold2<-apply(med.logfold,2,function(x) scale(x))
## rownames(med.logfold2)<-rownames(med.logfold)
## z.nom.med.logfold2<-twop2onez(p.nom,med.logfold2)
## z.left<-z.nom.med.logfold2$z.left ## z.left is enough

## don't scale, as results can be affected by outliers
z.nom.med.logfold<-twop2onez(p.nom,med.logfold)
z.left<-z.nom.med.logfold$z.left ## z.left is enough

## combine p-values

## sumz method ("Stouffer’s method", without weights)
sumz.left<-apply(z.left,1,sumz)

## convert one-sided z-scores back to two-sided p-values
combp.left<-onez2twop(sumz.left)
med.combp<-combp.left

## adjust p-values for multiple testing
med.combp.BH<-p.adjust(med.combp,method="BH")


## :plot:
## volcano for each replicate
pdf(file=paste0(myref,"/",mystage,"_ASE_counts.pdf"),height=3.25,width=9)
par(mfrow=c(1,ncol(p.nom)),mar=c(3.2,3.7,1.5,0.5))

for(i in 1:ncol(p.nom)){
    p.BH<-p.adjust(p.nom[,i],method="BH")

    volcols<-rep("gray",length(p.BH))
    volcols[med.logfold[,i]<= -1 &
                p.BH<=myAlpha]<-"dodgerblue3"
    volcols[med.logfold[,i]>=1 &
                p.BH<=myAlpha]<-"firebrick3"

    plot(med.logfold[,i],-log10(p.BH),cex=0.6,pch=16,
         col=alpha(volcols,0.5),axes=FALSE,xlab="",ylab="")
    abline(v=0,col="darkgray",lty=2)
    mtext(colnames(p.nom)[i],font=2)
    axis(1,hadj=0.5,padj=-0.6,cex.axis=1.1,las=1)
    axis(2,padj=0.5,hadj=0.8,cex.axis=1.1,las=2)
    mtext("log2(alt/ref) (ASE)",side=1,line=2.0)
    mtext("-log10(BH p-value) (ASE)",side=2,line=2.3)
}
dev.off()






## (3) read MBASED output
## ======================

mbase<-vector("list",length(hyb))

for(i in 1:length(hyb)){
    mbase[[i]]<-read.table(paste0(myref,"/ASE/",hyb[i],"/","mbased.out"),
                           header=TRUE,as.is=TRUE)
    names(mbase)[i]<-hyb[i]
}
## these should all have the same dimensions
##identical(mbase[[1]]$gene,mbase[[2]]$gene)


## !!! check that mbased and geneiase files are in same format
## ----
sum(gase[[1]]$feat==mbase[[1]]$gene)==nrow(mbase[[1]])


## only use heterogeneity p-values from mbased
## - - - - - - - - - - - - - - - - - - - - - -
p.het<-mbase[[1]]$pValueHeterogeneity
for(i in 2:length(hyb)){
    p.het<-cbind(p.het,mbase[[i]]$pValueHeterogeneity)
}
rownames(p.het)<-mbase[[1]]$gene
colnames(p.het)<-hyb

## !!! replace "0" in p.het by small value
## (e.g., by half of smallest non-zero value)
## [0 means that p<Nsimulations]
## [note that some p.het can be NA]
replzero<-mean(c(0,min(p.het[!is.na(p.het) & p.het!=0])))
p.het[!is.na(p.het) & p.het==0]<-replzero

## assuming that these are one-sided p-values,
##  directly compute z-values
z.het<-apply(p.het,2,function(x) qnorm(x))
rownames(z.het)<-rownames(p.het)

## sumz method ("Stouffer’s method", without weights)
sumz.het<-apply(z.het,1,sumz)

## convert back to one-sided p-values
combp.het<-pnorm(sumz.het,lower.tail=TRUE)

## adjust p-values for multiple testing
combp.het.BH<-p.adjust(combp.het,method="BH")





## ====
## DE
## ====


## (4) re-compute DE only for genes that are in ASE data set
## =========================================================

## featureCounts output
ftcts<-read.table(paste0(myref,"/DE/featureCounts_results.txt"),header=TRUE)
rownames(ftcts)<-ftcts$Geneid
colnames(ftcts)<-gsub("^.*bam\\.","",colnames(ftcts))

anno<-ftcts[,2:6] ## annotation info
ftcts<-ftcts[,-c(1:6)] ## count data

## info for samples
ftsaminfo<-data.frame(condition=gsub("_[0-9]+","",colnames(ftcts)),
                      species=gsub("(LS|MS|SS).*$","",colnames(ftcts)),
                      stage=gsub("^.*(LS|MS[0-9]|SS).*$","\\1",colnames(ftcts)),
                      replicate=gsub("^.*_([0-9]+)$","\\1",colnames(ftcts)),
                      row.names=colnames(ftcts)) ## sample info
##all(colnames(ftcts)==rownames(ftsaminfo))


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


## subsample parentals of correct stage
mysub<-match(c(sp1,sp2),colnames(ftcts))
if(sum(is.na(mysub))>0){
    print("Not all samples have data in featureCounts!")
}
if(sum(!is.na(unique(ftsaminfo$species[mysub])))!=2){
    stop("require two conditions")
}
##ftsaminfo[mysub,]

## subsample genes that have ASE data
mysubg<-match(gase[[1]]$feat,rownames(ftcts))
if(sum(is.na(mysubg))>0){
    print("Not all genes have data in featureCounts!")
}

## make DESeq data set
## -------------------

ftdds<-DESeqDataSetFromMatrix(countData=ftcts[mysubg,mysub],
                              colData=ftsaminfo[mysub,],
                              design= ~ species)

## add meta data
mcols(ftdds)<-DataFrame(mcols(ftdds),anno[mysubg,])

## important to drop factor levels that become irrelevant
ftdds$condition<-droplevels(ftdds$condition)
ftdds$species<-droplevels(ftdds$species)
ftdds$stage<-droplevels(ftdds$stage)
ftdds$replicate<-droplevels(ftdds$replicate)


## can be useful to always have specific factor first
if(myref=="Pax"){
    ftdds$species<-relevel(ftdds$species,"Ax")
}else if(myref=="Pex"){
    ftdds$species<-relevel(ftdds$species,"Ex")
}

## estimate size factors
ftdds<-estimateSizeFactors(ftdds)

## DE analysis
## -----------

ftdds <- DESeq(ftdds)
ftres <- results(ftdds, alpha = myAlpha)

## Log fold change shrinkage for visualization and ranking
##resultsNames(ftdds)
if(myref=="Pax"){
    ftresLFC <- lfcShrink(dds=ftdds, res=ftres,
                          coef="species_Ex_vs_Ax", type="apeglm")
}else if(myref=="Pex"){
    ftresLFC <- lfcShrink(dds=ftdds, res=ftres,
                          coef="species_Ax_vs_Ex", type="apeglm")
}
## other shrinkage estimators available, "apeglm" superior
## (cite 10.1093/bioinformatics/bty895 if used)

##plot(ftresLFC$log2FoldChange,med.resLFC.log2fold,cex=0.5)
cat("quantiles ftresLFC$log2FoldChange:\n")
cat(paste(c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1),
          quantile(ftresLFC$log2FoldChange,
                   probs=c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1),
                   na.rm=TRUE),
          sep=":\t",collapse="\n"),"\n")
cat("------\n")


## (5) Gene length correction of counts with NOISeq package
## ========================================================

## e.g., compute CPM (counts per million)
##       or RPKM (reads per kilobase of exon model per million mapped reads)
##       or FPKM (fragments per kilobase of exon model per million mapped reads)
##       (RPKM and FPKM are equivalent for SE reads)

## make NOISeq data set
## -------------------

mylength<-anno$gene.length
names(mylength)<-rownames(anno)

mychroms<-data.frame(Chr=anno$chr,
                     GeneStart=anno$gene.start,
                     GeneEnd=anno$gene.end,
                     stringsAsFactors=FALSE)
rownames(mychroms)<-rownames(anno)

mydata<-readData(data=ftcts[mysubg,mysub],length=mylength[mysubg],
                 chromosome=mychroms[mysubg,],
                 factors=droplevels(ftsaminfo[mysub,]))

## Normalization
## -------------

## correction for length and sequencing depth
myRPKM<-rpkm(assayData(mydata)$exprs,long=mylength[mysubg],k=0,lc=1)
## if long=1000 in rpkm function, CPM values (counts per million)
##  are returned (which only correct for sequencing depth, not length)

## !!! RPKM values should never be used for DE analysis or other
##     between-sample comparisons!

## compute median over samples
medRPKM<-cbind(
    apply(myRPKM[,1:length(sp1)],1,function(x) median(x)),
    apply(myRPKM[,-c(1:length(sp1))],1,function(x) median(x)))
colnames(medRPKM)<-c("sp1","sp2")

cat("quantiles medRPKM:\n")
cat(paste(c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1),
          quantile(medRPKM,probs=c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)),
          sep=":\t",collapse="\n"),"\n")
cat("------\n")
##sum(medRPKM[,1]<1 | medRPKM[,2]<1)



## =============================
## ALL CHANGES IN SAME DIRECTION
## =============================

## !!! Note: using med.logfold here, not med.logfold2, as the
##           former is more consistent with med.resLFC.log2fold
## - med.logfold ## (for each hybrid)
## - med.resLFC.log2fold ## (combined for all hybrids)
## - ftresLFC$log2FoldChange ## (parental species)

identical(rownames(med.logfold),names(med.resLFC.log2fold))
identical(rownames(med.logfold),names(ftresLFC$log2FoldChange))

samedir<-apply(cbind(med.logfold,
                     med.resLFC.log2fold,
                     ftresLFC$log2FoldChange),
               1,function(x) sum(x<0)==5 | sum(x>0)==5)

##sum(is.na(samedir)) ## (a couple in ftresLFC$log2FoldChange)

## head(cbind(med.logfold,
##            med.resLFC.log2fold,
##            ftresLFC$log2FoldChange)[is.na(samedir),],
##      n=100L)



## ========================
## READ CORRESPONDENCE OKAY
## ========================

corres<-vector("list",length(hyb))

for(i in 1:length(hyb)){
    corres[[i]]<-read.table(paste0(myref,"/ASE/",hyb[i],"/","mapReport.txt"),
                            header=TRUE,as.is=TRUE)
    names(corres)[i]<-hyb[i]
    corres[[i]]$gene<-gsub("\\.[[:digit:]]+","",corres[[i]]$gene)
}
## these should all have the same dimensions
##identical(corres[[1]]$gene,corres[[2]]$gene)

## !!! adjust order of mapReport outfiles in corres
## ----
for(i in 1:length(hyb)){
    neword<-match(names(med.resLFC.log2fold),corres[[i]]$gene)
    corres[[i]]<-corres[[i]][neword,]
}
sum(gase[[1]]$feat==corres[[1]]$gene)==nrow(corres[[1]])

## strict
goodmap<-matrix(1,nrow=nrow(corres[[1]]),ncol=length(hyb))
for(i in 1:length(hyb)){
    goodmap[,i]<-corres[[i]]$p.nonmap<=0.05 &
        corres[[i]]$p.unimap>=0.95 &
            !is.na(corres[[i]]$p.major) &
                corres[[i]]$p.major>=0.95 &
                    !is.na(corres[[i]]$p.major.bin) &
                        corres[[i]]$p.major.bin>=0.95
}
## loose
okaymap<-matrix(1,nrow=nrow(corres[[1]]),ncol=length(hyb))
for(i in 1:length(hyb)){
    okaymap[,i]<-corres[[i]]$p.nonmap<=0.1 &
        corres[[i]]$p.unimap>=0.9 &
            !is.na(corres[[i]]$p.major) &
                corres[[i]]$p.major>=0.9 &
                    !is.na(corres[[i]]$p.major.bin) &
                        corres[[i]]$p.major.bin>=0.9
}
## combine (allow one sample to fail goodmap when it still is
##          part of okaymap)

passmap<-rowSums(goodmap)>=length(hyb)-1 & rowSums(okaymap)==length(hyb)

##head(corres[[1]][!passmap,])

## define breaks for histograms
de.breaks<-seq(floor(min(ftresLFC$log2FoldChange,na.rm=TRUE)),
               ceiling(max(ftresLFC$log2FoldChange,na.rm=TRUE)),
               0.1)
ase.breaks<-seq(floor(min(med.resLFC.log2fold,na.rm=TRUE)),
                ceiling(max(med.resLFC.log2fold,na.rm=TRUE)),
                0.1)


pdf(file=paste0(myref,"/",mystage,"_readcorresp.pdf"),height=5,width=5)
## :plot:
hist(ftresLFC$log2FoldChange,breaks=de.breaks,
     xlab=paste("Log2_lfc",contrast,"(DE)"),main="")
hist(ftresLFC$log2FoldChange[!passmap],breaks=de.breaks,
     add=TRUE,col="red")

## :plot:
hist(med.resLFC.log2fold,breaks=ase.breaks,
     xlab="Log2_lfc Alt/Ref (ASE)",main="")
hist(med.resLFC.log2fold[!passmap],breaks=ase.breaks,
     add=TRUE,col="red")
dev.off()


## ===============
## COMBINE RESULTS
## ===============

##identical(gase[[1]]$n.vars,gase[[2]]$n.vars)

allres<-data.frame(gene=gase[[1]]$feat,
                   chr=anno$chr[mysubg],
                   start=anno$gene.start[mysubg],
                   end=anno$gene.end[mysubg],
                   strand=anno$gene.strand[mysubg],
                   length=anno$gene.length[mysubg],
                   nvars=gase[[1]]$n.vars,
                   log2ase=med.resLFC.log2fold,
                   BHase=med.combp.BH,
                   BHhet=combp.het.BH,
                   log2de=ftresLFC$log2FoldChange,
                   BHde=ftresLFC$padj,
                   rpkm1=medRPKM[,1],
                   rpkm2=medRPKM[,2],
                   samedirection=samedir,
                   passmapping=passmap,
                   stringsAsFactors=FALSE)

write.table(allres,file=paste0(myref,"/",mystage,"_results.csv"),
            col.names=TRUE,row.names=FALSE,quote=FALSE,sep=",")

## ---------------------------------------------------------------------

