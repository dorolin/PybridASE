## ---------------------------------------
## Correlation gene expression - phenotype
## ---------------------------------------

setwd("/Users/doro/Documents/Style/analysis/Pax/correlation/")

library("DESeq2")
library("apeglm") ## for lfcShrink
library("scales") ## for transparent colors

## table with trait measurements (Yarahmadov et al.)
## (rows: average values per stage x species or F1; columns: traits)
pheno<-read.csv("/Users/doro/Documents/Style/data/Growth_Sarah/2020-03-25_Supplemental.csv")

## previously identified pollen-specific genes
pollenSpec<-read.table("/Users/doro/Documents/Style/analysis/Pax/correlation/PollenGenes_Pax.txt", as.is=TRUE)[,1]

## featureCounts output
countsfile<-"/Users/doro/Documents/Style/analysis/Pax/DE/featureCounts_results.txt"


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


## make DESeq data set
## -------------------

## select three out of the AxLS, but AxLS_5, AxLS_6
##sample(paste0("AxLS_",1:4),3) ## "AxLS_3" "AxLS_2" "AxLS_1"
## also exclude ExLS_1, ExLS_3, ExLS_6,

myset<-c("Ax","AxEx","Ex")
mysub<-grepl(paste0("^(",myset[1],"|",myset[2],"|",myset[3],")$"),
             saminfo$species) &
    !is.element(rownames(saminfo),paste0("ExLS_",c(1,3,6))) &
        !is.element(rownames(saminfo),paste0("AxLS_",c(4,5,6)))

## sum(colnames(cts[,mysub])==rownames(saminfo[mysub,]))==sum(mysub)
## saminfo[mysub,]

dds<-DESeqDataSetFromMatrix(countData=cts[,mysub],
                            colData=saminfo[mysub,],
                            design= ~ stage + species)
## stage:species not possible with missing data for SS in F1s


dds

## add meta data

mcols(dds)<-DataFrame(mcols(dds),anno)

## important to drop factor levels that become irrelevant
dds$condition<-droplevels(dds$condition)
dds$species<-droplevels(dds$species)
dds$stage<-droplevels(dds$stage)
dds$replicate<-droplevels(dds$replicate)

## can be useful to always have specific factor first
dds$species<-relevel(dds$species,"Ax")

## estimate size factors
dds<-estimateSizeFactors(dds)
sizeFactors(dds)

## "rlog" transformation (blind = FALSE should be appropriate)
dds.rlog <- rlog(dds, blind = FALSE)
rlog.norm.counts <- assay(dds.rlog)

## "vst" transformation (blind = FALSE should be appropriate)
dds.vst <- varianceStabilizingTransformation(dds, blind = FALSE)
vst.norm.counts <- assay(dds.vst)


## CORRELATION ANALYSIS
## --------------------

## filter pollen genes
## -------------------

## correlating genes (remove with positive coefficient of >= 0.95)
rlog.pollen<-rlog.norm.counts[match(pollenSpec,rownames(rlog.norm.counts)),]
## -> all high correlation
## cts[match(pollenSpec,rownames(rlog.norm.counts)),mysub]

corPollenGenes<-character()
M<-rlog.norm.counts[apply(rlog.norm.counts,1,function(x) sum(x!=0)>=3),]
## note: more are removed than when using untransformed cts[,mysub]
for(i in 1:nrow(rlog.pollen)){
    tmp<-apply(M,1,function(x) cor(x,rlog.pollen[i,]))
    corPollenGenes<-c(corPollenGenes,names(tmp[tmp>=0.95 & !is.na(tmp)]))
}
rm(M)
corPollenGenes<-unique(corPollenGenes)


P<-rlog.norm.counts[match(corPollenGenes,rownames(rlog.norm.counts)),]
M<-rlog.norm.counts[-match(corPollenGenes,rownames(rlog.norm.counts)),]
colSums(P)/colSums(M)
## -> proportion higher in LS and MS2 stages
## -> exclude genes first and re-compute dds

## --------------------------
## exclude pollen genes first
## --------------------------
## (as for correlation analysis, different stages are compared, where
##  this could matter)
corPollenCTSrow<-match(corPollenGenes,rownames(cts))


dds<-DESeqDataSetFromMatrix(countData=cts[-corPollenCTSrow,mysub],
                            colData=saminfo[mysub,],
                            design= ~ stage + species)
## stage:species not possible with missing data for SS in F1s


dds

## add meta data

mcols(dds)<-DataFrame(mcols(dds),anno[-corPollenCTSrow,])

## important to drop factor levels that become irrelevant
dds$condition<-droplevels(dds$condition)
dds$species<-droplevels(dds$species)
dds$stage<-droplevels(dds$stage)
dds$replicate<-droplevels(dds$replicate)

## can be useful to always have specific factor first
dds$species<-relevel(dds$species,"Ax")

## estimate size factors
dds<-estimateSizeFactors(dds)
sizeFactors(dds)

## "rlog" transformation (blind = FALSE should be appropriate)
dds.rlog <- rlog(dds, blind = FALSE)
rlog.norm.counts <- assay(dds.rlog)

## "vst" transformation (blind = FALSE should be appropriate)
dds.vst <- varianceStabilizingTransformation(dds, blind = FALSE)
vst.norm.counts <- assay(dds.vst)



## compute mean gene expression per condition
## ------------------------------------------

##head(rlog.norm.counts)
##sum(colnames(rlog.norm.counts)==rownames(saminfo[mysub,]))==sum(mysub)

condition<-unique(saminfo$condition[mysub])
rlog.norm.counts.means<-matrix(NA,nrow=nrow(rlog.norm.counts),
                            ncol=length(condition))

for(i in 1:length(condition)){
    set<-which(saminfo$condition[mysub]==condition[i])
    rlog.norm.counts.means[,i]<-rowMeans(rlog.norm.counts[,set,drop=FALSE])
}
colnames(rlog.norm.counts.means)<-condition
rownames(rlog.norm.counts.means)<-rownames(rlog.norm.counts)
rlog.norm.counts.means<-as.data.frame(rlog.norm.counts.means)


## prepare phenotypes
## ------------------

## traits: length, cell wall elasticity, growth rate, cell division rate

pheno$species<-gsub("SS|MS1|MS2|LS","",pheno$speciesStage)
pheno$stage<-gsub("Ax|Ex|Pa|AxPa|ExPa|AxEx","",pheno$speciesStage)
pheno$phase<-factor(pheno$stage,levels=c("SS","MS1","MS2","LS"))
pheno$stage<-as.numeric(pheno$phase)
pheno$cellNumber<-NA
a<-pheno$cellDivisionRate; a[grepl("^(AxPa|ExPa|AxEx)$",pheno$species)&
                                 grepl("^SS$",pheno$phase)]<-NA ## !!!
pheno$cellDivisionRate<-a
pheno.sub<-pheno[match(condition,pheno$speciesStage),]

##sum(colnames(rlog.norm.counts.means)==pheno.sub$speciesStage)==length(condition)


## gene expression cellDivsionRate
## -------------------------------
## need to compute "logfold change gene expression between
##   S and M, M and M2, and M2 and L stages"

## this will not work on rlog-transformed values (as they can be <0)

cdr.stages<-cbind(c("AxExLS","AxExMS1","AxExMS2",
                    "AxLS","AxMS1","AxMS2","AxSS",
                    "ExLS","ExMS1","ExMS2","ExSS"),
                  c(NA,"AxExMS1","AxExMS2",
                    NA,"AxMS1","AxMS2","AxSS",
                    NA,"ExMS1","ExMS2","ExSS"),
                  c(NA,"AxExMS2","AxExLS",
                    NA,"AxMS2","AxLS","AxMS1",
                    NA,"ExMS2","ExLS","ExMS1"))
rownames(cdr.stages)<-cdr.stages[,1]
colnames(cdr.stages)<-c("speciesStage","nominator","denominator")
cdr.stages<-as.data.frame(cdr.stages,stringsAsFactors=FALSE)

## using DESeq to compute logfold change
## - - - - - - - - - - - - - - - - - - -
lfc.stages<-matrix(NA,nrow=nrow(cts[-corPollenCTSrow,]),ncol=nrow(cdr.stages))
colnames(lfc.stages)<-cdr.stages$speciesStage
rownames(lfc.stages)<-rownames(cts[-corPollenCTSrow,])

myAlpha<-0.05 ## significance level

## LS samples not to use
toexcl<-match(c(paste0("ExLS_",c(1,3,6)),paste0("AxLS_",c(4,5,6))),
              rownames(saminfo))

for(i in 1:nrow(cdr.stages)){
    if(is.na(cdr.stages[i,2]) | is.na(cdr.stages[i,3])){
        next
    }
    s1<-c(grep(paste0("^",cdr.stages[i,2],"$"),saminfo$condition),
          grep(paste0("^",cdr.stages[i,3],"$"),saminfo$condition))
    stagesub<-setdiff(s1,toexcl)

    ## make DESeq data set
    stagedds<-DESeqDataSetFromMatrix(countData=cts[-corPollenCTSrow,stagesub],
                                     colData=saminfo[stagesub,],
                                     design= ~ stage)

    ## add meta data
    mcols(stagedds)<-DataFrame(mcols(stagedds),anno[-corPollenCTSrow,])

    ## important to drop factor levels that become irrelevant
    stagedds$condition<-droplevels(stagedds$condition)
    stagedds$species<-droplevels(stagedds$species)
    stagedds$stage<-droplevels(stagedds$stage)
    stagedds$replicate<-droplevels(stagedds$replicate)


    ## can be useful to always have specific factor first
    stagedds$stage<-relevel(stagedds$stage,gsub("Ax|Ex|AxEx","",
                                                cdr.stages[i,3]))

    ## estimate size factors
    stagedds<-estimateSizeFactors(stagedds)

    ## DE analysis
    stagedds <- DESeq(stagedds)
    stageres <- results(stagedds, alpha = myAlpha)

    ## Log fold change shrinkage for visualization and ranking
    ##resultsNames(stagedds)
    stageresLFC <- lfcShrink(dds=stagedds, res=stageres,
                             coef=resultsNames(stagedds)[2], type="apeglm")
    ## other shrinkage estimators available, "apeglm" superior
    ## (cite 10.1093/bioinformatics/bty895 if used)

    ## check that direction of log2FoldChange is correct
    s1<-setdiff(grep(paste0("^",cdr.stages[i,2],"$"),saminfo$condition),toexcl)
    s2<-setdiff(grep(paste0("^",cdr.stages[i,3],"$"),saminfo$condition),toexcl)
    crude.lfc<-log2(rowSums(cts[-corPollenCTSrow,s1])/
                        rowSums(cts[-corPollenCTSrow,s2]))
    crude.lfc[crude.lfc==Inf]<-max(crude.lfc[crude.lfc!=Inf],na.rm=TRUE)+1
    crude.lfc[crude.lfc== -Inf]<-min(crude.lfc[crude.lfc!= -Inf],na.rm=TRUE)-1

    if(cor(crude.lfc,stageresLFC$log2FoldChange,use="complete.obs")>0){
        lfc.stages[,i]<-stageresLFC$log2FoldChange
    }else{
        warning(paste("switching log2FoldChange for",cdr.stages[i,1]))
        lfc.stages[,i]<-stageresLFC$log2FoldChange * -1
    }
}


## ALTERNATIVE / ADJUSTMENT, as NA's in lfc.stages
## - - - - - - - - - - - - - - - - - - - - - - - -
lfc.stages.cts<-matrix(NA,nrow=nrow(cts[-corPollenCTSrow,]),
                       ncol=nrow(cdr.stages))
colnames(lfc.stages.cts)<-cdr.stages$speciesStage
rownames(lfc.stages.cts)<-rownames(cts[-corPollenCTSrow,])

for(i in 1:nrow(cdr.stages)){
    if(is.na(cdr.stages[i,2]) | is.na(cdr.stages[i,3])){
        next
    }
    s1<-grep(paste0("^",cdr.stages[i,2],"$"),saminfo$condition)
    stagesub1<-setdiff(s1,toexcl)
    s2<-grep(paste0("^",cdr.stages[i,3],"$"),saminfo$condition)
    stagesub2<-setdiff(s2,toexcl)

    ## get counts
    stagects1<-cts[-corPollenCTSrow,stagesub1]
    stagects2<-cts[-corPollenCTSrow,stagesub2]

    lfc.stages.cts[,i]<-log2((rowSums(stagects1)+1)/(rowSums(stagects2)+1))
}

## for(i in 1:nrow(cdr.stages)){
##     if(is.na(cdr.stages[i,2]) | is.na(cdr.stages[i,3])){
##         next
##     }
##     cat(paste0(cdr.stages[i,1],": "))
##     cat(paste0(sum(lfc.stages.cts[is.na(lfc.stages[,i]),i]),"\n"))
## }
## ## -> all have indeed log2(counts1/counts2)=0

## replace NA's in lfc.stages by values in lfc.stages.cts
lfc.stages.adj<-lfc.stages
for(i in 1:nrow(cdr.stages)){
    if(is.na(cdr.stages[i,2]) | is.na(cdr.stages[i,3])){
        next
    }
    na<-is.na(lfc.stages[,i])
    lfc.stages.adj[na,i]<-lfc.stages.cts[na,i]
}

##sum(lfc.stages[!is.na(lfc.stages)]==0)


## get number of samples with counts per species and stage
## -------------------------------------------------------
speciesStages.dataPts<-matrix(NA,nrow=nrow(cts[-corPollenCTSrow,]),
                       ncol=nrow(cdr.stages))
colnames(speciesStages.dataPts)<-cdr.stages$speciesStage
rownames(speciesStages.dataPts)<-rownames(cts[-corPollenCTSrow,])

for(i in 1:nrow(cdr.stages)){
    s1<-grep(paste0("^",cdr.stages[i,1],"$"),saminfo$condition)
    stagesub1<-setdiff(s1,toexcl)

    ## get counts
    stagects1<-cts[-corPollenCTSrow,stagesub1]

    speciesStages.dataPts[,i]<-apply(stagects1,1,function(x) sum(x!=0))
}
##sum(colnames(rlog.norm.counts.means)==colnames(speciesStages.dataPts))==nrow(cdr.stages)



## ------------------------------------------------------
## compute trait correlations, original and permuted data
## ------------------------------------------------------

nPerm <- 1000
nHist <- 10
mybreaks<-seq(0,1,0.01)

## growthRate
## ----------
y<-pheno.sub$growthRate
names(y)<-pheno.sub$speciesStage
X<-rlog.norm.counts.means[,!is.na(y),drop=FALSE]
## potentially exclude genes with few data points
dp<-speciesStages.dataPts[,!is.na(y),drop=FALSE]
X<-X[apply(dp,1,function(x) sum(x)>=3),,drop=FALSE]
y<-y[!is.na(y)]
sum(names(y)==colnames(X))==length(y)

## simple correlation function cor (instead of cor.test)
cor.growthRate<-apply(X,1,function(x) cor(x,y,use="complete.obs"))

##quantile(abs(cor.growthRate),probs=c(0.95,0.99))

## calculate quantiles (5 % of absolute value) over genes per permutation
permQuant.growthRate<-matrix(NA,nrow=nPerm,ncol=2)
colnames(permQuant.growthRate)<-c("q0.95abs","q0.99abs")
dens.growthRate<-matrix(NA,nrow=length(mybreaks)-1,ncol=nHist)

set.seed(2048)
for(i in 1:nPerm){
    Xperm<-t(apply(X,1,function(x) sample(x)))
    cor.perm<-apply(Xperm,1,function(x) cor(x,y,use="complete.obs"))
    permQuant.growthRate[i,]<-quantile(abs(cor.perm),probs=c(0.95,0.99))
    if(i<=nHist){
        dens.growthRate[,i]<-hist(abs(cor.perm),breaks=seq(0,1,0.01),
                                  plot=FALSE)$density
    }
}


## length
## ------
y<-pheno.sub$length
names(y)<-pheno.sub$speciesStage
X<-rlog.norm.counts.means[,!is.na(y),drop=FALSE]
## potentially exclude genes with few data points
dp<-speciesStages.dataPts[,!is.na(y),drop=FALSE]
X<-X[apply(dp,1,function(x) sum(x)>=3),,drop=FALSE]
y<-y[!is.na(y)]
sum(names(y)==colnames(X))==length(y)

## simple correlation function cor (instead of cor.test)
cor.length<-apply(X,1,function(x) cor(x,y,use="complete.obs"))

##quantile(abs(cor.length),probs=c(0.95,0.99))

## calculate quantiles (5 % of absolute value) over genes per permutation
permQuant.length<-matrix(NA,nrow=nPerm,ncol=2)
colnames(permQuant.length)<-c("q0.95abs","q0.99abs")
dens.length<-matrix(NA,nrow=length(mybreaks)-1,ncol=nHist)

set.seed(2048)
for(i in 1:nPerm){
    Xperm<-t(apply(X,1,function(x) sample(x)))
    cor.perm<-apply(Xperm,1,function(x) cor(x,y,use="complete.obs"))
    permQuant.length[i,]<-quantile(abs(cor.perm),probs=c(0.95,0.99))
    if(i<=nHist){
        dens.length[,i]<-hist(abs(cor.perm),breaks=seq(0,1,0.01),
                              plot=FALSE)$density
    }
}


## cellWallElasticity
## ------------------
y<-pheno.sub$cellWallElasticity
names(y)<-pheno.sub$speciesStage
X<-rlog.norm.counts.means[,!is.na(y),drop=FALSE]
## potentially exclude genes with few data points
dp<-speciesStages.dataPts[,!is.na(y),drop=FALSE]
X<-X[apply(dp,1,function(x) sum(x)>=3),,drop=FALSE]
y<-y[!is.na(y)]
sum(names(y)==colnames(X))==length(y)

## simple correlation function cor (instead of cor.test)
cor.cellWallElasticity<-apply(X,1,function(x) cor(x,y,use="complete.obs"))

##quantile(abs(cor.cellWallElasticity),probs=c(0.95,0.99))

## calculate quantiles (5 % of absolute value) over genes per permutation
permQuant.cellWallElasticity<-matrix(NA,nrow=nPerm,ncol=2)
colnames(permQuant.cellWallElasticity)<-c("q0.95abs","q0.99abs")
dens.cellWallElasticity<-matrix(NA,nrow=length(mybreaks)-1,ncol=nHist)

set.seed(2048)
for(i in 1:nPerm){
    Xperm<-t(apply(X,1,function(x) sample(x)))
    cor.perm<-apply(Xperm,1,function(x) cor(x,y,use="complete.obs"))
    permQuant.cellWallElasticity[i,]<-quantile(abs(cor.perm),probs=c(0.95,0.99))
    if(i<=nHist){
        dens.cellWallElasticity[,i]<-hist(abs(cor.perm),breaks=seq(0,1,0.01),
                                          plot=FALSE)$density
    }
}


## cellDivisionRate
## ----------------
y<-pheno.sub$cellDivisionRate
names(y)<-pheno.sub$speciesStage
X<-lfc.stages.adj[,!is.na(y),drop=FALSE] ## !!!
## potentially exclude genes with few data points
dp<-speciesStages.dataPts[,!is.na(y),drop=FALSE]
X<-X[apply(dp,1,function(x) sum(x)>=3),,drop=FALSE]
y<-y[!is.na(y)]
sum(names(y)==colnames(X))==length(y)

##table(apply(X,1,function(x) sum(is.na(x))))

## simple correlation function cor (instead of cor.test)
cor.cellDivisionRate.adj<-apply(X,1,function(x) cor(x,y,use="complete.obs"))

##quantile(abs(cor.cellDivisionRate.adj),probs=c(0.95,0.99))

## calculate quantiles (5 % of absolute value) over genes per permutation
permQuant.cellDivisionRate.adj<-matrix(NA,nrow=nPerm,ncol=2)
colnames(permQuant.cellDivisionRate.adj)<-c("q0.95abs","q0.99abs")
dens.cellDivisionRate.adj<-matrix(NA,nrow=length(mybreaks)-1,ncol=nHist)

set.seed(2048)
for(i in 1:nPerm){
    Xperm<-t(apply(X,1,function(x) sample(x)))
    cor.perm<-apply(Xperm,1,function(x) cor(x,y,use="complete.obs"))
    permQuant.cellDivisionRate.adj[i,]<-quantile(abs(cor.perm),
                                                 probs=c(0.95,0.99))
    if(i<=nHist){
        dens.cellDivisionRate.adj[,i]<-hist(abs(cor.perm),breaks=seq(0,1,0.01),
                                            plot=FALSE)$density
    }
}




## ---------------
## combine results
## ---------------

cor.combined<-data.frame(growthRate=rep(NA,nrow(cts)),
                         length=rep(NA,nrow(cts)),
                         cellWallElasticity=rep(NA,nrow(cts)),
                         cellDivisionRate=rep(NA,nrow(cts)),
                         stringsAsFactors=FALSE)
rownames(cor.combined)<-rownames(cts)

cor.combined$growthRate[match(names(cor.growthRate),rownames(cor.combined))]<-cor.growthRate
cor.combined$length[match(names(cor.length),rownames(cor.combined))]<-cor.length
cor.combined$cellWallElasticity[match(names(cor.cellWallElasticity),rownames(cor.combined))]<-cor.cellWallElasticity
cor.combined$cellDivisionRate[match(names(cor.cellDivisionRate.adj),rownames(cor.combined))]<-cor.cellDivisionRate.adj

## coefficient of variation (based on all rlog count data)
cv<-apply(rlog.norm.counts,1,function(x) sd(x)/abs(mean(x)))
cor.combined$cv[match(names(cv),rownames(cor.combined))]<-cv
## mean of rlog count data
rlog.cts<-rowMeans(rlog.norm.counts)
cor.combined$rlog.cts[match(names(rlog.cts),rownames(cor.combined))]<-rlog.cts



## -------------------------------------------------------------------
## SAVINGS
## -------------------------------------------------------------------
save.image("correlation_analysis.RData")

tmp<-cbind(rownames(cor.combined),cor.combined)
colnames(tmp)<-c("geneid",colnames(cor.combined)[1:4],
                 "cv.rlog.cts","mean.rlog.cts")
write.table(tmp,file="cor_results.txt",col.names=TRUE,
            row.names=FALSE,quote=FALSE)


permQuant0.95Means<-c(mean(permQuant.growthRate[,1]),
                      mean(permQuant.length[,1]),
                      mean(permQuant.cellWallElasticity[,1]),
                      mean(permQuant.cellDivisionRate.adj[,1]))

permQuant0.99Means<-c(mean(permQuant.growthRate[,2]),
                      mean(permQuant.length[,2]),
                      mean(permQuant.cellWallElasticity[,2]),
                      mean(permQuant.cellDivisionRate.adj[,2]))

permQuant<-cbind(c("growthRate","length","cellWallElasticity",
                   "cellDivisionRate"),permQuant0.95Means,
                 permQuant0.99Means)
colnames(permQuant)[1]<-"trait"
permQuant<-as.data.frame(permQuant,stringsAsFactors=FALSE)
write.table(permQuant,file="cor_permQuantiles.txt",col.names=TRUE,
            row.names=FALSE,quote=FALSE)
## -------------------------------------------------------------------



## plot histograms with real data and permutation data (both thresholds)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## re-format permQuant
permQuant<-as.data.frame(t(read.table(file="cor_permQuantiles.txt",
                                      header=TRUE,row.names=1)))

## mybreaks<-seq(0,1,0.01)
mids<-rowMeans(cbind(mybreaks[-length(mybreaks)],mybreaks[-1]))

pdf("correlation_histograms.pdf",width=9,height=9,pointsize=14)

par(mfrow=c(2,2),mar=c(4,4,2,1))

## growthRate
m<-rowMeans(dens.growthRate)
sd<-apply(dens.growthRate,1,sd)
hist(abs(cor.combined$growthRate),breaks=seq(0,1,0.01),
     freq=FALSE,col=alpha("firebrick3",0.4),
     xlim=c(0,1),ylim=c(0,3),
     main="growthRate",xlab="cor(rlog.counts, growthRate)")
barplot(height=m,width=1/(length(mybreaks)-1),
        space=0,col=alpha("yellow",0.4),add=TRUE)
segments(x0=mids,x1=mids,y0=m-sd,y1=m+sd)
abline(v=permQuant$growthRate,lwd=1,lty=2,col="goldenrod3")
text(x=permQuant$growthRate,y=c(2.5,2.5),
     labels=c("0.95 quantile","0.99 quantile"),pos=2,
     srt=90,col="goldenrod3")

## length
m<-rowMeans(dens.length)
sd<-apply(dens.length,1,sd)
hist(abs(cor.combined$length),breaks=seq(0,1,0.01),
     freq=FALSE,col=alpha("firebrick3",0.4),
     xlim=c(0,1),ylim=c(0,3),
     main="length",xlab="cor(rlog.counts, length)")
barplot(height=m,width=1/(length(mybreaks)-1),
        space=0,col=alpha("yellow",0.4),add=TRUE)
segments(x0=mids,x1=mids,y0=m-sd,y1=m+sd)
abline(v=permQuant$length,lwd=1,lty=2,col="goldenrod3")
text(x=permQuant$length,y=c(2.5,2.5),
     labels=c("0.95 quantile","0.99 quantile"),pos=2,
     srt=90,col="goldenrod3")

## cellWallElasticity
m<-rowMeans(dens.cellWallElasticity)
sd<-apply(dens.cellWallElasticity,1,sd)
hist(abs(cor.combined$cellWallElasticity),breaks=seq(0,1,0.01),
     freq=FALSE,col=alpha("firebrick3",0.4),
     xlim=c(0,1),ylim=c(0,3),
     main="cellWallElasticity",xlab="cor(rlog.counts, cellWallElasticity)")
barplot(height=m,width=1/(length(mybreaks)-1),
        space=0,col=alpha("yellow",0.4),add=TRUE)
segments(x0=mids,x1=mids,y0=m-sd,y1=m+sd)
abline(v=permQuant$cellWallElasticity,lwd=1,lty=2,col="goldenrod3")
text(x=permQuant$cellWallElasticity,y=c(2.5,2.5),
     labels=c("0.95 quantile","0.99 quantile"),pos=2,
     srt=90,col="goldenrod3")

## cellDivisionRate
m<-rowMeans(dens.cellDivisionRate.adj)
sd<-apply(dens.cellDivisionRate.adj,1,sd)
hist(abs(cor.combined$cellDivisionRate),breaks=seq(0,1,0.01),
     freq=FALSE,col=alpha("firebrick3",0.4),
     xlim=c(0,1),ylim=c(0,3),
     main="cellDivisionRate",xlab="cor(rlog.counts, cellDivisionRate)")
barplot(height=m,width=1/(length(mybreaks)-1),
        space=0,col=alpha("yellow",0.4),add=TRUE)
segments(x0=mids,x1=mids,y0=m-sd,y1=m+sd)
abline(v=permQuant$cellDivisionRate,lwd=1,lty=2,col="goldenrod3")
text(x=permQuant$cellDivisionRate,y=c(2.5,2.5),
     labels=c("0.95 quantile","0.99 quantile"),pos=2,
     srt=90,col="goldenrod3")

dev.off()



## ## test effect of coefficient of variation on correlation
## ## - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## cv.breaks<-seq(min(log(cor.combined$cv),na.rm=TRUE),
##                max(log(cor.combined$cv),na.rm=TRUE),
##                length.out=11)
## cv.breaks[11]<-cv.breaks[11]*1.001
## cv.bins<-factor(nrow(cor.combined),levels=1:11)
## for(i in 1:10){
##     cv.bins[!is.na(log(cor.combined$cv)) &
##                 log(cor.combined$cv)>=cv.breaks[i] &
##                     log(cor.combined$cv)<cv.breaks[i+1]]<-i
## }
## cv.bins[is.na(log(cor.combined$cv))]<-11
## ##table(cv.bins)
## ## combine bins 7:10
## cv.bins[as.numeric(cv.bins)>=7 & as.numeric(cv.bins)<=10]<-7
## cv.bins<-droplevels(cv.bins)
## levels(cv.bins)[levels(cv.bins)==11]<-"NA"

## par(mfrow=c(2,2),mar=c(3,3,1,1))
## boxplot(abs(cor.combined$growthRate) ~ cv.bins,ylab="correlation")
## mtext("growthRate")
## boxplot(abs(cor.combined$length) ~ cv.bins,ylab="correlation")
## mtext("length")
## boxplot(abs(cor.combined$cellWallElasticity) ~ cv.bins,ylab="correlation")
## mtext("cellWallElasticity")
## boxplot(abs(cor.combined$cellDivisionRate) ~ cv.bins,ylab="correlation")
## mtext("cellDivisionRate")
## ## -> slight trend: intermediate cv -> higher correlation


## ## test effect of mean rlog counts on correlation
## ## - - - - - - - - - - - - - - - - - - - - - - - -
## cts.breaks<-seq(min(cor.combined$rlog.cts,na.rm=TRUE),
##                 max(cor.combined$rlog.cts,na.rm=TRUE),
##                 length.out=11)
## cts.breaks[11]<-cts.breaks[11]*1.001
## cts.bins<-factor(nrow(cor.combined),levels=1:11)
## for(i in 1:10){
##     cts.bins[!is.na(cor.combined$rlog.cts) &
##                  cor.combined$rlog.cts>=cts.breaks[i] &
##                      cor.combined$rlog.cts<cts.breaks[i+1]]<-i
## }
## cts.bins[is.na(cor.combined$rlog.cts)]<-11
## ##table(cts.bins)
## ## combine bins 7:10
## cts.bins[as.numeric(cts.bins)>=7 & as.numeric(cts.bins)<=10]<-7
## cts.bins<-droplevels(cts.bins)
## levels(cts.bins)[levels(cts.bins)==11]<-"NA"

## par(mfrow=c(2,2),mar=c(3,3,1,1))
## boxplot(abs(cor.combined$growthRate) ~ cts.bins,ylab="correlation")
## mtext("growthRate")
## boxplot(abs(cor.combined$length) ~ cts.bins,ylab="correlation")
## mtext("length")
## boxplot(abs(cor.combined$cellWallElasticity) ~ cts.bins,ylab="correlation")
## mtext("cellWallElasticity")
## boxplot(abs(cor.combined$cellDivisionRate) ~ cts.bins,ylab="correlation")
## mtext("cellDivisionRate")
## ## -> slight trend: more counts -> higher correlation


## plot(cv.bins,cor.combined$rlog.cts,xlab="cv bin",ylab="rlog.cts")
## plot(cts.bins,log(cor.combined$cv),xlab="rlog.cts bin",ylab="log(cv)")
## ## -> shows that slight trend: intermediate cv -> higher correlation
## ##    can be explained by intermediate cv = more counts


