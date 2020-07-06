#!/usr/bin/env Rscript

## requires library(MBASED), installed for personal use with
## module load vital-it/7
## module load R/3.6.1
## in R: BiocManager::install("MBASED")

## call with
## Rscript --vanilla runMBASED.R <annosnpsfile> <indir> <outdir> <outpfx> <nsims> <rho>

## ## [Arguments can be hardcoded below (in- and output files)]
args<-commandArgs(trailingOnly=TRUE)

if(length(args)!=6) {
    stop("Require six arguments: annotatedVars indir outdir sample nsims rho")
}


annosnpsfile<-args[1] ## annotated snps from ANNOVAR
indir<-args[2] ## GeneiASE in/out dir
outdir<-args[3]
outpfx<-args[4]
outpfx<-sub("^/","",outpfx)
nsims<-as.numeric(args[5])
rho<-as.numeric(args[6])



##annosnpsfile<-"$HOME/Style/Pax_genome/annovar/snps_filtered_sel_anno.variant_function"
##indir<-"$HOME/Style/step04_AxEx/Pax_genome/ASE/AxExMS1_1"
##outdir<-"$HOME/Style/step04_AxEx/Pax_genome/ASE/AxExMS1_1"
##outpfx<-"/"
##outpfx<-sub("^/","",outpfx)
##nsims = 10^5
##rho = 0.012

geninfile<-paste0(indir,"/input.tab")
biasfile<-paste0(indir,"/bias.txt")

mboutfile<-paste0(outdir,"/",outpfx,"mbased.out")
mbRData<-paste0(outdir,"/",outpfx,"mbased.RData")


print("=========================================================",quote=FALSE)
print(paste("infile:",geninfile,collapse="\t"),quote=FALSE)
print(paste("nsims:",nsims,collapse="\t"),quote=FALSE)
print(paste("rho:",rho,collapse="\t"),quote=FALSE)
print(paste("outfile:",mboutfile,collapse="\t"),quote=FALSE)
print("=========================================================",quote=FALSE)


## -----------------------------------------------------------------


library(MBASED)



## prepare annotations
## -------------------
snps<-read.table(annosnpsfile,as.is=TRUE)
colnames(snps)<-c("annotation","gene","chr","start","end","ref","alt")

## filter annotations

## table(snps$annotation)
## exclude snps that are intergenic, upstream, downstream, or both
tmp<-grep("(intergenic)|(downstream)|(upstream)",snps$annotation)
## table(snps$annotation[-tmp])
## exclude snps that belong to multiple genes
## or that are assigned to exonic;splicing or UTR5;UTR3
##  (then they can have different genes assigned too)
tmp2<-grep("(,)|(;)",snps$gene)
genes<-snps$gene
genes[c(tmp,tmp2)]<-NA
## remove info in parentheses
filtgene<-gsub("\\(.*\\)","",genes)

## append to snps
snps$filtgene<-filtgene
## original order
snps$id<-1:nrow(snps)
## ---------


## read input tab from geneiase
## ----------------------------
tmp<-read.table(geninfile,header=TRUE,as.is=TRUE)
cts<-tmp[,3:4]
rownames(cts)<-paste0(tmp$gene,".",tmp$snp.id)


## read bias input for geneiase
## ----------------------------
bias<-read.table(biasfile)[1,1]
## set pre-existing allelic bias for reference allele
mu = 1-bias



## prepare data
## ------------
snpid<-as.numeric(matrix(unlist(strsplit(rownames(cts),".",fixed=TRUE)),
                         ncol=2,byrow=TRUE)[,2])
snps.sub<-snps[match(snpid,snps$id),]
if(sum(snps.sub$id==snpid)!=length(snpid)){
    stop("snps in annotations file do not match the ones in input file")
}

mb<-data.frame(gene=snps.sub$filtgene,
               id=snps.sub$id,
               chr=snps.sub$chr,
               start=snps.sub$start,
               end=snps.sub$end,
               refAllele=snps.sub$ref,
               altAllele=snps.sub$alt,
               refCount=cts[,2],
               altCount=cts[,1])

mySNVs <- GRanges(seqnames=mb$chr,
                  ranges=IRanges(start=mb$start, width=1),
                  aseID=as.vector(mb$gene),
                  allele1=mb$refAllele,
                  allele2=mb$altAllele)
names(mySNVs) <- paste(mb$gene,mb$id,sep=".")


mySample <- SummarizedExperiment(
    assays=list(
        lociAllele1Counts=matrix(
            mb$refCount,
            ncol=1,
            dimnames=list(
                names(mySNVs),
                'mySample'
            )
        ),
        lociAllele2Counts=matrix(
            mb$altCount,
            ncol=1,
            dimnames=list(
                names(mySNVs),
                'mySample'
            )
        ),
        lociAllele1CountsNoASEProbs=matrix(
            rep(mu, length(mb$refCount)),
            ncol= 1,
            dimnames= list(names(mySNVs), 'mySample')
        ),
        lociCountsDispersions=matrix(
            rep(rho, length(mb$refCount)),
            ncol= 1,
            dimnames= list(names(mySNVs), 'mySample')
        )
    ),
    rowRanges=mySNVs
)

## run computations
## ----------------
ASEresults <- runMBASED(
    ASESummarizedExperiment=mySample,
    isPhased=TRUE,
    numSim=nsims,
    BPPARAM = SerialParam()
)

## summarize results
## -----------------
summarizeASEResults_1s <- function(MBASEDOutput) {
    geneOutputDF <- data.frame(
        majorAlleleFrequency=assays(MBASEDOutput)$majorAlleleFrequency[,1],
        pValueASE=assays(MBASEDOutput)$pValueASE[,1],
        pValueHeterogeneity=assays(MBASEDOutput)$pValueHeterogeneity[,1]
    )
    lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
    lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
    lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
    lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID)))
    return(
        list(
            geneOutput=geneOutputDF,
            locusOutput=lociOutputList
        )
    )
}

res.mb<-summarizeASEResults_1s(ASEresults)

res.gene<-data.frame(gene=rownames(res.mb$geneOutput),
                     res.mb$geneOutput,stringsAsFactors=FALSE)

## save results
## ------------
write.table(res.gene,file=mboutfile,quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep='\t')

save.image(file=mbRData)

