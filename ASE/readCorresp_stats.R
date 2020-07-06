#!/usr/bin/env Rscript

## note: the second part of the script is extremely slow and
##       there is definitely some room for optimization

## call with
## Rscript --vanilla readCorresp_stats.R <refgenefile> <correspfile> <indir> <outdir> <outpfx>

## ## [Arguments can be hardcoded below (in- and output files)]
args<-commandArgs(trailingOnly=TRUE)

if(length(args)!=5) {
    stop("Require five arguments: annotatedVars readCorresp indir outdir sample")
}


refgenefile<-args[1] ## annotated variants from ANNOVAR
correspfile<-args[2] ## read correspondence
indir<-args[3] ## GeneiASE in/out dir
outdir<-args[4]
outpfx<-args[5]
outpfx<-sub("^/","",outpfx)

genoutfile<-paste0(indir,"/geneiase.out")
geninfile<-paste0(indir,"/input.tab")

combaseoutfile<-paste0(outdir,"/",outpfx,"ASEsum.txt")
mapreportfile<-paste0(outdir,"/",outpfx,"mapReport.txt")

## ## -----------------------------------------------------------------
## ## require files:
## ## --------------
## ## ## annotated variants from ANNOVAR ($HOME/Style/Pax_genome/annovar/)
## refgenefile<-"$HOME/Style/Pax_genome/annovar/Peax304db/Peax304_refGene.txt"
## ## ## see 'ASE_analysis.R' on how to generate this file

## ## ## read correspondence (huge file)
## correspfile<-"$HOME/Style/step02_AxEx/corresp/AxExMS1_2_aligninfo_Pax.txt"
## ## ## see 'readCorresp.sh' and 'readCorresp.pl' on how to generate this file

## ## ## GeneiASE in/out
## genoutfile<-"$HOME/Style/step04_AxEx/Pax_genome/ASE/AxExMS1_2/geneiase.out"
## geninfile<-"$HOME/Style/step04_AxEx/Pax_genome/ASE/AxExMS1_2/input.tab"

## ## output files:
## ## -------------
## combaseoutfile<-"$HOME/Style/step04_AxEx/Pax_genome/ASE/AxExMS1_2/ASEsum.txt"
## mapreportfile<-"$HOME/Style/step04_AxEx/Pax_genome/ASE/AxExMS1_2/mapReport.txt"
## ## -----------------------------------------------------------------

## define when a read counts as overlapping a gene
## (approximation, as deviations in aligned read positions from
##  trimmed reads or gaps in alignment are ignored)
readlen<-101 ## read length
overl<-1 ## number of bases overlapping gene



## combine ASE stats and gene position info
## ----------------------------------------

## read annotated variants
refgene<-read.table(refgenefile,as.is=TRUE)
colnames(refgene)<-c("gene","chr","strand","start","end",
                     "cds.start","cds.end","n.cds","exon.starts","exon.ends")
refgene<-refgene[,1:10]
## compute mean gene positions
genpos<-data.frame(feat=gsub("\\.[0-9]","",refgene$gene),
                   chr=refgene$chr,
                   position=rowMeans(refgene[,4:5]),
                   stringsAsFactors=FALSE)

## GeneiASE output
genout<-read.table(genoutfile,header=TRUE,as.is=1)
## GeneiASE input
genin<-read.table(geninfile,header=TRUE,as.is=1)


## the effect size of static ASE for a variant, ASERNA-seq,
## taking (or not taking) the absolute value
## ASERNA-seq = | p − 0.5 | (Eq. 1), where p = calt/(calt + cref ), and c is
## the read count and alt = alternative allele, ref = reference allele
## (Edsgard et al. 2016; https://doi.org/10.1038/srep21134)


bias<-round(median(genin$alt.dp/(genin$alt.dp+genin$ref.dp)),digits=3)

ase<-(genin$alt.dp/(genin$alt.dp+genin$ref.dp)) - bias

logfold<-log2((genin$alt.dp + 1)/(genin$ref.dp + 1))


## take median for each gene (vals are still for each snp)
med.ase<-ave(ase,as.factor(genin$gene),FUN=median)
names(med.ase)<-genin$gene
med.logfold<-ave(logfold,as.factor(genin$gene),FUN=median)
names(med.logfold)<-genin$gene
## retain unique genes
med.ase<-med.ase[!duplicated(names(med.ase))]
med.logfold<-med.logfold[!duplicated(names(med.logfold))]

tmp<-merge(data.frame(feat=names(med.ase),
                      med.ase=med.ase,
                      stringsAsFactors=FALSE),genout)
tmp<-merge(data.frame(feat=names(med.logfold),
                      med.logfold=med.logfold,
                      stringsAsFactors=FALSE),tmp)

tmp$minuslog10<- -log10(tmp$fdr)


## make combined file
refgenesub<-gsub("\\.[0-9]","",refgene$gene)
comb<-cbind(refgene[match(tmp$feat,refgenesub),1:8],tmp[,c(4,2,12,3,5:11)])

## ## ----
## write.table(comb,file=combaseoutfile,quote=FALSE,
##             row.names=FALSE,col.names=TRUE,sep='\t')
## ## ----
## ## (don't write this file anymore as better ways to compute a measure
## ##  of logfold change, e.g., using DESeq2, as in 'ASE_DE_analysis.R')


## check where reads map in the other genome
## -----------------------------------------
## (this is a heavy loop)


corresp<-read.table(correspfile,as.is=TRUE)


report<-data.frame(gene=character(nrow(comb)),
                   n.reads=numeric(nrow(comb)),
                   p.nonmap=numeric(nrow(comb)),
                   p.unimap=numeric(nrow(comb)),
                   p.major=numeric(nrow(comb)),
                   p.major.bin=numeric(nrow(comb)),
                   stringsAsFactors=FALSE)

scaffs<-unique(comb$chr)

## this loop takes a while to run...
for(s in 1:length(scaffs)){
    subcombidx<-which(comb[,2]==scaffs[s])
    subcomb<-comb[subcombidx,,drop=FALSE]
    subcorresp<-corresp[which(corresp[,3]==scaffs[s]),,drop=FALSE]
    for(j in 1:nrow(subcomb)){

        ## find reads mapping to the region (ignoring multimappers)
        subcorrespps<-subcorresp[,4]
        ## excl. multimappers
        subcorrespps[grepl(",",subcorrespps,fixed=TRUE)]<-NA
        ## excl. nonmappers
        subcorrespps[grepl(".",subcorrespps,fixed=TRUE)]<-NA
        subcorrespps<-as.numeric(subcorrespps)
        copos<-subcorresp[which(!is.na(subcorrespps) &
                                    subcorrespps+readlen-1>=subcomb[j,4]+overl-1 &
                                        subcorrespps<=subcomb[j,5]-overl+1),,drop=FALSE]

        ## check whether they include multimappers, and whether they
        ##  reads map all (most) to the same region in other genome
        ##tmp2<-unique(unlist(strsplit(copos[,12],",",fixed=TRUE)))
        ## 255 for unique mappers

        ## if there are multimappers, want to see different locations they
        ##  are mapping to (ignoring for now how common this is)

        if(nrow(copos)>0){

            allocs<-data.frame(scaf=character(),pos=numeric(),
                               stringsAsFactors=FALSE)
            nmaps<-numeric(nrow(copos))

            for(k in 1:nrow(copos)){
                if(copos[k,10]=="."){
                    newlocs<-data.frame(scaf=NA,pos=NA,stringsAsFactors=FALSE)
                    nmaps[k]<-0
                }else if(grepl(",",copos[k,10],fixed=TRUE)){
                    sc<-unlist(strsplit(copos[k,10],",",fixed=TRUE))
                    ps<-as.numeric(unlist(strsplit(copos[k,11],",",fixed=TRUE)))
                    newlocs<-data.frame(scaf=sc,pos=ps,stringsAsFactors=FALSE)
                    nmaps[k]<-length(sc)
                }else{
                    newlocs<-data.frame(scaf=copos[k,10],
                                        pos=as.numeric(copos[k,11]),
                                        stringsAsFactors=FALSE)
                    nmaps[k]<-1
                }
                allocs<-rbind(allocs,newlocs)
            }

            ##nrow(allocs)
            ##nrow(copos)

            ## prop nonmappers
            pnon<-sum(nmaps==0)/length(nmaps)
            ## prop unique mappers
            puni<-sum(nmaps==1)/length(nmaps)

            ## prop major scaff
            tallocs<-sort(table(allocs[,1]),decreasing=TRUE)
            pmaj<-tallocs[1]/nrow(allocs)

            ## bin positions of major scaffold
            majps<-allocs[allocs[,1]==names(tallocs)[1] & !is.na(allocs[,1]),2]
            if(length(majps)>0){
                mybreaks<-seq(min(majps),max(majps),10000)
                mybreaks<-c(mybreaks,max(majps))
                histdat<-hist(majps,mybreaks,plot=FALSE)
                ## $breaks and $counts
                majord<-order(histdat$counts,decreasing=TRUE)
                majcts<-histdat$counts[majord][1]
                ## sum counts in two adjacent cells, if any
                if(majord[1]!=1){ ## preceding
                    majcts<-majcts+histdat$counts[(majord[1]-1)]
                }
                if(majord[1]!=length(majord)){ ## following
                    majcts<-majcts+histdat$counts[(majord[1]+1)]
                }
                ## prop major scaff and major bin(s)
                pmajb<-majcts/sum(histdat$counts)
            }else{
                pmajb<-NA
                ## case for PeaxChr5g0039570.1      1       0       NA      NA
            }

            ## report results
            report[subcombidx[j],1]<-subcomb[j,1]
            report[subcombidx[j],2]<-length(nmaps)
            report[subcombidx[j],3]<-pnon
            report[subcombidx[j],4]<-puni
            report[subcombidx[j],5]<-pmaj
            report[subcombidx[j],6]<-pmajb

        }else{

            ## report results
            report[subcombidx[j],1]<-subcomb[j,1]
            report[subcombidx[j],2]<-nrow(copos)
            report[subcombidx[j],3]<-NA
            report[subcombidx[j],4]<-NA
            report[subcombidx[j],5]<-NA
            report[subcombidx[j],6]<-NA

        }
    }
}

## ----
write.table(report,file=mapreportfile,quote=FALSE,
            row.names=FALSE,col.names=TRUE,sep='\t')
## ----
