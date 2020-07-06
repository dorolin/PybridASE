
## ## on ubelix:
## ## start interactive session:
## salloc --nodes=1 --ntasks-per-node=1 --mem-per-cpu=10G --time=01:00:00
## module load vital-it/7
## module load R/latest
## R

## filters
## -------

## SNPs covered by at least minrcov reads in at least minscov samples per set
## for species, counts are to minfix for correct allele (sp1: ref; sp2: alt)
minrcov<-10
minscov<-2
minfix<-0.98
mindist<-101 ## minimum distance between SNPs


## outfiles
## --------

outfile<-"input.tab"
outdir<-"." ## this is only prefix for input directories above


## ------------
read.data<-function(snps,sampledirs,samplefiles){
    ## !!! Note: no checking that files are correctly formatted

    ## initialize data frames to store allele counts
    ref.tab<-setNames(data.frame(matrix(ncol=length(sampledirs),
                                        nrow=nrow(snps))),sampledirs)
    alt.tab<-ref.tab
    cnts.tab<-list(ref = ref.tab, alt = alt.tab)

    for(i in 1:length(sampledirs)){

        infile<-paste0(sampledirs[i],"/",samplefiles)

        cnts<-read.table(infile, header=TRUE, as.is=TRUE)

        snpcnts<-merge(snps,cnts,sort=FALSE,all.x=TRUE,
                       by.x=c("chr","start","ref","alt"),
                       by.y=c("contig","position","refAllele","altAllele"))
        ## sort back to original order
        snpcnts<-snpcnts[order(snpcnts$id),]
        rownames(snpcnts)<-snpcnts$id

        ## add counts to big table
        cnts.tab$ref[,i]<-snpcnts$refCount
        cnts.tab$alt[,i]<-snpcnts$altCount

    }

    return(cnts.tab)
}
## ------------


## --------------------------------------------------------------------

for(s in 1:6){ ## loop over genomes (Pax/Pex) and stages (MS1/MS2/LS)

    ## infile names
    ## ------------

    if(s %in% 1:3){ ## Pax genome

        ## annotated snps from ANNOVAR
        annosnpsfile<-"$HOME/Style/Pax_genome/annovar/snps_filtered_sel_anno.variant_function"
        ## contains directories per sample with ASEReadCounter data
        indir<-"$HOME/Style/step04_AxEx/Pax_genome/ASE"

        ## directory names for ASEReadCounter data
        ## (sp1 must correspond to reference genome)
        if(s==1){
            stage<-"MS1"
            sp1<-c("AxMS1_1","AxMS1_2","AxMS1_3")
            sp2<-c("ExMS1_1","ExMS1_2","ExMS1_3")
            hyb<-c("AxExMS1_1","AxExMS1_2","AxExMS1_3")
        }else if(s==2){
            stage<-"MS2"
            sp1<-c("AxMS2_1","AxMS2_2","AxMS2_3")
            sp2<-c("ExMS2_1","ExMS2_2","ExMS2_3")
            hyb<-c("AxExMS2_1","AxExMS2_2","AxExMS2_3")
        }else if(s==3){
            stage<-"LS"
            sp1<-c("AxLS_1","AxLS_2","AxLS_3")
            sp2<-c("ExLS_2","ExLS_4","ExLS_5")
            hyb<-c("AxExLS_1","AxExLS_2","AxExLS_3")
            ## for LS, exclude: "ExLS_1", "ExLS_3", "ExLS_6"
            ##                  "AxLS_4", "AxLS_5", "AxLS_6"
        }

    }else if(s %in% 4:6){ ## Pex genome

        ## annotated snps from ANNOVAR
        annosnpsfile<-"$HOME/Style/Pex_genome/annovar/snps_filtered_sel_anno.variant_function"
        ## contains directories per sample with ASEReadCounter data
        indir<-"$HOME/Style/step04_AxEx/Pex_genome/ASE"

        ## directory names for ASEReadCounter data
        ## (sp1 must correspond to reference genome)
        if(s==4){
            stage<-"MS1"
            sp1<-c("ExMS1_1","ExMS1_2","ExMS1_3")
            sp2<-c("AxMS1_1","AxMS1_2","AxMS1_3")
            hyb<-c("AxExMS1_1","AxExMS1_2","AxExMS1_3")
        }else if(s==5){
            stage<-"MS2"
            sp1<-c("ExMS2_1","ExMS2_2","ExMS2_3")
            sp2<-c("AxMS2_1","AxMS2_2","AxMS2_3")
            hyb<-c("AxExMS2_1","AxExMS2_2","AxExMS2_3")
        }else if(s==6){
            stage<-"LS"
            sp1<-c("ExLS_2","ExLS_4","ExLS_5")
            sp2<-c("AxLS_1","AxLS_2","AxLS_3")
            hyb<-c("AxExLS_1","AxExLS_2","AxExLS_3")
            ## for LS, exclude: "ExLS_1", "ExLS_3", "ExLS_6"
            ##                  "AxLS_4", "AxLS_5", "AxLS_6"
        }

    }

    setwd(indir)



    ## read data
    ## ---------

    ## read snps annotation file

    snps<-read.table(annosnpsfile,as.is=TRUE)
    colnames(snps)<-c("annotation","gene","chr","start","end","ref","alt")

    ## filter annotations
    ## ------------------
    ## exclude snps that are intergenic, upstream, downstream, or both
    tmp<-grep("(intergenic)|(downstream)|(upstream)",snps$annotation)
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


    if(s %in% c(1,4)){
        write.table(snps,file=paste0(outdir,"/","SNPs.txt"),
                    col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
    }


    ## read ASEReadCounter output
    ## --------------------------
    cnts.sp1<-read.data(snps=snps,sampledirs=sp1,samplefiles="counts.txt")
    cnts.sp2<-read.data(snps=snps,sampledirs=sp2,samplefiles="counts.txt")
    cnts.hyb<-read.data(snps=snps,sampledirs=hyb,samplefiles="counts.txt")
    ## -> makes list of two matrices: $ref and $alt, each has samples in
    ##    columns and counts for SNPs in rows


    ## filter data
    ## -----------

    ## sp1
    tmp<-apply(cnts.sp1$ref + cnts.sp1$alt, 2, function(x) !is.na(x) & x>=minrcov)
    keep.sp1<-rowSums(tmp)>=minscov
    tmp<-rep(TRUE,nrow(cnts.sp1$ref))
    for(i in 1:ncol(cnts.sp1$ref)){
        mf<-cnts.sp1$ref[,i]/(cnts.sp1$ref[,i] + cnts.sp1$alt[,i]) >= minfix
        tmp[mf==FALSE | is.na(mf)]<-FALSE
    }
    keep.sp1<-keep.sp1 & tmp

    ## sp2
    tmp<-apply(cnts.sp2$ref + cnts.sp2$alt, 2, function(x) !is.na(x) & x>=minrcov)
    keep.sp2<-rowSums(tmp)>=minscov
    tmp<-rep(TRUE,nrow(cnts.sp2$ref))
    for(i in 1:ncol(cnts.sp2$ref)){
        mf<-cnts.sp2$alt[,i]/(cnts.sp2$ref[,i] + cnts.sp2$alt[,i]) >= minfix
        tmp[mf==FALSE | is.na(mf)]<-FALSE
    }
    keep.sp2<-keep.sp2 & tmp

    ## hyb
    tmp<-apply(cnts.hyb$ref + cnts.hyb$alt, 2, function(x) !is.na(x) & x>=minrcov)
    keep.hyb<-rowSums(tmp)>=minscov
    ##keep.hyb<-rowSums(tmp)==ncol(tmp)
    tmp<-rep(TRUE,nrow(cnts.hyb$ref))
    for(i in 1:ncol(cnts.hyb$ref)){
        ## exclude NAs for all, not only minscov
        tmp[is.na(cnts.hyb$ref[,i] + cnts.hyb$alt[,i])]<-FALSE
    }
    keep.hyb<-keep.hyb & tmp

    ## combine filters (incl. annotation to gene)
    keep<-keep.sp1 & keep.sp2 & keep.hyb & !is.na(snps$filtgene)


    ## check distance among SNPs (mindist)
    keptgene<-snps$filtgene[keep]
    keptpos<-rowMeans(cbind(snps$start,snps$end))[keep]
    keptid<-snps$id[keep]

    ## snp coverage across all samples
    snpcov<-cbind((cnts.sp1$ref-cnts.sp1$alt),(cnts.sp2$alt-cnts.sp2$ref),
                  (cnts.hyb$ref+cnts.hyb$alt))[keep,]
    snpcov<-rowSums(snpcov,na.rm=TRUE)

    ## in case file wasn't sorted
    tmpord<-data.frame(new=order(keptgene,keptpos),old=1:length(keptgene))
    keptgene<-keptgene[tmpord$new]
    keptpos<-keptpos[tmpord$new]
    keptid<-keptid[tmpord$new]
    snpcov<-snpcov[tmpord$new]


    tokeep<-rep(TRUE,length(keptpos))
    cntr<-0
    curgene<-keptgene[1]
    curpos<- -(mindist+1)
    while(cntr<length(keptpos)){
        idx<-cntr+1
        if(keptgene[idx]!=curgene){
            ## new gene
            curgene<-keptgene[idx]
            curpos<-keptpos[idx]
        }else if(keptpos[idx]>curpos+mindist){
            ## sufficient distance
            curpos<-keptpos[idx]
        }else{
            ## SNPs are too close
            collect<-1
            while(collect>0 & idx<length(keptpos)){
                if((keptpos[idx+1]>keptpos[idx]+mindist) |
                   (keptgene[idx+1]!=keptgene[idx])){
                    collect<-0
                }else{
                    idx<-idx+1
                }
            }
            ## set of SNPs that are close: cntr:idx
            testpos<-cntr:idx
            repr<-sample(x=testpos,size=1,prob=snpcov[testpos])
            ## set others to FALSE in tokeep
            tokeep[setdiff(testpos,repr)]<-FALSE
            curgene<-keptgene[idx]
            curpos<-keptpos[idx]
        }
        cntr<-idx
    }

    ## order back
    tokeep<-tokeep[order(tmpord$old[tmpord$new])]
    keptgene<-keptgene[order(tmpord$old[tmpord$new])]
    keptpos<-keptpos[order(tmpord$old[tmpord$new])]
    keptid<-keptid[order(tmpord$old[tmpord$new])]
    snpcov<-snpcov[order(tmpord$old[tmpord$new])]
    ## ## test that ordering is alright
    ## sum(snps$filtgene[keep]==keptgene)==length(keptgene)
    ## sum(rowMeans(cbind(snps$start,snps$end))[keep]==keptpos)==length(keptpos)
    ## sum(snps$id[keep]==keptid)==length(keptid)


    ## add the new filter
    keep2<-rep(FALSE,length(keep))
    keep2[keep]<-tokeep



    hyb.dat2<-data.frame(gene=snps$filtgene[keep2],
                         chr=snps$chr[keep2],
                         snp.id=snps$id[keep2],
                         snp.pos=snps$start[keep2],
                         ref.dp1=cnts.hyb$ref[keep2,1],
                         ref.dp2=cnts.hyb$ref[keep2,2],
                         ref.dp3=cnts.hyb$ref[keep2,3],
                         alt.dp1=cnts.hyb$alt[keep2,1],
                         alt.dp2=cnts.hyb$alt[keep2,2],
                         alt.dp3=cnts.hyb$alt[keep2,3])


    ## make output files for hybrids
    ## -----------------------------

    for(i in 1:length(hyb)){
        tmp<-data.frame(gene=snps$filtgene[keep2],
                        snp.id=snps$id[keep2],
                        alt.dp=cnts.hyb$alt[keep2,i],
                        ref.dp=cnts.hyb$ref[keep2,i])

        ## counts for geneiase
        write.table(tmp,file=paste0(outdir,"/",hyb[i],"/",outfile),
                    quote=FALSE,sep='\t',row.names=FALSE)

        ## bias (typically estimated from DNA data)
        write.table(round(median(tmp$alt.dp/(tmp$alt.dp+tmp$ref.dp)),digits=3),
                    file=paste0(outdir,"/",hyb[i],"/","bias.txt"),
                    quote=FALSE,col.names=FALSE,row.names=FALSE)

        ## hist
        pdf(paste0(outdir,"/",hyb[i],"/","bias.pdf"),height=4,width=4,
            pointsize=12)
        hist(tmp$ref.dp/(tmp$alt.dp+tmp$ref.dp),breaks=500,
             main=hyb[i],xlab="ref.dp/(alt.dp+ref.dp)")
        abline(v=median(tmp$ref.dp/(tmp$alt.dp+tmp$ref.dp)),
               col="firebrick",lty=1,lwd=2)
        abline(v=0.5,col="forestgreen",lty=3,lwd=2)
        dev.off()

        ## *bias.txt is used for -p parameter in geneiase: Betabinomial
        ## probability of success for the null dist
        ## (p = alternative allele count
        ## / total allele count). Typically estimated from DNA data and
        ## <0.5 due to mapping bias. Default: 0.49
        ## (-> it thus is the opposite of the red line in *bias.pdf)

    }

    ## in addition, save counts for parentals
    ## --------------------------------------

    nsamp<-length(sp1)+length(hyb)+length(sp2)

    allcts<-matrix(0,ncol=nsamp*2,nrow=sum(keep2))
    colnames(allcts)<-paste(rep(c(sp1,hyb,sp2),each=2),
                            rep(c("alt","ref"),times=nsamp),
                            sep=".")

    for(i in 1:length(sp1)){
        ##print((i-1)*2+1)
        ##print(i*2)
        allcts[,(i-1)*2+1]<-cnts.sp1$alt[keep2,i]
        allcts[,i*2]<-cnts.sp1$ref[keep2,i]
    }
    for(i in 1:length(hyb)){
        ##print(length(sp1)*2+(i-1)*2+1)
        ##print(length(sp1)*2+i*2)
        allcts[,length(sp1)*2+(i-1)*2+1]<-cnts.hyb$alt[keep2,i]
        allcts[,length(sp1)*2+i*2]<-cnts.hyb$ref[keep2,i]
    }
    for(i in 1:length(sp2)){
        ##print(length(c(sp1,hyb))*2+(i-1)*2+1)
        ##print(length(c(sp1,hyb))*2+i*2)
        allcts[,length(c(sp1,hyb))*2+(i-1)*2+1]<-cnts.sp2$alt[keep2,i]
        allcts[,length(c(sp1,hyb))*2+i*2]<-cnts.sp2$ref[keep2,i]
    }

    tmp<-data.frame(gene=snps$filtgene[keep2],
                    snp.id=snps$id[keep2])
    tmp<-cbind(tmp,allcts)

    write.table(tmp,file=paste0(outdir,"/","allcounts_",stage,".txt"),
                quote=FALSE,sep='\t',row.names=FALSE)


}

## --------------------------------------------------------------------
