## -----------------------------------------------------------------
## "bi-directional" analysis of ASE/DE, based on
## output tables from ASE_DE_analysis.R and OMA orthologs
## (i.e., finding genes that are interesting no matter to
##  which reference genome they were mapped)
## -----------------------------------------------------------------

library("Rmpfr") ## for working with small numbers
library("scales") ## for transparent colors
library("pheatmap") ## for 'pheatmap'

setwd("/Users/doro/Documents/Style/analysis/")


## =================================================================
## FUNCTIONS
## =================================================================

## PLOTTING FUNCTIONS

## compute x-axis positions
makeX<-function(lgs,pos,sname="Scf",gap=50,sgap=0.2,scale=0.000001){

    if (length(lgs)!=length(pos)){
        stop("Vectors lgs and pos are not of equal length.")
    }

    ## order data
    tmpord<-data.frame(new=order(lgs,pos),old=1:length(pos))
    pos<-pos[tmpord$new]
    lgs<-lgs[tmpord$new]
    chr<-lgs[1:(grep("Scf",lgs)[1]-1)]
    scaf<-lgs[(grep("Scf",lgs)[1]):length(lgs)]
    ## assuming that there are only 'Chr' and 'Scf', so they will
    ##  be ordered correctly by name
    ## these will still include NAs in scaf and pos

    fchr<-as.numeric(as.factor(chr))
    fchr<-c(fchr,rep(length(unique(chr))+1,length(scaf)))
    xpos<-pos*scale
    ## remove NA
    tokeep<-complete.cases(c(chr,scaf),xpos)
    allgs<-c(chr,scaf)[tokeep]
    fchr<-fchr[tokeep]
    xpos<-xpos[tokeep]
    nloc<-length(xpos)
    chrnames<-c(unique(chr[complete.cases(chr,pos[1:length(chr)])]),"Scf")

    ## adjust precision (Rmpfr package)
    xcum<-Rmpfr::mpfr(xpos,50)
    cstart<-as.numeric(xcum[1])
    cend<-numeric()
    for (i in 2:nloc){
        prev<-fchr[i-1]
        if (fchr[i]!=prev){
            xcum[i:nloc]<-xcum[i:nloc]+gap+xpos[i-1]
            cstart<-c(cstart,as.numeric(xcum[i]))
            cend<-c(cend,as.numeric(xcum[i-1]))
        }
    }
    ## also add small gap for new scaffolds
    for (i in (length(chr)+2):nloc){
        prev<-c(chr,scaf)[i-1]
        if (c(chr,scaf)[i]!=prev){
            xcum[i:nloc]<-xcum[i:nloc]+sgap+xpos[i-1]
        }
    }
    cend<-c(cend,as.numeric(xcum[length(xcum)]))

    ## re-add NA
    rpos<-rep(NA,length(pos))
    rpos[tokeep]<-as.numeric(xcum)

    ## order back
    rpos<-rpos[order(tmpord$old[tmpord$new])]

    return(list(xpos=rpos,lgstart=cstart,lgend=cend,lgnames=chrnames))
}


## compute gene density
compGeneDens<-function(myxvals,chrnames,step=10){

    if(length(myxvals$xpos)!=length(chrnames)){
        stop("myxvals$xpos and chrnames not of same length")
    }
    densbreaks<-numeric()
    denscounts<-numeric()
    idx<-data.frame(start=numeric(length(myxvals$lgnames)),
                    end=numeric(length(myxvals$lgnames)))

    for(i in 1:length(myxvals$lgnames)){

        if(myxvals$lgnames[i]=="Scf"){
            sub<-grep("Scf",chrnames)
        }else{
            sub<-which(chrnames==myxvals$lgnames[i])
        }

        ho<-hist(myxvals$xpos[sub],
                 breaks=seq(myxvals$lgstart[i],myxvals$lgend[i]+step,step),
                 plot=FALSE)
        densbreaks<-c(densbreaks,rowMeans(cbind(ho$breaks[-(length(ho$breaks))],
                                                ho$breaks[-1])))
        idx$start[i]<-length(denscounts)+1
        denscounts<-c(denscounts,ho$counts)
        idx$end[i]<-length(denscounts)
    }

    return(list(breaks=densbreaks,counts=denscounts,idx=idx))
}


## plot results along genome
makeGenomePlot<-function(myxvals,ostat,xtxt,ytxt,mtxt,yaxlim,
                         ptc,c.ic,col.c.ic,p.np,col.p.np,sa,gd){

    ## !!! NO ERROR CHECKING !!!

    plot(myxvals$xpos,ostat,cex=ptc,
         axes=F,type="n",xlab="",ylab="",ylim=yaxlim)
    ## complete/incomplete
    x<-col.c.ic==c.ic[1]
    points(myxvals$xpos[x],ostat[x],col=alpha(col.c.ic[x],0.4),
           pch=16,cex=ptc[x])
    x<-col.c.ic==c.ic[2]
    points(myxvals$xpos[x],ostat[x],col=alpha(col.c.ic[x],0.4),
           pch=16,cex=ptc[x])
    ## filter pass
    x<-col.p.np==p.np[1]
    points(myxvals$xpos[x],ostat[x],col=alpha(col.p.np[x],sa),
           pch=16,cex=ptc[x])
    x<-col.p.np==p.np[2]
    points(myxvals$xpos[x],ostat[x],col=alpha(col.p.np[x],sa),
           pch=16,cex=ptc[x])
    x<-col.p.np==p.np[3]
    points(myxvals$xpos[x],ostat[x],col=alpha(col.p.np[x],sa),
           pch=16,cex=ptc[x])
    ## gene density
    for(i in 1:nrow(gd$idx)){
        sub<-(gd$idx[i,1]):(gd$idx[i,2])
        lines(gd$breaks[sub],
              (gd$counts/max(gd$counts))[sub],
              col="yellow",lwd=1.1)
    }

    axis(2,padj=0.5,hadj=0.6,cex.axis=0.8,las=2)
    mtext(ytxt,side=2,line=1.8)

    ## plot horizontal lines for LGs
    par(xpd=NA)
    myy<-yaxlim[1]-(yaxlim[2]-yaxlim[1])*0.05
    segments(x0=myxvals$lgstart, x1=myxvals$lgend, y0=myy, y1=myy, lwd=3)
    par(xpd=FALSE)

    ## axis(1,padj=0.5,hadj=1,cex.axis=0.8,las=2,
    ##      at=rowMeans(cbind(myxvals$lgstart,myxvals$lgend)),
    ##      labels=myxvals$lgnames,tick=F)
    ## mtext(xtxt,side=1,line=2.8)
    axis(1,hadj=0.5,padj=1,cex.axis=0.8,las=1,
         at=rowMeans(cbind(myxvals$lgstart,myxvals$lgend)),
         labels=myxvals$lgnames,tick=F,line= -1.0)
    mtext(xtxt,font=3,side=1,line=1.5)

    mtext(mtxt,font=2)

    legend("bottomleft",
           legend=c("complete cases","incomplete cases",
               "pass filters focal","pass filters both",
               "pass filters non-focal"),
           pch=16,col=alpha(c(c.ic[2],c.ic[1],p.np[2],p.np[3],
                      p.np[1]),0.4),bty="n",cex=0.8,
           y.intersp=0.9)
}

## =================================================================



## =================================================================
## PER-STAGE ANALYSIS
## =================================================================

mystage<-"MS1" ## do this for each of ("MS1","MS2","LS")

paxres<-read.csv(paste0("Pax/",mystage,"_results.csv"),as.is=TRUE)
pexres<-read.csv(paste0("Pex/",mystage,"_results.csv"),as.is=TRUE)

## prepare OMA output
ids<-read.table("/Users/doro/Documents/Style/data/Orthologs/Pax304-Pex304/Output/Map-SeqNum-ID.txt",as.is=TRUE)
om<-read.table("/Users/doro/Documents/Style/data/Orthologs/Pax304-Pex304/Output/OrthologousMatrix.txt",header=TRUE)

paxids<-ids[ids[,1]=="PAX",]
pexids<-ids[ids[,1]=="PEX",]

og<-data.frame(pax=paxids[match(om$PAX,paxids[,2]),3],
               pex=pexids[match(om$PEX,pexids[,2]),3],
               stringsAsFactors=FALSE)

og$pax<-gsub("\\.[[:digit:]]+","",og$pax)
og$pex<-gsub("\\.[[:digit:]]+","",og$pex)

## merge all dataframes
comb<-merge(og,paxres,by.x="pax",by.y="gene",sort=FALSE,all=TRUE)
comb<-merge(comb,pexres,by.x="pex",by.y="gene",sort=FALSE,all=TRUE)


## ---------------------
## set filter thresholds
fil.logfold<-1
fil.signif<-0.05 ## should be equal 'myAlpha' in 'ASE_DE_analysis.R'
fil.rpkm<-1 ## >=1 is generally used to call a gene "expressed"
## ---------------------
## NOTE: genes with 1 SNP will get excluded here due to filtering
##       for BHhet>fil.signif, which is NA with 1 SNP. See below in
##       all-stages analysis on how to change this


pass.x<-which(abs(comb$log2ase.x)>=fil.logfold &
                  comb$BHase.x<=fil.signif &
                      comb$BHhet.x>fil.signif &
                          abs(comb$log2de.x)>=fil.logfold &
                              comb$BHde.x<=fil.signif &
                                  comb$rpkm1.x>=fil.rpkm &
                                      comb$rpkm2.x>=fil.rpkm &
                                          comb$samedirection.x==TRUE &
                                              comb$passmapping.x==TRUE)

pass.y<-which(abs(comb$log2ase.y)>=fil.logfold &
                  comb$BHase.y<=fil.signif &
                      comb$BHhet.y>fil.signif &
                          abs(comb$log2de.y)>=fil.logfold &
                              comb$BHde.y<=fil.signif &
                                  comb$rpkm1.y>=fil.rpkm &
                                      comb$rpkm2.y>=fil.rpkm &
                                          comb$samedirection.y==TRUE &
                                              comb$passmapping.y==TRUE)

##table(comb$chr.x[intersect(pass.x,pass.y)])





## ------------------------
## GENOME PLOT P. AXILLARIS
## ------------------------

lgs<-comb$chr.x
pos<-rowMeans(cbind(comb$start.x,comb$end.x))

myxvals<-makeX(lgs=lgs,pos=pos,
               gap=50,sgap=0.2,scale=0.000001)
## returns $xpos, $lgstart, $lgend, $lgnames

## gene density
gd<-compGeneDens(myxvals,comb$chr.x)

ytxt<-"Log2_lfc Alt/Ref (ASE)"
ostat<-comb$log2ase.x
xtxt<-"P. axillaris"
mtxt<-mystage
yaxlim<-range(ostat,na.rm=T)
##ptc<-rep(0.5,length(ostat)) ## point cex
ptc<- -log(comb$BHase.x,base=1e+10)+0.5
sa<-0.4 ## alpha value for colors of genes passing filters

c.ic<-c("gray","black")
p.np<-c("purple","orange","red")

col.c.ic<-rep(c.ic[1],length(ostat))
col.c.ic[complete.cases(comb)]<-c.ic[2]
col.p.np<-rep("NA",length(ostat))
col.p.np[pass.y]<-p.np[1]
col.p.np[pass.x]<-p.np[2]
col.p.np[intersect(pass.x,pass.y)]<-p.np[3]

pdf(paste0("Pax/",mystage,"_genome_ase.pdf"))
makeGenomePlot(myxvals,ostat,xtxt,ytxt,mtxt,yaxlim,
               ptc,c.ic,col.c.ic,p.np,col.p.np,sa,gd)
dev.off()


## ----------------------
## GENOME PLOT P. EXSERTA
## ----------------------

lgs<-comb$chr.y
pos<-rowMeans(cbind(comb$start.y,comb$end.y))

myxvals<-makeX(lgs=lgs,pos=pos,
               gap=50,sgap=0.2,scale=0.000001)
## returns $xpos, $lgstart, $lgend, $lgnames

## gene density
gd<-compGeneDens(myxvals,comb$chr.y)

ytxt<-"Log2_lfc Alt/Ref (ASE)"
ostat<-comb$log2ase.y
xtxt<-"P. exserta"
mtxt<-mystage
yaxlim<-range(ostat,na.rm=T)
##ptc<-rep(0.5,length(ostat)) ## point cex
ptc<- -log(comb$BHase.y,base=1e+10)+0.5
sa<-0.4 ## alpha value for colors of genes passing filters

c.ic<-c("gray","black")
p.np<-c("purple","orange","red")

col.c.ic<-rep(c.ic[1],length(ostat))
col.c.ic[complete.cases(comb)]<-c.ic[2]
col.p.np<-rep("NA",length(ostat))
col.p.np[pass.x]<-p.np[1]
col.p.np[pass.y]<-p.np[2]
col.p.np[intersect(pass.y,pass.x)]<-p.np[3]

pdf(paste0("Pex/",mystage,"_genome_ase.pdf"))
makeGenomePlot(myxvals,ostat,xtxt,ytxt,mtxt,yaxlim,
               ptc,c.ic,col.c.ic,p.np,col.p.np,sa,gd)
dev.off()


## ------------------------------------------------------------------------


## Exploring
## ---------


## correlations p-values
pdf(paste0(mystage,"_pval_cor.pdf"),height=4,width=8)
par(mfrow=c(1,2),mar=c(3,3,1,1))

plot(-log10(comb$BHase.x),-log10(comb$BHase.y),cex=0.5,
     axes=FALSE,xlab="",ylab="")
abline(a=0,b=1,lty=2,col="darkgray")
mtext(mystage,font=2)
axis(1,hadj=0.5,padj=-0.6,cex.axis=0.8,las=1)
axis(2,padj=0.5,hadj=0.8,cex.axis=0.8,las=2)
mtext("-log10(BH p-value) (ASE), P.ax reference",side=1,line=1.8)
mtext("-log10(BH p-value) (ASE), P.ex reference",side=2,line=2.0)
## -> ASE: substantial divergence

plot(-log10(comb$BHde.x),-log10(comb$BHde.y),cex=0.5,
     axes=FALSE,xlab="",ylab="")
abline(a=0,b=1,lty=2,col="darkgray")
mtext(mystage,font=2)
axis(1,hadj=0.5,padj=-0.6,cex.axis=0.8,las=1)
axis(2,padj=0.5,hadj=0.8,cex.axis=0.8,las=2)
mtext("-log10(BH p-value) (DE), P.ax reference",side=1,line=1.8)
mtext("-log10(BH p-value) (DE), P.ex reference",side=2,line=2.0)
## -> DE: much better correlated
dev.off()



## correlations log-fold changes
pdf(paste0(mystage,"_log2_cor.pdf"),height=4,width=8)
par(mfrow=c(1,2),mar=c(3,3,1,1))

plot(comb$log2ase.x,comb$log2ase.y,cex=0.5,
     axes=FALSE,xlab="",ylab="")
abline(a=0,b= -1,lty=2,col="darkgray")
grid()
abline(h=0,col="lightgray")
abline(v=0,col="lightgray")
mtext(mystage,font=2)
axis(1,hadj=0.5,padj=-0.6,cex.axis=0.8,las=1)
axis(2,padj=0.5,hadj=0.8,cex.axis=0.8,las=2)
mtext("log2(alt/ref) (ASE), P.ax reference",side=1,line=1.8)
mtext("log2(alt/ref) (ASE), P.ex reference",side=2,line=2.0)

plot(comb$log2de.x,comb$log2de.y,cex=0.5,
     axes=FALSE,xlab="",ylab="")
abline(a=0,b= -1,lty=2,col="darkgray")
grid()
abline(h=0,col="lightgray")
abline(v=0,col="lightgray")
mtext(mystage,font=2)
axis(1,hadj=0.5,padj=-0.6,cex.axis=0.8,las=1)
axis(2,padj=0.5,hadj=0.8,cex.axis=0.8,las=2)
mtext("log2(P.ex/P.ax) (DE), P.ax reference",side=1,line=1.8)
mtext("log2(P.ax/P.ex) (DE), P.ex reference",side=2,line=2.0)
dev.off()
## => ASE seems to be more susceptible to bias



## =================================================================
## ALL-STAGES ANALYSIS
## =================================================================


## prepare OMA output
ids<-read.table("/Users/doro/Documents/Style/data/Orthologs/Pax304-Pex304/Output/Map-SeqNum-ID.txt",as.is=TRUE)
om<-read.table("/Users/doro/Documents/Style/data/Orthologs/Pax304-Pex304/Output/OrthologousMatrix.txt",header=TRUE)

paxids<-ids[ids[,1]=="PAX",]
pexids<-ids[ids[,1]=="PEX",]

og<-data.frame(pax=paxids[match(om$PAX,paxids[,2]),3],
               pex=pexids[match(om$PEX,pexids[,2]),3],
               stringsAsFactors=FALSE)

og$pax<-gsub("\\.[[:digit:]]+","",og$pax)
og$pex<-gsub("\\.[[:digit:]]+","",og$pex)



## ---------------------
## set filter thresholds
fil.logfold<-1
fil.signif<-0.05 ## should be equal 'myAlpha' in 'ASE_DE_analysis.R'
fil.rpkm<-1 ## >=1 is generally used to call a gene "expressed"
fil.nvars<-1 ## when filtering for $BHhet, all NA get excluded, i.e.,
##              the genes with only 1 SNP. To include these genes, set
##              fil.nvars<-1; to exclude, set fil.nvars<-0
## ---------------------


## combine data per stage
## ----------------------
mystages<-c("MS1","MS2","LS")

paxres<-vector("list",length(mystages))
names(paxres)<-mystages
pexres<-paxres
comb<-paxres
pass.x<-paxres
pass.y<-paxres


for(i in 1:length(mystages)){
    paxres[[i]]<-read.csv(paste0("Pax/",mystages[i],"_results.csv"),as.is=TRUE)
    pexres[[i]]<-read.csv(paste0("Pex/",mystages[i],"_results.csv"),as.is=TRUE)
    ## merge all dataframes
    comb[[i]]<-merge(og,paxres[[i]],by.x="pax",by.y="gene",sort=FALSE,all=TRUE)
    comb[[i]]<-merge(comb[[i]],pexres[[i]],by.x="pex",by.y="gene",
                     sort=FALSE,all=TRUE)

    ## filter
    pass.x[[i]]<-which(abs(comb[[i]]$log2ase.x)>=fil.logfold &
                           comb[[i]]$BHase.x<=fil.signif &
                               (comb[[i]]$BHhet.x>fil.signif |
                                    comb[[i]]$nvars.x==fil.nvars) &
                                        abs(comb[[i]]$log2de.x)>=fil.logfold &
                                            comb[[i]]$BHde.x<=fil.signif &
                                                comb[[i]]$rpkm1.x>=fil.rpkm &
                                                    comb[[i]]$rpkm2.x>=fil.rpkm &
                                                        comb[[i]]$samedirection.x==TRUE &
                                                            comb[[i]]$passmapping.x==TRUE)

    pass.y[[i]]<-which(abs(comb[[i]]$log2ase.y)>=fil.logfold &
                           comb[[i]]$BHase.y<=fil.signif &
                               (comb[[i]]$BHhet.y>fil.signif |
                                    comb[[i]]$nvars.y==fil.nvars) &
                                        abs(comb[[i]]$log2de.y)>=fil.logfold &
                                            comb[[i]]$BHde.y<=fil.signif &
                                                comb[[i]]$rpkm1.y>=fil.rpkm &
                                                    comb[[i]]$rpkm2.y>=fil.rpkm &
                                                        comb[[i]]$samedirection.y==TRUE &
                                                            comb[[i]]$passmapping.y==TRUE)
}


## merge data across stages for genes passing filter
## -------------------------------------------------
paxnames<-character()
pexnames<-character()

for(i in 1:length(mystages)){
    paxnames<-c(paxnames,comb[[i]]$pax[union(pass.x[[i]],pass.y[[i]])])
    pexnames<-c(pexnames,comb[[i]]$pex[union(pass.x[[i]],pass.y[[i]])])
}
paxpexnames<-data.frame(pax=paxnames,pex=pexnames,stringsAsFactors=FALSE)
paxpexnames<-paxpexnames[!duplicated(paxpexnames),]
## the ones with NA in one column don't have orthologs

paxpexnamesNA<-paxpexnames
paxpexnamesNA[is.na(paxpexnames)]<-"NA"
paxpexnamesNA<-paste(paxpexnamesNA$pax,paxpexnamesNA$pex,sep=":")

all.ase<-paxpexnames


for(i in 1:length(mystages)){

    tmp.x<-comb[[i]]$pax
    tmp.x[is.na(comb[[i]]$pax)]<-"NA"
    tmp.y<-comb[[i]]$pex
    tmp.y[is.na(comb[[i]]$pex)]<-"NA"
    tmp.paxpex<-paste(tmp.x,tmp.y,sep=":")


    idx<-match(paxpexnamesNA,tmp.paxpex)
    ## includes NAs for non-matches
    idx<-idx[!is.na(idx)]
    tmp.name<-paxpexnamesNA[!is.na(match(paxpexnamesNA,tmp.paxpex))]

    tmp.ase<-data.frame(name=tmp.name,
                        pax=comb[[i]]$pax[idx],pex=comb[[i]]$pex[idx],
                        log2ase.x=comb[[i]]$log2ase.x[idx],
                        log2ase.y=comb[[i]]$log2ase.y[idx],
                        stringsAsFactors=FALSE)

    tmp.ase<-merge(data.frame(name=paxpexnamesNA,stringsAsFactors=FALSE),
                   tmp.ase,by="name",sort=FALSE,all.x=TRUE)

    idx<-match(paxpexnamesNA,tmp.ase$name)
    if(!identical(paxpexnamesNA,tmp.ase$name[idx])){
        stop("gene names do not match")
    }
    ## bind together, and inverse sign for P.ex reference results
    all.ase<-cbind(all.ase,tmp.ase[idx,4], -1 * tmp.ase[idx,5])

}
colnames(all.ase)[-c(1:2)]<-paste0(rep(mystages,each=2),
                                   rep("_ase",2*length(mystages)),
                                   rep(c(".pax",".pex"),times=2))




## -------
## HEATMAP
## -------

## pheatmap calls hclust, which fails with too many NAs
## -> require >50% observations per row
sub<-apply(all.ase[,-c(1:2)],1,function(x) sum(!is.na(x))>sum(is.na(x)))
sub.ase<-all.ase[sub,]

M<-as.matrix(sub.ase[,-c(1:2)])

## define meaningful breaks to not only have pale colors in plot
mymax<-ceiling(max(abs(quantile(M,probs=c(0.025,0.975),na.rm=TRUE))))
mybreaks<-seq(-mymax,mymax,0.2)
##mycolors<-grDevices::cm.colors(n=length(mybreaks)-1)
cl<-colorRampPalette(c("dodgerblue3","white"))(ceiling(length(mybreaks)/2))
cr<-colorRampPalette(c("white","firebrick3"))(ceiling(length(mybreaks)/2))
mycolors<-c(cl[-length(cl)],cr[-1])



annotation_col = data.frame(
    Reference = factor(rep(c("PAX","PEX"), times=length(mystages)),
        ordered=TRUE),
    Stage = factor(rep(mystages, each=2),ordered=TRUE)
)
rownames(annotation_col) = colnames(M)

ann_colors = list(
    Reference = c(PAX="skyblue1", PEX="tomato1"),
    Stage = c(MS1="palegreen1", MS2="palegreen3", LS="palegreen4")
)


pdf("ASE_heatmap.pdf")
pheatmap(M,na_col="gray90",color=mycolors,breaks=mybreaks,
         show_rownames=FALSE,show_colnames=FALSE,
         annotation_col=annotation_col,cluster_col = FALSE,
         annotation_colors = ann_colors)
dev.off()
## NOTE: if plotting is not working, just try multiple times...
## e.g., when Error in UseMethod("depth") :
##  no applicable method for 'depth' applied to an object of class "NULL"




## ---------------------------------
## OVERREPRESENTATION ON CHROMOSOMES
## ---------------------------------

propchr.x<-vector("list",length(mystages))
names(propchr.x)<-mystages
propchr.y<-propchr.x

for(i in 1:length(mystages)){

    ## Pax
    degchr.x<-comb[[i]]$chr.x[pass.x[[i]]]
    compl.x<-cbind(comb[[i]]$log2ase.x,
                   comb[[i]]$BHase.x,
                   comb[[i]]$BHhet.x,
                   comb[[i]]$log2de.x,
                   comb[[i]]$BHde.x,
                   comb[[i]]$rpkm1.x,
                   comb[[i]]$rpkm2.x,
                   comb[[i]]$samedirection.x,
                   comb[[i]]$passmapping.x)
    allchr.x<-comb[[i]]$chr.x[complete.cases(compl.x)]

    degchr.x<-gsub("^Scf[[:digit:]]+","Scf",degchr.x)
    allchr.x<-gsub("^Scf[[:digit:]]+","Scf",allchr.x)

    degchr.x<-table(degchr.x)
    allchr.x<-table(allchr.x)

    propchr.x[[i]]<-rbind(degchr.x/sum(degchr.x),allchr.x/sum(allchr.x),
                          degchr.x,allchr.x)
    rownames(propchr.x[[i]])<-c("observed","expected","obs.cnts","all.cnts")

    ## Pex
    degchr.y<-comb[[i]]$chr.y[pass.y[[i]]]
    compl.y<-cbind(comb[[i]]$log2ase.y,
                   comb[[i]]$BHase.y,
                   comb[[i]]$BHhet.y,
                   comb[[i]]$log2de.y,
                   comb[[i]]$BHde.y,
                   comb[[i]]$rpkm1.y,
                   comb[[i]]$rpkm2.y,
                   comb[[i]]$samedirection.y,
                   comb[[i]]$passmapping.y)
    allchr.y<-comb[[i]]$chr.y[complete.cases(compl.y)]

    degchr.y<-gsub("^Scf[[:digit:]]+","Scf",degchr.y)
    allchr.y<-gsub("^Scf[[:digit:]]+","Scf",allchr.y)

    degchr.y<-table(degchr.y)
    allchr.y<-table(allchr.y)

    propchr.y[[i]]<-rbind(degchr.y/sum(degchr.y),allchr.y/sum(allchr.y),
                          degchr.y,allchr.y)
    rownames(propchr.y[[i]])<-c("observed","expected","obs.cnts","all.cnts")
}

## >>>> NOTE: to make this properly, do a chisq.test or similar
##            relative to the total number of genes per chromosome

pdf("ASE_chrom.pdf",height=4.5,width=7)
par(mfrow=c(3,2),mar=c(2.5,4,2,0.5))
for(i in 1:length(mystages)){
    ## Pax ref
    barplot(propchr.x[[i]][1:2,],beside=TRUE,main=paste(mystages[i],"P.ax"),
            ylab="Proportion")
    ## show data
    text(x=seq(1.5,1.5*3*8,3),y=propchr.x[[i]][1,],adj=c(0.5,0.98),
         labels=propchr.x[[i]][3,],col="white")

    ## Pex ref
    barplot(propchr.y[[i]][1:2,],beside=TRUE,main=paste(mystages[i],"P.ex"),
            ylab="Proportion")
    ## show data
    text(x=seq(1.5,1.5*3*8,3),y=propchr.y[[i]][1,],adj=c(0.5,0.98),
         labels=propchr.y[[i]][3,],col="white")
}
dev.off()



## =================================================================
## SOME MORE DIAGNOSTIC AND OTHER PLOTS
## =================================================================

## DE volcano
pdf(paste0("Pax/","DE_volcano.pdf"),height=5,width=5)

for(i in 1:length(mystages)){

    volcols<-rep("gray",nrow(comb[[i]]))
    volcols[comb[[i]]$log2de.x<= -fil.logfold &
                comb[[i]]$BHde.x<=fil.signif]<-"dodgerblue3"
    volcols[comb[[i]]$log2de.x>=fil.logfold &
                comb[[i]]$BHde.x<=fil.signif]<-"firebrick3"

    plot(comb[[i]]$log2de.x,-log10(comb[[i]]$BHde.x),cex=0.5,pch=16,
         col=alpha(volcols,0.5),axes=FALSE,xlab="",ylab="")
    abline(v=0,col="darkgray",lty=2)
    mtext(mystages[i],font=2)
    axis(1,hadj=0.5,padj=-0.6,cex.axis=0.8,las=1)
    axis(2,padj=0.5,hadj=0.8,cex.axis=0.8,las=2)
    mtext("log2(P.ex/P.ax) (DE)",side=1,line=1.8)
    mtext("-log10(BH p-value) (DE)",side=2,line=2.0)
}

dev.off()

## ASE volcano
pdf(paste0("Pax/","ASE_volcano.pdf"),height=5,width=5)

for(i in 1:length(mystages)){

    volcols<-rep("gray",nrow(comb[[i]]))
    volcols[comb[[i]]$log2ase.x<= -fil.logfold &
                comb[[i]]$BHase.x<=fil.signif]<-"dodgerblue3"
    volcols[comb[[i]]$log2ase.x>=fil.logfold &
                comb[[i]]$BHase.x<=fil.signif]<-"firebrick3"

    plot(comb[[i]]$log2ase.x,-log10(comb[[i]]$BHase.x),cex=0.5,pch=16,
         col=alpha(volcols,0.5),axes=FALSE,xlab="",ylab="")
    abline(v=0,col="darkgray",lty=2)
    mtext(mystages[i],font=2)
    axis(1,hadj=0.5,padj=-0.6,cex.axis=0.8,las=1)
    axis(2,padj=0.5,hadj=0.8,cex.axis=0.8,las=2)
    mtext("log2(alt/ref) (ASE)",side=1,line=1.8)
    mtext("-log10(BH p-value) (ASE)",side=2,line=2.0)
}

dev.off()



## DE vs. ASE
pdf(paste0("Pax/","DE_ASE_scatter.pdf"),height=5,width=5)

for(i in 1:length(mystages)){

    volcols<-rep("black",nrow(comb[[i]]))

    plot(comb[[i]]$log2de.x,comb[[i]]$log2ase.x,cex=0.4,pch=16,
         col=alpha(volcols,0.4),axes=FALSE,xlab="",ylab="",type="n")
    abline(v=0,col=alpha("orange",0.6),lty=2,lwd=1.8)
    abline(h=0,col=alpha("dodgerblue3",0.3),lty=1,lwd=1.8)
    abline(a=0,b=1,col=alpha("firebrick3",0.3),lty=1,lwd=1.8)
    points(comb[[i]]$log2de.x,comb[[i]]$log2ase.x,cex=0.4,pch=16,
           col=alpha(volcols,0.4))
    mtext(mystages[i],font=2)
    axis(1,hadj=0.5,padj=-0.6,cex.axis=0.8,las=1)
    axis(2,padj=0.5,hadj=0.8,cex.axis=0.8,las=2)
    mtext("log2(P.ex/P.ax) (DE)",side=1,line=1.8)
    mtext("log2(alt/ref) (ASE)",side=2,line=2.0)
}

dev.off()





## ## ==================================================
## save(list=c("all.ase","comb","pass.x","pass.y","mystages",
##          "fil.logfold","fil.nvars","fil.rpkm","fil.signif"),
##      file="ASE_DE_analysis_bidir_fil.nvars1.RData")
## ## ==================================================
## ==================================================
##load(file="ASE_DE_analysis_bidir_fil.nvars1.RData")
## ==================================================

## look at a specific gene:
##t(comb[[1]][grepl("PeaxChr1g0013530",comb[[1]]$pax),])

