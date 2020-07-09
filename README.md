# Allele-specific expression analysis

Steps and scripts to analyse RNA-seq data for allele-specific expression (ASE) and differential expression (DE) analysis, with emphasis on the reduction of reference bias

## Main steps
 
1. [prepare raw reads](#1-prepare-raw-reads)
2. [prepare genomes](#2-prepare-genomes)
3. [align reads](#3-align-reads)
4. [process bam files](#4-process-bam-files)
5. [call SNPs](#5-call-snps)
6. [filter variants](#6-filter-variants)
7. [annotate SNPs](#7-annotate-snps)
8. [ASE and DE analysis](#8-ase-and-de-analysis)

* most computations were performed on [UBELIX](http://www.id.unibe.ch/hpc), the HPC cluster at the University of Bern
* some less intense computations were performed locally on macOS

## 1. Prepare raw reads

* raw RNA-seq data are on NCBI SRA [PRJNA533335](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA533335/)
* these are 72 RNA-seq accessions: `DataPrep/SRA_metadata.csv`, with sample names coded as follows:
  - species names: Ax, *Petunia axillaris*; Ex, *P. exserta*; Pa, *P. parodii*; AxEx/AxPa/ExPa, F1 hybrids
  - developmental stages: SS, small stage; MS1, medium 1 stage; MS2, medium 2 stage; LS, long stage
  - 3-6 biological replicates per condition
* 101bp single-end Illumina sequencing with stranded data
* trimming using [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) version 0.36 and quality checks using [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) version 0.11.7: `DataPrep/qualityControl.sh`

  &rarr; biased sequence content for the first couple of bases, which is supposed not to be a problem<br/>
  &rarr; high level of sequence duplication, which is expected for RNA-seq data

* only sequences from *P. axillaris*, *P. exserta*, and their F1 were used for further analysis: `files_AxEx.txt`


## 2. Prepare genomes

* working with the *Petunia axillaris* and the *P. exserta* genome in parallel
* *P. axillaris* v3.0.4 assembly (FASTA) from [CoGe](https://genomevolution.org/coge/) (not public) and annotations (GFF) from UBELIX `$HOME/../shared/304/Maker/Peax304_final_edited.gff` (not public)
* *P. exserta* v3.0.3 assembly (FASTA) from CoGe (not public) and v3.0.4 annotations (GFF) from UBELIX `$HOME/../shared/304/Maker/Peex304_final_edited.gff` (not public)  
  NOTE: as *P. exserta* v3.0.3 and v3.0.4 FASTA are identical except scaffold names, these were simply renamed in the GFF
* convert GFF to GTF with [gffread](https://github.com/gpertea/gffread) version 0.11.4 `gffread file.gff -T -E -o file.gtf`
* convert GFF to BED with [ea-utils](https://expressionanalysis.github.io/ea-utils/) (May 1, 2014) `gtf2bed file.gtf >file.bed`
* create dictionary and index FASTA with `picard-tools CreateSequenceDictionary file.fa` and `samtools faidx file.fa` 

### 2.1. Generate vcf with variant positions

* identify positions that differ between *P. axillaris* and *P. exserta*, and use them during alignment to reduce reference mapping bias
* required to [align reads](#3-align-reads) with STAR in WASP mode, and to [generate genome with IUPAC ambiguity codes](#22-generate-genome-with-iupac-ambiguity-codes) for GATK SplitNCigarReads
* as no genomic DNA reads for the used accessions were available, Illumina reads from the assemblies of the *P. exserta* and the *P. axillaris* genomes were used to identify SNPs
* *P. exserta* DNA reads were mapped to the *P. axillaris* genome, and vice versa, to produce bam files (done previously)
* pre-process bam files with mapped DNA reads using [picard-tools](https://broadinstitute.github.io/picard/) version 2.18.11 (i.e., SortSam, AddOrReplaceReadGroups, MarkDuplicates, and ReorderSam): `GenomePrep/genomeSNPs_preprocessing.sh`
* to reduce runtime for SNP calling, the *P. axillaris* genome was split into eight intervals, but this was not done for the *P. exserta* genome (see below for how the genomes were split into many more [intervals](#23-split-genome-into-intervals) for SNP calling with RNA-seq reads)
* call SNPs for DNA reads: `GenomePrep/genomeSNPs_call.sh` (picard-tools CreateSequenceDictionary and [samtools](http://www.htslib.org) faidx to make genome .dict and .fai files; samtools index to prepare bam files; GATK HaplotypeCaller in single sample mode to call SNPs)
* filter SNPs with GATK: `GenomePrep/genomeSNPs_filter.sh` (first round of SelectVariants and VariantFiltration for quality control, second round to find SNPs that are ALT/ALT, i.e., *P. exserta* is homozygous for different allele)  
  &rarr; generates *Exreads_on_Pax304_snps_filtered.tsv* (quality control) and *Exreads_on_Pax304_snps_filtered_2_sel.vcf* (ALT/ALT positions)

* check filters with GATK VariantsToTable (included in `GenomePrep/genomeSNPs_filter.sh`) and plot results with R

     ```R
     st<-read.table("Exreads_on_Pax304_snps_filtered.tsv",header=TRUE)
     sts<-st$DP[st$DP<400]

     pdf("Pex_on_Pax_DP.pdf")
     plot(density(sts,na.rm=TRUE),xlim=c(0,400),
          xlab="Read depth",main="P. exserta reads on P. axillaris genome")
     dev.off()

     pdf("Pex_on_Pax.pdf")
     for(i in 8:13){
        plot(density(st[st[,i]>=quantile(st[,i],probs=0.05,na.rm=TRUE) & 
                     st[,i]<=quantile(st[,i],probs=0.95,na.rm=TRUE),i],na.rm=TRUE), 
                     xlim=c(quantile(st[,i],probs=c(0.05,0.95),na.rm=TRUE)),
                     xlab=colnames(st)[i],main="P. ex on P. ax")
     }
     dev.off()
     ```

### 2.2. Generate genome with IUPAC ambiguity codes

* re-code ALT/ALT homozygous positions to REF/ALT heterozygous positions and index with bcftools: `GenomePrep/processVcfFiles_index.sh`
(calls `GenomePrep/manipulate_vcf.pl`: NOTE: this script is not overly sophisticated and assumes certain input format)
* re-code these REF/ALT positions with IUPAC ambiguity codes with GATK FastaAlternateReferenceMaker: `GenomePrep/manipulateReference.sh` (additionally calls `GenomePrep/manipulate_fasta.pl` and includes picard-tools CreateSequenceDictionary and samtools faidx to make new genome .dict and .fai files)  
  &rarr; generates *Petunia_axillaris_304_PexIUPAC.fa*
  
### 2.3. Split genome into intervals

* split genome into intervals to make SNP calling manageable in terms of memory and run time
* make (initially 200) genome intervals by splitting chromosomes/scaffolds by Ns: `GenomePrep/splitIntervals.sh` (picard-tools ScatterIntervalsByNs and GATK SplitIntervals)

* count number of intervals in each list

     ```bash
     IDIR=$HOME/Style/Pax_genome/intervals/intervalLists/
     for INTVL in {0..199}; do
        IFILE=$(ls ${IDIR}/*scattered.interval_list | sed -n "$(( ${INTVL} + 1 ))p")
        cnt=`grep -v "^@" ${IFILE} | grep -c "^Peax"`
        echo "${INTVL} ${cnt}" >>intvlcnts.txt
     done
     ````

  &rarr; a large number of intervals is contained in lists from 0172 to 0199 (because of a bug?)

* make second round of interval splitting for lists that contain a large number of intervals: `GenomePrep/splitIntervals_2nd.sh`

* count number of intervals in new lists

     ```bash
     ## run in $HOME/Style/Pax_genome/intervals/intervalLists/split_2nd/
     grep -c "^Peax" ./* >intvlcnt.txt
     cut -d ':' -f2 intvlcnt.txt | sort | uniq
     ```

  &rarr; again, a large number of intervals is contained in list 0582

* make third round of interval splitting for interval 0582: `GenomePrep/splitIntervals_3rd.sh`

* count number of new intervals

     ```bash
     ## run in $HOME/Style/Pax_genome/intervals/intervalLists/split_3rd/
     grep -c "^Peax" ./* >intvlcnt.txt
     cut -d ':' -f2 intvlcnt.txt | sort | uniq
     ```
  &rarr; done; created a total of 588 intervals, which are used below to [call SNPs](#5-call-snps) in RNA-seq data


### 2.4. Identify orthologs

* using [OMA StandAlone](https://genome.cshlp.org/content/29/7/1152.abstract) version 2.3.1 to identify orthologs between *P. axillaris* and *P. exserta*
* prepare OMA data base: `orthologs/oma_prep.sh` (calls `orthologs/fastaAdjust.pl`, `orthologs/fastaFilter.pl`, and `orthologs/fastaInfo.pl`; the perl scripts `orthologs/fastaAdjust.pl` and `orthologs/fastaFilter.pl` are based on code from the orthoMCL orthomclAdjustFasta tool from the [OrthoMCL software](https://orthomcl.org/common/downloads/software/v2.0/) ([Li et al. 2003](https://genome.cshlp.org/content/13/9/2178.abstract)))
* run OMA in tree steps (copy and paste contents of scripts): 
* step 1: `orthologs/oma_ube_0.sh`
* step 2: `orthologs/oma_ube_1st.sh`
* step 3: `orthologs/oma_ube_2nd.sh`

  &rarr; generates *Map-SeqNum-ID.txt* and *OrthologousMatrix.txt*


## 3. Align reads

* using [STAR](https://github.com/alexdobin/STAR) version 2.6.0c
* working with the *P. axillaris* and the *P. exserta* genome in parallel
* STAR in WASP mode requires vcf file with variant positions (e.g., *Exreads_on_Pax304_snps_filtered_2_sel.vcf*, [generated above](#21-generate-vcf-with-variant-positions))
* step 1: generate genome indexes: `Alignment/prepareGenome.sh`  
  NOTE: there should be nothing else in the genome directory as this step may fail with an unspecific message otherwise
* step 2: mapping in two-pass mode: `Alignment/mapReads_firstPass.sh` and `Alignment/mapReads_secondPass.sh`  
  NOTE: a potentially better way to deal with a large number of reads not aligned because they are "too short" is to modify the penalty for indels during mapping instead of setting --outFilterScoreMinOverLread 0.3 and --outFilterMatchNminOverLread 0.3 (needs to be tested)
  
* check mapping performance

     ```bash
     ## get basic mapping statistics for each sample
     dir=`sed -n '1p' $HOME/Style/step02_AxEx/files_AxEx.txt`
     echo "sample" >mapping_stats.txt
     grep "%" ${dir}/Log.final.out | cut -f1 | sed 's/ |$//g' | sed -E 's/^\s+//g' >>mapping_stats.txt
     while read dir; do
        echo "$dir" >tmp
        grep "%" ${dir}/Log.final.out | cut -f2 | sed 's/%$//g' >>tmp
        paste mapping_stats.txt <(cat tmp) >tmp2
        mv -f tmp2 mapping_stats.txt
        rm -f tmp
     done <$HOME/Style/step02_AxEx/files_AxEx.txt
     ```
     
     ```R
     ## make plots with R
     st<-read.table("mapping_stats.txt",header=TRUE,sep='\t',as.is=TRUE)

     pax<-which(grepl("Ax",colnames(st)) & !grepl("Ex",colnames(st)))
     pex<-which(grepl("Ex",colnames(st)) & !grepl("Ax",colnames(st)))
     f1<-which(grepl("AxEx",colnames(st)))
     
     grp<-character(ncol(st)-1)
     grp[pax-1]<-"1"
     grp[pex-1]<-"3"
     grp[f1-1]<-"2"
     grp<-as.factor(grp)
     
     pdf("mapping_stats.pdf",pointsize=16)
     par(mfrow=c(2,2),mar=c(3,4,1,1))
     for(i in c(1,2,5,8)){
        boxplot(as.numeric(st[i,-1]) ~ grp, ylab=st[i,1], names=c("Pax","F1","Pex"), xlab="")
     }
     dev.off()
     
     pdf("mapping_stats_indel.pdf",pointsize=14,height=3.5)
     par(mfrow=c(1,2),mar=c(3,4,1,1))
     for(i in c(3,4)){
        boxplot(as.numeric(st[i,-1]) ~ grp, ylab=st[i,1], names=c("Pax","F1","Pex"), xlab="")
     }
     dev.off()
     
     pdf("mapping_stats_other.pdf",pointsize=16)
     par(mfrow=c(2,2),mar=c(3,4,1,1))
     for(i in c(6,7,9,10)){
        boxplot(as.numeric(st[i,-1]) ~ grp, ylab=st[i,1], names=c("Pax","F1","Pex"), xlab="")
     }
     dev.off()
     ```

* quality control of aligned reads with [MultiQC](https://github.com/ewels/MultiQC) version 1.8: `Alignment/mapReads_stats_multiQC.sh`


## 4. Process bam files

* index original STAR output bam files with samtools index: `SNPs/processBamFiles_index.sh`
* mark duplicates (and index) with picard-tools MarkDuplicates: `SNPs/processBamFiles.sh`
* split reads that contain Ns in their cigar string with GATK SplitNCigarReads: `SNPs/splitCigar_md.sh`, using manipulated reference genome with IUPAC ambiguity codes (e.g., *Petunia_axillaris_304_PexIUPAC.fa*, [generated above](#22-generate-genome-with-iupac-ambiguity-codes))  
  NOTE: setting maximum number of mismatches allowed in the overhang to 2 (STAR --outFilterMismatchNoverReadLmax 0.04 used above means in a 50:50 overhang could be on average 2 mismatches: --max-mismatches-in-overhang 2)

  &rarr; this also adjusts mapping quality from 255 to 60 and only keeps HI, nM, AS flags (index of multimapping read, number of mismatches, alignment score)


## 5. Call SNPs

* step 1: `SNPs/SNPs_call_md.sh`, removing read duplicates (picard-tools AddOrReplaceReadGroups and GATK HaplotypeCaller)
* step 2: `SNPs/SNPs_callCohortDB_md_array.sh` (GATK GenomicsDBImport and GATK GenotypeGVCFs)  
  NOTE: requires a huge amount of memory; splitting the genome into a few hundred [intervals](#23-split-genome-into-intervals) makes it manageable
* check whether all runs finished by grepping log files for "OutOfMemoryError" and "IllegalStateException"
* step 3: `SNPs/SNPs_callCohortDB_md_gather.sh` (picard-tools MergeVcfs)
* after finishing, remove all the interval_* directories, as they contain thousands of files


## 6. Filter variants

* filter variants using hard filters with GATK: `SNPs/SNPs_filter_md.sh` (SelectVariants, VariantFiltration, and VariantsToTable)  
  &rarr; generates *cohort_md_snps_filtered.tsv* and *cohort_md_snps_filtered_sel.vcf*

* make plot with R

     ```R
     st<-read.table("cohort_md_snps_filtered.tsv",header=TRUE)
     sts<-st$DP[st$DP<400]
     ##plot(density(sts,na.rm=TRUE),xlim=c(0,400),xlab="Read depth",main="RNA on P. ax")

     pdf("RNA_on_Pax.pdf")
     for(i in 8:13){
        plot(density(st[st[,i]>=quantile(st[,i],probs=0.05,na.rm=TRUE) & 
                     st[,i]<=quantile(st[,i],probs=0.95,na.rm=TRUE),i],na.rm=TRUE),
             xlim=c(quantile(st[,i],probs=c(0.05,0.95),na.rm=TRUE)),xlab=colnames(st)[i],
             main="RNA on P. ax")
     }
     dev.off()
     ```


## 7. Annotate SNPs

* using [ANNOVAR](https://academic.oup.com/nar/article/38/16/e164/1749458) (2019Oct24), following steps as described in the [user guide](https://annovar.openbioinformatics.org/en/latest/user-guide/gene/#what-about-gff3-file-for-new-species)

     ```bash
     ## use gtfToGenePred tool to convert GTF file to GenePred file
     $HOME/bin/gtfToGenePred -genePredExt \
                             $HOME/Style/Pax_genome/Peax304_final_edited.gtf \
                             Peax304_refGene.txt
     ## generate a transcript FASTA file with script provided by ANNOVAR
     perl $HOME/bin/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile \
                             $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
                             Peax304_refGene.txt --out Peax304_refGeneMrna.fa

     mkdir Peax304db
     mv Peax304_* Peax304db/

     ## modify file with variants
     grep "^Peax" $HOME/Style/step04_AxEx/Pax_genome/SNPs/cohort_md_snps_filtered_sel.vcf >tmp
     paste <(cut -f1,2 tmp) <(cut -f2 tmp) <(cut -f4,5 tmp) > snps_filtered_sel.txt
     rm -f tmp

     ## issue gene-based annotation
     perl $HOME/bin/annovar/annotate_variation.pl --outfile snps_filtered_sel_anno \
                                                  --buildver Peax304 \
                                                  snps_filtered_sel.txt Peax304db/
     ```

  &rarr; generates *snps_filtered_sel_anno.variant_function*


## 8. ASE and DE analysis

### 8.1. Data preparation and initial analysis

* calculate read counts per REF and ALT allele for each sample with GATK ASEReadCounter, not removing read duplicates: `ASE/countReadsASE.sh`
* filter SNPs and make infiles for GeneiASE from info for each triplet of samples per stage (two parental species and F1): `ASE/ASE_analysis_infiles.R`  
  &rarr; generates *input.tab*, *bias.txt*, and *SNPs.txt*
* run [GeneiASE](https://github.com/edsgard/geneiase) version 1.0.1 separately for each F1 hybrid: `ASE/runGeneiASE.sh`  
  &rarr; generates *geneiase.out*
* run [MBASED](http://bioconductor.org/packages/release/bioc/html/MBASED.html) version 1.20.0 separately for each F1 hybrid: `ASE/runMBASED.sh` (calls `ASE/runMBASED.R`)  
  &rarr; generates *mbased.out*
* check where reads that map to a *P. axillaris* gene map in the *P. exserta* genome, separately for each F1 hybrid: `ASE/readCorresp.sh` (calls `ASE/readCorresp.pl`) and `ASE/readCorresp_stats.sh` (calls `ASE/readCorresp_stats.R`)  
  &rarr; generates *mapReport.txt*
* merge bam files for all samples using picard-tools MergeSamFiles: `DE/mapReads_merge.sh`
* run featureCounts from the [Subread package](http://subread.sourceforge.net) version 1.6.4 to count reads per gene, not removing read duplicates: `DE/countReads.sh`  
  &rarr; generates *featureCounts_results.txt*
* quality control of counted reads with [MultiQC](https://github.com/ewels/MultiQC) version 1.8: `DE/countReads_stats_multiQC.sh`


### 8.2. In-depth analysis

* analyse results per stage with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) version 1.24.0 (using [apeglm](https://doi.org/10.1093/bioinformatics/bty895) shrinkage), [NOISeq](https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html) version 2.28.0, and custom functions: `ASE/ASE_DE_analysis.R`  
  &rarr; generates *\<stage\>_results.csv*
* combine results from parallel runs that either used the *P. axillaris* or the *P. exserta* genome: `ASE/ASE_DE_analysis_bidir.R`

* additionally, perform DE analysis per stage for all genes with DESeq2 (using [apeglm](https://doi.org/10.1093/bioinformatics/bty895) shrinkage): `DE/DE_analysis_quick.R`  
  &rarr; generates *results_\<stage\>.csv*
* additionally, compute correlation between gene expression and phenotypic traits across stages for each gene: `DE/correlation_analysis.R`






