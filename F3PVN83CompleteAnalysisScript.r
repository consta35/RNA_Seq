setwd ("E:\\PVN_Sequencing\\PVN_AllGen83\\F3_PVN_TOPHAT_BAM83")
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)

#Read in bam files
readsF3Br1=readGAlignments("F3Beta_564_SGM10_accepted_hits.bam")
readsF3Br2=readGAlignments("F3Beta_604_SGM18_accepted_hits.bam")
readsF3Br3=readGAlignments("F3Beta_614_SGM2_accepted_hits.bam")
readsF3Br4=readGAlignments("F3Beta_639_SGM6_accepted_hits.bam")
readsF3Cr1=readGAlignments("F3Control_579_SGM13_accepted_hits.bam")
readsF3Cr2=readGAlignments("F3Control_584_SGM14_accepted_hits.bam")
readsF3Cr3=readGAlignments("F3Control_599_SGM17_accepted_hits.bam")
readsF3Cr4=readGAlignments("F3Control_634_SGM5_accepted_hits.bam")
readsF3Cr5=readGAlignments("F3Control_654_SGM9_accepted_hits.bam")
readsF3Cr6=readGAlignments("F3Control_559_SGM1_accepted_hits.bam")

##New reads
readsF3Br5=readGAlignments("F3Beta_32_SGM21_3accepted_hits.bam")
readsF3Br6=readGAlignments("F3Beta_37_SGM22_1accepted_hits.bam")
readsF3Cr41=readGAlignments("F3Control_22_SGM5_1accepted_hits.bam")
readsF3Cr51=readGAlignments("F3Control_27_SGM9_1accepted_hits.bam")

#Make exon file
txdb=makeTxDbFromGFF("U_AllGenCuffmerge83_grep.gtf", format = "gtf")
ex_by_gene<- exonsBy(txdb,'gene')


#Get Count data 
countsF3Br1= countOverlaps(ex_by_gene,readsF3Br1)
countsF3Br2= countOverlaps(ex_by_gene,readsF3Br2)
countsF3Br3= countOverlaps(ex_by_gene,readsF3Br3)
countsF3Br4= countOverlaps(ex_by_gene,readsF3Br4)
countsF3Br5= countOverlaps(ex_by_gene,readsF3Br5)
countsF3Br6= countOverlaps(ex_by_gene,readsF3Br6)
countsF3Cr1= countOverlaps(ex_by_gene,readsF3Cr1)
countsF3Cr2= countOverlaps(ex_by_gene,readsF3Cr2)
countsF3Cr3= countOverlaps(ex_by_gene,readsF3Cr3)
countsF3Cr4= countOverlaps(ex_by_gene,readsF3Cr4)
countsF3Cr5= countOverlaps(ex_by_gene,readsF3Cr5)
countsF3Cr6= countOverlaps(ex_by_gene,readsF3Cr6)
countsF3Cr41= countOverlaps(ex_by_gene,readsF3Cr41)
countsF3Cr51= countOverlaps(ex_by_gene,readsF3Cr51)


##Make count table

countTable<- data.frame(Control1=countsF3Cr1,Control2=countsF3Cr2, Control3=countsF3Cr3, Control4=countsF3Cr4,Control5=countsF3Cr5, Control6=countsF3Cr6, Control41=countsF3Cr41, Control51=countsF3Cr51, Beta1=countsF3Br1,Beta2=countsF3Br2,Beta3=countsF3Br3,Beta4=countsF3Br4,Beta5=countsF3Br5, Beta6=countsF3Br6,  stringsAsFactors=FALSE)

### Remove rows with zero counts- all analysis techniques cannot handle genes with zero counts
x <- rowSums(countTable<=0)!=ncol(countTable)
newCountTable <- countTable[x,]

write.table(newCountTable, file="F3PVN_Rawcount_zero_elim.csv", sep= ",")

### replace XLOC codes with gene names by using join function in linux shell
cd /scratch/m/matthew7/consta35/F3_PVN_TOPHAT_BAM83

dos2unix F3PVN_geneIDs.txt F3PVN_Rawcount_zero_elim.txt

join <(sort F3PVN_geneIDs.txt) <(sort F3PVN_Rawcount_zero_elim.txt) > F3PVN_Rawcount_zero_elim_NAMEs.txt

### sum duplicated gene name values
library(reshape)
library(reshape2)
counts<- read.csv("F3PVN_Rawcount_zero_elim_NAMEs.csv", header=TRUE)
md<-melt(counts, id= "Gene")
cd<-cast(md,formula=Gene~variable,sum)

write.table(cd,file="cd_F3PVN_Rawcount_zero_elim_NAMEs.csv",sep=",")

setwd("E:\\PVN_Sequencing\\PVN_AllGen83\\F3_PVN_TOPHAT_BAM83")

## re-enter data after separating out original reads and resequenced reads into 2 separate files
dataori<-read.csv("cd_F3PVN_Rawcount_zero_elim_lncRNA_EnsRM_NAMEs_Original.csv", header=TRUE, row.names=1)
datareseq<-read.csv("cd_F3PVN_Rawcount_zero_elim_lncRNA_EnsRM_NAMEs_reseq.csv", header=TRUE, row.names=1)

### NORMALIZATION FOR ORIGNAIL
### pull out genes that have at least 5 counts in at least 5 samples
filter <- apply(dataori, 1, function(x) length(x[x>5])>=5)
filtered <- dataori[filter,]

###RUV Upperquantile normalization. Combat requires 'clean and normalized data' in order to work properly. So, original and resequenced data were normalized separately before being joined and corrected with ComBat.
library(RUVSeq)
x <- as.factor(c(1,1,1,1,1,1,2,2,2,2))
set <- newSeqExpressionSet(as.matrix(filtered),
			phenoData = data.frame(x, row.names=colnames(filtered)))
genes <- rownames(filtered)		
libSizes <- as.vector(colSums(filtered))	

## Upperquantile normalization
set <- betweenLaneNormalization(set, which="upper")

abc<- normCounts(set)
write.table(abc, file="F3PVN_NormCounts_ORIGINAL_UQ.csv", sep=",",quote=F)

### NORMALIZATION FOR RESEQUENCED
## filtering here is for minimum of 5 counts in 2 samples, as there are 2 controls and 2 betas. 
filter <- apply(datareseq, 1, function(x) length(x[x>5])>=2)
filtered <- datareseq[filter,]

###RUV Upperquantile normalization. Combat requires 'clean and normalized data' in order to work properly. So, original and resequenced data were normalized separately before being joined and corrected with ComBat.
library(RUVSeq)
x <- as.factor(c(1,1,1,1,1,1,2,2,2,2))
set <- newSeqExpressionSet(as.matrix(filtered),
			phenoData = data.frame(x, row.names=colnames(filtered)))
genes <- rownames(filtered)		
libSizes <- as.vector(colSums(filtered))	

## Upperquantile normalization
set <- betweenLaneNormalization(set, which="upper")

### RUV visualization 
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")

### Plots
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)	

abc<- normCounts(set)
write.table(abc, file="F3PVN_NormCounts_RESEQ_UQ.csv", sep=",",quote=F)

### JOIN norm counts this was done in linux terminal
cd /scratch/m/matthew7/consta35/F3_PVN_TOPHAT_BAM83

dos2unix F3PVN_NormCounts_ORIGINAL_UQ.txt F3PVN_NormCounts_RESEQ_UQ.txt

join <(sort -t'|' F3PVN_NormCounts_RESEQ_UQ.txt) <(sort -t'|' F3PVN_NormCounts_ORIGINAL_UQ.txt)> F3PVN_NormCounts_JOINED_UQ.txt

## COMBAT on JOINED COUNTS
datajoined<-read.csv("F3PVN_NormCounts_JOINED_UQ1.csv", header=TRUE, row.names=1)

pheno<- read.csv("pheno.csv", header=TRUE, row.names=1)
library(sva)
batch = pheno$Batch
modcombat = model.matrix(~Treatment, data=pheno)
combat_data = round(ComBat(dat=datajoined, batch=batch, mod=modcombat, par.prior=TRUE, prior.plot=TRUE), 0)

## combat makes it so that 4 graphs are shown in a single instance, instead of one at a time. This code [par(mfrow=c(1,1))] corrects this and makes it 1 image per instance again
par(mfrow=c(1,1))

## Get PCA of ComBat corrected Data
pca<-prcomp(t(as.matrix(combat_data)))
plot(pca$x)
text(pca$x[,1], pca$x[,2], colnames(combat_data), pos=2)

### ComBat will give some genes negative count values (unclear as to why) possibly due to low gene count before correction. Further analyses cannot handle negative gene counts and therefore data must be filtered to have at least 1 count per sample. 
filter <- apply(combat_data, 1, function(x) length(x[x>1])>=14)
data <- combat_data[filter,]

write.table(data, file="F3PVN_NormCounts_JOINED_Combat.csv", sep=",",quote=F)

### Remove Technical replicates. since Technical replicates in control are still different enough from the original sequencing (1 year old RNA, completely different library prep, sequenced separately, longer sequencing length, greater sequencing depth) it does not make sense to average these data into the original sequencing data. Therefore, the technical replicates were removed after batch correction. 
data<-read.csv("F3PVN_NormCounts_JOINED_Combat_TR_RM.csv", header=TRUE, row.names=1)

## RUVSeq
library(RUVSeq)
x <- as.factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
set <- newSeqExpressionSet(as.matrix(data),
			phenoData = data.frame(x, row.names=colnames(data)))
genes <- rownames(data)		
libSizes <- as.vector(colSums(data))

### Plots
plotRLE(set, outline=FALSE, ylim=c(-6, 6), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)	

### RUVr
design<- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y,method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
set1<- RUVr(set, genes, k=2, res)
## k=2 is used in F3 to account for the fact that 2 samples in Beta are sequenced more in depth and a full year later- adding variability to overall results
plotPCA(set1, col=colors[x], cex=1.2)		
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])

## EDGER analysis
design<- model.matrix(~x + W_1+ W_2, data=pData(set1))
d <- DGEList(counts=counts(set1),group=x)
d <- calcNormFactors(d, method="upperquartile")
d <- estimateGLMCommonDisp(d, design,verbose=T)
#Disp = 0.00476 , BCV = 0.069 
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

summary(decideTestsDGE(lrt, adjust.method= "fdr", p.value=0.05))
   [,1] 
-1    31
0  12026
1     11

## Make processed data file
abc<-as.data.frame(normCounts(set1))
abc$twd <- d$tagwise.dispersion
abc<- cbind(abc, lrt$table)
abc$PValue_fdr <- p.adjust(method="fdr",p=abc$PValue)
write.table(abc, file="F3EDGER_COUNTS_JOINED_COMBAT_RUVr_TechRepRM__NORMCOUNTS_k2.csv", sep=",",quote=F)

# Make a basic volcano plot
with(abc, plot(logFC, -log10(PValue), pch=20, main="F3 PVN Gene Expression k2", xlim=c(-4,4), ylim=c(0,7)))

# Add colored points: red if PValue_fdr<0.05)
with(subset(abc, PValue_fdr<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))

###heatmap
vplot_sig<-subset(abc, PValue_fdr<.05 )
vvplotsig<-vplot_sig[,1:12]
vvplotsig<-as.matrix(vvplotsig)

lmno <- cpm(vvplotsig, prior.count=2, log=TRUE)


library(heatmap3)
library(ggplot2)
library(gplots)
heatmap3(lmno, Rowv = NULL, Colv = NA,
  distfun = function(x) as.dist(1 - cor(t(x), use = "pa")),
  balanceColor = F,showColDendro = F,
  showRowDendro = F, col = greenred(1050), legendfun, method = "complete", ColAxisColors = 0,
  RowAxisColors = 0, hclustfun = hclust, reorderfun = function(d, w)
  reorder(d, w), add.expr, symm = FALSE,
  scale = "row", na.rm = TRUE, 
  ColSideWidth = 0.4, 
  file = "heatmap3.pdf", topN = NA, filterFun = sd, margins = c(7, 7), lasRow = 2, lasCol = 2,
  labRow =NA, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL,
  keep.dendro = FALSE, verbose = getOption("verbose"))

































