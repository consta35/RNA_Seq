
setwd("E:\\PVN_Sequencing\\PVN_AllGen83\\F2_PVN_Tophat_Bam83")
library(GenomicFeatures)
library(GenomicRanges)

txdb=makeTxDbFromGFF("U_AllGenCuffmerge83_grep.gtf", format = "gtf")

# Use the function transcriptsBy(txdb,'gene') for the whole genic region instead of just exons 
ex_by_gene<- exonsBy(txdb,'gene')

# load the samtools library for R
library(Rsamtools)
library(GenomicAlignments)

## Increase memory limit to allow processes to occur 
memory.limit(size=800000000)
#read the F2 sequencing read alignment into R

readsF3Br1=readGAlignments("Tophat45_SGM21.bam")
readsF3Br2=readGAlignments("Tophat70_SGM26.bam")
readsF3Br3=readGAlignments("Tophat90_SGM30.bam")
readsF3Br4=readGAlignments("Tophat60_SGM24.bam")
readsF3Br5=readGAlignments("Tophat65_SGM25.bam")
readsF3Cr1=readGAlignments("Tophat50_SGM22.bam")
readsF3Cr2=readGAlignments("Tophat55_SGM23.bam")
readsF3Cr3=readGAlignments("Tophat75_SGM27.bam")
readsF3Cr4=readGAlignments("Tophat80_SGM28.bam")
readsF3Cr5=readGAlignments("Tophat85_SGM29.bam")

#count reads overlapping the exons
countsF3Br1= countOverlaps(ex_by_gene,readsF3Br1)
countsF3Br2= countOverlaps(ex_by_gene,readsF3Br2)
countsF3Br3= countOverlaps(ex_by_gene,readsF3Br3)
countsF3Br4= countOverlaps(ex_by_gene,readsF3Br4)
countsF3Br5= countOverlaps(ex_by_gene,readsF3Br5)
countsF3Cr1= countOverlaps(ex_by_gene,readsF3Cr1)
countsF3Cr2= countOverlaps(ex_by_gene,readsF3Cr2)
countsF3Cr3= countOverlaps(ex_by_gene,readsF3Cr3)
countsF3Cr4= countOverlaps(ex_by_gene,readsF3Cr4)
countsF3Cr5= countOverlaps(ex_by_gene,readsF3Cr5)


# create count table
countTable<- data.frame(Control1=countsF3Cr1,Control2=countsF3Cr2, Control3=countsF3Cr3, Control4=countsF3Cr4,Control5=countsF3Cr5, Beta1=countsF3Br1,Beta2=countsF3Br2,Beta3=countsF3Br3,Beta4=countsF3Br4,Beta5=countsF3Br5, stringsAsFactors=FALSE)

x <- rowSums(countTable<=0)!=ncol(countTable)
newCountTable <- countTable[x,]

write.table(newCountTable, file="F2PVN_Rawcount_zero_elim.csv", sep= ",")

data<-newCountTable

## Get PCA of data
pca<-prcomp(t(as.matrix(data)))
plot(pca$x)
text(pca$x[,1], pca$x[,2], colnames(data), pos=2)

plot(pca$rotation)
text(pca$rotation[,1], pca$rotation[,2], rownames(data), pos=2)

## run DESEQ2 to find outliers
library(DESEQ2)
data<-as.matrix(data)
pheno<- read.csv("pheno.csv", header=TRUE, row.names=1)


dds <- DESeqDataSetFromMatrix(countData = data,
colData = pheno,
design = ~ Treatment)

dds <- DESeq(dds)
res <- results(dds)
## find out how many outliers
summary(res)

#out of 20567 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 224, 1.1% 
#low counts [2]   : 0, 0% 
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
out of 12680 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 0, 0% 
LFC < 0 (down)   : 0, 0% 
outliers [1]     : 148, 1.2% 
low counts [2]   : 0, 0% 
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


## pull out Cooks Scores
maxCooks <- apply(assays(dds)[["cooks"]], 1, max)

write.table(maxCooks, "maxcooksf2pvn_round2.csv", sep=",")

## re-enter data after outliers elim
data<-read.csv("F3_PVN_GSEA_Ready_outliersrm.csv", header=TRUE, row.names=1)

### pull out genes that have at least 5 counts in at least 5 samples
filter <- apply(data, 1, function(x) length(x[x>5])>=5)
filtered <- data[filter,]

write.table(filtered, "F2_PVN_filtered.csv", sep=",")

### Use Join to get gene names this is done in linux terminal
cd /scratch/m/matthew7/consta35/PVN_Sequencing/F2_PVNSeq

dos2unix 83_GeneIDsF2PVN.txt F2_PVN_filteredAC.txt

join <(sort 83_GeneIDsF2PVN.txt) <(sort F2_PVN_filteredAC.txt) > F2PVN_filtered_wnames.txt


### Read output from Join back into R. In this data, we will have multiple rows with the same gene name- for ex. we will have multiple rows with counts for Y_RNA, and Sno. Melt and Cast allow us to sum all the counts that have the same gene name.
library(reshape)
library(reshape2)
data<- read.csv("F2PVN_filtered_wnames1.csv", header=TRUE)
md<-melt(data, id= "Gene")
cd<-cast(md,formula=Gene~variable,sum)

write.table(cd,file="cd_F2PVNraw_zeroelimnames.csv",sep=",")
### After we sum all the counts of duplicated gene names, I then manually (in Excel)remove lncRNA and Ensemble genes (these are genes that are not yet fully annotated in the Guinea pig Genome)

## re-enter data after outliers elim, and lncRNA and Ensemble genes are removed
data<-read.csv("cd_F2PVNraw_zero_lncRna_ens_Outliers_elimnames.csv", header=TRUE, row.names=1)

###RUV setup
library(RUVSeq)
x <- as.factor(c(1,1,1,1,1,2,2,2,2,2))
set <- newSeqExpressionSet(as.matrix(filtered),
			phenoData = data.frame(x, row.names=colnames(filtered)))
genes <- rownames(filtered)		
libSizes <- as.vector(colSums(filtered))	

### RUV visualization 
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")

### Plots
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
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

plotPCA(set1, col=colors[x], cex=1.2)		
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])

### EdgeR analysis on RUVr data
design<- model.matrix(~x + W_1, data=pData(set1))
d <- DGEList(counts=counts(set1),group=x)
d <- calcNormFactors(d, method="upperquartile")
d <- estimateGLMCommonDisp(d, design,verbose=T)
d <- estimateGLMTagwiseDisp(d, design)
#Disp = 0.01872 , BCV = 0.1368 
fit <- glmFit(d, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

summary(decideTestsDGE(lrt, adjust.method= "fdr", p.value=0.05))
    [,1] 
-1    84
0  12008
1    115

## Make processed data file
abc<-as.data.frame(normCounts(set1))
abc$twd <- d$tagwise.dispersion
abc<- cbind(abc, lrt$table)
abc$PValue_fdr <- p.adjust(method="fdr",p=abc$PValue)
write.table(abc, file="F2EDGER_RUVr__LRT_OutliersElim_NORMcounts_PVN_k1.csv", sep=",",quote=F)

# Make a basic volcano plot
with(abc, plot(logFC, -log10(PValue), pch=20, main="F2 PVN Gene Expression k2", xlim=c(-4,4), ylim=c(0,15)))

# Add colored points: red if PValue_fdr<0.05)
with(subset(abc, PValue_fdr<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))

###heatmap
vplot_sig<-subset(abc, PValue_fdr<.05 )
vvplotsig<-vplot_sig[,1:10]
vvplotsig<-as.matrix(vvplotsig)

lmno <- cpm(vvplotsig, prior.count=2, log=TRUE)

heatmap3(lmno, Rowv = NULL, Colv = NA,
  distfun = function(x) as.dist(1 - cor(t(x), use = "pa")),
  balanceColor = F,showColDendro = F,
  showRowDendro = F, col = greenred(1050), legendfun, method = "complete", ColAxisColors = 0,
  RowAxisColors = 0, hclustfun = hclust, reorderfun = function(d, w)
  reorder(d, w), add.expr, symm = FALSE,
  scale = "row", na.rm = TRUE, 
  ColSideWidth = 0.4, 
  file = "heatmap3.pdf", topN = NA, filterFun = sd, margins = c(5, 5), lasRow = 2, lasCol = 2,
  labRow =NA, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL,
  keep.dendro = FALSE, verbose = getOption("verbose"))







