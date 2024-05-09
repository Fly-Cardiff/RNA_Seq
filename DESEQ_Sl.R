workDir <- "SET ME"
setwd(workDir)
devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")

# Load packages 

library("DESeq2")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("SARTools")
library("RColorBrewer")
library("pheatmap")

# Determining the files for markdup (not removed)

sample_ID <- c("B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4")
RNAi <- factor(c("KD", "KD", "KD", "KD", "Con", "Con", "Con", "Con"),levels=c("Con","KD"))
filenames <- c("B1_S5_merge.markdup.featurecount","B2_S6_merge.markdup.featurecount","B3_S7_merge.markdup.featurecount","B4_S8_merge.markdup.featurecount","C1_S9_merge.markdup.featurecount","C2_S10_merge.markdup.featurecount","C3_S11_merge.markdup.featurecount","C4_S12_merge.markdup.featurecount")
sampleTable = data.frame(sample_ID,filenames,RNAi)
colData_Exp = data.frame(RNAi)

#These are for the analysis without C4 - as outlier in the PCA's

sample_ID <- c("B1", "B2", "B3", "B4", "C1", "C2", "C3")
RNAi <- factor(c("KD", "KD", "KD", "KD", "Con", "Con", "Con"),levels=c("Con","KD"))
filenames <- c("B1_S5_merge.markdup.featurecount","B2_S6_merge.markdup.featurecount","B3_S7_merge.markdup.featurecount","B4_S8_merge.markdup.featurecount","C1_S9_merge.markdup.featurecount","C2_S10_merge.markdup.featurecount","C3_S11_merge.markdup.featurecount")
sampleTable = data.frame(sample_ID,filenames,RNAi)
colData_Exp = data.frame(RNAi)

#removing ncRNA and transposables to enrich for mRNA
featuresToRemove <- as.character(read.table("PATH/ncRNAs_insertions_transposable.txt", stringsAsFactors = FALSE))
counts <- loadCountData(target=sampleTable, rawDir=workDir, featuresToRemove=featuresToRemove)

#Assign dds
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData_Exp, design = ~RNAi, tidy=FALSE, ignoreRank = FALSE)

# Run DESeq
dds <-DESeq(dds, test="Wald", fitType = "parametric")

#View sizefactors used 
sizeFactors(dds)

######################### RESULTS #################
resultsNames(dds)

#W/o complicated design 
dds$RNAi
res = results(dds,independentFiltering = TRUE, alpha  = 0.05, pAdjustMethod = "BH")
summary(res)

#Finding those that are significant using a cut off of 1.5 fold which is 0.58 log2fold change  
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.05)
Up <- subset(resSig, log2FoldChange>0.58)
Down <- subset(resSig, log2FoldChange< -0.58)

#PCA and dendrograms 
VSTdds <-vst(dds, blind = FALSE) #Variance stabilisation transformation 
sampleDists <- dist(t(assay(VSTdds)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(VSTdds))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#PCA
rv <- rowVars(assay(VSTdds))
select <- order(rv, decreasing=TRUE)
pca_results <- prcomp(t(assay(VSTdds)[select,]))
PC1=pca_results$x[,1]
PC2=pca_results$x[,2]
PC3=pca_results$x[,3]
Genotype <- dds$RNAi
PCA_three <- data.frame(PC1,PC2,PC3,Genotype)

#Getting eigenvalues
eigenvalues = pca_results$sdev^2
eigenvalues = eigenvalues/sum(eigenvalues)
percentVar = round(eigenvalues*100)

#Scree plot
barplot(percentVar,xlab="Principal Component",ylab="Percent variance (%)",names.arg=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15"),ylim=c(0,50),las=2)


#PCA

p1<- ggplot(PCA_three, aes(x=PC1, y=PC2, col=Genotype), label=rownames(PCA_three)) + xlim(-40,40) + ylim(-40,40) + geom_point(size=2.5) + geom_text_repel(label=rownames(PCA_three),size = 3.5) +
  xlab(paste0("PC1: ",percentVar[1],"% Variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% Variance")) 
p2 <- p1 + theme_light() + theme(axis.text= element_text(family="Trebuchet MS", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p3 <- p2 + theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=11)) + theme(panel.grid = element_line(size = 0.75, linetype = 2))
p3 + scale_color_manual(values=c("black","#006600")) 

p1<- ggplot(PCA_three, aes(x=PC1, y=PC3, col=Genotype), label=rownames(PCA_three)) + xlim(-40,40) + ylim(-40,40) + geom_point(size=2.5) + geom_text_repel(label=rownames(PCA_three),size = 3.5) +
  xlab(paste0("PC1: ",percentVar[1],"% Variance")) +
  ylab(paste0("PC3: ",percentVar[2],"% Variance")) 
p2 <- p1 + theme_light() + theme(axis.text= element_text(family="Trebuchet MS", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p3 <- p2 + theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=11)) + theme(panel.grid = element_line(size = 0.75, linetype = 2))
p3 + scale_color_manual(values=c("black","#006600")) 

p1<- ggplot(PCA_three, aes(x=PC2, y=PC3, col=Genotype), label=rownames(PCA_three)) + xlim(-40,40) + ylim(-40,40) + geom_point(size=2.5) + geom_text_repel(label=rownames(PCA_three),size = 3.5) +
  xlab(paste0("PC2: ",percentVar[1],"% Variance")) +
  ylab(paste0("PC3: ",percentVar[2],"% Variance")) 
p2 <- p1 + theme_light() + theme(axis.text= element_text(family="Trebuchet MS", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())
p3 <- p2 + theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=11)) + theme(panel.grid = element_line(size = 0.75, linetype = 2))
p3 + scale_color_manual(values=c("black","#006600")) 


#Volcano plot

#For scaling
max(Up$log2FoldChange, na.rm = T)
min(Down$log2FoldChange, na.rm = T)
max(-log10(Up$padj), na.rm = T)
max(-log10(Down$padj), na.rm = T)

#Dispersion estimates
plotDispEsts(dds)

#Plotting Volcano
plotdf <- data.frame(res$log2FoldChange, res$padj)
rownames(plotdf) = rownames(res) #Row names are gene names 
plotdf$diffex<-"Non-Sig"
plotdf$diffex[res$log2FoldChange > 0.58 & res$padj < 0.05] <- "Up"
plotdf$diffex[res$log2FoldChange < -0.58 & res$padj < 0.05] <- "Down"

p <- ggplot(plotdf, aes(x=res$log2FoldChange, y=-log10(res$padj), col=diffex)) + geom_point(size=0.6) + theme_classic() + xlim(-10,10) + ylim(0,50)
p2<- p + geom_vline(xintercept=c(-0.58,0.58), col="grey", lty="longdash") + geom_hline(yintercept=-log10(0.05), col="grey",lty="longdash")
p3 <- p2 + scale_color_manual(values=c("red", "black", "green"))
p4 <- p3 + labs(y="-Log10 Pqdj",x="Log2 Fold Change") + theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=11)) 
p4 + theme(axis.text= element_text(family="Trebuchet MS", face = "bold", colour = "black", size = 10)) + theme(legend.title =element_blank())

                     