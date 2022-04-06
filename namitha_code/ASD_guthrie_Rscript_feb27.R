#Autism Spectrum Disorder Methylation analysis - Guthrie cards - at birth (Within twin pair analysis - 11 twin pairs)

#Setting base directory and loading libraries required for analysis
baseDir=("/group/canc2/puumba/Data/InfiniumData/JeffCraig/ASD_EPIC_DATA")
library(missMethyl)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(ruv)
library(RColorBrewer)
library(matrixStats)
library(gplots)
library(WGCNA)
library(FlowSorted.Blood.450k)

#Reading in annotation for EPIC methylation arrays
ann = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

#Reading in the sample summary sheet (metadata)
targets_gen = read.metharray.sheet(baseDir, pattern = "ASD_guthrie_summarysheet.csv")

#Reading array raw data and assigning sample names with array data
#targets$ID = paste(targets$Sample_Group,targets_gen$Sample_Name,sep=".")
rgSet = read.metharray.exp(targets = targets_gen)
sampleNames(rgSet) = targets_gen$SampleID

#Testing for MZ status
snpBetas = getSnpBeta(rgSet)
d = dist(t(snpBetas))
hr = hclust(d, method = "complete", members=NULL)
plot(hr)

#Quality control
detP = detectionP(rgSet)
qcReport(rgSet, sampNames = targets_gen$Sample_Name, pdf="qc-report_ASD_guthrie_Apr8.pdf")
pdf("mean_detection_ASDguthrie_Apr8.pdf",width=14)
par(mfrow=c(1,2))
cols=brewer.pal(4,"Set1")
barplot(apply(detP,2,mean),col=as.numeric(factor(targets_gen$Sample_Name)),las=2,cex.names= 0.5, cex.axis=0.75,main="Mean detection p-values of probe signals",ylab="Mean detection p-value")
barplot(apply(detP,2,mean),col=as.numeric(factor(targets_gen$Sample_Name)),las=2,cex.names= 0.5, cex.axis=0.75,ylim=c(0,0.010),main="Mean detection p-values of probe signals",ylab="Mean detection p-value")
dev.off()

#Preprocessing
mset.raw = preprocessRaw(rgSet)

#Data exploration
#Multi-dimensional scaling (MDS) plots before filtering
pdf("mds_plots_ASDguthrie_Apr8.pdf",width=14)
par(mfrow=c(1,2))
mdsPlot(mset.raw, sampGroups = targets_gen$Sample_Name, sampNames=targets_gen$Social_interaction_on_ADOS,legendPos="bottom")
mdsPlot(mset.raw, sampGroups = targets_gen$Sex, sampNames=targets_gen$SampleID,legendPos="bottom")
dev.off()

#Normalisation
mSetSw = preprocessSWAN(rgSet)
densityPlot(getBeta(mSetSw),main="SWAN")

#SQN method
mSetSq = preprocessQuantile(rgSet)

par(mfrow=c(1,3))
densityPlot(rgSet, sampGroups = targets_gen$Social_interaction_on_ADOS,main="Raw", legend = FALSE)
densityPlot(getBeta(mSetSw), sampGroups = targets_gen$Social_interaction_on_ADOS,main="SWAN", legend = FALSE)
densityPlot(getBeta(mSetSq), sampGroups = targets_gen$Social_interaction_on_ADOS,main="SQN", legend = FALSE)

#Cell type composition analysis
library(FlowSorted.Blood.450k)
rgSet$Slide <- as.numeric(rgSet$Slide)
rgSet$Sex<- as.character(rgSet$Sex)
rgSet$Sample_Name<- as.character(rgSet$Sample_Name)
cellCounts_new <- estimateCellCounts(rgSet, compositeCellType = "Blood", processMethod = "auto", probeSelect = "auto", cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"), referencePlatform = c("IlluminaHumanMethylation450k"),returnAll = FALSE, meanPlot = TRUE)
#plot cell type composition by sample group
par(mfrow=c(1,1))
a = cellCounts_new[targets_gen$Diagnosis..L.0..M.1..H.2. == "0",]
b = cellCounts_new[targets_gen$Diagnosis..L.0..M.1..H.2. == "1",]
c = cellCounts_new[targets_gen$Diagnosis..L.0..M.1..H.2. == "2",]
age.pal <- brewer.pal(8,"Set1")
boxplot(a, at=0:5*3 + 1, xlim=c(0, 18), ylim=range(a,c), xaxt="n",
        col=age.pal[1], main="", ylab="Cell type proportion")
boxplot(b, at=0:5*3 + 2, xaxt="n", add=TRUE, col=age.pal[2])
axis(1, at=0:5*3 + 1.5, labels=colnames(a), tick=TRUE) 
boxplot(c, at=0:5*3 + 3, xaxt="n", add=TRUE, col=age.pal[3])
axis(1, at=0:5*3 + 1.5, labels=colnames(a), tick=TRUE) 
legend("topleft", legend=c("Low","Moderate", "High"), fill=age.pal, cex=.7)

wilcox.test(b[,"Gran"],c[,"Gran"], paired=FALSE)
#None of the p-values are less than 0.05 - separate word document in folder

#Filtering low quality probes
keepProbes = rowSums(detP < 0.01) == ncol(detP)
mSetSqFlt = mSetSq[keepProbes,]
gmSetSqFlt = mapToGenome(mSetSqFlt)

#remove SNPs
gmSetSqFlt = dropLociWithSnps(gmSetSqFlt, snps = c("CpG", "SBE"))

#Removing cross-reactive probes
Xreact = read.csv(file="/group/canc2/puumba/Data/InfiniumData/NamithaData/Rprojects/Autism/Guthrie_analysis/EPIC_850k_crossreactiveProbes.csv", stringsAsFactors=FALSE)
#Xreact = read.csv(file="~/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)
noXreact = !(featureNames(gmSetSqFlt) %in% Xreact$TargetID)
gmSetSqFlt = gmSetSqFlt[noXreact,]

#Removing probes on X and Y chromosomes
autosomes = !(featureNames(gmSetSqFlt) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
gmSetSqFlt = gmSetSqFlt[autosomes,]

#Relative log expression (RLE plot)
mValsSw = getM(gmSetSqFlt)
medSw = apply(mValsSw, 1, median)
YSw = mValsSw - medSw
par(mfrow=c(1,2))
pdf("rle_plot_Apr8.pdf")
boxplot(YSw,outline=FALSE,ylim=c(-1.5,1.5), ylab="Relative Log Methylation Value", cols=as.character(factor(targets_gen$Social_interaction_on_ADOS,)),xaxt="none")
title(xlab="Samples",cex=2, line=1)
dev.off()

#MDS plots generation after filtering
pal = brewer.pal(8, "Dark2")
mds1Sw = plotMDS(mValsSw, top=1000, gene.selection="common",dim.plot=c(1,2))
mds2Sw = plotMDS(mValsSw, top=1000, gene.selection="common",dim.plot=c(1,3))
mds3Sw = plotMDS(mValsSw, top=1000, gene.selection="common",dim.plot=c(2,3))
mds4Sw = plotMDS(mValsSw, top=1000, gene.selection="common",dim.plot=c(3,4))
pdf("MDS_plot_normalised_ASDguthrie_Apr8_language.pdf",height=14,width=14)
par(mfrow=c(2,2))
plotMDS(mds1Sw, xlab="Dimension 1", ylab="Dimension 2",col=pal[as.factor(targets_gen$language)],dim=c(1,2), labels=targets_gen$SampleID)
legend("bottomright",bg="white",col=pal,cex=.7,pch=1,legend=0:1)
plotMDS(mds2Sw, xlab="Dimension 1", ylab="Dimension 3",col=pal[as.factor(targets_gen$language)],dim=c(1,3), labels=targets_gen$SampleID)
legend("bottomright",bg="white",col=pal,cex=.7,pch=1,legend=0:1)
plotMDS(mds3Sw, xlab="Dimension 2", ylab="Dimension 3",col=pal[as.factor(targets_gen$language)],dim=c(2,3), labels=targets_gen$SampleID)
legend("bottomright",bg="white",col=pal,cex=.7,pch=1,legend=0:1)
plotMDS(mds4Sw, xlab="Dimension 3", ylab="Dimension 4",col=pal[as.factor(targets_gen$language)],dim=c(3,4), labels=targets_gen$SampleID)
legend("bottomright",bg="white",col=pal,cex=.7,pch=1,legend=0:1)
dev.off()

#Principal Component Analysis (PCA)
fit <- prcomp(t(mValsSw),center = TRUE, scale = TRUE,retx=TRUE)
loadings = fit$x
plot(fit,type="lines")
nGenes = nrow(mValsSw)
nSamples = ncol(mValsSw)
datTraits = targets_gen[,7:15]
moduleTraitCor = cor(loadings[,1:6], datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
par(cex=0.75, mar = c(6, 8.5, 3, 3))
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("PCA_heatmap_ASDblood_Apr8.pdf",height=8,width=8)
labeledHeatmap(Matrix = t(moduleTraitCor), xLabels = colnames(loadings), yLabels = names(datTraits), colorLabels = FALSE, colors = blueWhiteRed(6), textMatrix = t(textMatrix), setStdMargins = FALSE, cex.text = 0.5, cex.lab.y = 0.6, zlim = c(-1,1), main = paste("PCA-trait relationships: Top principal components"))
dev.off()

#-------------------------------------------------------------#
#Regression analysis - within twin pair on ADOS
mDat = mValsSw
Twin_Group <- factor(targets_gen$Family_ID)
ADOS=targets_gen$Social_interaction_on_ADOS
design_tw <- model.matrix(~ADOS+Twin_Group)
#design <- model.matrix(~targets_gen$Social_interaction_on_ADOS)
fit_tw <- lmFit(mDat, design_tw)
fit2_tw <- eBayes(fit_tw)
summary(decideTests(fit2_tw,p.value = 0.05))
top = topTable(fit2_tw,coef="ADOS",num=Inf, sort.by = "P")
head(top)
nsig <- sum(top$adj.P.Val < 0.05)
sum(top$P.Value< 0.05)
ann_sub = 
  ann[,c("chr","pos","strand","Name","Islands_Name","Relation_to_Island","UCSC_RefGene_Name","UCSC_RefGene_Group")]
output = merge(ann_sub,top,by.x="Name",by.y="row.names")
write.csv(output[order(output$P.Value),],file="ASD_guthrie_top_dmps_onADOS_withintw_Apr8_limma.csv",row.names=FALSE)

top1000 = topTable(fit2_tw,coef="ADOS",num=1000, sort.by = "P")
output_1000 = merge(ann_sub,top1000,by.x="Name",by.y="row.names")
write.csv(output_1000[order(output_1000$P.Value),],file="ASD_guthrie_top_dmps_onADOS_Apr8_withintw_limma_top1k.csv",row.names=FALSE)

#Converting to beta from M values
bDat = ilogit2(mDat)
bDat_new = getBeta(gmSetSqFlt)
#View(bDat)
write.csv(bDat,file="ASD_guthrie_beta_onADOS_withintw_Apr8.csv",row.names=TRUE)

#RUV
#extract illumina negative control
INCs <- getINCs(rgSet)
head(INCs)
#add negative control data to M-values
Mc <- rbind(mValsSw,INCs)
ctl <- rownames(Mc) %in% rownames(INCs)
table(ctl)

rfit1 <- RUVfit(data=Mc, design=design_tw, coef=2, ctl=ctl) # Stage 1 analysis
rfit2 <- RUVadj(rfit1)
top1 <- topRUV(rfit2, num=Inf)
head(top1)
#ctl <- rownames(mValsSw) %in% rownames(top1[top1$p.ebayes > 0.5,])
#bottom 80% used as negative controls
ctl <- rownames(mValsSw) %in% rownames(top1[ceiling(0.2*nrow(top1)):nrow(top1),])
table(ctl)
# Perform RUV adjustment and fit
rfit1 <- RUVfit(data=mValsSw, design=design_tw, coef=2, ctl=ctl) # Stage 2 analysis 
rfit2 <- RUVadj(rfit1)
top_ruv = topRUV(rfit2, number = Inf)
top_ruv_1000 = topRUV(rfit2, number = 1000)
head(top_ruv_1000)
#dat_subset_ruv = as.matrix(subset(mValsSw, rownames(mValsSw) %in% rownames(top_ruv_1000)))
#heatmap.2(dat_subset_ruv, cexRow=0.8, cexCol=0.8, col=colorRampPalette(c("Blue","Yellow")), density.info="none", trace="none", scale="row", main="diff meth probes due to epilepsy")
output_ruv = merge(ann,top_ruv,by.x="Name",by.y="row.names")
output_ruv_1000 = merge(ann,top_ruv_1000,by.x="Name",by.y="row.names")
#output_full_epicanno = merge(ann_epic,top,by.x="Name",by.y="row.names")
write.csv(output_ruv[order(output_ruv$p.ebayes),],file="ASD_guthrie_topruv80%_dmps_onADOS_withintw_Apr8.csv",row.names=FALSE)
write.csv(output_ruv_1000[order(output_ruv_1000$p.ebayes),],file="ASD_guthrie_topruv80%_1kdmps_onADOS_withintw_Apr8.csv",row.names=FALSE)

#To visualise the effect of RUV adjustment
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("missMethyl", version = "3.8")

##to update bioconductor - from Jovana (jane's email)
install.packages("BiocManager")
BiocManager::version()
BiocManager::valid(lib.loc = .libPaths()[1])
BiocManager::install(lib.loc = .libPaths()[1])

#Madj <- getAdjusted(mValsSw, rfit2)
Madj <- missMethyl::getAdjusted(mValsSw, rfit2)

par(mfrow=c(1,2))
plotMDS(mValsSw, labels=targets_gen$Sample_Name, col=as.integer(factor(targets_gen$Social_interaction_on_ADOS)),
        main="Unadjusted")
#legend(legend=c("1","0"),pch=7,cex=0.7,col=1:2)

plotMDS(Madj, labels=targets_gen$Sample_Name, col=as.integer(factor(targets_gen$Social_interaction_on_ADOS)),
        main="Adjusted: RUV-inverse")
#legend(legend=c("1","0"),pch=7,cex=0.7,col=1:2)

#Trying different k values
rfit3 <- RUVfit(data=mValsSw, design=design_tw, coef=2, ctl=ctl, method = "ruv4", k=1) # Stage 2 analysis 
rfit4 <- RUVadj(rfit1)

rfit5 <- RUVfit(data=mValsSw, design=design_tw, coef=2, ctl=ctl, method = "ruv4", k=2) # Stage 2 analysis 
rfit6 <- RUVadj(rfit1)

rfit7 <- RUVfit(data=mValsSw, design=design_tw, coef=2, ctl=ctl, method = "ruv4", k=8) # Stage 2 analysis 
rfit8 <- RUVadj(rfit1)

#get adjusted values
Madj1 <- missMethyl::getAdjusted(mValsSw, rfit4)
Madj2 <- missMethyl::getAdjusted(mValsSw, rfit6)
Madj3 <- missMethyl::getAdjusted(mValsSw, rfit8)

#Visualise plots
par(mfrow=c(1,2))
plotMDS(mValsSw, labels=targets_gen$Sample_Name, col=as.integer(factor(targets_gen$Diagnosis..L.0..M.1..H.2.)),dim=c(14,15),
        main="Unadjusted")
#legend("topleft",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)

plotMDS(Madj, labels=targets_gen$Sample_Name, col=as.integer(factor(targets_gen$Diagnosis..L.0..M.1..H.2.)),dim=c(14,15),
        main="Adjusted: RUV-inverse")
#legend("topleft",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)

plotMDS(Madj1, labels=targets_gen$Sample_Name, col=as.integer(factor(targets_gen$Diagnosis..L.0..M.1..H.2.)),dim=c(14,15),
        main="Adjusted: RUV-4, k=1")
#legend("bottomleft",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)

plotMDS(Madj2, labels=targets_gen$Sample_Name, col=as.integer(factor(targets_gen$Diagnosis..L.0..M.1..H.2.)),dim=c(14,15),
        main="Adjusted: RUV-4, k=2")
#legend("bottomright",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)

plotMDS(Madj3, labels=targets_gen$Sample_Name, col=as.integer(factor(targets_gen$Diagnosis..L.0..M.1..H.2.)),dim=c(14,15),
        main="Adjusted: RUV-4, k=8")
#legend
#dev.off()

## Gene ontology of top DMPs ##
#Gene Ontology analysis (gometh): top 1000 probes
res <- read.csv("ASD_guthrie_topruv80%_dmps_onADOS_withintw_Apr8.csv", header = TRUE)
sigCpGs_1k = res$Name[1:1000]
#total number
length(sigCpGs_1k)
sigCpGs_1k = as.character(sigCpGs_1k)
all = res$Name
length(all)
#gometh()
par(mfrow=c(1,1))
gometh_1k <- gometh(sig.cpg=sigCpGs_1k, all.cpg=all, plot.bias=TRUE)
gometh_kegg <- gometh(sig.cpg = sigCpGs_1k, all.cpg = all, collection = "KEGG", prior.prob=TRUE)
topKEGG(gometh_kegg)
write.csv(gometh_1k, file = "GOterms_1k_Apr8.csv", row.names = TRUE)
write.csv(gometh_kegg, file = "GOmeth_kegg_Apr8.csv", row.names = TRUE)
#top GO terms
topGO(gometh_1k, ontology = "BP", number = 20L)
## Total number of significant GO categories at 5% FDR
sum(gometh_1k$FDR<0.5) 

#Plotting effect sizes of top DMPs
plot(ilogit2(mDat[rownames(mDat) == "cg02344606",]), ADOS, main = "IFRD1: cg02344606", ylab = "ADOS z-score", xlab = "Beta Value")
regl5<-lm(ADOS ~ ilogit2(mDat[rownames(mDat) == "cg02344606",]))
abline(regl5)
cor.test(ADOS, ilogit2(mDat[rownames(mDat) == "cg02344606",]), method = "pearson")

plot(ilogit2(mDat[rownames(mDat) == "cg03003115",]), ADOS, main = "ZNF720: cg03003115", ylab = "ADOS z-score", xlab = "Beta Value")
regl5<-lm(ADOS ~ ilogit2(mDat[rownames(mDat) == "cg03003115",]))
abline(regl5)
cor.test(ADOS, ilogit2(mDat[rownames(mDat) == "cg03003115",]), method = "pearson")

plot(ilogit2(mDat[rownames(mDat) == "cg15740041",]), ADOS, main = "KLHDC4: cg15740041", ylab = "ADOS z-score", xlab = "Beta Value")
regl5<-lm(ADOS ~ ilogit2(mDat[rownames(mDat) == "cg15740041",]))
abline(regl5)
cor.test(ADOS, ilogit2(mDat[rownames(mDat) == "cg15740041",]), method = "pearson")

plot(ilogit2(mDat[rownames(mDat) == "cg12902426",]), ADOS, main = "TEAD4: cg12902426", ylab = "ADOS z-score", xlab = "Beta Value")
regl5<-lm(ADOS ~ ilogit2(mDat[rownames(mDat) == "cg12902426",]))
abline(regl5)
cor.test(ADOS, ilogit2(mDat[rownames(mDat) == "cg12902426",]), method = "pearson")

#Effect sizes of differences in ADOS score versus difference in beta value wthin twin pairs

diff_ADOS_1 = abs(ADOS[1] - ADOS[2])
diff_ADOS_2 = abs(ADOS[3] - ADOS[4])
diff_ADOS_3 = abs(ADOS[5] - ADOS[6])
diff_ADOS_4 = abs(ADOS[7] - ADOS[8])
diff_ADOS_5 = abs(ADOS[9] - ADOS[10])
diff_ADOS_6 = abs(ADOS[11] - ADOS[12])
diff_ADOS_7 = abs(ADOS[13] - ADOS[14])
diff_ADOS_8 = abs(ADOS[15] - ADOS[16])
diff_ADOS_9 = abs(ADOS[17] - ADOS[18])
diff_ADOS_10 = abs(ADOS[19] - ADOS[20])
diff_ADOS_11 = abs(ADOS[21] - ADOS[22])
diff_ADOS_score = c(9,4,3,7,11,5,8,2,2,5,12)

diff_effectsize_1 = ilogit2(mDat[rownames(mDat) == "cg02344606"][1]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][2])
diff_effectsize_2 = ilogit2(mDat[rownames(mDat) == "cg02344606"][3]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][4])
diff_effectsize_3 = ilogit2(mDat[rownames(mDat) == "cg02344606"][5]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][6])
diff_effectsize_4 = ilogit2(mDat[rownames(mDat) == "cg02344606"][7]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][8])
diff_effectsize_5 = ilogit2(mDat[rownames(mDat) == "cg02344606"][9]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][10])
diff_effectsize_6 = ilogit2(mDat[rownames(mDat) == "cg02344606"][11]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][12])
diff_effectsize_7 = ilogit2(mDat[rownames(mDat) == "cg02344606"][13]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][14])
diff_effectsize_8 = ilogit2(mDat[rownames(mDat) == "cg02344606"][15]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][16])
diff_effectsize_9 = ilogit2(mDat[rownames(mDat) == "cg02344606"][17]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][18])
diff_effectsize_10 = ilogit2(mDat[rownames(mDat) == "cg02344606"][19]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][20])
diff_effectsize_11 = ilogit2(mDat[rownames(mDat) == "cg02344606"][21]) - ilogit2(mDat[rownames(mDat) == "cg02344606"][22])
diff_effectsize_tot = c(diff_effectsize_1,diff_effectsize_2,diff_effectsize_3,diff_effectsize_4,diff_effectsize_5,diff_effectsize_6,diff_effectsize_7,diff_effectsize_8,diff_effectsize_9,diff_effectsize_10,diff_effectsize_11)
diff_effectsize_tot

#IFRD1: cg02344606
plot(diff_effectsize_tot, diff_ADOS_score, main = "IFRD1: cg02344606", ylab = "ADOS difference z-score", xlab = "Delta Beta Value")
regl5<-lm(diff_ADOS_score ~ diff_effectsize_tot)
abline(regl5)
cor.test(diff_ADOS_score, diff_effectsize_tot, method = "pearson")

#ZNF720: cg03003115
diff_effectsize_1 = ilogit2(mDat[rownames(mDat) == "cg03003115"][1]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][2])
diff_effectsize_2 = ilogit2(mDat[rownames(mDat) == "cg03003115"][3]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][4])
diff_effectsize_3 = ilogit2(mDat[rownames(mDat) == "cg03003115"][5]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][6])
diff_effectsize_4 = ilogit2(mDat[rownames(mDat) == "cg03003115"][7]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][8])
diff_effectsize_5 = ilogit2(mDat[rownames(mDat) == "cg03003115"][9]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][10])
diff_effectsize_6 = ilogit2(mDat[rownames(mDat) == "cg03003115"][11]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][12])
diff_effectsize_7 = ilogit2(mDat[rownames(mDat) == "cg03003115"][13]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][14])
diff_effectsize_8 = ilogit2(mDat[rownames(mDat) == "cg03003115"][15]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][16])
diff_effectsize_9 = ilogit2(mDat[rownames(mDat) == "cg03003115"][17]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][18])
diff_effectsize_10 = ilogit2(mDat[rownames(mDat) == "cg03003115"][19]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][20])
diff_effectsize_11 = ilogit2(mDat[rownames(mDat) == "cg03003115"][21]) - ilogit2(mDat[rownames(mDat) == "cg03003115"][22])
diff_effectsize_tot_cg03003115 = c(diff_effectsize_1,diff_effectsize_2,diff_effectsize_3,diff_effectsize_4,diff_effectsize_5,diff_effectsize_6,diff_effectsize_7,diff_effectsize_8,diff_effectsize_9,diff_effectsize_10,diff_effectsize_11)
diff_effectsize_tot_cg03003115

plot(diff_effectsize_tot_cg03003115, diff_ADOS_score, main = "ZNF720: cg03003115", ylab = "ADOS difference z-score", xlab = "Delta Beta Value")
regl5<-lm(diff_ADOS_score ~ diff_effectsize_tot_cg03003115)
abline(regl5)
cor.test(diff_ADOS_score, diff_effectsize_tot_cg03003115, method = "pearson")

#KLHDC4: cg15740041
diff_effectsize_1 = ilogit2(mDat[rownames(mDat) == "cg15740041"][1]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][2])
diff_effectsize_2 = ilogit2(mDat[rownames(mDat) == "cg15740041"][3]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][4])
diff_effectsize_3 = ilogit2(mDat[rownames(mDat) == "cg15740041"][5]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][6])
diff_effectsize_4 = ilogit2(mDat[rownames(mDat) == "cg15740041"][7]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][8])
diff_effectsize_5 = ilogit2(mDat[rownames(mDat) == "cg15740041"][9]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][10])
diff_effectsize_6 = ilogit2(mDat[rownames(mDat) == "cg15740041"][11]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][12])
diff_effectsize_7 = ilogit2(mDat[rownames(mDat) == "cg15740041"][13]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][14])
diff_effectsize_8 = ilogit2(mDat[rownames(mDat) == "cg15740041"][15]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][16])
diff_effectsize_9 = ilogit2(mDat[rownames(mDat) == "cg15740041"][17]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][18])
diff_effectsize_10 = ilogit2(mDat[rownames(mDat) == "cg15740041"][19]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][20])
diff_effectsize_11 = ilogit2(mDat[rownames(mDat) == "cg15740041"][21]) - ilogit2(mDat[rownames(mDat) == "cg15740041"][22])
diff_effectsize_tot_cg15740041 = c(diff_effectsize_1,diff_effectsize_2,diff_effectsize_3,diff_effectsize_4,diff_effectsize_5,diff_effectsize_6,diff_effectsize_7,diff_effectsize_8,diff_effectsize_9,diff_effectsize_10,diff_effectsize_11)
diff_effectsize_tot_cg15740041

plot(diff_effectsize_tot_cg15740041, diff_ADOS_score, main = "KLHDC4: cg15740041", ylab = "ADOS difference z-score", xlab = "Delta Beta Value")
regl5<-lm(diff_ADOS_score ~ diff_effectsize_tot_cg15740041)
abline(regl5)
cor.test(diff_ADOS_score, diff_effectsize_tot_cg15740041, method = "pearson")

#TEAD4: cg12902426
diff_effectsize_1 = ilogit2(mDat[rownames(mDat) == "cg12902426"][1]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][2])
diff_effectsize_2 = ilogit2(mDat[rownames(mDat) == "cg12902426"][3]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][4])
diff_effectsize_3 = ilogit2(mDat[rownames(mDat) == "cg12902426"][5]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][6])
diff_effectsize_4 = ilogit2(mDat[rownames(mDat) == "cg12902426"][7]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][8])
diff_effectsize_5 = ilogit2(mDat[rownames(mDat) == "cg12902426"][9]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][10])
diff_effectsize_6 = ilogit2(mDat[rownames(mDat) == "cg12902426"][11]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][12])
diff_effectsize_7 = ilogit2(mDat[rownames(mDat) == "cg12902426"][13]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][14])
diff_effectsize_8 = ilogit2(mDat[rownames(mDat) == "cg12902426"][15]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][16])
diff_effectsize_9 = ilogit2(mDat[rownames(mDat) == "cg12902426"][17]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][18])
diff_effectsize_10 = ilogit2(mDat[rownames(mDat) == "cg12902426"][19]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][20])
diff_effectsize_11 = ilogit2(mDat[rownames(mDat) == "cg12902426"][21]) - ilogit2(mDat[rownames(mDat) == "cg12902426"][22])
diff_effectsize_tot_cg12902426 = c(diff_effectsize_1,diff_effectsize_2,diff_effectsize_3,diff_effectsize_4,diff_effectsize_5,diff_effectsize_6,diff_effectsize_7,diff_effectsize_8,diff_effectsize_9,diff_effectsize_10,diff_effectsize_11)
diff_effectsize_tot_cg12902426

plot(diff_effectsize_tot_cg12902426, diff_ADOS_score, main = "TEAD4: cg12902426", ylab = "ADOS difference z-score", xlab = "Delta Beta Value")
regl5<-lm(diff_ADOS_score ~ diff_effectsize_tot_cg12902426)
abline(regl5)
cor.test(diff_ADOS_score, diff_effectsize_tot_cg12902426, method = "pearson")

#DMRCate - differentially methylated region analysis
source("http://bioconductor.org/biocLite.R")
biocLite("DMRcate")
library("DMRcate")
#design matrix in regression 
design_tw <- model.matrix(~ADOS+Twin_Group)
#design_tw_dmrc <- model.matrix(~ADOS)
#myannotation <- cpg.annotate("array", mDat, analysis.type="differential", design=design_tw, coef=2, fdr=1)
myannotation <- cpg.annotate("array", mDat, what = "M", arraytype = "EPIC", analysis.type="differential", design_tw, coef=2)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, pcutoff=0.001)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19") 
head(dmrcoutput$results)
length(dmrcoutput$results)
write.csv(dmrcoutput$results, file = "dmrcoutput_Apr8.csv", row.names = TRUE)
