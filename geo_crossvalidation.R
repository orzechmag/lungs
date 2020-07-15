setwd("~/Desktop/validation")

library(WGCNA)

options(stringsAsFactors = F)

myData = read.table("valid_exprs.txt", sep = "\t", dec = ".", header = T)
rownames(myData) = myData$ID_REF
datExpr0 = as.data.frame(t(myData[, -1]))
names(datExpr0) = myData$ID_REF
rownames(datExpr0) = names(myData)[-1]

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

traitData = read.table("valid_annot.txt", sep = "\t", dec = ".", header = T)
lungs = rownames(datExpr0)
traitRows = match(lungs, traitData$id)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]

sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits$type, signed = F)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits[2]), main = "Sample dendrogram and trait heatmap")

#enableWGCNAThreads()

# STEP-BY-STEP #
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

k = softConnectivity(datE = datExpr0, power = 20)
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(k)
scaleFreePlot(k)

softPower = 20
adjacency = adjacency(datExpr0, power = softPower)

TOM = TOMsimilarityFromExpr(datExpr0, power = 20)
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.05
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "geneDendro-merged.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEsmoduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

plotTOM = dissTOM^20
diag(plotTOM) = NA
#jpeg(file = "ERBB_TOMplot.jpeg", width = 1569, height = 803)
#TOMplot(plotTOM, geneTree, as.character(moduleColors), main = "TOM heatmap plot, module genes", terrainColors = T)
#dev.off()

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits$type, use = "p")
moduleTraitPValue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTrait = as.data.frame(cbind(moduleTraitCor, moduleTraitPValue))
write.table(moduleTrait, "moduletrait_table.txt", sep = "\t", dec = ".", row.names = T, col.names = T)

sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPValue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits[2]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

type = as.data.frame(datTraits$type)
names(type) = "type"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr0, type, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(type), sep="")
names(GSPvalue) = paste("p.GS.", names(type), sep="")

plotMEpairs(MEs, y = datTraits$type)

sizeGrWindow(8,7)
which.module="turquoise"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
turquoise = datExpr0[,moduleColors==which.module ]
write.table(turquoise, "ERBB_turquoise.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "turquoise"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

sizeGrWindow(8,7)
which.module="brown"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
brown = datExpr0[,moduleColors==which.module ]
write.table(brown, "brown.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "brown"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

sizeGrWindow(8,7)
which.module="blue"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
blue = datExpr0[,moduleColors==which.module ]
write.table(blue, "blue.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

sizeGrWindow(8,7)
which.module="black"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
black = datExpr0[,moduleColors==which.module ]
write.table(black, "black.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "black"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

sizeGrWindow(8,7)
which.module="greenyellow"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
greenyellow = datExpr0[,moduleColors==which.module ]
write.table(greenyellow, "greenyellow.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "greenyellow"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

sizeGrWindow(8,7)
which.module="lightgreen"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
lightgreen = datExpr0[,moduleColors==which.module ]
write.table(lightgreen, "lightgreen.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "lightgreen"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

sizeGrWindow(8,7)
which.module="magenta"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
magenta = datExpr0[,moduleColors==which.module ]
write.table(magenta, "magenta.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "magenta"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

sizeGrWindow(8,7)
which.module="red"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=250,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
red = datExpr0[,moduleColors==which.module ]
write.table(red, "red.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

module = "red"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for LUSC vs LUAD",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

GS1 = as.numeric(cor(datTraits$type, datExpr0, use = "p"))
GeneSignificance = abs(GS1)
ModuleSignificance = tapply(GeneSignificance, moduleColors, mean, na.rm = T)

sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,moduleColors, las = 2)

status = as.data.frame(datTraits$type)
names(status) = "status"
MET = orderMEs(cbind(MEs, status))
sizeGrWindow(5, 7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8,
                      xLabelsAngle = 90)

library(gplots)

lighten <- function(color, factor = 5){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue = 255)
  col
}

col = redgreen(50)
col = sapply(col, lighten)

eset_turquoise = as.matrix(log2(t(turquoise)))
heatmap.2(eset_turquoise, Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", labRow = F, labCol = F, xlab = "patients", ylab = "genes", ColSideColors = c(rep("blue", 130), rep("yellow", 75)))
coord = locator(1)
legend(coord, legend = c("LUSC", "LUAD"), col = c("blue", "yellow"), lty= 1, lwd = 10, xpd = T)

eset_blue = as.matrix(log2(t(blue)))
heatmap.2(eset_blue, Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", labRow = F, labCol = F, xlab = "patients", ylab = "genes", ColSideColors = c(rep("blue", 130), rep("yellow", 75)))
legend(coord, legend = c("LUSC", "LUAD"), col = c("blue", "yellow"), lty= 1, lwd = 10, xpd = T)

eset_brown = as.matrix(log2(t(brown)))
heatmap.2(eset_brown, Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", labRow = F, labCol = F, xlab = "patients", ylab = "genes", ColSideColors = c(rep("blue", 130), rep("yellow", 75)))
legend(coord, legend = c("LUSC", "LUAD"), col = c("blue", "yellow"), lty= 1, lwd = 10, xpd = T)

eset_black = as.matrix(log2(t(black)))
heatmap.2(eset_black, Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", labRow = F, labCol = F, xlab = "patients", ylab = "genes", ColSideColors = c(rep("blue", 130), rep("yellow", 75)))
legend(coord, legend = c("LUSC", "LUAD"), col = c("blue", "yellow"), lty= 1, lwd = 10, xpd = T)

eset_greenyellow = as.matrix(log2(t(greenyellow)))
heatmap.2(eset_greenyellow, Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", labRow = F, labCol = F, xlab = "patients", ylab = "genes", ColSideColors = c(rep("blue", 130), rep("yellow", 75)))
legend(coord, legend = c("LUSC", "LUAD"), col = c("blue", "yellow"), lty= 1, lwd = 10, xpd = T)

eset_lightgreen = as.matrix(log2(t(lightgreen)))
heatmap.2(eset_lightgreen, Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", labRow = F, labCol = F, xlab = "patients", ylab = "genes", ColSideColors = c(rep("blue", 130), rep("yellow", 75)))
legend(coord, legend = c("LUSC", "LUAD"), col = c("blue", "yellow"), lty= 1, lwd = 10, xpd = T)

eset_magenta = as.matrix(log2(t(magenta)))
heatmap.2(eset_magenta, Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", labRow = F, labCol = F, xlab = "patients", ylab = "genes", ColSideColors = c(rep("blue", 130), rep("yellow", 75)))
legend(coord, legend = c("LUSC", "LUAD"), col = c("blue", "yellow"), lty= 1, lwd = 10, xpd = T)

eset_red = as.matrix(log2(t(red)))
heatmap.2(eset_red, Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", dendrogram = "row", labRow = F, labCol = F, xlab = "patients", ylab = "genes", ColSideColors = c(rep("blue", 130), rep("yellow", 75)))
legend(coord, legend = c("LUSC", "LUAD"), col = c("blue", "yellow"), lty= 1, lwd = 10, xpd = T)

install.packages("~/Downloads/anRichmentMethods_0.91-94.tar.gz", repos = NULL, type = "source"); 
install.packages("~/Downloads/anRichment_1.10-1.tar.gz", repos = NULL, type = "source")
library(anRichmentMethods)
library(anRichment)

GOcollection = MSigDBCollection("msigdb_v7.1.xml", 
                      MSDBVersion = "7.1", organism = "human", excludeCategories = c("C1", "C3", "C4", "C6", "C7"))
knownGroups(GOcollection)
MSigDB_coll = subsetCollection(GOcollection, tags = c("MSigDB C5: GO gene sets - BP", "MSigDB C2: curated gene sets - CP:KEGG"))

probes = names(datExpr0)
annot_entrez = read.table("annot_entrez.txt", sep = "\t", dec = ".", header = T)
probes2annot = match(probes, annot_entrez$Affy.ID)# Get the corresponding Locuis Link IDs
allLLIDs = annot_entrez$Gene[probes2annot]
list = c("turquoise", "blue", "brown", "black", "greenyellow", "lightgreen", "magenta", "red")
GOenr = enrichmentAnalysis(moduleColors, allLLIDs, 
                           refCollection = MSigDB_coll,
                           useBackground = "allOrgGenes",
                           threshold = 0.001,
                           thresholdType = "FDR",
                           getOverlapEntrez = T,
                           getOverlapSymbols = F,
                           maxReportedOverlapGenes = 500,
                           ignoreLabels = c("grey"))
tab = GOenr$enrichmentTable
write.table(tab, "valid_GOenr.txt", sep = "\t", dec = ".", row.names = F, col.names = T)

library(FactoMineR)
library(factoextra)

dat_notch = read.table("notch_pca_valid.txt", sep = "\t", dec = ".", header = T)
res_notch = PCA(dat_notch[,3:74])
hab = as.factor(dat_notch$type)
fviz_pca_ind(res_notch, geom = "point", habillage = hab)

dat_wnt = read.table("wnt_pca_valid.txt", sep = "\t", dec = ".", header = T)
res_wnt = PCA(dat_wnt[,3:952])
fviz_pca_ind(res_wnt, geom = "point", habillage = hab)

dat_erbb = read.table("erbb_pca_valid.txt", sep = "\t", dec = ".", header = T)
res_erbb = PCA(dat_erbb[,3:490])
fviz_pca_ind(res_erbb, geom = "point", habillage = hab)

dat_hh = read.table("hh_pca_valid.txt", sep = "\t", dec = ".", header = T)
res_hh = PCA(dat_hh[,3:145])
fviz_pca_ind(res_hh, geom = "point", habillage = hab)
