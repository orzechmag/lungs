library(WGCNA)

options(stringsAsFactors = F)

myData = read.table("lungs.txt", sep = "\t", dec = ".", header = T) # separate files according to the signaling pathway
rownames(myData) = myData$NAME
datExpr0 = as.data.frame(t(myData[, -1]))
names(datExpr0) = myData$NAME
rownames(datExpr0) = names(myData)[-1]

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

traitData = read.table("clinical_annot.txt", sep = "\t", dec = ".", header = T)
lungs = rownames(datExpr0)
traitRows = match(lungs, traitData$id)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]

sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits$type, signed = F)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits[2]), main = "Sample dendrogram and trait heatmap")

enableWGCNAThreads()

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

k = softConnectivity(datE = datExpr0, power = 6) # power = 6 regarding Notch, Hh & Wnt; power = 8 regarding ErbB
sizeGrWindow(10, 5)
par(mfrow = c(1, 2))
hist(k)
scaleFreePlot(k)

softPower = 6 # softPower according to the chosen power (6 or 8)
adjacency = adjacency(datExpr0, power = softPower)

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

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

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
pdf(file = "mergedModuleDendrogram.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEsmoduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

plotTOM = dissTOM^6 # ^ according to the chosen power (6 or 8)
diag(plotTOM) = NA
TOMplot(plotTOM, geneTree, as.character(moduleColors), main = "TOM heatmap plot, module genes", terrainColors = T)

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits$type, use = "p")
moduleTraitPValue = corPvalueStudent(moduleTraitCor, nSamples)

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

sizeGrWindow(8, 9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="blue"
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),nrgcols=300,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )

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

MMblue = data.frame(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]))
write.table(MMblue, "MMblue.txt", row.names = T, col.names = T, sep = "\t", dec = ".")

# repeat lines 161-195 for each module color

GS1 = as.numeric(cor(datTraits$type, datExpr0, use = "p"))
GeneSignificance = abs(GS1)
ModuleSignificance = tapply(GeneSignificance, moduleColors, mean, na.rm = T)

sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,moduleColors)

status = as.data.frame(datTraits$type)
names(status) = "status"
MET = orderMEs(cbind(MEs, status))
sizeGrWindow(5, 7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8,
                      xLabelsAngle = 90)

hubs = chooseTopHubInEachModule(datExpr0, colorh = moduleColors, omitColors = "grey", power = 6)

###
modules = c("blue", "green")
probes = names(datExpr0)
inModule = is.finite(match(moduleColors, modules))
modProbes = (probes[inModule])
modGenes = myData$NAME[match(modProbes, myData$NAME)]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

exportNetworkToCytoscape(modTOM,
                         edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                         nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                         weighted = TRUE,
                         threshold = 0.01,
                         nodeNames = modProbes,
                         altNodeNames = modGenes,
                         nodeAttr = moduleColors[inModule])

# repeat lines  218-233 for each module or modules to be exported to Cytoscape

annot = read.table("annot.txt", sep = "\t", dec = ".")
probes = (names(datExpr0))
probes2annot = match(probes, annot$V1)
allIDs = annot$V2[probes2annot]

library(anRichment)
library(MSigDBCollection)

GOcollection = buildGOcollection(organism = "human") 
GO.BPcollection = subsetCollection(GOcollection, tags = "GO.BP") 
MSigDB = MSigDBCollection(organism = "human")
knownGroups(MSigDB)
MSigDB_coll = subsetCollection(MSigDB, tags = c("MSigDB C2: curated gene sets - CP:KEGG", "MSigDB C5: GO gene sets - BP"))

GOenr = enrichmentAnalysis(moduleColors, allIDs, 
                           refCollection = MSigDB_coll,
                           useBackground = "allOrgGenes",
                           threshold = 0.001,
                           thresholdType = "FDR",
                           getOverlapEntrez = F,
                           getOverlapSymbols = T,
                           maxReportedOverlapGenes = 500,
                           ignoreLabels = c("grey"))
tab = GOenr$enrichmentTable

write.table(tab, "GOenr.txt", sep = "\t", dec = ".", row.names = F, col.names = T)
