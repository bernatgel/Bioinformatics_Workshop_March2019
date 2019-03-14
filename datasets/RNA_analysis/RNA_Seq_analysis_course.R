####################################################
##
####################################################

library('NOISeq')
library('cqn')
library('biomaRt')
library('sva')   #need to find mgcv nlme genefilter BiocParallel
library('limma')
library('maSigPro')
library('mdgsa')

rna.table <- read.table('GSE112221_RNA_seq.txt', header = 1)
rna.table <- rna.table[!duplicated(rna.table$geneID),]
rownames(rna.table) <- rna.table$geneID
rna.table <- rna.table[3:12]
rna.table <- rna.table[which(rowSums(rna.table) != 0),]
rna.table <- rna.table[1:8]
colnames(rna.table) <- c("NL1", "HCC1", "NL2", "HCC2", "NL3", "HCC3", "NL4", "HCC4")

mydesign <- data.frame("cellType" = c("NL", "HCC", "NL", "HCC", "NL", "HCC", "NL", "HCC"), 
                       "patient" = c(1, 1, 2, 2, 3, 3, 4, 4))
rownames(mydesign) <- c("NL1", "HCC1", "NL2", "HCC2", "NL3", "HCC3", "NL4", "HCC4")


# Gene annotation from Biomart --------------------------------------------

biomartHuman = useMart("ensembl", dataset="hsapiens_gene_ensembl")  
atributes = listAttributes(biomartHuman)
atributes[grep("gc", atributes$description, ignore.case = TRUE),]


myannot = getBM(attributes = c("external_gene_name", "percentage_gene_gc_content", "gene_biotype", "transcript_length"),
                         filters = "external_gene_name", values=rownames(rna.table), mart=biomartHuman)
myannot <- myannot[!duplicated(myannot$external_gene_name),]
rna.table <- rna.table[which(rownames(rna.table)%in%myannot$external_gene_name),]

table_design <- list('rna.table' = rna.table, 'mydesign' = mydesign, 'myannot' = myannot)

#save('rna.table' = rna.table, 'mydesign' = mydesign, 'myannot' = myannot, file = "table_design.RData")
#load("table_design.RData")

# NOISeq Quality Control --------------------------------------------------

noiseqData = readData(data = rna.table, factors = mydesign, 
                      gc = myannot[,1:2], biotype = myannot[,c(1,3)]) 

### BIODETECTION PLOT:
mybiodetection = NOISeq::dat(noiseqData, type = "biodetection", logtransf = T)

pdf("noiseqPlot_biodetection.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
for (i in 1:4) {
  NOISeq::explo.plot(mybiodetection, plottype = "persample", samples = c(i,i+2))
}
dev.off()


### BIOTYPE DISTRIBUTION PLOT: 
## D36 - Task 3B  ---->>>  TO DO

mycountsbio = NOISeq::dat(noiseqData, type = "countsbio")

pdf("noiseqPlot_biotype_distribution_1sample.pdf", width = 7, height = 7)
for (i in 1:8) {
  NOISeq::explo.plot(mycountsbio, plottype = "boxplot", toplot = 1, samples=i)
}
dev.off()


### SATURATION PLOT: 

mysaturation = NOISeq::dat(noiseqData, type = "saturation")

pdf("noiseqPlot_biotype_saturation.pdf", width = 7, height = 7)
NOISeq::explo.plot(mysaturation, toplot = 1, samples = c(1:8), yleftlim = NULL, yrighlim = NULL)
dev.off()


### COUNT DISTRIBUTION PLOT
## D36 - Task 3D

mycountsbio = NOISeq::dat(noiseqData, type = "countsbio", factor = NULL)

pdf("noiseqPlot_count_distribution_global.pdf", width = 10, height = 5)
NOISeq::explo.plot(mycountsbio, plottype = "boxplot", samples = NULL)
dev.off()


# LOW COUNTS
myCounts = NOISeq::dat(noiseqData, type = "countsbio")
pdf("noiseqPlot_distribution_lowcounts.pdf", width = 7, height = 7)
NOISeq::explo.plot(myCounts, plottype = "barplot", samples = NULL)
dev.off()


# Filtering low counts ----------------------------------------------------
## D44 - Task 1

rnaseq = NOISeq::filtered.data(rna.table, factor = mydesign$cellType, 
                               norm = FALSE, method = 1, cpm = 2, cv.cutoff = 500)

noiseqData = readData(data = rnaseq, factors = mydesign)

myCounts = NOISeq::dat(noiseqData, type = "countsbio")
pdf("noiseqPlot_distribution_lowcounts.pdf", width = 7, height = 7)
NOISeq::explo.plot(myCounts, plottype = "barplot", samples = NULL)
dev.off()


# Length bias

noiseqData = NOISeq::readData(data = rnaseq, 
                              factors = mydesign,
                              gc = myannot[,c(1,2)],
                              length = myannot[, c(1,4)])

mylengthbias = NOISeq::dat(noiseqData, type = "lengthbias")
pdf("noiseqPlot_lengthbias.pdf", width = 14, height = 7)
x11()
par(mfrow = c(2, 2))
for(i in 1:8) {
  NOISeq::explo.plot(mylengthbias, samples = i, toplot = "global")
}
dev.off()


myGCcontent <- NOISeq::dat(noiseqData, k = 0, type = "GCbias")

pdf("noiseqPlot_GCcontent.pdf", width = 6*3, height = 4*3)
x11()
par(mfrow = c(2,2))
for(i in 1:8) {
  NOISeq::explo.plot(myGCcontent, samples = i)
}
dev.off()


# PCA

myPCA = NOISeq::dat(noiseqData, type = "PCA", norm = FALSE, logtransf = T)
pdf("noiseqPlot_PCA_before_normalization.pdf", width = 5*2, height = 5)
x11()
par(mfrow = c(1,2))
NOISeq::explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "cellType")
NOISeq::explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "patient")


#TMM normalization

rnaseq <- rnaseq[!is.na(rownames(rnaseq[myannot$external_gene_name[myannot$gene_biotype=="protein_coding"],])),]

myTMM = tmm(rnaseq, long = 1000, lc = 0, k = 0, refColumn = 1,
            logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)

save("myTMM" = myTMM, "mydesign" = mydesign, "myannot" = myannot, file = 'myTMM.RData')
#load(myTMM.RData)

noiseqData = NOISeq::readData(data = myTMM,
                              factors = mydesign,
                              gc = myannot[,c(1,2)],
                              length = myannot[,c(1,4)])

myPCA = dat(noiseqData, type = "PCA", norm = TRUE, logtransf = TRUE)
X11()
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "cellType")
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "patient")

dev.off()


### Differential expression
cellType <- mydesign$cellType
DE.design = model.matrix(~0 + cellType)

fit = lmFit(myTMM, DE.design)
fit = eBayes(fit)
topTable(fit)

#################

myContrMatrix = makeContrasts(
  HCC.NL = (cellTypeHCC - cellTypeNL),
  levels = DE.design)
myContrMatrix
fitTIME = contrasts.fit(fit, myContrMatrix)
fitTIME = eBayes(fitTIME)
head(fitTIME$cov.coefficients)
topTable(fitTIME, coef = "HCC.NL")

DEGlimma = decideTests(fitTIME, method="separate", adjust.method="BH", p.value=0.01, lfc=0)
head(DEGlimma)

DE.genes <- rownames(fitTIME$p.value)[which(fitTIME$p.value<0.01)]

save("DE.genes" = DE.genes, file = "DE_genes.RData")
#load(DE_genes.RData)

# Functional Enrichment Analysis ------------------------------------------

go.annot = getBM(attributes = c("external_gene_name","go_id", "name_1006"),
                                   filters = "external_gene_name", values=DE.genes, mart=biomartHuman)

hsapiens.genes <- getBM(attributes = c("external_gene_name", "gene_biotype"), mart=biomartHuman)
hsapiens.genes <- hsapiens.genes[which(hsapiens.genes$gene_biotype=="protein_coding"),1]

EnrichALLterms = function (test, notTest, annotation, p.adjust.method = "fdr") {
  
  annot2test = unique(annotation[,3])
  
  resultat = t(sapply(annot2test, Enrich1term, test = test, notTest = notTest, annotation = annotation))
  
  return (data.frame(resultat, 
                     "adjPval" = p.adjust(as.numeric(resultat[,"pval"]), method = p.adjust.method), 
                     stringsAsFactors = F))
  
}

Enrich1term = function (term, test, notTest, annotation) {
  
  annotTest = length(intersect(test, annotation[annotation[,3] == term,1]))
  
  if ((annotTest) > 0) {
    annotNOTtest = length(intersect(notTest, annotation[annotation[,3] == term,1]))
    mytest = matrix(c(annotTest, length(test)-annotTest, annotNOTtest, length(notTest)-annotNOTtest), ncol = 2)
    resultat = c(term, annotTest, length(test), annotNOTtest, length(notTest), 
                 fisher.test(mytest, alternative = "greater")$p.value)
    names(resultat) = c("term", "annotTest", "test", "annotNotTest", "notTest", "pval")
  } else {
    resultat = c(term, 0, 0, 0, 0, 100)
    names(resultat) = c("term", "annotTest", "test", "annotNotTest", "notTest", "pval")
  }
  
  return(resultat)
  
}

myEnrichResults <- EnrichALLterms(DE.genes, hsapiens.genes, go.annot)
myEnrichResults$pval <- as.numeric(myEnrichResults$pval)