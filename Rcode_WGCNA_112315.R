############## Install WGCNA ####################
source("http://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
biocLite("impute") 
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg"); 
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6)); 
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep=""); 
biocLite(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages("WGCNA")
install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 

## Analysis by Guanjing Hu from November 23th, 2015
## Modified from a previous version (July 3th, 2015) to remove a D5 20dpa sample that turns out to be from polyploid cotton
## This the final (hopefully) version of the integrated network analysis for the cotton seed transcriptomic data.
## R code was adpated from from chapter 12 of the book " Horvath S (2011) Weighted Network Analysis. Applications in Genomics and Systems Biology. Springer Book. ISBN: 978-1-4419-8818-8"
## The analysis was coducted in following steps
## 1. Basic data processing and cleaning - input DESeq2 rlog table and trait table for sample clustering and detecting outliers.
## 2. Choosing the soft-thresholding power - default or powers making good fit of scale free topology.
## 3. Network Construction - single block, corType = "pearson" (not "bicor", see discussion below), networkType = "signed"
## 4. General network topology analysis - produce a few plots for exploration





############### Step 1. Use raw count and DESeq2  ###############
## unix command: R CMD BATCH s1.dataProcessing.R
################

date()
getwd()
library(WGCNA);
library(RColorBrewer)
library(flashClust);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# enableWGCNAThreads();

# load expression data
data = read.table("All ADT not normalized.txt",header=TRUE, sep="\t");
# Take a quick look at what is in the data set:
dim(data);  #37223   157
names(data);
names(data)<-gsub(".bam","",names(data))

# here is the trouble: D5.20_3
colSums(data[,-1])
# D5.10_2.T  D5.10_3.D  D5.10_3.T  D5.20_1.D  D5.20_1.T  D5.20_2.D  D5.20_2.T 
#  11798417   10504123   14636706    7456999   10126531    7399781   10004142 
# D5.20_3.D  D5.20_3.T  D5.30_1.D  D5.30_1.T  D5.30_2.D  D5.30_2.T  D5.30_3.D 
#   4986424   14520980    9865550   13304659   12960104   17623242    8482872 


#  Make each row corresponds to a gene and column to a sample or auxiliary information.
datExprT = as.data.frame(t(data[, -1]));
names(datExprT) = data[,1];
dim(datExprT<-datExprT[-grep(".A$|.D$",rownames(datExprT),perl=TRUE ), ] )
rownames(datExprT)<-gsub(".T","",rownames(datExprT))

write.table(t(datExprT),file="s1.count.raw.txt", sep="\t")
pdf("s1.raw.boxplot.pdf")
boxplot(log2(t(datExprT) ), las=2)
dev.off()


# normalization using DESeq2 rlog
count<-as.data.frame( t(datExprT) )
coldata<-data.frame( sample = gsub("_.","", names(count)), genome = gsub("[.].*","", names(count)), dpa = gsub(".*[.]|_.","", names(count)), rep = gsub(".*_","", names(count)) )
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
# rlog transformation, note that with defalt blind=TRUE, design is actually ~1
rld <- rlog(dds)
head(assay(rld))
pdf("s1.rld_pca.pdf")
#PCA plot
#Another way to visualize sample-to-sample distances is a principal-components analysis (PCA). In this ordination method, the data points (i.e., here, the samples) are projected onto the 2D plane such that they spread out in the two directions which explain most of the differences in the data. The x-axis is the direction (or principal component) which separates the data points the most. The amount of the total variance which is contained in the direction is printed in the axis label.
plotPCA(rld, intgroup = c("genome", "dpa"))
# Here, we have used the function plotPCA which comes with DESeq2. The two terms specified by intgroup are the interesting groups for labeling the samples; they tell the function to use them to choose colors.

# We can also build the PCA plot from scratch using ggplot2. This is done by asking the plotPCA function to return the data used for plotting rather than building the plot. See the ggplot2 documentation for more details on using ggplot.
library(genefilter)
sumPCA<-
function (x, intgroup = "condition", ntop = 500)
{
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]),
    1, paste, collapse = " : "))
    return(pca)
}
data<-sumPCA(rld, intgroup = c("genome", "dpa"))
summary(data)
### Importance of components:
### PC1     PC2     PC3      PC4      PC5     PC6     PC7     PC8    PC9    PC10    PC11    PC12    PC13    PC14    PC15    PC16
### Standard deviation     50.4710 21.0492 15.2206 13.39625 11.68574 9.38357 7.38384 6.89884 6.7829 6.11396 5.69427 4.70881 4.59109 4.24778 4.01172 3.87475
### Proportion of Variance  0.6256  0.1088  0.0569  0.04408  0.03354 0.02163 0.01339 0.01169 0.0113 0.00918 0.00796 0.00545 0.00518 0.00443 0.00395 0.00369
### Cumulative Proportion   0.6256  0.7345  0.7913  0.83543  0.86897 0.89059 0.90398 0.91567 0.9270 0.93615 0.94412 0.94956 0.95474 0.95917 0.96312 0.96681
library(ggplot2)
qplot(PC1, PC2, main="top 500", color=coldata$genome, shape=coldata$dpa, data=as.data.frame(data$x)) +
xlab(paste0("PC1: 62.6% variance")) +
ylab(paste0("PC2: 10.9% variance"))
# From both PCA visualizations, we see that the differences between genomes are not stronger than the differences due to dpa.

data<-sumPCA(rld, intgroup = c("genome", "dpa"), ntop=37223)
summary(data)
### Importance of components:
### PC1     PC2     PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10     PC11     PC12
### Standard deviation     94.5976 63.7282 54.8723 52.64956 40.00305 32.31420 28.16078 26.19804 24.86588 22.63824 19.82023 18.47903
### Proportion of Variance  0.3079  0.1398  0.1036  0.09539  0.05507  0.03593  0.02729  0.02362  0.02128  0.01764  0.01352  0.01175
### Cumulative Proportion   0.3079  0.4477  0.5513  0.64669  0.70175  0.73768  0.76497  0.78859  0.80987  0.82750  0.84102  0.85277
library(ggplot2)
qplot(PC1, PC2, main="all 37223", color=coldata$genome, shape=coldata$dpa, data=as.data.frame(data$x)) +
xlab(paste0("PC1: 30.8% variance")) +
ylab(paste0("PC2: 14.0% variance"))
# From both PCA visualizations, we see that the differences between genomes are not stronger than the differences due to dpa.

dev.off()

expr<-as.data.frame(assay(rld) )
names(expr)<-names(count)
expr<-expr[,-which(names(expr) %in% "D5.20_3")]
write.table(expr,"s1.count.rlog.txt", sep="\t")

# make subsets for each genome
datExprT<-t(expr)
Adata<-datExprT[grep("A2",rownames(datExprT) ),]
Ddata<-datExprT[grep("D5",rownames(datExprT) ),]
TM1data<-datExprT[grep("TM1",rownames(datExprT) ),]
Yucdata<-datExprT[grep("Yuc",rownames(datExprT) ),]
AD3data<-datExprT[grep("AD3",rownames(datExprT) ),]
Fdata<-(Adata[-which(rownames(Adata) %in% "A2.20_3"),]+Ddata)/2
rownames(Fdata) <- gsub("A2","ADs",rownames(Fdata))
# construct expr
dataT<- rbind(Adata, Ddata, AD3data, Yucdata, TM1data)
dataF<- rbind(dataT, Fdata)   # including ADs
dataP<- rbind(Adata,TM1data, Ddata, Fdata)  # as proteomic set


# Now we read in the physiological trait data
oil<-read.table("oil_content073115.txt", header=TRUE, sep="\t")
oil<-oil[!is.na(oil$percentage_oil_content),]
oil$sample<-paste(oil$genome, oil$dpa,sep=".")
oil_content<-aggregate(oil$percentage_oil_content,list(oil$sample),mean)

weight<-read.table("seed_weight.txt", header=TRUE, sep="\t")
weight<-weight[!is.na(weight$weight),]
weight$sample<-paste(weight$genome, weight$DPA,sep=".")
seed_weight<-aggregate(weight$weight,list(weight$sample),mean)
traitData = merge(seed_weight, oil_content, by="Group.1",all.y=TRUE)
names(traitData)<-c("sample", "seed_weight","oil_content")
traitData$dpa<-as.numeric(gsub(".*[.]","",traitData$sample) )
traitData$sample<-gsub("YUC","Yuc",gsub("AD1","TM1",traitData$sample))
traitData
#    sample seed_weight oil_content dpa
# 1   A2.10 0.003076923    4.353808  10
# 2   A2.20 0.077333333    4.834227  20
# 3   A2.30 0.108055556    8.320123  30
# 4   A2.40 0.107333333   13.716703  40
# 5  TM1.10 0.015320911    5.050659  10
# 6  TM1.20 0.166055556    8.681095  20
# 7  TM1.30 0.233333333   20.902414  30
# 8  TM1.40 0.211000000   22.499568  40
# 9  AD3.10          NA    4.491244  10
# 10 AD3.20          NA    4.739112  20
# 11 AD3.30          NA   13.470681  30
# 12 AD3.40          NA   15.596463  40
# 13  D5.10 0.005007407    4.206810  10
# 14  D5.20 0.040119048    4.495644  20
# 15  D5.30 0.047500000   17.461650  30
# 16  D5.40 0.049333333   21.295438  40
# 17 Yuc.10          NA    5.653722  10
# 18 Yuc.20          NA    3.129453  20
# 19 Yuc.30          NA   10.135646  30
# 20 Yuc.40          NA   19.571135  40


# give the same trait value to each replicated sample, match order
sampleR<-data.frame(sampleR=rownames(dataT) )
sampleR$sample<-gsub("_.","",sampleR$sampleR)
traitData<-merge(sampleR, traitData, by="sample")
traitData<-traitData[,-1]
# Order the rows of Traits so that they match those of dataT
traitRows = match(rownames(dataT), traitData$sampleR)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
# show that row names agree
table(rownames(datTraits)==rownames(dataT))
# TRUE 59
# SO the traits and expression data have been aligned correctly.


# sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(t(dataT),type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)

# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor1=ifelse(Z.k<thresholdZ.k,"red","black")

thresholdZ.k= -2.5
outlierColor2=ifelse(Z.k<thresholdZ.k ,"red","black")


pdf("s1.sample_dendrogram_and_trait_heatmap.rlog.pdf")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation:
# where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits),"C",sep="")
datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2, traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
colors=datColors,main="Sample dendrogram and trait heatmap")
dev.off()
## "D5.40_2" and "D5_40.3"are outliers, but it seems real to me biologically
## compared the dendrogram to previous analysis, D5.20_3 is more likely to be TM1.20


# Next we make a multi-set data, considering 3 different combinations
## dataT<- rbind(Adata, Ddata, AD3data, Yucdata, TM1data)
## dataF<- rbind(dataT, Fdata)   # including ADs
## dataP<- rbind(Adata,TM1data, Ddata, Fdata)  # as proteomic set
nSets = 3
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Five genomes", "Additional synthetic" ,"Proteomic set")
shortLabels = c("setT", "setF" ,"setP")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = dataT);
multiExpr[[2]] = list(data = dataF);
multiExpr[[3]] = list(data = dataP);

# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)
# $nSets 3
# $nGenes 37223
# $nSamples 59 70 46
# $structureOK  TRUE

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding ZERO genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
# Update exprSize
checkSets(multiExpr)
# $nSets 3
# $nGenes 36560
# $nSamples 59 70 46
# $structureOK  TRUE


pdf(file = "s1.SampleClusteringS.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;

save(multiExpr,nSets, setLabels, shortLabels, datTraits, file = "R-01-dataInput.RData")


############### Step 2.  Choosing the soft-thresholding power: analysis of network topology  ###############
## R CMD BATCH s2.choosePower.R &
################

date()
library(WGCNA)
library(RColorBrewer)
library(ggplot2);
options(stringsAsFactors = FALSE)
# enableWGCNAThreads(nThreads=10)

lnames = load(file = "R-01-dataInput.RData")
lnames
nSets
nGenes<-checkSets(multiExpr)$nGenes #36560
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels
shortLabels


# Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
types<-c("unsigned", "signed", "signed hybrid")

for (type in types)
{
    powers = c(c(1:10), seq(from = 12, to=40, by=2))
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type)[[2]])      }
    collectGarbage()
    
    # Plot the results:
    colors=brewer.pal(5,"Set1")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
        for (col in 1:length(plotCols))
        {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        }
    }
    
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    sizeGrWindow(8, 6)
    pdf(paste("s2.ChooseSoftThresholdPower_",gsub(".* ","", type), ".pdf", sep="") )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
        if (set==1)
        {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
            addGrid()
        }
    if (col==1)
    {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
        labels=powers,cex=cex1,col=colors[set]);
    } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
    if (col==1)
    {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
    } else
    legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
    
    dev.off()
    assign(paste("powerTables.",gsub(".* ","", type),sep=""),powerTables)
# repeat above for powerTables.signed, and powerTables.hybrid
}

save(powerTables.unsigned,powerTables.signed,powerTables.hybrid , file = "R-02-choosePower.RData")


# Inspect "s2.ChooseSoftThresholdPower_signed.pdf". Default Power=12 is reasonable enough.
# Basically, the default power is 6 for unsigned network, 12 for signed network; if the observed power choice is smaller than default, I will choose the observed power, otherwise use default.





############### Step 3.  Network Construction  ###############
## nohup R CMD BATCH s3.buildNetwork.R &
################

date()
library(WGCNA)
options(stringsAsFactors = FALSE);
# enableWGCNAThreads(nThreads=10)

library(RColorBrewer)
library(ggplot2);

lnames = load(file = "R-01-dataInput.RData")
lnames     #  "multiExpr"   "nSets"       "setLabels"   "shortLabels" "datTraits"
nSets      # 3
setLabels  #  "Five genomes"         "Additional synthetic" "Proteomic set"
shortL#  "setT" "setF" "setP"
nGenes<-checkSets(multiExpr)$nGenes
nGenes     #  36560

powers<-c(12)

# work with individual genome, then work with different soft threshold
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    for (j in powers )
    {
        softPower = j
        print(paste("Start building network for ",shortLabels[set]," using soft threshold ",j,"......",sep=""))
        # Network construction
        
        net = blockwiseModules(
             # Input data
             subDat,
             # Data checking options
             checkMissingData = TRUE,
             
             # Options for splitting data into blocks
             blocks = NULL,
             randomSeed = 12345,
             maxBlockSize = nGenes,  # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
             
             # Network construction arguments: correlation options, use bicor instead of default pearson
             corType = "pearson",
             # Adjacency and topology overlap function options
             power = j, networkType = "signed", TOMType = "signed",
             
             # Saving or returning TOM
             saveTOMs = TRUE,
             saveTOMFileBase = paste(shortLabels[set],"_power",j,"_TOM",sep=""),
             
             # Basic tree cut options
             deepSplit = 2,  #default, known to reasonable
             minModuleSize = min(30, ncol(subDat)/2 ), #default 20, use 30 for transcriptome
             pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
             
             # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
             mergeCutHeight = 0.25,
             
             # others
             reassignThreshold = 0,
             numericLabels = TRUE,
             verbose = 3)
             
        assign(paste(shortLabels[set],"net",j,sep=""), net)
        }
}
save(list=grep(".+net.+",ls(), value=TRUE), file = "R-03-buildNetwork.RData")

# OK this is pretty tricky here. The use of "bicor" correlation is supposed to be more robust than "pearson" correlation, because it is good at removing outliers from correlation calculation. BUT we have a different problem here: "pearson" network results in over 30 modules, while "bicor" gives only 2 modules; "pearson" resulted modules are much better correlated with gene significance than "bicor" modules. This is because some genome&dpa specific expressions (represented by only 3 samples) were treated as outliers, and their effects were removed. Below code and pdf output illustrated this point. SO we should stick to "Pearson" correlation.
##########################
pdf("s3.CompareBicor&Pearson.pdf")
netB = blockwiseModules( multiExpr[[1]]$data, power = 12, maxBlockSize = nGenes, corType = "bicor", networkType = "signed", TOMType = "signed")
netP = blockwiseModules( multiExpr[[1]]$data, power = 12, maxBlockSize = nGenes, corType = "pearson", networkType = "signed", TOMType = "signed")
# Use seed weight to define a gene significance variable
GS.weight=as.numeric(cor(multiExpr[[1]]$data,datTraits$seed_weight,use="p"))
# This translates the numeric values into colors
GS.weightColor=numbers2colors(GS.weight,signed=T)
# Use oil to define a gene significance variable
GS.oil=as.numeric(cor(multiExpr[[1]]$data,datTraits$oil_content,use="p"))
# This translates the numeric values into colors
GS.oilColor=numbers2colors(GS.oil,signed=T)
# Use dpa to define a gene significance variable
GS.dpa=as.numeric(cor(multiExpr[[1]]$data,datTraits$dpa,use="p"))
# This translates the numeric values into colors
GS.dpaColor=numbers2colors(GS.dpa,signed=T)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netB$dendrograms[[1]],colors=data.frame(labels2colors(netB$colors),labels2colors(netP$colors), GS.weightColor, GS.oilColor, GS.dpaColor), groupLabels=c("Bicor","Pearson", "Seed_weight", "oil_content", "dpa"),dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05, main="Cluster dendrogram, bicor")
plotDendroAndColors(netP$dendrograms[[1]],colors=data.frame(labels2colors(netP$colors),labels2colors(netB$colors), GS.weightColor, GS.oilColor, GS.dpaColor), groupLabels=c("Pearson","Bicor", "Seed_weight", "oil_content", "dpa"),dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05, main="Cluster dendrogram, pearson")
dev.off()
##########################



############### Step 4.  General network topology analysis  ###############
## nohup R CMD BATCH s4.networkTopology.R &
################


library(WGCNA);
library(RColorBrewer);
library(scatterplot3d);
library(flashClust);
library(ggplot2);
options(stringsAsFactors = FALSE);

remove(list=ls())
source('multiplot.r', chdir = TRUE)
source('multiscalefreeplot.r', chdir = TRUE)
source('summarySE.r', chdir = TRUE)
source('WGCNA_missingFun.r', chdir = TRUE)

lnames = load(file = "R-01-dataInput.RData")
lnames     #  "multiExpr"   "nSets"       "setLabels"   "shortLabels", "datTraits"
nSets      # 3
setLabels  #  "Five genomes"         "Additional synthetic" "Proteomic set"
shortLabels# "setT" "setF" "setP"
nGenes<-checkSets(multiExpr)$nGenes #36560

load("R-03-buildNetwork.RData")



powers<-c(12)
pname<-paste("power=", powers, sep="")


# work with individual genome, then work with different soft threshold
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    genome <-  shortLabels[set]
    net1<-get(paste(genome,"net",powers[1],sep="") )
#   net2<-get(paste(genome,"net",powers[2],sep="") )
#   net3<-get(paste(genome,"net",powers[3],sep="") )
#   net4<-get(paste(genome,"net",powers[4],sep="") )
    gg<-net1$goodGenes    #get the good genes descision
    
    pdf(paste("s4.",genome, "_connectivity.pdf", sep=""))

    # Network topology and other concepts need to be explored more!!!!!!!!!!!!!!!!!!
    # The following computes the network connectivity (Connectivity)
    # pdf("A2_network_connectivity.pdf")
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=length(powers)))
    names(degree)=pname
    for(j in 1:length(powers) ){
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        # Creating Scale Free Topology Plots (SUPPLEMENTARY FIGURE S1 in our article)
        scaleFreePlot(k, truncated=T,main= paste("beta=",powers[j]))
        degree[,j]<-k
    }
# multiScaleFreePlot(degree, paste("Connectivity distribution - ", genome, sep=""))
    
    # Plot dendrograms all three clustering results
    plotDendroAndColors( net1$dendrograms[[1]], main = paste( pname[1], "dendrogram" ), labels2colors(net1$colors[gg]),pname,dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
#   plotDendroAndColors( net2$dendrograms[[1]], main = paste( pname[2], "dendrogram" ), cbind(labels2colors(net1$colors[gg]),labels2colors(net2$colors[gg])),pname,dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
#   plotDendroAndColors( net3$dendrograms[[1]], main = paste( pname[3], "dendrogram" ),  cbind(labels2colors(net1$colors[gg]),labels2colors(net2$colors[gg]),labels2colors(net3$colors[gg]) ),pname,dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
#   plotDendroAndColors( net4$dendrograms[[1]], main = paste( pname[4], "dendrogram" ), cbind(labels2colors(net1$colors[gg]),labels2colors(net2$colors[gg]),labels2colors(net3$colors[gg]), labels2colors(net4$colors[gg])),pname,dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

    dev.off()
}


# plot all powers for all genomes
pdf("s4.Allsets_Connectivity.pdf")
for(j in 1:length(powers))
{
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=4))
    names(degree)<- shortLabels
    for(set in 1:nSets)
    {
        subDat <-  multiExpr[[set]]$data
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        degree[,set]<-k
    }
    multiScaleFreePlot(degree, paste("Connectivity distribution, Power=",powers[j],sep=""))
}
dev.off()


# According to above results, use power=12
# Next we make a multi-set data, considering 3 different combinations
## dataT<- rbind(Adata, Ddata, AD3data, Yucdata, TM1data)
## dataF<- rbind(dataT, Fdata)   # including ADs
## dataP<- rbind(Adata,TM1data, Ddata, Fdata)  # as proteomic set

samples = rep( c("A10", "A20", "A30" , "A40", "D10", "D20", "D30" , "D40",  "Tom10", "Tom20", "Tom30",   "Tom40", "Yuc10", "Yuc20", "Yuc30", "Yuc40", "TM10", "TM20", "TM30", "TM40", "syn10", "syn20", "syn30", "syn40" ), each=3)
samples<-samples[c(1:17,19:63,65:72)]

sample= list(samples[1:59], samples, samples[c(1:12,48:59, 13:23, 60:70)])

for (set in 1:3 ) # all and trio
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    genome <-  shortLabels[set]
    net<-get(paste(genome,"net12",sep="") )
    gg<-net$goodGenes    #get the good genes descision
#    adjacency = adjacency(subDat, power = 12, type = "signed")
    Nmodules= dim(net$MEs)[2]
    assign(paste("net", genome, sep=""),net)
    
#   load(net$TOMFiles)
    print(net$TOMFiles)
    print(paste("Number of modules in ",genome," network is ",Nmodules,sep=""))
    colorsa<- labels2colors(net$colors)
    
    pdf(paste("s4.",genome,"_modules.pdf",sep=""))
    
    # Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
    # Displaying module heatmap and the eigengene
    # sizeGrWindow(8,7);
    MEs<-net$MEs
    plots <- list()  # new empty list
    # dpa=as.factor(rep(seq(10,40,by=10),each=3))
    ss = as.factor(sample[[set]])
    for(me in 0:(Nmodules-1)) {
       which.module=paste("ME",me,sep="")
       module.color=labels2colors(me)
       #heatmap
       par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
       plotMat(t(scale(subDat[,net$colors==me ]) ),
               nrgcols=30,rlabels=T,rcols=module.color,
               main=paste(which.module, module.color, sep=": "), cex.main=2)
       #barplot
       par(mar=c(5, 4.2, 0, 0.7))
       barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
       ylab="eigengene expression",xlab="seed development (dpa)", names.arg=as.character(ss))
       #line, anova
       df<-data.frame(ME=MEs[,which.module], ss, module = which.module )
       fit<-aov(ME~ss,df)
       dfc<-summarySE(df, measurevar="ME", groupvars=c("ss", "module"))
       dfc$genome <- gsub(".0","",dfc$ss)
       dfc$dpa <-gsub("A|D|syn|TM|Tom|Yuc","",dfc$ss)
       plots[[me+1]]<- ggplot(dfc, aes(x=genome, y=ME, fill = dpa)) +
                       geom_bar(stat="identity",position=position_dodge(), color="black", size=0.3) +
       geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.3,position=position_dodge(0.9)) +
       ggtitle(paste(which.module," ",module.color,", P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
       theme_bw() +
       theme(plot.title=element_text( size=11),legend.position = "none")
    }
    for(page in 1:ceiling(Nmodules/9))
    {
        if(Nmodules>(9*page))
        {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
        else
        {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    }

   plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)

   dev.off()# end the big plot
}

# "all_power12_TOM-block.1.RData"
# "Number of modules in all network is 36"
# "tri_power12_TOM-block.1.RData"
# "Number of modules in tri network is 338"
# "di_power12_TOM-block.1.RData"
# "Number of modules in di network is 27"
# "poly_power12_TOM-block.1.RData"
# "Number of modules in poly network is 31"

# "setT_power12_TOM-block.1.RData"
# "Number of modules in setT network is 57"
# "setF_power12_TOM-block.1.RData"
# "Number of modules in setF network is 56"
# "setP_power12_TOM-block.1.RData"
# "Number of modules in setP network is 62"

save(sample, netsetT, netsetF, netsetP, file = "R-04-networkTopology.RData")







############### Step 5.  Relate modules to phenotype and functional gene sets  ###############
## nohup R CMD BATCH s5.ModuleInterpretation.R &
################


library(WGCNA);
library(flashClust);
library(RColorBrewer);
options(stringsAsFactors = FALSE);

remove(list=ls())

lnames = load(file = "R-01-dataInput.RData")
lnames     #  "multiExpr"   "nSets"       "setLabels"   "shortLabels", "datTraits"
nSets      # 3
setLabels  #  "Five genomes"         "Additional synthetic" "Proteomic set"
shortLabels# "setT" "setF" "setP"
nGenes<-checkSets(multiExpr)$nGenes #36560
nGenes


load('R-04-networkTopology.RData')

# work on only the network with all and synthetic
net<-netsetF
# get eigengenes
MEs<-net$MEs
# eigengene~sample, anova
samples <- sample[[2]]
ss<-as.factor(samples)
pval<-apply(MEs,2,function(x){round(anova(aov(x~ss) )$"Pr(>F)"[1],4)})
pval<-as.data.frame(pval)
pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
pval$numeric<-as.numeric(substring(rownames(pval),3) )
pval<-pval[order(pval$numeric),]
pval$symbol[1]<-" "  # ME0 always meaningless
pval




# assign mediate values of A and D to syn
syn<-(datTraits[c(1:5,7:12),] + datTraits[13:23,] )/2
rownames(syn)<-gsub("A2","syn",rownames(syn) )
datTraits4<-rbind(datTraits,syn)

# define a gene significance variable. We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait.
GS = cor(multiExpr[[2]]$data,datTraits4,use="p")
# This translates the numeric values into colors
GS.Color=numbers2colors(GS,signed=T)

# conduct DE analysis
## source('deseq2.r', chdir = TRUE)

# inport differential expression results
load('DEid.Rdata')

pdf("s5.relateTraits.pdf")
sigDEs2colors<-function(x,sigC) { ifelse( colnames(multiExpr[[1]]$data) %in% x, sigC, "white")}
DE_dev.Color   <- as.data.frame( sapply(DEid_dev,    function(x) sigDEs2colors(x, sigC="purple")) )
DE_con1.Color<- as.data.frame( sapply(DEid_con[1:3], function(x) sigDEs2colors(x, sigC="darkgreen")) )
DE_con2.Color<- as.data.frame( sapply(DEid_con[4:10], function(x) sigDEs2colors(x, sigC="brown")) )


#  Hierarchical cluster tree (average linkage, dissTOM) of the 37222 proteins. The first color bands provide a simple visual comparison of module assignments (branch cuttings) based on the dynamic hybrid branch cutting method. Other color bands visualize the gene significance measure: "red" indicates a high positive correlation with seed weight, oil content and dpa. Note that the brown, blue and red module contain many genes that have high positive correlations with traits.
plotColor <- cbind( labels2colors(net$colors), GS.Color, DE_dev.Color, DE_con1.Color, DE_con2.Color)
plotLabel <- c("module",paste("GS.",colnames(GS),sep=""), paste("DE.dev",names(DE_dev.Color),sep=""),  paste("DE.",names(DEid_con),sep="") )

plotDendroAndColors(net$dendrograms[[1]],colors=plotColor, groupLabels= plotLabel ,dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05, main="Cluster dendrogram with gene significance measure")


# Relate eigengenes to external traits or sample conditions
# Add the weight to existing module eigengenes
MET=orderMEs(cbind(MEs,datTraits4))
#Visualization of the eigengene network representing the relationships among the modules and sample traits. The top panel shows a hierarchical clustering dendrogram of the eigengenes based on the dissimilarity diss(q_1,q_2)=1-cor(E^{(q_1)},E^{(q_2)}). The bottom panel shows the shows the eigengene adjacency A_{q1,q2}=0.5+0.5 cor(E^{(q_1)},E^{(q_2)}).
plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)

# plot eigengene network for only the significant modules
module.sig<-rownames(pval[pval$symbol=="*",])
MET.sig<-MET[,module.sig]
plotEigengeneNetworks(MET.sig,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
# plot it again with color names
names(MET.sig)<-paste("ME",labels2colors(as.numeric(substring(names(MET.sig),3) ) ),sep="")
plotEigengeneNetworks(MET.sig,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)

######### I want to plot a color bar for the sig only dendrogram ##
sigModule<-as.numeric(gsub("ME","",rownames(pval[pval$symbol=="*",]) ) )
ch.col <- labels2colors(sigModule)
plot(1:length(sigModule), 1:length(sigModule), type="n", yaxt="n", ylab="")
for (k in 1:length(sigModule)) {   rect(k, 1, k+1,  2, col = ch.col[k], border=NA) }
#########

# graphical representation for corelation with modules
moduleTraitCor = cor(MEs, datTraits4, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples=48);
MEcolors<-paste("ME",labels2colors(as.numeric(gsub("ME","",names(MEs)) ) ), sep="")
# sizeGrWindow(5,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1), mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
# Table of module-trait correlations and p-values. Each cell reports the correlation (and p-value) resulting from  correlating module eigengenes (rows) to traits (columns). The table is color-coded by correlation according to the color legend.
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(GS), yLabels = MEcolors, ySymbols = names(MEs), colorLabels = TRUE, colors = blueWhiteRed(50), textMatrix = as.matrix(textMatrix), setStdMargins = FALSE, cex.text = 0.7,zlim = c(-1,1), main = paste("Module-trait relationships"))
# plot it for only significant modules
where.sig<-sort(match(module.sig, rownames(moduleTraitCor)) )
moduleTraitCor.sig <- moduleTraitCor[where.sig,]
textMatrix.sig <- textMatrix[where.sig,]
labeledHeatmap(Matrix = moduleTraitCor.sig, xLabels = colnames(GS), yLabels = MEcolors[where.sig], ySymbols = rownames(moduleTraitCor.sig), colorLabels = TRUE, colors = blueWhiteRed(50), textMatrix = as.matrix(textMatrix.sig), setStdMargins = FALSE, cex.text = 0.7,zlim = c(-1,1), main = paste("Module-trait relationships: sig only"))


# For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.
# calculate the module membership values
# (aka. module eigengene based connectivity kME):
MM=signedKME(multiExpr[[2]]$data, MEs)
rownames(MM)<-names(multiExpr[[2]]$data)
# Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules. As an example, we look at the brown module that a high correlation with body weight. We plot a scatterplot of Gene Significance vs. Module Membership in select modules...
colorOfColumn=substring(names(MM),4)
par(mfrow = c(3,1))
for (module in sigModule) {
    column = match(module,colorOfColumn)
    restModule=net$colors==module
    for (trait in colnames(GS)){
        verboseScatterplot(MM[restModule,column], GS[restModule,trait], xlab=paste("Module Membership of ME",module),ylab=paste("GS.",trait),    main=paste("kME.",module,"vs. GS"), col=labels2colors(module))
    }
}

dev.off()


# write gene with corresponding module assignment and annotation
aa<-load('D5annotation.Rdata')
aa  # "annotation221" "annot"
aa<-annotation221[,c("transcript","tair10.defline", "pfam","panther","kog", "kegg.ec", "kegg.orthology", "go", "tair10", "tair10.symbol")]
dim(aa<-aa[grep("[.]1$",aa$transcript),])  #37505
aa$gene = gsub("[.]1$","", aa$transcript)
me<-data.frame(gene = colnames(multiExpr[[2]]$data), ME = net$colors)
dim(me<-merge(me, aa, all.x=TRUE, by="gene"))
write.table(me, file="s5.module&annotation.txt", row.names=FALSE,sep="\t")



## Functional enrichment analysis
# add Gorai id to module group
MM$ID<-colnames(multiExpr[[2]]$data)
for(me in module.sig)
{
    me<-paste("k",me,sep="")
    write.table(MM[,c("ID",me)], file=paste("GSEA/MM_rnk/",me,".rnk",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
# Run GseaPreranked analysis
MM<-MM[,order(names(MM))]
write.table(MM,file="s5.moduleMembership.txt",sep="\t",row.names=FALSE,quote=FALSE)
save(net, datTraits4, pval, GS, MM, file = "R-05-all12_GS&MM.RData")


me0<-me
names(MM)[1]<-"gene"
dim(me0<-merge(me0,MM,by="gene",all.x=TRUE, all.y=TRUE) )
# me0[me0$kME9>0.9 & me0$ME==9, c("gene","tair10.defline")]
write.table(me0,file="s5.moduleMembership&annotation.txt",sep="\t",row.names=FALSE)

## GseaPreranked analysis
# http://www.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html?_GSEAPreranked_Page
#
# 1. launch javaGSEA desktop Application (http://www.broadinstitute.org/gsea/downloads.jsp)
# 2. Click “Load data” tab on the left “Steps in GSEA analysis” tool bar
# 3. Load ranked gene lists and Gene sets files
#    - Gene Set database: “gorai_051515.gmt”, combined from GO, mapman (similar to KEGG pathways), and Ran’s FA related gene sets. GMT format, each row represents a gene set
#    - Ranked gene list: “MM_rnk/kME1”, lists of module representative Gorai and module membership kME.
# 4. Go to top bar “Tools” and click “GseaPreranked”
#    Required fields
#    - Gene set database: gorai_051515.gmt
#    - Number of permutations: 1000
#    - Ranked List: kME1,etc.
#    - Collapse dataset to gene symbols: FALSE
#    Basic fields
#    - Analysis name: “kME1”
#    - Enrichment statistic: weighted (default): p=1
#    - Max size: 200
#    - Min size: 5 After filtering from the gene sets any gene not in the expression dataset, gene sets smaller than this are excluded from the analysis.
#    - Save results in this folder: #####
# 5. Use Cytoscape EnrichmentMap to visualize results, p<0.05, q<0.05


table(net$colors)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
# 2420 5714 4167 2195 1937 1852 1737 1704 1595 1575 1489 1442 1162  984  749  736
# 16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31
# 653  476  402  385  370  304  272  134  134  133  128  127  107  104  101   89
# 32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47
# 86   84   79   69   59   59   57   56   50   47   47   44   44   44   43   43
# 48   49   50   51   52   53   54   55
# 40   36   35   34   33   33   31   30
####### NOTE that the highest rank of significant gene sets soon becomes bigger than module size in small modules, which makes very weak association between module and gene sets.

# write GseaPreranked results to ouput files
library(xlsx)
# folders
folders<-grep("kME",list.files("GSEA/") , value=TRUE)
gseaF<-data.frame()
for(ff in folders )
{
    path<-paste("GSEA/",ff,"/",sep="")
    file<-grep("^gsea.*pos.*xls",list.files(path), value=TRUE)
    pos<-read.delim(paste(path,file,sep=""))
    print(ff)
    pos<-pos[pos$NOM.p.val<0.05 & pos$FDR.q.val<0.05,]
    print( dim(pos) )
    if(dim(pos)[1]>0){
        pos$module<-gsub("k|[.].*","",ff)
        gseaF<-rbind(gseaF, pos)
    }

}
gseaF<-gseaF[,c(13,1:12)]
write.table(gseaF, file="s5.gseaResults.txt", row.names=FALSE,sep="\t", quote=FALSE)



# Annotated moduls with topGO
library(topGO)
load('D5annotation.Rdata')
# "annotation221" "annot"
load('cottonGOenrich.RData')
# "geneID2GO"       "go4profile"  "goraiUniverse"   "gsc"             "level1_BP_terms" "level1_MF_terms"  "level2_BP_terms" "level2_MF_terms" "level3_BP_terms" "level3_MF_terms"   "level4_BP_terms" "level4_MF_terms" "level5_BP_terms" "level5_MF_terms"
# all included in individual network
universe<-colnames(multiExpr[[2]]$data)
# for each module containing genes of interest
GOresults<-data.frame()
for(module in sigModule)
{
    genes<-universe[net$colors==module]
    geneList <- factor(as.integer(universe %in% genes))
    names(geneList) <- universe
    
    pdf(file=paste("topGO/ME",module,".pdf", sep=""))
    
    # topGO analysis
    remove(enrich)
    for(on in c("MF","BP","CC"))
    {
        print(on)
        # Make topGO object
        GOdata <- new("topGOdata", ontology = on, allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        # fisher test
        result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        results.table <- GenTable(GOdata, result, topNodes = length(result@score))
        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
        results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
        # label ontology type
        results.table$ontology<-on
        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
        keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
        if(exists("enrich")) enrich<- rbind(enrich, keep)
        if(!exists("enrich")) enrich<- keep
        
        # draw figure for GO terms pval<=0.05 before FDR correction
        if(is.na(sigNo<-length(keep$ontology))){next}
        showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
        mtext(on, line=-1)
    }
    dev.off()
    if(dim(enrich)[1]>0)
    {
        enrichME<-enrich
        enrichME$ME=module
        GOresults<-rbind(GOresults,enrichME)   }
}
write.table(GOresults, file="s8.consensus_topGO/GOresults.txt", sep="\t", row.names=FALSE)









############### Step 6. export networks for cytoscape  ###############
nohup R CMD BATCH s6......R &
################
library(WGCNA);
library(flashClust);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE);

remove(list=ls())
load('R-05-all12_GS&MM.RData') # net, datTraits4, pval, GS, MM
load('R-01-dataInput.RData')   # multiExpr, nSets, setLabels, shortLabels, datTraits"


# Extract only FA genes
FAs<-read.table("FAs.txt",header=TRUE,sep="\t")
dim(FAs) #657
length(FAgenes<-unique(FAs$nodeName)) # only 657 unique genes, oops
probes = colnames(multiExpr[[2]]$data)   #36560
# locate FA genes
asFAs = is.finite(match(probes, FAgenes))
table(asFAs)
# FALSE  TRUE
# 35907   653


# locate FAs in modules
moduleFAs<-as.data.frame(table(net$colors[asFAs]))
names(moduleFAs)<-c("moduleLabels","FAs")
moduleAll <-as.data.frame(table(net$colors))
names(moduleAll)<-c("moduleLabels","All")
moduleFAs<-merge(moduleFAs, moduleAll, by="moduleLabels" ,all.y=TRUE)
moduleFAs$FAs[is.na(moduleFAs$FAs)]<-0
## calculate enrichment
tt<-colSums(moduleFAs[,2:3])
moduleFAs$fisherP<-round( apply(moduleFAs[,2:3],1,
function(x)fisher.test(matrix(as.numeric(c(x,tt-x)), nrow = 2, dimnames = list( c( "FAs","all"),c("inModule", "out"))) ,alternative="greater" )$p.value) ,3)
# FAs enriched in below modules
moduleFAs[moduleFAs$fisherP<0.05,]
# moduleLabels FAs  All fisherP
#             1 137 5714   0.000
#             2  97 4167   0.005
#             4  48 1937   0.016
#            11  43 1442   0.001
#            18  13  402   0.033
#            23   6  134   0.037

# list modules for each FA gene
modules<-data.frame(nodeName=probes, module=net$colors)
FAs<-merge(FAs, modules, by="nodeName", all.x=TRUE)
dim(FAs<-FAs[!is.na(FAs$module),])
unique(FAs$Biological.process)
#[1] "Plastid Fatty Acid Synthesis from Pyruvate"                     "β-Oxidation"
#[3] "Other Acyl Lipid Related"                                       "Eukaryotic Glycerolipid Synthesis"
#[5] "Oil Storage"                                                    "Sphingolipid Synthesis"
#[7] "Plastidial Glycerolipid, Galactolipid and Sulfolipid Synthesis" "Transcription Factors Associated with Lipid Synthesis"
#[9] "Mitochondrial Fatty Acid and Lipoic Acid Synthesis"             "TAG synthesis"
#[11] "Lipid Transfer Proteins"                                        "Miscellaneous"
#[13] "Fatty Acid Elongation"                                          "Fatty acid desaturation (out of plastid)"
#[15] "Plastid Lipid Trafficking"
xtabs(~Biological.process+module,data=FAs)
# now I want to do fishers test to show which module is each category enriched at
pdf("s6.FAs2modules.pdf",width=10,height=7)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
coln<-dim(pval)[1]                             # put modules in columns
rown<-length(unique(FAs$Biological.process) )  # put FAs categories in rows
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(NA, nrow = rown, ncol = coln);
CountTbl = matrix(0, nrow = rown, ncol = coln);
# color list of MEs in the color of decreasing numbers of memebers
colModule  <- as.numeric(gsub("ME","",rownames(pval)))
# category
rowCategory <- unique(FAs$Biological.process)
# colors for each gene
colColors  <- labels2colors(1:coln)
# anova significance sybol
colP  <- pval$symbol
# Execute all pairwaise comparisons
for (rmod in 1:rown)
    for (cmod in 1:coln)
    {
        rMembers = (FAs$Biological.process == rowCategory[rmod] );
        cMembers = (FAs$module == colModule[cmod] );
        CountTbl[rmod, cmod] = sum(FAs$Biological.process == rowCategory[rmod]  & FAs$module == colModule[cmod]   )
        if(CountTbl[rmod, cmod]>0) pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
    }
# display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
    # Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
rModTotals = apply(CountTbl, 1, sum)
cModTotals = apply(CountTbl, 2, sum)
select<-(cModTotals>0)
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap( Matrix = pTable[,select], colorLabels = TRUE,
    xLabels = paste(colP[select],rownames(pval)[select]), yLabels = paste(" ", rowCategory),
    textMatrix = CountTbl[,select], colors = blueWhiteRed(100)[50:100],
    main = "Correspondence of FA-related gene categories to modules ",
    cex.text = 0.5, cex.lab = 0.5, setStdMargins = FALSE      )
dev.off()

# check some FA families
xtabs(~Protein.Gene.Abbreviation+module, FAs[FAs$Biological.process=="Lipid Transfer Proteins",])
#                            module
#   Protein.Gene.Abbreviation  0  1  2  3  4  5  7  8  9 10 11 13 16 18 25 31 35 52
#                        ACT   0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#                        LTP1  2  5  5  0  0  1  0  2  0  1  1  1  0  1  1  0  0  0
#                        LTP2  1  1  3  0  0  0  0  0  1  0  1  0  1  0  0  0  0  0
#                        LTP3  1  3  1  1  0  0  0  1  1  0  0  2  1  0  0  0  1  1
#                        LTP4  1  0  1  2  0  1  0  0  0  0  1  0  0  0  0  1  0  0
#                        LTP5  1 11  9  1  1  2  3  1  0  0  3  1  1  0  0  0  0  0
xtabs(~Protein.Gene.Abbreviation+module, FAs[FAs$Biological.process=="Transcription Factors Associated with Lipid Synthesis",])
#                          module
#   Protein.Gene.Abbreviation 0 2 8 9 11 12 13 28 44
#                        ABI3 0 1 0 0  0  0  0  0  0
#                        FUS3 0 1 0 0  0  0  0  0  0
#                        HSL1 0 1 0 1  0  0  2  0  0
#                        HSL2 0 1 0 0  0  0  0  0  0
#                        LEC1 0 0 1 0  0  0  1  0  0
#                        LEC2 0 0 0 0  0  0  2  0  0
#                        PKL  0 0 0 0  0  1  0  1  0
#                        WRI1 1 1 0 0  2  0  0  0  1



# Export FA network to cytoscape
# connect to speedy, or whatever used to build network and store the TOM file (ISOLONE, biocrunch, etc.)
# Cytoscape [2] allows the user to input an edge file and a node file, allowing the user to specify for example the link weights and the node colors.
## ssh hugj2006@speedy.ent.iastate.edu
## R
# Get topological overlap, "TOM"
load("setF_power12_TOM-block.1.RData")   #local use
# Select the corresponding Topological Overlap
subTOM = as.matrix(TOM)[asFAs, asFAs];
str(subTOM)
subProbes = probes[asFAs];
dimnames(subTOM) = list(subProbes, subProbes)
subColors<-net$colors[asFAs]
quantile(TOM, 1:10/10)
aa<-read.table("s5.moduleMembership&annotation.txt", header=TRUE,sep="\t")
rownames(aa)<-aa$gene
aa_FAs<-aa[subProbes,]
aa_FAs<-aa_FAs[,-1]
table(rownames(aa_FAs) == FAs$nodeName )
nn<-cbind(FAs,aa_FAs)
nn$module.color<-labels2colors(nn$module)
library(gplots)
nn$module.hex<-col2hex(nn$module.color)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(subTOM,
    edgeFile = "CytoscapeInput-edges-FAs.txt",
    nodeFile = "CytoscapeInput-nodes-FAs.txt",
    weighted = TRUE,
    threshold = 0.02,   #
    nodeNames = nn$nodeName,
    nodeAttr = nn[,-1] )



############### Step 7. individual genome analysis  ###############
nohup R CMD BATCH s7......R &
################
library(WGCNA);
library(RColorBrewer)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# construct individual networks for each
load("R-01-dataInput.RData") # multiExpr, nSets, setLabels, shortLabels, datTraits"
datExprT<-multiExpr[[2]]$data    # setF
# make subsets for each genome
Adata<-datExprT[grep("A2",rownames(datExprT) ),]
Ddata<-datExprT[grep("D5",rownames(datExprT) ),]
TM1data<-datExprT[grep("TM1",rownames(datExprT) ),]
Yucdata<-datExprT[grep("Yuc",rownames(datExprT) ),]
AD3data<-datExprT[grep("AD3",rownames(datExprT) ),]
ADsdata<-datExprT[grep("ADs",rownames(datExprT) ),]

# We work with six sets:
nSets = 6
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development", "Syn seed devleopment")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3", "Syn")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = Adata);
multiExpr[[2]] = list(data = Ddata);
multiExpr[[3]] = list(data = TM1data);
multiExpr[[4]] = list(data = Yucdata);
multiExpr[[5]] = list(data = AD3data);
multiExpr[[6]] = list(data = ADsdata);
# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding 14670 genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
# Update exprSize
checkSets(multiExpr)
# $nSets
# [1] 6
# $nGenes
# [1] 36560
# $nSamples
# [1] 12 11 12 12 12 11
# $structureOK
# [1] TRUE


# Still too many genes to make networks, reduce the size by applying anova test
lanova<-data.frame(gene=colnames(multiExpr[[1]]$data))
# 12 samples
dpa=as.factor(rep(seq(10,40,by=10),each=3))
panova<-function(x){anova(aov(x~dpa))$"Pr(>F)"[1]}  # p-value of anova tests
for (set in c(1,3:5)){ lanova[,shortLabels[set]] <- apply(multiExpr[[set]]$data, 2, panova ) }
# 11 samples
dpa=dpa[c(1:5,7:12)]
panova<-function(x){anova(aov(x~dpa))$"Pr(>F)"[1]}  # p-value of anova tests
for (set in c(2,6)){ lanova[,shortLabels[set]] <- apply(multiExpr[[set]]$data, 2, panova ) }
# take min P
lanova$minP<-apply(lanova[,-1],1,min )
table(lanova$A2<0.05)   # FALSE  TRUE 21098 15462
table(lanova$D5<0.05)   # FALSE  TRUE 22598 13962
table(lanova$TM1<0.05)  # FALSE  TRUE 18460 18100
table(lanova$Yuc<0.05)  # FALSE  TRUE 24270 12290
table(lanova$AD3<0.05)  # FALSE  TRUE 16809 19751
table(lanova$minP<0.05) # FALSE  TRUE 5969 30591
multiExpr0 <- multiExpr  #save previous as 0
for (set in 1:nSets)
{
    multiExpr[[set]]$data = multiExpr[[set]]$data[, lanova$minP<0.05];
}
# Update exprSize
checkSets(multiExpr)
# $nSets 5
# $nGenes 30591
# $nSamples 12 11 12 12 11
# $structureOK  TRUE

# We now cluster the samples on their Euclidean distance, separately in each set.
pdf(file = "s7.SampleClustering0.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr0[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
pdf(file = "s7.SampleClusteringS.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;

save(multiExpr, lanova, file = "R-07-prep.RData")


# Choose a set of soft-thresholding powers, use "signed" network for this analysis
type<-"signed"
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for(set in 1:nSets){
powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type)[[2]])      }
collectGarbage()
    
# Plot the results:
colors=brewer.pal(nSets,"Set1")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
    {
        for (col in 1:length(plotCols))
        {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        }
}
    
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(paste("s7.ChooseSoftThresholdPower_",gsub(".* ","", type), ".pdf", sep="") )
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
        if (set==1)
        {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
            addGrid()
        }
        if (col==1)
        {
            text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
            labels=powers,cex=cex1,col=colors[set]);
        } else
        text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
        if (col==1)
        {
            legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
        } else
        legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
dev.off()
assign(paste("powerTables.",gsub(".* ","", type),sep=""),powerTables)
# examine power choice figure, default =12 is still reasonable
save(multiExpr, lanova, powerTables.signed, file = "R-07-prep.RData")

# Construct networks
nGenes<-checkSets(multiExpr)$nGenes
nGenes
powers<-12
# work with individual genome, then work with different soft threshold
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    for (j in powers )
    {
        softPower = j
        print(paste("Start building network for ",shortLabels[set]," using soft threshold ",j,"......",sep=""))
        # Network construction
        
        net = blockwiseModules(
             # Input data
             subDat,
             # Data checking options
             checkMissingData = TRUE,
             
             # Options for splitting data into blocks
             blocks = NULL,
             randomSeed = 12345,
             maxBlockSize = nGenes,  # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
             
             # Network construction arguments: correlation options, use bicor instead of default pearson
             corType = "pearson",
             # Adjacency and topology overlap function options
             power = j, networkType = "signed", TOMType = "signed",
             
             # Saving or returning TOM
             saveTOMs = TRUE,
             saveTOMFileBase = paste(shortLabels[set],"_power",j,"_TOM",sep=""),
             
             # Basic tree cut options
             deepSplit = 2,  #default, known to reasonable
             minModuleSize = min(30, ncol(subDat)/2 ), #default 20, use 30 for transcriptome
             pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
             
             # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
             mergeCutHeight = 0.25,
             
             # others
             reassignThreshold = 0,
             numericLabels = TRUE,
             verbose = 3)
             
        assign(paste(shortLabels[set],"net",j,sep=""), net)
        }
}
save(multiExpr, lanova, powerTables.signed, file = "R-07-prep.RData")
save(list=grep(".+net.+",ls(), value=TRUE), file = "R-07-individualNetworks.RData")


# Module level analysis
source('multiscalefreeplot.r', chdir = TRUE)
source('summarySE.r', chdir = TRUE)
source('multiplot.r', chdir = TRUE)

# 12 samples
dpa12=as.factor(rep(seq(10,40,by=10),each=3))
# 11 samples
dpa11=dpa12[c(1:5,7:12)]
# make dpa data frame
dpas<-list(A2=dpa12,D5=dpa11,TM1=dpa12,Yuc=dpa12,AD3=dpa12,Syn=dpa11)

library(ggplot2)
for (set in 1:nSets )
{
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    genome <-  shortLabels[set]
    net<-get(paste(genome,"net12",sep="") )
    MEs<-net$MEs
    nSamples = checkSets(multiExpr)$nSamples[set]
    dpa <- dpas[[set]]
    # gg<-net$goodGenes    #get the good genes descision
    # adjacency = adjacency(subDat, power = 20, type = "signed")
    Nmodules= dim(net$MEs)[2]
    print(paste("Number of modules in ",genome," network is ",Nmodules,sep=""))
    pdf(paste("s7.",genome,"_modules.pdf",sep=""))
    
    # Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
    # Displaying module heatmap and the eigengene
    # sizeGrWindow(8,7);

    plots <- list()  # new empty list
    for(me in 0:(Nmodules-1)) {
        which.module=paste("ME",me,sep="")
        module.color=labels2colors(me)
        #heatmap
        par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
        plotMat(t(scale(subDat[,net$colors==me ]) ),
        nrgcols=30,rlabels=T,rcols=module.color,
        main=paste(which.module, module.color, sep=": "), cex.main=2)
        #barplot
        par(mar=c(5, 4.2, 0, 0.7))
        barplot(MEs[,which.module], col=module.color, main="", cex.main=2, ylab="eigengene expression",xlab="seed development (dpa)", names.arg=dpa)
        #line, anova
        df<-data.frame(ME=MEs[,which.module], dpa, module = which.module )
        fit<-aov(ME~dpa,df)
        dfc<-summarySE(df, measurevar="ME", groupvars=c("dpa", "module"))
        plots[[me+1]]<- ggplot(dfc, aes(x=dpa, y=ME, group=module)) +
        geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.1) +
        geom_line(colour=module.color) + geom_point( ) +
        ggtitle(paste(which.module," ",module.color,", anova P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
        theme(plot.title=element_text( size=11))
    }
    for(page in 1:ceiling(Nmodules/9))
    {
    if(Nmodules>(9*page))
    {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    else
    {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    }
    
    plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
    
    # Relate eigengenes to external traits or sample conditions
    moduleTraitCor = cor(MEs, dpa, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
    # graphical representation for corelation
    # sizeGrWindow(5,6)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    par(mfrow=c(1,1), mar = c(6, 8.5, 3, 3))
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor, xLabels = "dpa", yLabels = names(MEs),
    ySymbols = names(MEs), colorLabels = TRUE, colors = blueWhiteRed(50),
    textMatrix = as.matrix(textMatrix), setStdMargins = FALSE,
    cex.text = 0.7,zlim = c(-1,1), main = paste("Module-development relationships"))
    # Another way to visualize the correlations with eigegene dendrogram and adjacency heatmap, BUT the heatmap color bar is NOT RIGHT, should be (-1, 1)
    MET = orderMEs(cbind(MEs, as.numeric(as.character(dpa))))
    # sizeGrWindow(5,7.5);
    par(cex = 0.9)
    plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
    
    dev.off()# end the big plot
}
# "Number of modules in A2 network is 39"
# "Number of modules in D5 network is 23"
# "Number of modules in TM1 network is 27"
# "Number of modules in Yuc network is 38"
# "Number of modules in AD3 network is 32"
# "Number of modules in Syn network is 30"




############### Step 8. cross-genome comparison  ###############
nohup R CMD BATCH s7......R &
################
library(WGCNA);
library(flashClust);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE);

remove(list=ls())
load("R-07-prep.RData") # "multiExpr"          "lanova"             "powerTables.signed"
load("R-07-individualNetworks.RData") # "A2net12"  "AD3net12" "D5net12"  "Synnet12" "TM1net12" "Yucnet12"
source('multiscalefreeplot.r', chdir = TRUE)
source('summarySE.r', chdir = TRUE)
source('multiplot.r', chdir = TRUE)

checkSets(multiExpr)
nSets<-checkSets(multiExpr)$nSets   #6
nGenes<-checkSets(multiExpr)$nGenes #30591
nSamples<-checkSets(multiExpr)$nSamples[1] #12 11 12 12 12 11
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development", "Syn seed devleopment")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3", "Syn")

softPower = 12

# put all net in one list in correct order
nets<-list(A2=A2net12,D5=D5net12,TM1=TM1net12,Yuc=Yucnet12,AD3=AD3net12, Syn=Synnet12)
names(nets)==shortLabels  # TRUE make sure correct order
all.colors<-cbind(A2net12$colors, D5net12$colors, TM1net12$colors, Yucnet12$colors, AD3net12$colors, Synnet12$colors)
# put all anova p together
# 12 samples
dpa12=as.factor(rep(seq(10,40,by=10),each=3))
# 11 samples
dpa11=dpa12[c(1:5,7:12)]
# make dpa data frame
dpas<-list(A2=dpa12,D5=dpa11,TM1=dpa12,Yuc=dpa12,AD3=dpa12,Syn=dpa11)
anovaP<-list()
for(i in 1:nSets)
{
    MEs<-nets[[i]]$MEs
    dpa<-dpas[[i]]
    pval<-apply(MEs,2,function(x){round(anova(aov(x~dpa) )$"Pr(>F)"[1],4)})
    pval<-as.data.frame(pval)
    pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
    pval$numeric<-as.numeric(substring(rownames(pval),3) )
    pval<-pval[order(pval$numeric),]
    pval$symbol[1]<-" "  # ME0 always meaningless
    anovaP[[i]]<-pval
}
names(anovaP)<-shortLabels


# for any consensus module analysis, define new CmultiExpr and TOMs
comparisons <- list( all =1:6          ## all five consensus
# diploid = c(1,2),  ## diploid consensus, A2 and D5
# polyploids =c(3,4,5,6),
# AD1 =c(3,4)      ## G. hirsutum
)

# RUN MY HUGE LOOP FOR ALL THESE COMPARISONS
for(li in 1:1){
   comp<-comparisons[[li]]
    print(paste("Consensus module analysis for ",names(comparisons)[li], ": ",sep=""))
    print(shortLabels[comp])
    
    
    ######### Start of my consensus analysis unit ###########
    nSets<-length(comp)
    CmultiExpr = vector(mode = "list", length = nSets)
    for(i in 1:nSets){ CmultiExpr[[i]] = multiExpr[[comp[i]]]}
    fileName<-paste(shortLabels[comp],sep="",collapse="_")
    
    print("Start to read in TOMs.")
    # Initialize an appropriate array to hold the TOMs
    TOMs = array(0, dim = c(nSets, nGenes, nGenes));
    # load and store TOMs from each individual data set
    for (set in 1:nSets)
    {
        #  load(paste("/Volumes/jfw-lab/home/jing/oilseedNetwork/concensusTotal/",shortLabels[comp[set]], "_power20_TOM-block.1.RData", sep="") )   #local use
        load(paste(shortLabels[comp[set]], "_power12_TOM-block.1.RData", sep="") )   #biocrunch
        TOMs[set, , ]<-as.matrix(TOM)
    }
    
    
    # Step-by-step consensus network construction
    # consensus is defined as the component-wise minimum of the multiple TOMs.
    
    # first, TOMs need to be scaled to be comparable, using the 95th percentile
    # BUT BE CAREFUL, for each individual network, original TOMs and scaled TOMs might generate different numbers of modules!!!!
    print("Scale TOMs using 95th percentile: values, scaling powers")
    scaleP = 0.95
    # Set RNG seed for reproducibility of sampling
    set.seed(12345)
    # Sample sufficiently large number of TOM entries
    nSamples = as.integer(1/(1-scaleP) * 1000)
    # Choose the sampled TOM entries
    scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
    TOMScalingSamples = list();
    # These are TOM values at reference percentile
    scaleQuant = rep(1, nSets)
    # Scaling powers to equalize reference TOM values
    scalePowers = rep(1, nSets)
    # prepare the TOM scaled
    TOMscaled <-TOMs
    # Loop over sets
    for (set in 1:nSets)
    {
        # Select the sampled TOM entries
        TOMScalingSamples[[set]] = as.dist(TOMs[set, , ])[scaleSample]
        # Calculate the 95th percentile
        scaleQuant[set] = quantile(TOMScalingSamples[[set]], probs = scaleP, type = 8)
        # Scale the other TOMs
        if (set>1)
        {
            scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
            TOMscaled[set, ,] = TOMs[set, ,]^scalePowers[set];
        }
    }
    # check the scaling achieved using a quantile-quantile plot
    scaledTOMSamples = list();
    for (set in 1:nSets)
    scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
    # Open a suitably sized graphics window
    pdf(paste("s8.checkTOMscaling-",fileName,".pdf",sep=""))
    pwset<-combn(nSets,2)
    for(i in 1:choose(nSets,2))
    {
        # qq plot of the unscaled samples
        qqUnscaled = qqplot(TOMScalingSamples[[pwset[1,i]]], TOMScalingSamples[[pwset[2,i]]], plot.it = TRUE, cex = 0.6, xlab = paste("TOM in", setLabels[comp[pwset[1,i]]]), ylab = paste("TOM in", setLabels[comp[pwset[2,i]]]), main = "Q-Q plot of TOM", pch = 20)
        # qq plot of the scaled samples
        qqScaled = qqplot(scaledTOMSamples[[pwset[1,i]]], scaledTOMSamples[[pwset[2,i]]], plot.it = FALSE)
        points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
        abline(a=0, b=1, col = "blue")
        legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
    }
    dev.off()
    # the red scaled quantile plot should be closer than black unscaled plot to the theoretic blue line
    print(scaleQuant)   # 0.2002059 0.2124997 0.2477379 0.2039054 0.2607832
    print(scalePowers)  # 1.000000 1.038477 1.152664 1.011515 1.196674
    
    # Second, calculate the consensus Topological Overlap by taking the component-wise ("parallel") minimum of the TOMs in individual sets.
    # Between diploids - diploid consensus network vs A2 and D5 each
    if(nSets==2) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ])
    if(nSets==3) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ])
    if(nSets==4) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ], TOMscaled[4, , ])
    if(nSets==5) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ], TOMscaled[4, , ], TOMscaled[5, , ])
    if(nSets==6) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ], TOMscaled[4, , ], TOMscaled[5, , ], TOMscaled[6, , ])
    # Clustering
    consTree = flashClust(as.dist(1-consensusTOM), method = "average")
    # We like large modules, so we set the minimum module size relatively high:
    minModuleSize = 30;
    # Module identification using dynamic tree cut:
    unmergedLabels = cutreeDynamic(dendro = consTree, distM = as.matrix(1-consensusTOM), deepSplit = 2, cutHeight = 0.995, minClusterSize = minModuleSize, pamRespectsDendro = TRUE )
    unmergedColors = labels2colors(unmergedLabels)
    # a quick summary of the module detection,
    print(paste("Before merging, ", length(unique(unmergedLabels)), " modules were resulted.",sep=""))
    # Calculate module eigengenes
    unmergedMEs = multiSetMEs(CmultiExpr, colors = NULL, universalColors = unmergedColors)
    # Calculate consensus dissimilarity of consensus module eigengenes
    consMEDiss = consensusMEDissimilarity(unmergedMEs);
    # Cluster consensus modules
    consMETree = flashClust(as.dist(consMEDiss), method = "average");
    # merge modules
    merge = mergeCloseModules(CmultiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)
    # Numeric module labels
    moduleLabels = merge$colors;
    # the corresponding colors for large module ID need to be adjusted, otherwise they cannot be plotted
    sortedModules = sort(unique(moduleLabels))
    moduleLabels.adjusted = match(moduleLabels,sortedModules)-1
    # Convert labels to colors
    moduleColors = labels2colors(moduleLabels.adjusted)
    # Eigengenes of the new merged modules:
    consMEs = merge$newMEs;
    # Calculate new module eigengenes
    mergedMEs = multiSetMEs(CmultiExpr, colors = NULL, universalColors = moduleLabels)
    # Calculate consensus dissimilarity of consensus module eigengenes, Cluster consensus modules
    consMETree.merged = flashClust(as.dist(consensusMEDissimilarity(mergedMEs) ), method = "average");
    
    # Save useful info for the consensus modules
    save(consMEs, unmergedLabels, moduleColors, moduleLabels, moduleLabels.adjusted, consTree, file = paste("s8.moduleConsensus.", fileName, ".RData", sep="") )
    print("Save consensus modules information and start to draw figures.")
    
    
    # Plot the consensus module results
    pdf(paste("s8.consensus-",fileName,".pdf",sep=""))
    par(mfrow = c(1,1))
    plot(consMETree, main = "Consensus clustering of consensus module eigengenes before merging",xlab = "", sub = "")
    abline(h=0.25, col = "red")
    plot(consMETree.merged, main = "Consensus clustering of consensus module eigengenes after merging",xlab = "", sub = "")
    plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors), c("Unmerged", "Merged"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    plotDendroAndColors(consTree, cbind(moduleColors, labels2colors(all.colors) ), c("Consensus",shortLabels[comp]), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Consensus gene dendrogram and module colors")
    # Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
    # Displaying module heatmap and the eigengene
    # sizeGrWindow(8,7);
    for(i in comp)
    {   MEs<-consMEs[[i]]$data
        Nmodules<-dim(MEs)[2]
        module.names<-names(MEs)
        dpa<-dpas[[i]]
        plots <- list()  # new empty list
        for(me in 1:Nmodules)
        {
            which.module=module.names[me]
            module.color=labels2colors(match(as.numeric(substring(which.module,3)),sortedModules)-1)
            #heatmap
            par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
            plotMat(t(scale(CmultiExpr[[i]]$data[,moduleColors==module.color ]) ), nrgcols=30,rlabels=T,rcols=module.color,
            main=paste(shortLabels[comp[i]],which.module, module.color, sep=": "), cex.main=2)
            #barplot
            par(mar=c(5, 4.2, 0, 0.7))
            barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
            ylab="eigengene expression",xlab="seed development (dpa)", names.arg=as.character(dpa) )
            #line, anova
            df<-data.frame(ME=MEs[,me], dpa, module = which.module )
            fit<-aov(ME~dpa,df)
            dfc<-summarySE(df, measurevar="ME", groupvars=c("dpa", "module"))
            plots[[me]]<- ggplot(dfc, aes(x=dpa, y=ME, group=module)) +
            geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.1) +
            geom_line(colour=module.color) + geom_point( ) +
            ggtitle(paste(which.module," ",module.color,", anova P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
            theme(plot.title=element_text( size=8))
        }
        for(page in 1:ceiling(Nmodules/9))
        {
            if(Nmodules>(9*page))
            {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
            else
            {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
        }
    }
    dev.off()
    
    
    # Calcultae and Plot consensus module preservation results
    pdf(paste("s8.modulePreservasion-",fileName,".pdf",sep=""),width=8, height=10)
    # Recalculate consMEs to give them color names, or
    # consMEs = multiSetMEs(multiExpr, universalColors = merge$colors);
    #sizeGrWindow(8,10);
    par(cex = 0.4)  #very import for cex
    plotEigengeneNetworks(consMEs, setLabels[comp], marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
    # Characterizing consensus modules by differential expression of their corresponding eigengenes in the various time points. Red means over-expression, green under-expression; numbers in each cell give the corresponding t-test p-value. Each column corresponds to an eigengene and each row corresponds to a time point.
    # setCorrelationPreservation(consMEs, setLabels[comp], excludeGrey = TRUE, greyLabel = "grey")
    # pairwise plot
    pwset<-combn(nSets,2)
    for(i in 1:choose(nSets,2) )
    {par(cex=0.4);
        plotEigengeneNetworks(list(consMEs[[pwset[1,i]]],consMEs[[pwset[2,i]]]), setLabels[c(pwset[1,i],pwset[2,i])], marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
    }
    # with GS and MM
    par(mfrow=c(4,1) )
    for(i in 1:choose(nSets,2) )
    {
        # the pair for comparison
        MEs1<-consMEs[[pwset[1,i]]]$data
        MEs2<-consMEs[[pwset[2,i]]]$data
        dpa1<-dpas[[pwset[1,i]]]
        dpa2<-dpas[[pwset[2,i]]]
        # diffTable <- apply(MEs1-MEs2,2, function(x){unlist(lapply(split(x,dpa),mean ))} )
        diffTable <- apply(MEs1,2, function(x){unlist(lapply(split(x,dpa1),mean ))} ) - apply(MEs2,2, function(x){unlist(lapply(split(x,dpa2),mean ))} )
        pvalTable <- matrix(0, nrow = 4, ncol = dim(MEs1)[2]);
        for(cc in 1:dim(MEs1)[2])
        {
            for(rr in 1:4)
            { pvalTable[rr,cc] = t.test(split(MEs1[,cc],dpa1)[[rr]],split(MEs2[,cc],dpa2)[[rr]])$p.value }
        }
        labeledHeatmap( Matrix = diffTable, xLabels = names(MEs1), yLabels = c(10,20,30,40),
        colorLabels = TRUE, colors = blueWhiteRed(50),
        textMatrix = round(as.matrix(pvalTable),2),
        cex.text = 0.7,  zlim = c(-0.8,0.8),setStdMargins = FALSE,
        main = paste("Consensus MEs differential expression: ", paste(shortLabels[comp[pwset[,i] ]],collapse=" vs " ), sep=""  ) )
    }
    dev.off()
    
    
    # Note that preservation measures can also be generaged through (setCorrelationPreservation)
    # do a permutation test for significance
    resultP<-data.frame(matrix(ncol=choose(nSets,2)))
    names(resultP)<-c()
    pwset<-combn(nSets,2)
    for(i in 1:choose(nSets,2) )
        {
        names(resultP)[i]<-paste(shortLabels[pwset[1,i]],"vs",shortLabels[pwset[2,i]],sep="")
        }
    # make 1000 times permutation
    nP = 1000
    for(i in 1:nP) {
    # permutation
    colorsP<-sample(moduleLabels, size=length(moduleLabels), replace=FALSE)
    # permutate module color assignment
    consMEsP = multiSetMEs(multiExpr, universalColors = colorsP);
    # given module colors, calculate conMEs
    dd<-setCorrelationPreservation(consMEsP, setLabels, excludeGrey = TRUE, greyLabel = "grey")
    resultP[i,]<-as.matrix(dd)[c(2:6,9:12,16:18,23:24,30)]
    }
    # observed top 2.5% values should be larger than
    apply(resultP, 2, function(x){sort(x, decreasing=TRUE)[25]})
    # A2vsD5   A2vsTM1   A2vsYuc   A2vsAD3   A2vsSyn   D5vsTM1   D5vsYuc   D5vsAD3
    # 0.6323412 0.7954535 0.7701976 0.7919692 0.7893501 0.6331537 0.7026363 0.6745426
    # D5vsSyn  TM1vsYuc  TM1vsAD3  TM1vsSyn  YucvsAD3  YucvsSyn  AD3vsSyn
    # 0.7951399 0.8111917 0.8373168 0.7414110 0.8384389 0.7931451 0.7507333
    apply(resultP, 2, function(x){sort(x, decreasing=FALSE)[25]})
    # A2vsD5   A2vsTM1   A2vsYuc   A2vsAD3   A2vsSyn   D5vsTM1   D5vsYuc   D5vsAD3
    # 0.5081714 0.5940392 0.5694629 0.5855468 0.5831047 0.5074908 0.5310748 0.5194917
    # D5vsSyn  TM1vsYuc  TM1vsAD3  TM1vsSyn  YucvsAD3  YucvsSyn  AD3vsSyn
    # 0.5886688 0.5928129 0.6167381 0.5482962 0.6171534 0.5840223 0.5574039
    ############### actual values are
    D<-setCorrelationPreservation(consMEs, setLabels[comp], excludeGrey = TRUE, greyLabel = "grey")
    #                           A2 seed development D5 seed development
    # A2 seed development            0.0000000           0.6787252
    # D5 seed development            0.6787252           0.0000000
    # TM1 seed development           0.7377839           0.6740817
    # Yuc seed development           0.7550656           0.7795420
    # AD3 seed development           0.7470093           0.7593282
    # Syn seed devleopment           0.7999661           0.8483083
    #                           TM1 seed development Yuc seed development
    # A2 seed development             0.7377839            0.7550656
    # D5 seed development             0.6740817            0.7795420
    # TM1 seed development            0.0000000            0.7779310
    # Yuc seed development            0.7779310            0.0000000
    # AD3 seed development            0.7689547            0.8661618
    # Syn seed devleopment            0.7246077            0.8138153
    #                           AD3 seed development Syn seed devleopment
    # A2 seed development             0.7470093            0.7999661
    # D5 seed development             0.7593282            0.8483083
    # TM1 seed development            0.7689547            0.7246077
    # Yuc seed development            0.8661618            0.8138153
    # AD3 seed development            0.0000000            0.7887938
    # Syn seed devleopment            0.7887938            0.0000000
    save(D, resultP, file = "s8.preservationD.RData" )

    
    # Plot marginal analysis between consensus modules and each individual dataset
    pdf(paste("s8.Consensus&marginal-",fileName,".pdf",sep=""),width=10,height=7)
    par(mfrow=c(1,1));
    par(cex = 1.0);
    par(mar=c(8, 10.4, 2.7, 1)+0.3);
    consMEs.no<-length(unique(moduleLabels))
    consMEs.module<-labels2colors(match(as.numeric(substring(names(consMEs[[1]]$data ),3) ), sortedModules)-1 )  #colors in order
    # loop pairwise comparison
    for(i in comp)  # compare to all 5 , no matter which consensus dealt with, but only look at those being compared!!!!!!!
    {
        coln<-consMEs.no                  #put consensus modules in columns
        rown<-ncol(nets[[i]]$MEs )  # put each individual set modules in rows
        # Initialize tables of p-values and of the corresponding counts
        pTable = matrix(0, nrow = rown, ncol = coln);
        CountTbl = matrix(0, nrow = rown, ncol = coln);
        # color list of MEs in the color of decreasing numbers of memebers
        colModule  <- consMEs.module
        rowModule  <- labels2colors(as.numeric(names(table(nets[[i]]$colors)) ))
        # colors for each gene
        colColors  <- moduleColors
        rowColors  <- labels2colors(nets[[i]]$colors )
        # anova significance sybol
        rowP  <- anovaP[[i]]$symbol
        # Initialize tables of p-values and of the corresponding counts
        pTable = matrix(0, nrow = rown, ncol = coln);
        CountTbl = matrix(0, nrow = rown, ncol = coln);
        # Execute all pairwaise comparisons
        for (rmod in 1:rown)
        for (cmod in 1:coln)
        {
            rMembers = (rowColors == rowModule[rmod] );
            cMembers = (colColors == colModule[cmod] );
            pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
            CountTbl[rmod, cmod] = sum(rowColors == rowModule[rmod]  & colColors == colModule[cmod]  )
        }
        
        # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
        # Truncate p values smaller than 10^{-50} to 10^{-50}
        pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
        pTable[pTable>50 ] = 50 ;
        # Marginal counts (really module sizes)
        rModTotals = apply(CountTbl, 1, sum)
        cModTotals = apply(CountTbl, 2, sum)
        # Use function labeledHeatmap to produce the color-coded table with all the trimmings
        labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
        xLabels = paste(" ", colModule), yLabels = paste(" ", rowModule),
        xSymbols = paste(names(consMEs[[1]]$data ),"-", colModule, ": ", cModTotals, " ", sep=""),
        ySymbols = paste(names(nets)[i], rowModule, ": ", rModTotals, rowP, sep=""),
        textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
        main = "Correspondence of dataset-specific to consensus modules ",
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
    }
    dev.off()
    
    # Plot marginal analysis between consensus modules and each individual dataset
    # comp<-comparisons[[li]]
    # nSets<-length(comp)
    # fileName<-paste(shortLabels[comp],sep="",collapse="_")
    # sortedModules = sort(unique(moduleLabels))
    pdf(paste("s8.Consensus&marginal-",fileName,".sig.pdf",sep=""),width=10,height=7)
    par(mfrow=c(1,1));
    par(cex = 1.0);
    par(mar=c(8, 10.4, 2.7, 1)+0.3);
    consMEs.no<-length(unique(moduleLabels))
    consMEs.module<-labels2colors(match(as.numeric(substring(names(consMEs[[1]]$data ),3) ), sortedModules)-1 )  #colors in order
    # loop pairwise comparison
    for(i in comp)  # compare to all 5 , no matter which consensus dealt with, but only look at those being compared!!!!!!!
    {
        coln<-consMEs.no                  #put consensus modules in columns
        rown<-ncol(nets[[i]]$MEs )  # put each individual set modules in rows
        # Initialize tables of p-values and of the corresponding counts
        pTable = matrix(0, nrow = rown, ncol = coln);
        CountTbl = matrix(0, nrow = rown, ncol = coln);
        # color list of MEs in the color of decreasing numbers of memebers
        colModule  <- consMEs.module
        rowModule  <- labels2colors(as.numeric(names(table(nets[[i]]$colors)) ))
        # colors for each gene
        colColors  <- moduleColors
        rowColors  <- labels2colors(nets[[i]]$colors )
        # anova significance sybol
        rowP  <- anovaP[[i]]$symbol
        # Initialize tables of p-values and of the corresponding counts
        pTable = matrix(0, nrow = rown, ncol = coln);
        CountTbl = matrix(0, nrow = rown, ncol = coln);
        # Execute all pairwaise comparisons
        for (rmod in 1:rown)
        for (cmod in 1:coln)
        {
            rMembers = (rowColors == rowModule[rmod] );
            cMembers = (colColors == colModule[cmod] );
            pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
            CountTbl[rmod, cmod] = sum(rowColors == rowModule[rmod]  & colColors == colModule[cmod]  )
        }
        
        # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
        # Truncate p values smaller than 10^{-50} to 10^{-50}
        pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
        pTable[pTable>50 ] = 50 ;
        # Marginal counts (really module sizes)
        rModTotals = apply(CountTbl, 1, sum)
        cModTotals = apply(CountTbl, 2, sum)
        # Use function labeledHeatmap to produce the color-coded table with all the trimmings
        labeledHeatmap( Matrix = pTable[rowP=="*",], colorLabels = TRUE,
        xLabels = paste(" ", colModule), yLabels = paste(" ", rowModule[rowP=="*"]),
        xSymbols = paste(names(consMEs[[1]]$data ),"-", colModule, ": ", cModTotals, " ", sep=""),
        ySymbols = paste(names(nets)[i], rowModule[rowP=="*"], ": ", rModTotals[rowP=="*"], rowP[rowP=="*"], sep=""),
        textMatrix = CountTbl[rowP=="*",], colors = blueWhiteRed(100)[50:100],
        main = "Correspondence of dataset-specific to consensus modules ",
        cex.text = 0.5, cex.lab = 0.5, setStdMargins = FALSE      )
    }
    dev.off()

    
    ######### End of my consensus analysis unit ###########
}

########################################################
###### Step 8 additional. module preservasion  #########
nohup R CMD BATCH s7......R &
################
library(WGCNA);
library(flashClust);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE);

remove(list=ls())
load("R-07-prep.RData") # "multiExpr"          "lanova"             "powerTables.signed"
load("R-07-individualNetworks.RData") # "A2net12"  "AD3net12" "D5net12"  "Synnet12" "TM1net12" "Yucnet12"
load('s8.moduleConsensus.A2_D5_TM1_Yuc_AD3_Syn.RData')

checkSets(multiExpr)
nSets<-checkSets(multiExpr)$nSets   #6
nGenes<-checkSets(multiExpr)$nGenes #30591
nSamples<-checkSets(multiExpr)$nSamples # 12 11 12 12 12 11
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development", "Syn seed devleopment")
shortLabels<-gsub(" .*","",setLabels)
pwset<-combn(nSets,2)

multiExpr.test<-multiExpr
names(multiExpr.test) <-shortLabels
multiColor1<-list(A2=moduleLabels, D5=moduleLabels, TM1=moduleLabels, Yuc=moduleLabels, AD3=moduleLabels, Syn=moduleLabels )

# The number of permutations drives the computation time of the module preservation function. For a publication use 200 permutations.
# But for brevity, let's use a small number
nPermutations1=200
# Set it to a low number (e.g. 3) if only the medianRank statistic and other observed statistics are needed.
# Permutations are only needed for calculating Zsummary and other permutation test statistics.
# set the random seed of the permutation test analysis
set.seed(1)
system.time({
    mp = modulePreservation(multiExpr, multiColor1, networkType="signed", referenceNetworks = pwset[1,], testNetworks=as.list(pwset[2,]), nPermutations = nPermutations1,
    randomSeed = 1, quickCor = 0, verbose = 3)
})
# Save the results of the module preservation analysis
save(mp, file = "s8.preservationZ&Medianrank.RData")
# view A2 vs D5 median rank
mp$preservation$observed[[1]][[2]] ->mr
mr[order(mr$medianRank.pres,decreasing=TRUE),1:2]

# If needed, reload the data:
load(file = "s8.preservationZ&Medianrank.RData")

pdf("s9.modulePreservation.pdf")
for(i in 1:choose(nSets,2) )
{
    # specify the reference and the test networks
    ref=pwset[1,i]; test = pwset[2,i]
    Obs.PreservationStats= mp$preservation$observed[[i]][[test]]
    Z.PreservationStats=mp$preservation$Z[[i]][[test]]
    # Look at the observed preservation statistics
    Obs.PreservationStats
    # Z statistics from the permutation test analysis
    Z.PreservationStats
    
    # Let us now visualize the data.
    modIDs = rownames(Obs.PreservationStats)
    modColors=labels2colors(order(as.numeric(modIDs) )-1 )
    moduleSize = Obs.PreservationStats$moduleSize
    # we will omit the grey module (background genes)
    # and the gold module (random sample of genes)
    selectModules = !(modColors %in% c("grey", "gold"))
    # Text labels for points
    point.label = modIDs[selectModules]
    # Composite preservation statistics
    medianRank=Obs.PreservationStats$medianRank.pres
    Zsummary=Z.PreservationStats$Zsummary.pres
    
    par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
    # plot medianRank versus module size
    plot(moduleSize[selectModules],medianRank[selectModules],col=1, bg=modColors[selectModules],
    pch = 21,main=paste("medianRank -",shortLabels[ref], "vs",shortLabels[test]),
    cex = 2, ylab ="medianRank",xlab="Module size", log="x")
    labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)
    
    # plot Zsummary versus module size
    plot(moduleSize[selectModules],Zsummary[selectModules], col = 1, bg=modColors[selectModules],pch = 21,
    main=paste("Zsummary -",shortLabels[ref], "vs",shortLabels[test]),
    cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
    labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
    # Add threshold lines for Zsummary
    abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
}
dev.off()



#################################################################
###### Step 9. functional analysis of consensus modules #########
nohup R CMD BATCH s7......R &
################
library(WGCNA);
library(flashClust);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE);

load("R-07-prep.RData") # "multiExpr"          "lanova"             "powerTables.signed"
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development", "Syn seed devleopment")
shortLabels<-gsub(" .*","",setLabels)

load("s8.moduleConsensus.A2_D5_TM1_Yuc_AD3_Syn.RData")

# Extract only FA genes
FAs<-read.table("FAs.txt",header=TRUE,sep="\t")
dim(FAs) #657
length(FAgenes<-unique(FAs$nodeName))
probes = colnames(multiExpr[[2]]$data)   #30591
# locate FA genes
asFAs = is.finite(match(probes, FAgenes))
table(asFAs)
# FALSE  TRUE
# 30008   583  only 583 FAs were included for consensus analysi

# locate FAs in modules
moduleFAs<-as.data.frame(table(moduleLabels[asFAs]))
names(moduleFAs)<-c("moduleLabels","FAs")
moduleAll <-as.data.frame(table(moduleLabels))
names(moduleAll)<-c("moduleLabels","All")
moduleFAs<-merge(moduleFAs, moduleAll, by="moduleLabels" ,all.y=TRUE)
moduleFAs$FAs[is.na(moduleFAs$FAs)]<-0
## calculate enrichment
tt<-colSums(moduleFAs[,2:3])
moduleFAs$fisherP<-round( apply(moduleFAs[,2:3],1,
function(x)fisher.test(matrix(as.numeric(c(x,tt-x)), nrow = 2, dimnames = list( c( "FAs","all"),c("inModule", "out"))) ,alternative="greater" )$p.value) ,3)
# FAs enriched in below modules
moduleFAs[moduleFAs$fisherP<0.05,]
# moduleLabels FAs  All fisherP
#            18  47 1558   0.002
#            40  20  669   0.036
#            47  99 3714   0.000
#            93   8   70   0.000
save(moduleFAs, file = "s9.consensusModuleFunctions.RData")

# list modules for each FA gene
modules<-data.frame(nodeName=probes, module=moduleLabels)
FAs<-merge(FAs, modules, by="nodeName", all.x=TRUE)
dim(FAs<-FAs[!is.na(FAs$module),])
unique(FAs$Biological.process)
#[1] "Plastid Fatty Acid Synthesis from Pyruvate"                     "β-Oxidation"
#[3] "Other Acyl Lipid Related"                                       "Eukaryotic Glycerolipid Synthesis"
#[5] "Oil Storage"                                                    "Sphingolipid Synthesis"
#[7] "Plastidial Glycerolipid, Galactolipid and Sulfolipid Synthesis" "Transcription Factors Associated with Lipid Synthesis"
#[9] "Mitochondrial Fatty Acid and Lipoic Acid Synthesis"             "TAG synthesis"
#[11] "Lipid Transfer Proteins"                                        "Miscellaneous"
#[13] "Fatty Acid Elongation"                                          "Fatty acid desaturation (out of plastid)"
#[15] "Plastid Lipid Trafficking"
xtabs(~Biological.process+module,data=FAs)
# now I want to do fishers test to show which module is each category enriched at
pdf("s9.FAs2consensusmodules.pdf",width=10,height=7)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
o_cm<-c("ME39","ME5","ME6","ME84","ME125","ME121","ME126","ME33","ME124","ME19","ME112","ME18","ME105","ME24","ME93","ME122","ME85","ME32","ME83","ME107","ME119","ME128","ME108","ME129","ME130","ME69","ME92","ME117","ME123","ME55","ME37","ME53","ME38","ME118","ME17","ME89","ME110","ME127","ME111","ME116","ME48","ME10","ME99","ME86","ME13","ME30","ME44","ME50","ME88","ME29","ME72","ME67","ME40","ME65","ME31","ME22","ME46","ME66","ME102","ME47","ME94","ME64","ME58","ME96","ME9","ME49","ME36","ME71","ME91")                             # put modules in columns
coln<-length(o_cm)
rown<-length(unique(FAs$Biological.process) )  # put FAs categories in rows
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(NA, nrow = rown, ncol = coln);
CountTbl = matrix(0, nrow = rown, ncol = coln);
# color list of MEs in the color of decreasing numbers of memebers
colModule  <- as.numeric(gsub("ME","",o_cm))
# category
rowCategory <- unique(FAs$Biological.process)
# Execute all pairwaise comparisons
for (rmod in 1:rown)
for (cmod in 1:coln)
{
    rMembers = (FAs$Biological.process == rowCategory[rmod] );
    cMembers = (FAs$module == colModule[cmod] );
    CountTbl[rmod, cmod] = sum(FAs$Biological.process == rowCategory[rmod]  & FAs$module == colModule[cmod]   )
    if(CountTbl[rmod, cmod]>0) pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
}
# display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
# make infinitely small p-value to the smallest finite values
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
rModTotals = apply(CountTbl, 1, sum)
cModTotals = apply(CountTbl, 2, sum)
select<-(cModTotals>0)
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap( Matrix = pTable[,select], colorLabels = TRUE,
xLabels = o_cm[select], yLabels = paste(" ", rowCategory),
textMatrix = CountTbl[,select], colors = blueWhiteRed(100)[48:100],
main = "Correspondence of FA-related gene categories to consensus modules ",
cex.text = 0.5, cex.lab = 0.5, setStdMargins = FALSE      )
dev.off()
write.table(FAs, file="s9.FAs2consensusmodules.txt", row.names=FALSE, sep="\t", quote=FALSE)


# Annotated consensus moduls with topGO
library(topGO)
load('D5annotation.Rdata')
# "annotation221" "annot"
load('cottonGOenrich.RData')
# all included in individual network
universe<-colnames(multiExpr[[1]]$data)
# genes of interest

GOresults<-data.frame()
for(module in unique(moduleLabels))
{
    genes<-universe[moduleLabels==module]
    geneList <- factor(as.integer(universe %in% genes))
    names(geneList) <- universe
    
    pdf(file=paste("topGO_consensus/cME",module,".pdf", sep=""))

    # topGO analysis
    if(exists("enrich")) {remove(enrich)}
    for(on in c("MF","BP","CC"))
    {
        print(on)
        # Make topGO object
        GOdata <- new("topGOdata", ontology = on, allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        # fisher test
        result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        results.table <- GenTable(GOdata, result, topNodes = length(result@score))
        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
        results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
        # label ontology type
        results.table$ontology<-on
        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
       keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
        if(exists("enrich")) enrich<- rbind(enrich, keep)
        if(!exists("enrich")) enrich<- keep
        
        # draw figure for GO terms pval<=0.05 before FDR correction
        if(is.na(sigNo<-length(keep$ontology))){next}
        showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
        mtext(on, line=-1)
    }
    dev.off()
    if(dim(enrich)[1]>0)
    {
        enrichME<-enrich
        enrichME$ME=module
        GOresults<-rbind(GOresults,enrichME)   }
}
write.table(GOresults, file="topGO_consensus/consensus_GOresults.txt", sep="\t", row.names=FALSE)

# generate mega-module for dissimiliarity =0.2
pdf("topGO_consensus/megaModules.pdf")
for(i in 1:6)
{
    A = adjacency(consMEs[[i]]$data, power = 1, type="signed");
    tre<-flashClust(as.dist(1-A),method="average")
    plot(tre, main=setLabels[i], cex=0.7 );
    abline(h=0.2,  col = "red")
    abline(h=0.25, col = "blue")
    cut<-cutree(tre,h=0.25)
    assign(paste("cut",shortLabels[i],sep=""),cut)
}
dev.off()

## MAke netwrok visualization for eigengenes
# plot eigengene network
# plotEigengeneNetworks(consMEs[[1]]$data,setLabels[1])
for(i in 1:6){
    A = adjacency(consMEs[[i]]$data, power = 1, type="signed")
    # below code does the same to plot the eigengene network dendrogram
    plot(flashClust(as.dist(1-A),method="average") )
    #abline(h=0.2, col = "red")
    # export A to cytoscape for netwrok visualization
    cyt = exportNetworkToCytoscape(A,
    edgeFile = paste("cytoscape_consensus/consMEs",shortLabels[i],"-edge.txt", sep=""),
    nodeFile = paste("cytoscape_consensus/consMEs",shortLabels[i],"-nodes.txt", sep=""),
    weighted = TRUE,
    threshold = 0,  #if higher than minValue, losing nodes
    nodeNames = colnames(consMEs[[i]]$data),
    altNodeNames = NA,
    nodeAttr = NA)
}
# read edge files
i=1
x<-read.table(file=paste("cytoscape_consensus/consMEs",shortLabels[i],"-edge.txt", sep=""),header=TRUE)
x$name<-paste(x$fromNode,x$toNode,sep=".")
x<-x[,c(7,3)]
names(x)[2]<-shortLabels[i]
edge<-x
for(i in 2:6)
{
    x<-read.table(file=paste("cytoscape_consensus/consMEs",shortLabels[i],"-edge.txt", sep=""),header=TRUE)
    x$name<-paste(x$fromNode,x$toNode,sep=".")
    x<-x[,c(7,3)]
    names(x)[2]<-shortLabels[i]
    edge<-merge(edge,x, by="name",all.x=TRUE,all.y=TRUE)
}
edge$fromNode <-gsub("[.].*","",edge$name)
edge$toNode <-gsub(".*[.]","",edge$name)
edge$direction<-"undirected"
edge<-edge[,c("fromNode","direction","toNode","A2","D5","TM1","Yuc","AD3","Syn")]
write.table(edge, "cytoscape_consensus/consMEs-edges.txt",sep="\t",quote=FALSE,row.names=FALSE)
# order MEs according to the perservation figure
o_cm # put modules in columns
o<-data.frame(nodeName=o_cm,order=1:length(o_cm))
oo<-as.data.frame(table(moduleLabels))
oo$nodeName<-paste("ME",oo$moduleLabels,sep="")
names(oo)[2]<-"size"
dim(o<-merge(o,oo[,2:3],by="nodeName") )
# incorperate mega Module info
o$megaM_A2<-cutA2[o$nodeName]
o$megaM_D5<-cutD5[o$nodeName]
o$megaM_TM1<-cutTM1[o$nodeName]
o$megaM_Syn<-cutSyn[o$nodeName]
o$megaM_AD3<-cutAD3[o$nodeName]
o$megaM_Yuc<-cutYuc[o$nodeName]
# give color
# col=brewer.pal(11,"Paired")
# m<-o[,4:9]; for(i in 1:11){m[m==i]<-col[i]}
# o[,4:9]<-m
# incorperate cor with dpa
# 12 samples
dpa12=rep(seq(10,40,by=10),each=3)
# 11 samples
dpa11=dpa12[c(1:5,7:12)]
dpaC <- data.frame( nodeName=names(consMEs[[1]]$data),
                    datC_A2 =cor(consMEs[[1]]$data,dpa12),
                    datC_D5 =cor(consMEs[[2]]$data,dpa11),
                    datC_TM1=cor(consMEs[[3]]$data,dpa12),
                    datC_Yuc=cor(consMEs[[4]]$data,dpa12),
                    datC_AD3=cor(consMEs[[5]]$data,dpa12),
                    datC_Syn=cor(consMEs[[6]]$data,dpa11)    )
dpaC[,2:7]<-numbers2colors(dpaC[,2:7])
dim(o<-merge(o,dpaC,by="nodeName") )

write.table(o, "cytoscape_consensus/consMEs-nodes.txt",sep="\t",quote=FALSE,row.names=FALSE)
##########################################
## Cytoscape
## 1. File => Import Network File => choose edge file
## 2. File => Import Table File   => choose node file
## 3. select all nodes but ME0, extract new network from selection
## 4. Layout => Attribute Circle Layout => all nodes by "order"; Layout => Bundle Edges => all nodes and edges => Number of handles =1
## 5. style: node label color BLACK, node shape Ellipse, check box "Lock node width and height", node size by continuous mapping of "size" 50 to 90,
## 6. style for each genome: transparency 20-150 for 0.9-1 (A2,D5,TM1, Yuc, AD3, Syn), Edge Stroke color (unselected) ("#F8766D" "#00BF7D" "#00B0F6" "#E76BF3" "#A3A500","619CFF"), Node color by datC
## paste in powerpoints, resize 4.44*6.1




############### Step 10. consider homoeolog ratios for consensus modules  ###############
nohup R CMD BATCH s10......R &
################

# get library factors from total counts
library(DESeq2)
total<-read.table("s1.count.raw.txt",header=TRUE,sep="\t")
total<-total[,names(total)!="D5.20_3"]
sampleT<-data.frame( sample = gsub("_.","", names(total)), genome = gsub("[.].*","", names(total)), dpa = gsub(".*[.]|_.","", names(total)), rep = gsub(".*_","", names(total)) )
dds <- DESeqDataSetFromMatrix( countData = total, colData = sampleT, design = ~ sample)
# rlog transformation, note that with defalt blind=TRUE, design is actually ~1
dds<-estimateSizeFactors(dds)
sizeFactors(dds)
totalSizeFactors<-data.frame(sampleS=names(total),sf=sizeFactors(dds) )

# Get homoeolog read counts
# load expression data
data = read.table("All ADT not normalized.txt",header=TRUE, sep="\t");
# Take a quick look at what is in the data set:
dim(data);  #37223   157
names(data);
names(data)<-gsub(".bam","",names(data))
rownames(data)<-data$gene
# remove columes of "D5.20_3.D" "D5.20_3.T"
dim(data<- data[,-grep("D5.20_3",names(data))] )
data<-data[,-1]
data<-data[,grep("[.]A|[.]D",names(data))]
sample<-data.frame(sample=names(data),sampleS=gsub("[.]A|[.]D","",names(data) ), homoeolog=gsub(".*[.]","",names(data) ) )
sample<-merge(sample,totalSizeFactors,by="sampleS")
sample$dpa<-gsub(".*[.]|_.","",sample$sampleS)
sample$genome<-gsub("[.].*","",sample$sampleS)

# pdf("s10.homoeologBinSize.pdf")
# par(mfrow=c(2,1))
# select<-which(sample$genome=="TM1" )
# boxplot(colSums(data[,select])~sample$homoeolog[select]+sample$dpa[select])
# dev.off()

## determine DE between diploids
######### 1. compare total counts from A2 and D5
select<-which( sampleT$genome %in% c("A2","D5") )
count<-total[,select]
coldata<-sampleT[select,]
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ genome + dpa)
dds <- DESeq(dds)
res <- results(dds, c("genome", "A2","D5"))
summary(res, alpha=0.05)
# log2 fold change (MAP): genome A2 vs D5
# out of 35774 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 8646, 24%
# LFC < 0 (down)   : 8760, 24%
assign("diploid_t",res)
######### 2. compare partition counts A2-A amd D5-D
select<-which( sample$genome %in% c("A2","D5") )
count<-data[,select]
coldata<-sample[select,]
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ homoeolog + dpa)
# instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
sizeFactors(dds)<-coldata$sf
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds, c("homoeolog", "A","D"))
summary(res, alpha=0.05)
# log2 fold change (MAP): homoeolog A vs D
# out of 34259 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 8314, 24%
# LFC < 0 (down)   : 8237, 24%
assign("diploid_h_tsf",res)   # based on homoeolog partition while using total   size factor
# use its own libarzy size
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ homoeolog + dpa)
dds <- DESeq(dds)
res <- results(dds, c("homoeolog", "A","D"))
summary(res, alpha=0.05)
# based on homoeolog partition while using partitioned bin size as factor
# out of 34259 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 8500, 25%
# LFC < 0 (down)   : 8035, 23%
assign("diploid_h",res)

######### let us now compare results diploid_t, diploid_h, diploid_h_tsf
At<- which(!is.na(diploid_t$padj)     & diploid_t$padj    <0.05 & diploid_t$log2FoldChange>0     )
Dt<- which(!is.na(diploid_t$padj)     & diploid_t$padj    <0.05 & diploid_t$log2FoldChange<0     )
At_h<- which(!is.na(diploid_h$padj)     & diploid_h$padj    <0.05 & diploid_h$log2FoldChange>0     )
Dt_h<- which(!is.na(diploid_h$padj)     & diploid_h$padj    <0.05 & diploid_h$log2FoldChange<0     )
At_h_tsf<- which(!is.na(diploid_h_tsf$padj)     & diploid_h_tsf$padj    <0.05 & diploid_h_tsf$log2FoldChange>0     )
Dt_h_tsf<- which(!is.na(diploid_h_tsf$padj)     & diploid_h_tsf$padj    <0.05 & diploid_h_tsf$log2FoldChange<0     )
pdf("s10.vennComparingDiploidDEs.pdf")
par(mfrow=c(2,1))
venn(list(At=At,Dt=Dt,At_h=At_h,Dt_h=Dt_h))
venn(list(At=At,Dt=Dt,At_h_tsf=At_h_tsf,Dt_h_tsf=Dt_h_tsf))
venn(list(At=At,At_h=At_h,At_h_tsf=At_h_tsf))
venn(list(Dt=Dt,Dt_h=Dt_h,Dt_h_tsf=Dt_h_tsf))
venn(list(At_h=At_h,Dt_h=Dt_h,At_h_tsf=At_h_tsf,Dt_h_tsf=Dt_h_tsf))
dev.off()
# hard to tell statistically whether diploid_h or diploid_h_tsf works better to recover TRUE diploid_t results
# diploid_t identided most DEs, due to the use of total counts
# diploid_hsf shows 24-24% as diploid_t, maybe slightly better than diploid_h
# still don't want in allopolyploids At and Dt were corrected by different library size, since they should be directly comparable
# anyway, use library size correction for homoeolog


## determine homoeolog expression bias: At!=Dt
for( genome in  c("Yuc", "AD3", "TM1") )
{
    select<-which(sample$genome==genome)
    dds <- DESeqDataSetFromMatrix( countData = data[,select], colData = sample[select,], design = ~ homoeolog+dpa)
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    sizeFactors(dds)<-sample$sf[select]
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds, c("homoeolog", "A","D"))
    print(genome)
    print(res)
    summary(res, alpha=0.05)
    assign(genome,res)
}
# Yuc: LFC > 0 (up)     : 12657, 38%;  LFC < 0 (down)   : 4908, 15%
# AD3: LFC > 0 (up)     : 13284, 39%;  LFC < 0 (down)   : 5330, 16%
# TM1: LFC > 0 (up)     : 14467, 43%;  LFC < 0 (down)   : 5940, 18%

# and repeat Ran's comparison at each stage: Ran used RPKM instead
for( dd in  c("10", "20", "30", "40") )
{
    select<-which(sample$genome=="TM1" & sample$dpa ==dd)
    dds <- DESeqDataSetFromMatrix( countData = data[,select], colData = sample[select,], design = ~ homoeolog)
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    sizeFactors(dds)<-sample$sf[select]
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    print(dd)
    print(res)
    summary(res, alpha=0.05)
}
# 10dpa: LFC > 0 (up)     : 2966, 9.5%; LFC < 0 (down)   : 7886, 25%
# 20dpa: LFC > 0 (up)     : 2983, 9.4%; LFC < 0 (down)   : 8053, 25%
# 30dpa: LFC > 0 (up)     : 3015, 9.7%; LFC < 0 (down)   : 8342, 27%
# 40dpa: LFC > 0 (up)     : 2938, 9.3%; LFC < 0 (down)   : 7838, 25%

# and let At and Dt each corrected by the partitioned library size
for( dd in  c("10", "20", "30", "40") )
{
    select<-which(sample$genome=="TM1" & sample$dpa ==dd)
    dds <- DESeqDataSetFromMatrix( countData = data[,select], colData = sample[select,], design = ~ homoeolog)
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    # sizeFactors(dds)<-sample$sf[select]
    dds <- DESeq(dds)
    res <- results(dds)
    print(dd)
    print(res)
    summary(res, alpha=0.05)
}
# 10dpa: LFC > 0 (up)     : 5196, 17%; LFC < 0 (down)   : 5502, 18%
# 20dpa: LFC > 0 (up)     : 5248, 16%; LFC < 0 (down)   : 5488, 17%
# 30dpa: LFC > 0 (up)     : 4569, 15%; LFC < 0 (down)   : 4604, 15%
# 40dpa: LFC > 0 (up)     : 4775, 15%; LFC < 0 (down)   : 4735, 15%

#### so overall, the At bin size is bigger than Dt bin size???
save(diploid_t, diploid_h, diploid_h_tsf, Yuc, TM1, AD3, file = "R-11-homoeologBias.RData")
# work with results: TM1, Yuc, AD3, get modules
load('R-10-homoeologBias.RData')
load('R-07-prep.RData')
load("s8.moduleConsensus.A2_D5_TM1_Yuc_AD3_Syn.RData")

# order MEs according to the perservation figure
o<-data.frame(nodeName=c("ME39","ME5","ME6","ME84","ME125","ME121","ME126","ME33","ME124","ME19","ME112","ME18","ME105","ME24","ME93","ME122","ME85","ME32","ME83","ME107","ME119","ME128","ME108","ME129","ME130","ME69","ME92","ME117","ME123","ME55","ME37","ME53","ME38","ME118","ME17","ME89","ME110","ME127","ME111","ME116","ME48","ME10","ME99","ME86","ME13","ME30","ME44","ME50","ME88","ME29","ME72","ME67","ME40","ME65","ME31","ME22","ME46","ME66","ME102","ME47","ME94","ME64","ME58","ME96","ME9","ME49","ME36","ME71","ME91"),order=1:69)
o$module<-gsub("ME","",o$nodeName)


for(i in c("AD3", "Yuc","TM1") )
{bias<-get(i)
    names(moduleLabels)<-colnames(multiExpr[[1]]$data)
    bias$module<-moduleLabels[rownames(bias)]
    cutoff<-log2(1.5)
    all<-as.data.frame(table(bias$module))
    At<-as.data.frame(table(bias$module[bias$padj<0.05 & !is.na(bias$padj) & bias$log2FoldChange>cutoff]))
    Dt<-as.data.frame(table(bias$module[bias$padj<0.05 & !is.na(bias$padj) & bias$log2FoldChange< -cutoff]))
    t<-merge(merge(all,At,by="Var1",all.x=TRUE,all.y=TRUE),Dt,by="Var1",all.x=TRUE,all.y=TRUE)
    names(t)<-c("module","all","At","Dt")
    t[is.na(t)]<-0
    t$towards<-ifelse(t$At==t$Dt ,"-",ifelse(t$At>t$Dt,"At","Dt"))
    t$pval<-apply(t[,3:4],1,function(x) {ee<-(x[1]+x[1])/2; tail<-ifelse(x[1]>x[2],"greater","less"); fisher.test(matrix(c(x[1],x[2],ee,ee), nrow = 2, dimnames = list( c( "Observed","Expected"),c("At", "Dt"))), alternative=tail )$p.value })
    table(t$padj<-p.adjust(t$pval,method="bonferroni")<0.05 )
    t<-merge(t,o, by="module")
    print(i)
    # print(t)
    print( xtabs(padj~towards, data=t) )
}
# basically all unbalanced towards At bias
#  "AD3"
# towards
# At Dt
# 27  0
# "Yuc"
# towards
# At Dt
3 31  0
# "TM1"
#towards
# At Dt
# 28  0

######################################################
######################################################
# try to get new read counts from 4.1 new mapping
setwd("~/jfw-lab/Projects/Eflen/eflen_seed")
grep("AD-",grep("noeflen",list.files(), value=TRUE) , value=TRUE))->files
seed<-matrix();for(i in files){temp<-read.table(i,header=TRUE,sep="\t"); seed<-cbind(seed,temp) }

grep("[.]A[.]|[.]D[.]", grep("AD-",grep("noeflen",list.files(), value=TRUE) , value=TRUE) , value=TRUE)->files

x<-as.data.frame( seed[,grep("bam",names(seed))] )
x<-apply(x,2,as.numeric)
dim( y<-aggregate(x,by=list(seed$gene), FUN=sum) )
rownames(y)<-y$Group.1
y<-y[,-1]
z<-y[,grep("sort.bam",names(y))]  #z as total
dds <- DESeqDataSetFromMatrix( countData = z, colData = data.frame(ss=gsub("sort.bam","",names(z))), design = ~ ss)
dds<-estimateSizeFactors(dds)
sf<-data.frame(ss=gsub(".sort.bam","",names(z)),sf= sizeFactors(dds))
# deal with homoeologs counts
y<-y[,grep("[.]A[.]|[.]D[.]",names(y))]
sample<-data.frame(sample=gsub("[.]sort|[.]bam","",names(y) ), dpa=gsub("AD.|.R.*","",names(y) ) )
sample$homoeolog <-gsub(".*[.]","",sample$sample)
sample$ss <-gsub("[.]A|[.]D","",sample$sample)
sample<-merge(sample,sf,by="ss")
#by dpa
for( dd in  c("10", "20", "30", "40") )
{
    select<-which( sample$dpa ==dd)
    dds <- DESeqDataSetFromMatrix( countData = y[,select], colData = sample[select,], design = ~ homoeolog)
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    # sizeFactors(dds)<-sample$sf[select]
    dds <- DESeq(dds)
    res <- results(dds)
    print(dd)
    print(res)
    summary(res, alpha=0.05)
}
# LFC > 0 (up)     : 5470, 17% ; LFC < 0 (down)   : 5574, 17%
# LFC > 0 (up)     : 2836, 8.4%; LFC < 0 (down)   : 3000, 8.9%
# LFC > 0 (up)     : 5429, 17%; LFC < 0 (down)   : 6009, 18%
# LFC > 0 (up)     : 2401, 7.3%; LFC < 0 (down)   : 2512, 7.6%


#### with correction
for( dd in  c("10", "20", "30", "40") )
{
    select<-which( sample$dpa ==dd)
    dds <- DESeqDataSetFromMatrix( countData = y[,select], colData = sample[select,], design = ~ homoeolog)
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    sizeFactors(dds)<-sample$sf[select]
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    print(dd)
    print(res)
    summary(res, alpha=0.05)
}
# LFC > 0 (up)     : 3675, 11%;  LFC < 0 (down)   : 5641, 17%
# LFC > 0 (up)     : 1106, 3.3%; LFC < 0 (down)   : 4552, 14%
# LFC > 0 (up)     : 4250, 13%;  LFC < 0 (down)   : 6115, 19%
# LFC > 0 (up)     : 2173, 6.6%; LFC < 0 (down)   : 2959, 9%


#### overall, with sf correction
dds <- DESeqDataSetFromMatrix( countData = y, colData = sample, design = ~ homoeolog+dpa)
# instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
sizeFactors(dds)<-sample$sf
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds, c("homoeolog", "A","D"))
res
summary(res, alpha=0.05)
# LFC > 0 (up)     : 9654, 27%, LFC < 0 (down)   : 4960, 14%
## So with more homoeolog reads mapped and partitioned, the magnitude of bias droped from 60% to 40%, while the direction towards A remained unchanged








############### Step 11. export individual networks for cytoscape  ###############
nohup R CMD BATCH s10......R &
################
library(WGCNA);
library(flashClust);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE);

remove(list=ls())

load("R-07-prep.RData") # "multiExpr"          "lanova"             "powerTables.signed"
load("R-07-individualNetworks.RData") # "A2net12"  "AD3net12" "D5net12"  "Synnet12" "TM1net12" "Yucnet12"
checkSets(multiExpr)
nSets<-checkSets(multiExpr)$nSets   #6
nGenes<-checkSets(multiExpr)$nGenes #30591
nSamples<-checkSets(multiExpr)$nSamples #12
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development", "Syn seed devleopment")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3", "Syn")


# Extract only FA genes
FAs<-read.table("s9.FAs2consensusmodules.txt",header=TRUE,sep="\t")
dim(FAs) #584
dim(FAs<-FAs[FAs$module!=0,]) #454
FAs$moduleC<-labels2colors(as.numeric(factor(FAs$module)))
library(gplots)
FAs$moduleC_hex<-col2hex(FAs$moduleC)

probes = colnames(multiExpr[[2]]$data)   #30591
# locate FA genes
asFAs = is.finite(match(probes, FAs$nodeName))
table(asFAs)
# FALSE  TRUE
# 30158   433

# and i also want to check homoeolog bias for those genes
load('R-11-homoeologBias.RData') #"diploid_t"     "diploid_h"     "diploid_h_tsf" "Yuc"           "TM1"           "AD3"
library(DESeq2)


# Export FA network to cytoscape
# connect to speedy, or whatever used to build network and store the TOM file (ISOLONE, biocrunch, etc.)
# Cytoscape [2] allows the user to input an edge file and a node file, allowing the user to specify for example the link weights and the node colors.
# Get topological overlap, "TOM"
for(i in 1:nSets)
{
    load(paste(shortLabels[i],"_power12_TOM-block.1.RData", sep=""))
    # Select the corresponding Topological Overlap
    subTOM = as.matrix(TOM)[asFAs, asFAs];
    # add bias data to FAs
    if(i %in% 3:5){
        bias<-get(shortLabels[i])
        dim(bias<-bias[FAs$nodeName,])  #454
        bias$log2FoldChange[bias$padj<0.05]<-0
        # table(FAs$nodeName==rownames(bias) )
        FAs$bias<-bias$log2FoldChange
    }
    # str(subTOM)
    remove(TOM)
    subProbes = probes[asFAs];
    dimnames(subTOM) = list(subProbes, subProbes)
    # Export the network into edge and node list files Cytoscape can read
    cyt = exportNetworkToCytoscape(subTOM,
    edgeFile = paste(shortLabels[i],"-edges-FAs.txt",sep=""),
    nodeFile = paste(shortLabels[i],"-nodes-FAs.txt",sep=""), weighted = TRUE, threshold = quantile(subTOM, 0.8),
    nodeNames = subProbes,
    nodeAttr = FAs[,-1])
}
##########################################
## Cytoscape
## 1. File => Import Network File => choose edge file
## 2. File => Import Table File   => choose node file
## 3. select all nodes but ME0, extract new network from selection
## 4. Customize Layout and Style as saved images in folder "cytoscape_individual/Setting/"
## 5. Tools => NetworkAnalyzer => Network Analysis. Note that this analysis considers network as binary network, not weighted network, so the results are for refernce only.





# I want to describe density and clustering coeffiency for TM1 and Yuc subnetworks
# processes<-c("Plastid Fatty Acid Synthesis from Pyruvate", "Other Acyl Lipid Related", "Oil Storage", "Plastidial Glycerolipid, Galactolipid and Sulfolipid Synthesis", "Mitochondrial Fatty Acid and Lipoic Acid Synthesis", "TAG synthesis", "Sphingolipid Synthesis", "Lipid Transfer Proteins", "Miscellaneous", "Fatty Acid Elongation", "Transcription Factors Associated with Lipid Synthesis", "β-Oxidation", "Fatty acid desaturation (out of plastid)", "Plastid Lipid Trafficking")
ttt<-as.data.frame(table(FAs$Biological.process) )
names(ttt)<-c("Biological.process","Gene.number")
processes<-ttt$Biological.process

str(subDat<-multiExpr[[3]]$data)  #TM1
adjacency = adjacency(subDat, power = 12, type = "signed")
table(colnames(adjacency)==probes)
ajFAs<-adjacency[asFAs,asFAs]
mean(ajFAs)                        # density 0.086531
mean(ccFAs<-clusterCoef(ajFAs)  )  # 0.325477
for(i in 1:length(processes))
{
    # must make sure the order match
    # table(rownames(ajFAs)==FAs$nodeName)
    as.p<-(FAs$Biological.process==processes[i] )
    p_adjacency<-ajFAs[as.p,as.p]
    # colnames(p_adjacency)==FAs$nodeName[FAs$Biological.process==processes[i]]
    ttt$TM1.density[i] = mean(p_adjacency)
    ttt$TM1.mean.clusterCoef[i] = mean( clusterCoef(p_adjacency) )
}

str(subDat<-multiExpr[[4]]$data)  #Yuc
adjacency = adjacency(subDat, power = 12, type = "signed")
table(colnames(adjacency)==probes)
ajFAs<-adjacency[asFAs,asFAs]
mean(ajFAs)         # density   0.06101973
mean(ccFAs<-clusterCoef(ajFAs)  )  # 0.2580747
for(i in 1:length(processes))
{
    # must make sure the order match
    # table(rownames(ajFAs)==FAs$nodeName)
    as.p<-(FAs$Biological.process==processes[i] )
    p_adjacency<-ajFAs[as.p,as.p]
    # colnames(p_adjacency)==FAs$nodeName[FAs$Biological.process==processes[i]]
    ttt$Yuc.density[i] = mean(p_adjacency)
    ttt$Yuc.mean.clusterCoef[i] = mean( clusterCoef(p_adjacency) )
}
write.table(ttt,"s11.FAs.network.parameters.txt", sep="\t")

#

#########################################################
#########################################################
#####################  2015/12  #########################
#####   If I need to re-run this analysis again    ######
#####  after the 2014, 2015 and current version    ######
#####  I am fXXking quitting science all together  ######
#########################################################
#########################################################
#########################################################
