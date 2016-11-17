################# 10/18/16 Differential Coexpression analysis

library(WGCNA);
library(flashClust);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE);
remove(list=ls())
load("R-07-prep.RData") # "multiExpr"          "lanova"             "powerTables.signed"
load("R-07-individualNetworks.RData") # "A2net12"  "AD3net12" "D5net12"  "Synnet12" "TM1net12" "Yucnet12"

for(i in 1:6)
{
    gg <- unique(gsub("[.].*","",rownames(multiExpr[[i]]$data) ) )
    print(gg)
    subDat    <-  multiExpr[[i]]$data
    subDat   <-  apply(subDat,2,as.numeric)
    adj = adjacency(subDat, power = 12, type = "signed")
    assign(paste(gg,"adj",sep=""),adj)
}

# percentage of rewired edges R/C
# C = n(n-1)/2
# R = sum(abs(Xadj-Yadj))/2

n = checkSets(multiExpr)$nGenes # 30591
C = n*(n-1)/2
rate = sum(abs(A2adj-D5adj))/2/C

sum(abs(A2adj-D5adj))/2/C
# 0.05304749
sum(abs(TM1adj-Yucadj))/2/C
# 0.03911227
sum(abs(AD3adj-Yucadj))/2/C
# 0.03560141
sum(abs(AD3adj-TM1adj))/2/C
# 0.04546805
sum(abs(A2adj-ADsadj))/2/C
# 0.0363874
sum(abs(D5adj-ADsadj))/2/C
# 0.03647
sum(abs(Yucadj-ADsadj))/2/C
# 0.04002498
sum(abs(AD3adj-ADsadj))/2/C
# 0.04709474
sum(abs(TM1adj-ADsadj))/2/C
# 0.04859758


####################################
## Differential correlation analysis
####################################
library(DiffCorr)
#  comp.2.cc.fdr(output.file="res.txt", t(multiExpr[[1]]$data), t(multiExpr[[2]]$data), threshold=0.05)

checkSets(multiExpr)
nSets<-checkSets(multiExpr)$nSets   #6
nGenes<-checkSets(multiExpr)$nGenes #30591
nSamples<-checkSets(multiExpr)$nSamples# #12 11 12 12 12 11
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development", "Syn seed devleopment")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3", "Syn")

softPower = 12

# Analysis of differntial coexpression gene pairs
pwset<-combn(nSets,2)
for(i in 1:choose(nSets,2))
{
    X<-t(multiExpr[[ pwset[1,i] ]]$data)
    Y<-t(multiExpr[[ pwset[2,i] ]]$data)
    
    outfile <- paste0(shortLabels[ pwset[1,i] ],"vs",shortLabels[ pwset[2,i] ],".res.txt" )
    
    comp.2.cc.fdr(output.file=outfile, X, Y, threshold=0.05)
    
}

# move all files into the "DC" folder
# check rewired pairs in consensus module
load("s8.moduleConsensus.A2_D5_TM1_Yuc_AD3_Syn.RData") ->lnames
lnames
# [1] "consMEs"               "unmergedLabels"        "moduleColors"
# [4] "moduleLabels"          "moduleLabels.adjusted" "consTree"
probes = colnames(multiExpr[[2]]$data)   #30591
names(moduleLabels) <- probes
mSizes<-as.data.frame(table(moduleLabels))




for(i in 1:choose(nSets,2))
{
    print(paste0("Comparing ", shortLabels[ pwset[1,i] ], " vs ", shortLabels[ pwset[2,i] ]))
    
    infile <- paste0("DC/",shortLabels[ pwset[1,i] ],"vs",shortLabels[ pwset[2,i] ],".res.txt" )
    x<-read.table(infile,header=TRUE, sep="\t")
    p= 2*nrow(x)/nGenes/(nGenes-1)
    
    x<-x[,1:2]
    x$X<-moduleLabels[x$molecule.X]
    x$Y<-moduleLabels[x$molecule.Y]
    # get rewired connnections considering consensus modules
    a <- as.matrix(xtabs(~X+Y, data=x))
    
    # sum up both directions
    new <- a
    new[lower.tri(new)]<-NA
    new[upper.tri(new)] <- a[upper.tri(a)] + t(a)[upper.tri(a)]
    # make result table
    n<- as.data.frame(new)
    n<-n[!is.na(n$Freq),]
    n$type <- ifelse(n$X==n$Y,"within","between")
    
    n$X.size <- mSizes$Freq[match(n$X, mSizes$moduleLabels)]
    n$Y.size <- mSizes$Freq[match(n$Y, mSizes$moduleLabels)]
    n$edge.total <- n$X.size*n$Y.size
    n$edge.total[n$type=="within"] <- (n$X.size[n$type=="within"] -1) * n$X.size[n$type=="within"] /2
    n$rewired.percentage<-n$Freq/n$edge.total
    
    #
    n$p = p
    # Probability of observed number >= n$Freq based on the overall p value, if <0.05, significantly high rewiring
    n$binom.pvalue <- pbinom(n$Freq-1, size=n$edge.total, prob=p, lower.tail=FALSE)
    n$qvalue <- p.adjust(n$binom.pvalue, "BH")
    
    outfile <- paste0("DC/",shortLabels[ pwset[1,i] ],"vs",shortLabels[ pwset[2,i] ],".CM.txt" )
    write.table(n, outfile, row.names=FALSE, sep="\t" )



}
# Are those edges being rewired between A2 and D5 continue to be changed later in polyploids, or new pairs involved? No, only a small fractions kept being rewired. I cannot interpret this.
# how are rewired edges located in consensus modules?
# a lot of rewired edges are connecting cME0, which makes sense though since cME0 genes are not in conserved modules.


###########################################
# Identify differential coexpression genes
###########################################
for(i in 1:choose(nSets,2))
{
    print(paste0("Comparing ", shortLabels[ pwset[1,i] ], " vs ", shortLabels[ pwset[2,i] ]))
    
    infile <- paste0(shortLabels[ pwset[1,i] ],"vs",shortLabels[ pwset[2,i] ],".res.txt" )
    x<-read.table(infile,header=TRUE, sep="\t")
    
    ## percentage of rewired pairs
    # C=n*(n-1)/2
    # p=nrow(x)/C
    p= 2*nrow(x)/nGenes/(nGenes-1)
    print(paste0("Differential coexpression pairs: ", nrow(x)," (", p,")" ) )
    
    ## number of rewired edges for each node
    ks <- as.data.frame(table(c(x$molecule.X, x$molecule.Y)) )
    ks$P <- apply(ks, 1, function(x) pbinom(as.numeric(x[2]), size=n, prob=p, lower.tail=FALSE) ) 
    ks$q.bh<-p.adjust(ks$P, "BH")
    ks.sig<-ks[ks$P<0.05 & ks$q.bh<0.05,]
    n.dcg <- nrow(ks.sig)
    p.dcg <- n.dcg/nGenes
    print(paste0("Differential coexpression geness: ", n.dcg," (", p.dcg,")" ) )

    outfile <- paste0(shortLabels[ pwset[1,i] ],"vs",shortLabels[ pwset[2,i] ],".gene.dcg.txt" )
    write.table(ks.sig, outfile, row.names=FALSE, sep="\t" )
}

# "Comparing A2 vs D5"
# "Differential coexpression pairs: 15708780 (0.0335737074756426)"
# "Differential coexpression geness: 10849 (0.354646791539995)"
# "Comparing A2 vs TM1"
# "Differential coexpression pairs: 3676001 (0.00785656061477527)"
# "Differential coexpression geness: 6804 (0.222418358340688)"
# "Comparing A2 vs Yuc"
# "Differential coexpression pairs: 1982313 (0.00423671327672572)"
# "Differential coexpression geness: 4479 (0.14641561243503)"
# "Comparing A2 vs AD3"
# "Differential coexpression pairs: 3636533 (0.00777220733675822)"
# "Differential coexpression geness: 6409 (0.209506063874996)"
# "Comparing A2 vs Syn"
# "Differential coexpression pairs: 568595 (0.00121523391390757)"
# "Differential coexpression geness: 1883 (0.061554051845314)"
# "Comparing D5 vs TM1"
# "Differential coexpression pairs: 19240674 (0.0411222743274908)"
# "Differential coexpression geness: 10600 (0.346507142623647)"
# "Comparing D5 vs Yuc"
# "Differential coexpression pairs: 4843495 (0.010351795893108)"
# "Differential coexpression geness: 6654 (0.217514955379033)"
# "Comparing D5 vs AD3"
# "Differential coexpression pairs: 11697128 (0.0249997742521792)"
# "Differential coexpression geness: 9222 (0.301461214082573)"
# "Comparing D5 vs Syn"
# "Differential coexpression pairs: 891335 (0.00190501239133796)"
# "Differential coexpression geness: 1998 (0.0653133274492498)"
# "Comparing TM1 vs Yuc"
# "Differential coexpression pairs: 1124591 (0.00240354052088961)"
# "Differential coexpression geness: 2738 (0.0895034487267497)"
# "Comparing TM1 vs AD3"
# "Differential coexpression pairs: 2744883 (0.00586652170931569)"
# "Differential coexpression geness: 4875 (0.1593605962538)"
# "Comparing TM1 vs Syn"
# "Differential coexpression pairs: 5174383 (0.0110589887444434)"
# "Differential coexpression geness: 7658 (0.250335065869046)"
# "Comparing Yuc vs AD3"
# "Differential coexpression pairs: 799532 (0.00170880574337507)"
# "Differential coexpression geness: 1999 (0.0653460168023275)"
# "Comparing Yuc vs Syn"
# "Differential coexpression pairs: 962655 (0.00205744159444366)"
# "Differential coexpression geness: 2267 (0.0741067634271518)"
# "Comparing AD3 vs Syn"
# "Differential coexpression pairs: 3286858 (0.00702486182924298)"
# "Differential coexpression geness: 5849 (0.191200026151482)"


# move all files into the "DC" folder
files<-grep("dcg.txt",list.files("DC"),value=TRUE)
DCid<-list()
for(f in files)
{
    x<-read.table(paste0("DC/",f), header=TRUE, row.names=1,sep="\t")
    DCid[[ gsub(".gene.dcg.txt","",f) ]]<-rownames(x)
}
sapply(DCid,length)
save(DCid, file="DCid.Rdata")

# compare DE with DC
load("DEid.Rdata")
tt<-as.data.frame( sapply(DEid_con,length) )
names(tt)<-"DE_No"
tt$DE_perc <- tt$DE_No/37223

for(i in names(DCid))
{
    tt[i,"DC_No"]<-length(DCid[[i]])
    tt[i,"overlap"]<-length(intersect(DEid_con[[i]], DCid[[i]]))
}
tt$overlap_perc<-tt$overlap/tt$DC_No


# Annotated DCGs with topGO
library(topGO)
load('D5annotation.Rdata')
# "annotation221" "annot"
load('cottonGOenrich.RData')
# "geneID2GO"       "go4profile"  "goraiUniverse"   "gsc"             "level1_BP_terms" "level1_MF_terms"  "level2_BP_terms" "level2_MF_terms" "level3_BP_terms" "level3_MF_terms"   "level4_BP_terms" "level4_MF_terms" "level5_BP_terms" "level5_MF_terms"
# all included in individual network
universe<-goraiUniverse
# for each module containing genes of interest
GOresults<-data.frame()
for(i in names(DCid))
{
    genes<-DCid[[i]]
    geneList <- factor(as.integer(universe %in% genes))
    names(geneList) <- universe
    
    pdf(file=paste("DC/topgo_",i,".pdf", sep=""))
    
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
        enrichCompr<-enrich
        enrichCompr$compr=i
        GOresults<-rbind(GOresults,enrichCompr)   }
}
write.table(GOresults, file="DC/GOresults.txt", sep="\t", row.names=FALSE)

# functional annotation of dcg and dc-gene-pairs
# are those being rewired between A2 and D5 continue to be changed later in polyploids, or new pairs involved
# how are these located in consensus, the simplified network?


# Annotated DCGs with topGO
library(topGO)
load('D5annotation.Rdata')
# "annotation221" "annot"
load('cottonGOenrich.RData')
# "geneID2GO"       "go4profile"  "goraiUniverse"   "gsc"             "level1_BP_terms" "level1_MF_terms"  "level2_BP_terms" "level2_MF_terms" "level3_BP_terms" "level3_MF_terms"   "level4_BP_terms" "level4_MF_terms" "level5_BP_terms" "level5_MF_terms"
# all included in individual network
universe<-goraiUniverse
# for each module containing genes of interest
GOresults<-data.frame()
# Annotation DCGs that were not found in DE
for(i in names(DCid))
{
    print(i)
    genes<-setdiff( DCid[[i]], DEid_con[[i]])
    geneList <- factor(as.integer(universe %in% genes))
    names(geneList) <- universe
    
    pdf(file=paste("DC/topgo_",i,".onlyDC.pdf", sep=""))
    
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
        enrichCompr<-enrich
        enrichCompr$compr=i
        GOresults<-rbind(GOresults,enrichCompr)   }
}
write.table(GOresults, file="DC/GOresults_onlyDC.txt", sep="\t", row.names=FALSE)

# functional annotation of dcg and dc-gene-pairs doesn't seem to generate interpratable results, because without partintioning into modular structure, the lists are too complicated.

######################################################################
## Check rewired edges for oil relaated network between TM1 and Yuc ##
######################################################################

# Extract only FA genes
FAs<-read.table("s9.FAs2consensusmodules.txt",header=TRUE,sep="\t")
dim(FAs) #583
dim(FAs<-FAs[FAs$module!=0,]) #433

# get DC results for TM1 vs Yuc
i=10
print(paste0("Comparing ", shortLabels[ pwset[1,i] ], " vs ", shortLabels[ pwset[2,i] ]))
    
infile <- paste0("DC/",shortLabels[ pwset[1,i] ],"vs",shortLabels[ pwset[2,i] ],".res.txt" )
dc<-read.table(infile,header=TRUE, sep="\t")
nrow(dc)  # 1124591 rewired edges
p= 2*nrow(dc)/nGenes/(nGenes-1)
p   # percentage of rewired edges 0.002403541


dim( dcFAs <- dc[dc$molecule.X %in% FAs$nodeName & dc$molecule.Y %in% FAs$nodeName, ] )  # 73 rewired edges among FA related genes
length(unique(c(dcFAs$molecule.X,dcFAs$molecule.Y))) # 83 genes with rewired edges



FAs$moduleC[FAs$module!=0]<-labels2colors(as.numeric(factor(FAs$module[FAs$module!=0])))
FAs$moduleC[FAs$module==0]<-"grey"
library(gplots)
FAs$moduleC_hex<-col2hex(FAs$moduleC)
write.table(FAs,file="DC/dcFAs.nodes.txt",row.names=FALSE, sep="\t")

# differential coexpression genes in FAs
dcg<-read.table("DC/TM1vsYuc.gene.dcg.txt",head=TRUE,sep="\t")


