#Package installation
install.packages("BiocManager") #only once
BiocManager::install() #update BiocManager			
BiocManager::install(c("SNPRelate")) #only once
BiocManager::install(c("gdsfmt")) #only once
install.packages("Rcpp") #only once   		
install.packages("ape") #only once	

#Load libraries
rm(list=ls()) #remove the objects in memory
library(Rcpp)			
library(gdsfmt)		
library(SNPRelate)		
library(tidyverse)

#Load the work directory and generate GDS
setwd("~/Users/Documents/Phylogenetics/")
snpgdsVCF2GDS("VCF.vcf", "VCF.gds", ignore.chr.prefix="Bchr")
genofile <- snpgdsOpen("VCF.gds")
snpgdsSummary("VCF.gds")
set.seed(1000) 	#Following https://benbowlab.github.io/Phylogeny.html

#Generate IBS, matrix, and tree
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))
dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL, autosome.only=FALSE, remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)	snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)	
cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL, col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL, label.H=FALSE, label.Z=TRUE, verbose=TRUE)		
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram,main="Tree")

#Export tree as NEWICK
library(ape)	
my_tree<-as.phylo(ibs.hc$hclust)
write.tree(phy=my_tree, file="VCF.newick")
snpgdsClose(genofile)	#Close the file#
