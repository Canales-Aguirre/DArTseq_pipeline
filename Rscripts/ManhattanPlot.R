#install packages | only necessary once
install.packages("BiocManager")
BiocManager::install("karyoploteR")
BiocManager::install("regioneR")

#load karyoploteR and regioneR
library(regioneR)
library(karyoploteR)

#load vcf, construct genome
set.seed(123456)
vcf <- read.table("VCF.vcf")
snps <- toGRanges(vcf[,c(1,2,2)])
snps$pval <- rnorm(n=NROW(snps), mean=0.5, sd=1)
snps$pval[snps$pval<0] <- -1*snps$pval[snps$pval<0]
snps$pval <- 10^(-1*snps$pval)
snps.ranges <- toGRanges(snps)
custom.genome <- toGRanges("ChromStructure.txt")

#construct plot
kp <- plotKaryotype(custom.genome, plot.type=4)
transf.pval <- -log10(snps.ranges$pval)
points.col <- colByValue(transf.pval, colors=c("#BBBBBB00", "grey"))
kp <- kpPlotManhattan(kp, data=snps.ranges, highlight="GENE.bed", highlight.col="red", points.col=points.col)
ymax <- kp$latest.plot$computed.values$ymax #y-axis 
ticks <- c(0, seq_len(floor(ymax)))
kpAxis(kp, ymin=0, ymax=ymax, tick.pos = ticks)
