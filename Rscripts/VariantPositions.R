#install packages | only necessary once
install.packages("BiocManager")
BiocManager::install("karyoploteR")

#load karyoploteR
library(karyoploteR)

#read vcf, construct genome
vars <- read.table("VCF.vcf")
custom.genome <- toGRanges("ChromStructure.txt")
vars <- toGRanges(vars[,c(1,2,2,3:length(vars))])

#plot variants as points
kp <- plotKaryotype(genome=custom.genome, plot.type=6)
kpPoints(kp, data=vars, y=0.5)

#plot variants as density
kp <- plotKaryotype(genome=custom.genome, plot.type=1)
kpPlotDensity(kp, data=vars)

#plotting genomic regions
regions.reads <- toGRanges("GENE.bed")
kp <- plotKaryotype(genome=custom.genome)
kpPlotRegions(kp, data=regions.reads)

#plotting gene of interest-related variants
regions.gene <- toGRanges("GENE.bed")
kp <- plotKaryotype(genome=custom.genome, chromosomes="chr*")
kpPlotRegions(kp, data=regions.gene)

#plot the per base coverage (of gene of interest)
regions.coverage <- toGRanges("GENE.bed")
kp <- plotKaryotype(genome=custom.genome, chromosomes="chr*")
kpPlotCoverage(kp, data=regions.coverage, border="blue", col="orchid")
