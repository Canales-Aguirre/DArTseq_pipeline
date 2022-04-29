#install vcfR | only necessary once
install.packages('vcfR')

#load vcfR
library(vcfR)

#load vcf, reference
vcf <- read.vcfR("VCF.vcf.gz", verbose = TRUE )
dna <- ape::read.dna("REF.fasta", format = "fasta")

#create chromR object
chrom <- create.chromR(name="ORGANISM_NAME", vcf=vcf, seq=dna, verbose=TRUE)
chrom <- masker(chrom, min_DP = 10, min_MQ = 20)
chrom <- proc.chromR(chrom, verbose = TRUE)
head(chrom)
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp) #number the variants (not required)
head(dp)
heatmap.bp(dp) #all variants

#generate PNG image of heatmap 
png(filename="VCF_SNPs.png", height=5000, width=5000, res=300)
heatmap.bp(dp[1001:1500,]) #only variant 1001 to 1500
dev.off()
