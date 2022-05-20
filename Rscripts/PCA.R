# create distance matrix with plink1.9
# plink --vcf VCF --double-id --allow-extra-chr --distance-matrix --out VCF_DM

# load distance matrix, pop IDs, and ind IDs
dist_populations <- read.table("VCF_DM.mdist",header=F)
pop <- data.frame(pop=read.table("VCF_DM.mdist.id")[,1])
popInd <- data.frame(popInd=read.table("VCF_DM.mdist.id")[,2])

# generate eigenvectors and PCA
mds_populations <- cmdscale(dist_populations,eig=T,5)
eigenvec_populations <- cbind(pop,popInd,mds_populations$points)
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)

# packages
install.packages("tidyverse", dependencies = TRUE) # only once
install.packages("ggforce", dependencies = TRUE) # only once
library(tidyverse)
library(ggforce)

# PCA plot
ggplot(data = eigenvec_populations, aes(x=`1`, y=`2`, label=popInd)) +
     geom_point(mapping=aes(x=`1`, y=`2`, color=pop), show.legend=FALSE) +
     geom_text(hjust="inward", vjust="inward", size=1.5) +
     geom_mark_ellipse(aes(color=pop, fill=pop), show.legend=FALSE) +
     geom_hline(yintercept=0, linetype="dotted") + 
     geom_vline(xintercept=0, linetype="dotted") +
     labs(title = "PCA of Vietnam_China",
          x = paste0("PC1 (",eigen_percent[1]," %)"),
          y = paste0("PC2 (",eigen_percent[2]," %)")) + 
     theme_minimal()
