# RTM-GWAS result plot

# Required Package
library(qqman)

# import RTM-GWAS results *.pval

rtm_gwas <- read.delim(choose.files(), header = T)

# imported GWAS results
str(rtm_gwas)

# Manhattan plot visualization

manhattan(
  x = rtm_gwas,
  chr = "Chromosome",
  bp = "Position",
  p = "Y",
  snp = "Locus"
)

# Customization  of plot

manhattan(
  x = rtm_gwas,
  chr = "Chromosome",
  bp = "Position",
  p = "Y",
  snp = "Locus", 
  suggestiveline = F, 
  genomewideline = F
)


# Q-Q plot

qq(rtm_gwas$Y)