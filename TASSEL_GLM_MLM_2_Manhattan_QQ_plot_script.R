library(qqman)
library(dplyr)

# import TASSEL results
# note 
TASSEL_MLM_Out <- read.table("Tassel_glm_out.txt", header = T, sep = "\t")

# Number of traits
head(unique(TASSEL_MLM_Out$Trait))

# note: for each plot trait name must be specificed  

# first trait as example (i.e., EarHT)
Trait1 <-  TASSEL_MLM_Out %>% filter(.$Trait == "EarHT")

# Bonferroni correction threshold
nmrk <- nrow(Trait1)
(GWAS_Bonn_corr_threshold <- -log10(0.05 / nmrk))

# Manhattan plot
(Mann_plot <- manhattan(
  TASSEL_MLM_Out,
  chr = "Chr",
  bp = "Pos",
  snp = "Marker",
  p = "p",
  col = c("red", "blue"),
  annotateTop = T,
  genomewideline = GWAS_Bonn_corr_threshold,
  suggestiveline = F
)
)

# QQ plot
QQ_plot <- qq(TASSEL_MLM_Out$p)

# Manhattan and Q-Q plot arranged in 1 rows and 2 columns 
old_par <- par()

par(mfrow=c(1,2))
(Mann_plot <- manhattan(
  TASSEL_MLM_Out,
  chr = "Chr",
  bp = "Pos",
  snp = "Marker",
  p = "p",
  col = c("red", "blue"),
  annotateTop = T,
  genomewideline = GWAS_Bonn_corr_threshold,
  suggestiveline = F,
  main = "EarHT" # trait name 
)
)

(QQ_plot <- qq(TASSEL_MLM_Out$p,  main = "EarHT" ))
