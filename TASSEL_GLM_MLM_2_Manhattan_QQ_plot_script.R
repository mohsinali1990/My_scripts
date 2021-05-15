library(qqman)
library(dplyr)

# import TASSEL results
# note 
GLM_results <- read.table("Tassel_glm_out.txt", header = T, sep = "\t")

#traits name
head(unique(GLM_results$Trait))

# note: trait name must be specific

# first trait as example (i.e., EarHT)
Trait1 <-  GLM_results %>% filter(.$Trait == "EarHT")

# Bonferroni correction threshold
nmrk <- nrow(Trait1)
(GWAS_Bonn_corr_threshold <- -log10(0.05 / nmrk))
 
# Manhattan plot
(Mann_plot <- manhattan(
  GLM_results,
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
QQ_plot <- qq(GLM_results$p)

# Manhattan and Q-Q plot arranged in 1 rows and 2 columns 
old_par <- par()

par(mfrow=c(1,2))
(Mann_plot <- manhattan(
  GLM_results,
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

(QQ_plot <- qq(GLM_results$p,  main = "EarHT" ))
