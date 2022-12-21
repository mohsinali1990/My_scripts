# Author: Dr.Mohsin Ali

# install.packages("pcadapt", dependencies = T)


library(pcadapt)
library(qvalue)

# pcadapt has been developed to detect genetic markers involved in biological adaptation. pcadapt provides statistical tools for outlier detection based on Principal Component Analysis (PCA).

# pcadapt package can perform genome scans for selection based on individual genotype data

myVCF <- read.pcadapt(
  "Repeat_Het0.2_494GBS_TaxaOrdered_Maf0.05_Miss0_QCout_WP.recode.vcf",
  type = "vcf",
  type.out = c("bed", "matrix"),
  allele.sep = c("/", "|")
)

# why pcadapt? => pcadapt is robust to admixture and does not assume prior knowledge of population structure.

# pcadapt performs proceeds in two steps: 1) pca is performed on centered and scaled genotypic data;
                                        # 2) second step involved computation of test statistics and p-vales on the correltions between SNPs and the first "k" PCs.
# note: an outlier test pcadapt (Luu et al., 2017) was based on allele frequency differentiation

# Luu, K., Bazin, E., & Blum, M. G. (2017). pcadapt: an R package to perform genome scans for selection based on principal component analysis. Molecular ecology resources, 17(1), 67-77.

#PCA only

myVCF_PCA_only <- pcadapt(
  myVCF,
  K = 10, # number of PCs
  method = "mahalanobis", #mahalanobis (default): the robust Mahalanobis distance is computed for each genetic marker using a robust estimate of both mean and covariance matrix between the K vectors of z-scores.
  min.maf = 0.05,
  ploidy = 2,
  LD.clumping = NULL,
  pca.only = FALSE,
  tol = 1e-04
)

PCA_Var_out<- myVCF_PCA_only$singular.values
print(PCA_Var_out)*100

barplot(PCA_Var_out)


# data visualization 
plot(myVCF_PCA_only, option = "screeplot")
plot(myVCF_PCA_only , option = "manhattan")
plot(myVCF_PCA_only , option = "qqplot")
plot(myVCF_PCA_only , option = "stat.distribution")



Popinfo <- read.delim("POPINFO.txt")

#2D pca plot
plot(
  myVCF_PCA_only,
  option = "scores",
  pch = 16,
  pop = Popinfo$Group..Full.number.,
  xlab = "PCA1",
  ylab = "PCA2"
)

#pairwise PCs plot with grouping
pairs(
  myVCF_PCA_only$scores,
  col = c(Popinfo$Color),
  pch=16,
  labels = c(
    "PCA1",
    "PCA2",
    "PCA3",
    "PCA4",
    "PCA5",
    "PCA6",
    "PCA7",
    "PCA8",
    "PCA9",
    "PCA10"
    )
)


plot(myVCF_PCA_only, option = "screeplot")
plot(myVCF_PCA_only, option = "scores")
plot(myVCF_PCA_only , option = "manhattan")


qval <- qvalue(myVCF_PCA_only$pvalues)$qvalues
alpha <- 0.1
outliers_pcadapt <- which(qval < alpha)
print(outliers_pcadapt)

length(outliers_pcadapt) # 1428 outliers


alpha <- 0.05 # use of a more stringent threshold to detect outliers
outliers <- which(qval < alpha)
print(outliers)
length(outliers_pcadapt) # 1428 outliers

