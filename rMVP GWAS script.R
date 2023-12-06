################################################################################
#                                 rMVP GWAS 
################################################################################

# devtools::install_github("xiaolei-lab/rMVP")
install.packages("rMVP")
library(rMVP)

################################################################################
# Step1: formatting of data 
################################################################################

rm(list = ls()) #clear R environment

MVP.Data(fileHMP="mdp_genotype.hmp.txt", #Genotypic data in HapMap format
         filePhe="mdp_traits.txt",
         sep.hmp="\t",
         sep.phe="\t",
         SNP.effect="Add",
         fileKin=T,
         out="mvp.hmp"
)

################################################################################
# Step2: import formatted data (MVP.Data) from working directory 
################################################################################

genotypic_dat <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype_dat <- read.table("mvp.hmp.phe",head=TRUE)
map_info <- read.table("mvp.hmp.geno.map" , head = TRUE)
popstr <- read.table("mdp_population_structure.txt", header = T, skip = 1)
Kinship <- MVP.K.VanRaden(genotypic_dat, verbose = T)

################################################################################
##Step2: GWAS run 
################################################################################

################################################################################
# PCA as covariate
################################################################################
GWAS_PCA_mvp <- MVP(
  phe = phenotype_dat, 
  geno = genotypic_dat, 
  map = map_info,
  method =  c("GLM", "MLM", "FarmCPU"),
  nPC.GLM = 5,
  nPC.MLM = 5, 
  nPC.FarmCPU = 5
)

################################################################################
# Kinship as covariate
################################################################################

dir.create("PC_K")
setwd("./PC_K/")
getwd()

GWAS_Kin_mvp <- MVP(
  phe = phenotype_dat, 
  geno = genotypic_dat, 
  map = map_info,
  method =  c("GLM", "MLM", "FarmCPU"),
  nPC.GLM = 0,
  nPC.MLM = 0, 
  nPC.FarmCPU = 0,
  K = Kinship
)

################################################################################
# PCA and Kinship as covariates
################################################################################

dir.create("../PcaKinship" )
setwd("../PcaKinship")

GWAS_PCA_Kin_mvp <- MVP(
  phe = phenotype_dat, 
  geno = genotypic_dat, 
  map = map_info,
  method =  c("GLM", "MLM", "FarmCPU"),
  nPC.GLM = 5,
  nPC.MLM = 5, 
  nPC.FarmCPU = 5,
  K = Kinship
)
