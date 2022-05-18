################################################################################
#                Genotypic data quality control and PCA analysis
################################################################################

# In this tutorial, we will use PLINK to pre-process (i.e. quality control) of rice SNP data. 
# then we will perform PCA analysis, and then we will make PCA plot using ggplot2 

####################
# Lets get started 
###################

# For this tutorial we will use genotypic data of 36 rice accessions. Genotypic data type is VCF format.
# Plink options
system("plink.exe --help") # to get more info on plink options

################################################################################
#  For detailed information about options visit link below
#                           https://www.cog-genomics.org/plink/1.9/data
################################################################################

# lets start with VCF and recode into plink format for down stream analysis

system("plink.exe --vcf rice36indi.vcf --recode --nonfounders --allow-no-sex --out rice")

list.files("./")

# lets make .bed formated plink file
system("plink.exe --file rice --make-bed --nonfounders --allow-no-sex --out rice_bed")

# QC options 

# 1. Missingness per SNP: --geno  
# 2. Missingness per individual: --mind 
# 3. Minor allele frequency: --maf

system("plink.exe --bfile rice_bed --geno 0.2 --mind 0.2 --maf 0.05 --make-bed --out rice_qc_out")

# lets make vcf file
system("plink.exe --bfile rice_qc_out --recode vcf-iid --out rice_qc_out_vcf")

# Now, we have a fully filtered VCF we can start some analyses with it.
# First of all we will investigate population genetic structure using principal component analysis (PCA).
#
# PCA analysis: a model-free method -------------------------------------------------------------

# one of the assumptions of PCA analysis is that SNP data we use is independent (i.e., there are no spurious correlation among measured variables)
# this is not true for most of SNP dataset as allele frequencies are correlated due to physical linkage and linkage disequilibrium .

# step 1 identify prune sites

system(
  "plink.exe --vcf rice_qc_out_vcf.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --recode vcf-iid --indep-pairwise 50 10 0.2 --out rice_PRUNED"
)


# step 2: PCA analysis
system(
  "plink.exe --vcf rice_qc_out_vcf.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract rice.prune.in --pca --make-bed --out rice_prune_pca"
)


# pca plot
library(tidyverse)

plinkPCA <- read_table2("rice_prune_pca.eigenvec", col_names = F)
plinkPCA <- plinkPCA[,c(-1,-2)] # remove first two columns
EigenValue <- scan("rice_prune_pca.eigenval")
view(EigenValue)

# set columns names
names(plinkPCA)[1:ncol(plinkPCA)] <- paste0("PC", 1:(ncol(plinkPCA)))


# percentage variance explained
pve <- data.frame(PC = 1:20, pve = EigenValue/sum(EigenValue)*100)
view(pve)

# PCA plot

# make plot
P <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
P + ylab("Percentage variance explained") + theme_light()


# plot pca
 ggplot(plinkPCA, aes(PC1, PC2)) + geom_point(size = 3)+ coord_equal() +
   theme_light()+
   coord_equal() +
   theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# lets divide population in two group
# (pop <- rep(c("A","B"), each=18))
#  pop <- as.data.frame(pop)
#  plinkPCA$pop <- pop$pop
 mypop = read_tsv("pop.txt", col_names = F)
 
 plinkPCA$pop <- mypop$X1
 
 ggplot(plinkPCA, aes(PC1, PC2, color = as.factor(pop))) + geom_point(size = 3) + coord_equal() +
   theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
   labs(color= "Group")
 