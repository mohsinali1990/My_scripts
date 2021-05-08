# Author: Dr.Mohsin Ali; Email: mohsinali@caas.cn 

#install rBLAST using devtools

devtools::install_github("mhahsler/rBLAST")

library(ggplot2)
library(rBLAST)


my_query = readDNAStringSet("90K_test.fasta")

# setting up BLAST
RefGenome = '161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta'

#build database


makeblastdb(RefGenome, dbtype= "nucl")

# run BLAST query
blast_run <- blast(db=RefGenome, type="blastn")
rBLAST_Results <- predict(blast_run, my_query)
# str(rBLAST_Results)
# summary(rBLAST_Results)


#save blast results
write.csv(rBLAST_Results, "Rblast_results.csv")

max_spos <- max(rBLAST_Results$S.start)

ggplot(rBLAST_Results) +
  geom_density(aes(x=S.start), kernel= "rectangular", n=3000) +
  xlim(0, max_spos) + 
  theme_linedraw()
