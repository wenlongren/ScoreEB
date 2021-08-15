library(BEDMatrix)
library(data.table)

dir1 <- "C:/Users/ThinkPad/Desktop/Lucky/ScoreEB/Frontiers_in_Genetics_Resubmit/Data/Sample_Size_Data/mrMLM_Format/2000/"

file_sim_bed <- paste0(dir1, "Sim_Sample_2000.bed")
file_sim_bim <- paste0(dir1, "Sim_Sample_2000.bim")
file_sim_fam <- paste0(dir1, "Sim_Sample_2000.fam")

sim_bed <- BEDMatrix(file_sim_bed, simple_names = TRUE)
sim_bim <- fread(file_sim_bim)
sim_fam <- fread(file_sim_fam)

len_sample <- dim(sim_fam)[1]
len_marker <- dim(sim_bim)[1]

sim_bed <- as.matrix(sim_bed)
for(i in 1:len_sample){
  sim_bed[i,which(is.na(sim_bed[i,])==TRUE)] <- 0
}

geno_for_code1 <- sim_bim[,5]
mrMLM_geno <- cbind(sim_bim[,c(2,1,4)], geno_for_code1, t(sim_bed[,]))

colnames(mrMLM_geno) <- c("rs#", "chrom", "pos", "genotype for code 1", as.matrix(sim_fam[,1]))

### compute kinship ###
kin <- sim_bed%*%t(sim_bed)
mean_diag <- mean(diag(kin))
kin <- kin/mean_diag

rownames(kin) <- c(as.matrix(sim_fam[,1]))

write.table(mrMLM_geno, paste0(dir1, "mrMLM_Genotype_2000.csv"), sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(kin, paste0(dir1, "mrMLM_Kinship_2000.csv"), sep = ",", quote = FALSE, col.names = FALSE, row.names = TRUE)

pheno <- sim_fam[,-c(2:5)]
colnames(new_pheno) <- c("<Phenotype>", paste("Trait",c(1:(dim(pheno)[2]-1))))

write.table(pheno,paste0(dir1, "mrMLM_Phenotype_2000.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)







