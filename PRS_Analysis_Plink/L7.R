#Step 1 (do not redo). In R, harmonize the 1KG dataset and the GWAS summary statistics (select the intersection of SNPs)
library(data.table)
bmi <- fread("PASS_BMI1.sumstats", header = T)
t2d <- fread("PASS_Type_2_Diabetes.sumstats", header = T)
eur_1kg <- fread("1000G_eur_exclude_dups_chr18.bim", header = F)

keep_snps_bmi <- which(bmi$SNP %in% eur_1kg$V2)
keep_snps_t2d <- which(t2d$SNP %in% eur_1kg$V2)

bmi <- bmi[keep_snps_bmi, ]
t2d <- t2d[keep_snps_t2d, ]

#Tiffany made fake beta values, because these GWAS summary statistics don't have them. 
#In real life, you need the real beta values from the GWAS. :) 
#beta = z*se, and se is a function of sample size. 
bmi$BETA <- bmi$Z*(1/bmi$N)*1e4
t2d$BETA <- t2d$Z*(1/t2d$N)*1e4

bmi$SE <- (1/bmi$N)*1e4
t2d$SE <- (1/t2d$N)*1e4

bmi$P <- pnorm(q = abs(bmi$Z), lower.tail = F)
t2d$P <- pnorm(q = abs(t2d$Z), lower.tail = F)

write.table(bmi, file = "PASS_BMI1.chr18.sumstats", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(t2d, file = "PASS_Type_2_Diabetes.chr18.sumstats", row.names = F, col.names = T, sep = "\t", quote = F)


bmi <- fread("D:/UCSD/DSC291_stats_genomics/HW1_Files/HW1_Files/PASS_BMI1.chr18.sumstats", header = TRUE)
t2d <- fread("D:/UCSD/DSC291_stats_genomics/HW1_Files/HW1_Files/PASS_Type_2_Diabetes.chr18.sumstats", header = TRUE)

### #Step 5. Prepare Score File in R (or awk if you are an awk wizard!)
for (trait in c("T2D","BMI")){
  y <- fread(paste0("D:/UCSD/DSC291_stats_genomics/HW1_Files/HW1_Files/PRS.SNPs.",trait), header = F)$V1
  w <- which(nchar(y) > 0)
  y <- y[w]
  
  m <- match(y, bmi$SNP)
  score <- cbind(bmi$SNP[m], bmi$A1[m], bmi$BETA[m]) #want to make sure BETA corresponds to A2 (not A1)
  write.table(score, file = paste0("D:/UCSD/DSC291_stats_genomics/HW1_Files/HW1_Files/", trait, "_score_file.txt"), 
              row.names = F, col.names = F, sep = "\t", quote = F)
}

