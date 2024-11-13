
#PRS analysis: Goal - predict the risk for high BMI and T2D for 1000 Genomes individuals (using only chr 18 data) 

#Step 1-2: Check out the SNPs in the GWAS data and the SNPs in the 1000 Genomes dataset (we are only interested in the intersection of these variants) 
#Step 3: Select which subset of SNPs are actually going into the risk score (PRS). Via pruning and thresholding. 
#Step 4: Identify which subset of SNPs survived pruning and thresholding - and extract these from 1000 Genomes files. 
#Step 5: Use --score to calculate the polygenic risk score 
#Step 6: Compare scores between BMI and T2D. 

#Step 1: In R, harmonize the 1KG dataset and the GWAS summary statistics (select the intersection of SNPs)
> library(data.table)
> bmi <- fread("PASS_BMI1.sumstats", header = T)
> t2d <- fread("PASS_Type_2_Diabetes.sumstats", header = T)
> eur_1kg <- fread("1000G_eur_exclude_dups_chr18.bim", header = F)

> keep_snps_bmi <- which(bmi$SNP %in% eur_1kg$V2)
> keep_snps_t2d <- which(t2d$SNP %in% eur_1kg$V2)

> bmi <- bmi[keep_snps_bmi, ]
> t2d <- t2d[keep_snps_t2d, ]

#Tiffany made fake beta values, because these GWAS summary statistics don't have them. 
#In real life, you need the real beta values from the GWAS. :) 
#beta = z*se, and se is a function of sample size. 
bmi$BETA <- bmi$Z*(1/bmi$N)*1e4
t2d$BETA <- t2d$Z*(1/t2d$N)*1e4

bmi$SE <- (1/bmi$N)*1e4
t2d$SE <- (1/t2d$N)*1e4

bmi$P <- pnorm(q = abs(bmi$Z), lower.tail = F)
t2d$P <- pnorm(q = abs(t2d$Z), lower.tail = F)

> write.table(bmi, file = "PASS_BMI1.chr18.sumstats", row.names = F, col.names = T, sep = "\t", quote = F)
> write.table(t2d, file = "PASS_Type_2_Diabetes.chr18.sumstats", row.names = F, col.names = T, sep = "\t", quote = F)

#Step 2. Extract the kept SNPs from the 1000G_eur_exclude_dups_chr18 plink files. 

awk 'NR!=1{print $1}' PASS_BMI1.chr18.sumstats > SNPs_BMI.txt #prints column 1 but removes column header 
awk 'NR!=1{print $1}' PASS_Type_2_Diabetes.chr18.sumstats > SNPs_T2D.txt

plink \
        --bfile 1000G_eur_exclude_dups_chr18 \
        --extract SNPs_BMI.txt \
        --make-bed \
        --out 1000G_eur_exclude_dups_chr18_BMI

plink \
        --bfile 1000G_eur_exclude_dups_chr18 \
        --extract SNPs_T2D.txt \
        --make-bed \
        --out 1000G_eur_exclude_dups_chr18_T2D


#Step 3. Select SNPs that we want to keep for the PRS: perform pruning and thresholding. 

trait=BMI
gwas=PASS_BMI1.chr18.sumstats
plink --bfile 1000G_eur_exclude_dups_chr18_${trait} --clump-p1 0.0001 --clump-r2 0.1 --clump-kb 250 --clump ${gwas} --clump-snp-field SNP --clump-field P --out 1000G_eur_exclude_dups_chr18_clumped_${trait}

trait=T2D
gwas=PASS_Type_2_Diabetes.chr18.sumstats
plink --bfile 1000G_eur_exclude_dups_chr18_${trait} --clump-p1 0.0001 --clump-r2 0.1 --clump-kb 250 --clump ${gwas} --clump-snp-field SNP --clump-field P --out 1000G_eur_exclude_dups_chr18_clumped_${trait}

#OR (for easier viewing of parameters, but my unix shell did not handle multiple arguments well)

plink \ 
	--bfile 1000G_eur_exclude_dups_chr18_${trait} \
	--clump-p1 0.0001 \
	--clump-r2 0.1 \
	--clump-kb 250 \
	--clump ${gwas} \
	--clump-snp-field SNP \
	--clump-field P \
	--out 1000G_eur_exclude_dups_chr18_clumped_${trait}

	
#notes: clump-p1: all snps included regardless of p-value. 
#clump-r2: snps with r2 of 0.1 or greater will be removed (pruned) from index variant 
#clump-kb snps within 250kb of the index variant will be pruned

#Step 4. Extract the SNPs that survived pruning. 

for trait in T2D BMI
do
awk 'NR!=1{print $3}' 1000G_eur_exclude_dups_chr18_clumped_${trait}.clumped > PRS.SNPs.${trait}
done

for trait in T2D BMI
do
plink --bfile 1000G_eur_exclude_dups_chr18_${trait} --extract PRS.SNPs.${trait} --make-bed --out 1000G_eur_exclude_dups_chr18_PRS_${trait}
done

#Step 5. Prepare Score File in R (or awk if you are an awk wizard!)

> for (trait in c("T2D","BMI")){
  y <- fread(paste0("PRS.SNPs.",trait), header = F)$V1
  w <- which(nchar(y) > 0)
  y <- y[w]

  m <- match(y, bmi$SNP)
  score <- cbind(bmi$SNP[m], bmi$A1[m], bmi$BETA[m]) #want to make sure BETA corresponds to A2 (not A1)
  write.table(score, file = paste0(trait, "_score_file.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
  }

for trait in T2D BMI
do
plink --bfile 1000G_eur_exclude_dups_chr18_PRS_${trait} --out PRS_chr18_${trait} --score ${trait}_score_file.txt 1 2 3 
done

#Step 6. Let's see what the distribution of these scores look like and if they are correlated. 

> prs_t2d <- fread("PRS_chr18_T2D.profile", header = T)
> prs_bmi <- fread("PRS_chr18_BMI.profile", header = T)

> hist(prs_t2d$SCORE)
> hist(prs_bmi$SCORE)

> plot(prs_t2d$SCORE, prs_bmi$SCORE, xlab = "T2D", ylab = "BMI", main = "Correlation of PRS for Two Related Phenotypes")
> abline(lm(prs_bmi$SCORE ~ prs_t2d$SCORE))
> cor.test(prs_t2d$SCORE, prs_bmi$SCORE) #cor = 0.53, p-val < 2.2e-16
> cor.test(prs_t2d$SCORE, prs_bmi$SCORE)$p.val #exact p-value 1.59e-37

