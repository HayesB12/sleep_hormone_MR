# Chronotype + Total Testosterone > BrCa
# Chronotype + Bioavailable Testosterone > BrCa

# Update script throughout depending on data being used #

#install.packages("devtools")
library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
library(MendelianRandomization)
#use Wes's MVMR package
library(devtools)
#install_github("WSpiller/MVMR")
library(MVMR)
library(tidyr)
library(dplyr)
library(RadialMR)
#install.packages("meta")
library(meta)
#install.packages("hablar")
library(hablar)
library(data.table)

##########################
### FEMALE
##########################

exposure_dat <- read_exposure_data("chronotype_female.txt", sep = "\t",
                                 snp_col = "SNP", beta_col = "beta",
                               se_col = "se", effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval")

exposure_chrono_F <- clump_data(exposure_dat)

exposure_chrono_F$samplesize.exposure <- 244207

###

exposure_dat <- read_exposure_data("totaltest_female_inv_fin.txt", sep = "\t",
                                 snp_col = "SNP", beta_col = "BETA",
                                se_col = "SE", effect_allele_col = "ALLELE1",
                               other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM_INF")

exposure_totaltest_F <- clump_data(exposure_dat)

exposure_totaltest_F$samplesize.exposure <- 199596

###

exposure_dat <- read_exposure_data("biotest_female_inv_fin.txt", sep = "\t",
                                  snp_col = "SNP", beta_col = "BETA",
                                 se_col = "SE", effect_allele_col = "ALLELE1",
                                other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM_INF")

exposure_biotest_F <- clump_data(exposure_dat)

exposure_biotest_F$samplesize.exposure <- 180386


################

chr_tot_snplist <- as.data.frame(c(exposure_chrono_F$SNP, exposure_totaltest_F$SNP))

chr_bio_snplist <- as.data.frame(c(exposure_chrono_F$SNP, exposure_biotest_F$SNP))

colnames (chr_tot_snplist)[1] <- "SNP"
colnames (chr_bio_snplist)[1] <- "SNP"

#####################################
# CHRONOTYPE AND TOTAL TESTOSTERONE #
#####################################

# READ IN EXPOSURE DATA #

# Chronotype & Total Testosterone SNPs extracted from chronotype GWAS #

exposure_dat1 <- read_exposure_data("mvmr_chrono_and_totaltest_females.txt", sep = "\t",snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")
exposure_dat1$exposure <- "chronotype"

exposure_dat1 <- exposure_dat1[exposure_dat1$SNP %in% chr_tot_snplist$SNP,]
dim(exposure_dat1)

# Chronotype & Total Testosterone SNPs extracted from Total Testosterone GWAS #

exposure_dat2 <- read_exposure_data("mvmr_totaltest_and_chrono_females.txt", sep = "\t",snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")
exposure_dat2$exposure <- "total_testosterone" 

exposure_dat2 <- exposure_dat2[exposure_dat2$SNP %in% chr_tot_snplist$SNP,]
dim(exposure_dat2)

exposure_dat1$samplesize.exposure <- 244207
exposure_dat2$samplesize.exposure <- 199596

# combine these datasets & clump data #

exposure_dat <- rbind(exposure_dat1, exposure_dat2)
exposure_dat <- clump_data(exposure_dat)

# READ IN OUTCOME DATA #

outcome_dat <- read_outcome_data("outcome_mvMR_chrono_totaltest_BrCa_overall.txt", sep = "\t",
                                 snp_col = "SNP", beta_col = "Beta.meta",
                                 se_col = "sdE.meta", effect_allele_col = "Effect.Meta",
                                 other_allele_col = "Baseline.Meta", eaf_col = "EAF", pval_col = "p.meta")

outcome_dat$samplesize.outcome <- 247173

# HARMONISE EXPOSURE AND OUTCOME DATA #

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

# PERFORM MVMR ANALYSIS #

res_mv <- mv_multiple(mvdat)

or_mv_chr_tot <- generate_odds_ratios(res_mv$result)

#write.csv(or_mv_chr_tot, "or_mvmr_chrono_tottest_BrCa.csv", row.names=F, quote=F)


# SENSITIVITY ANALYSES #

# Reformat data to allow assessment of instrument strength & pleiotropy #

XGs_betas <- mvdat$exposure_beta
XGs_se <- mvdat$exposure_se

YG_betas <- mvdat$outcome_beta 
YG_se <- mvdat$outcome_se

mvmr <- format_mvmr(XGs_betas[,c(1:2)], YG_betas, XGs_se[,c(1:2)], YG_se, row.names(XGs_betas)) 
mvmr_res <- mvmr(mvmr, 0, 1)

# MR-PRESSO for mvMR #
# mr_presso(BetaOutcome = "Y_effect", BetaExposure = c("E1_effect", "E2_effect"), SdOutcome = "Y_se", SdExposure = c("E1_se", "E2_se"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)

mvmr_pleiotropy <- as.data.frame(pleiotropy_mvmr(mvmr))
mvmr_strength <- as.data.frame(strength_mvmr(mvmr))

mvmr_results_chr_tot_F <- cbind.data.frame(mvmr_pleiotropy, mvmr_strength)

colnames (mvmr_results_chr_tot_F)[3] <- "F-stat Chronotype" 
colnames (mvmr_results_chr_tot_F)[4] <- "F-stat Testosterone"
rownames (mvmr_results_chr_tot_F)[1] <- "Chrono + Total Testosterone (F)"

################################
########################################
################################
########################################
################################

# READ IN EXPOSURE DATA #

# Chronotype & Bioavailable Testosterone SNPs extracted from chronotype GWAS #

exposure_dat1 <- read_exposure_data("mvmr_chrono_and_biotest_females.txt", sep = "\t",snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")
exposure_dat1$exposure <- "chronotype"

exposure_dat1 <- exposure_dat1[exposure_dat1$SNP %in% chr_bio_snplist$SNP,]
dim(exposure_dat1)

# Chronotype & Total Testosterone SNPs extracted from Bioavailable Testosterone GWAS #

exposure_dat2 <- read_exposure_data("mvmr_biotest_and_chrono_females.txt", sep = "\t",snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")
exposure_dat2$exposure <- "bioavailable_testosterone" 

exposure_dat2 <- exposure_dat2[exposure_dat2$SNP %in% chr_bio_snplist$SNP,]
dim(exposure_dat2)

exposure_dat1$samplesize.exposure <- 244207
exposure_dat2$samplesize.exposure <- 180386

# combine these datasets & clump data #

exposure_dat <- rbind(exposure_dat1, exposure_dat2)
exposure_dat <- clump_data(exposure_dat)

# READ IN OUTCOME DATA #

outcome_dat <- read_outcome_data("outcome_mvMR_chrono_biotest_BrCa_overall.txt", sep = "\t",
                                 snp_col = "SNP", beta_col = "Beta.meta",
                                 se_col = "sdE.meta", effect_allele_col = "Effect.Meta",
                                 other_allele_col = "Baseline.Meta", eaf_col = "EAF", pval_col = "p.meta")

outcome_dat$samplesize.outcome <- 247173

# HARMONISE EXPOSURE AND OUTCOME DATA #

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

# PERFORM MVMR ANALYSIS #

res_mv <- mv_multiple(mvdat)

or_mv_chr_bio <- generate_odds_ratios(res_mv$result)

#write.csv(or_mv_chr_bio, "or_mvmr_chrono_biotest_BrCa.csv", row.names=F, quote=F)


# SENSITIVITY ANALYSES #

# Reformat data to allow assessment of instrument strength & pleiotropy #

XGs_betas <- mvdat$exposure_beta
XGs_se <- mvdat$exposure_se

YG_betas <- mvdat$outcome_beta 
YG_se <- mvdat$outcome_se

mvmr <- format_mvmr(XGs_betas[,c(1:2)], YG_betas, XGs_se[,c(1:2)], YG_se, row.names(XGs_betas)) 
mvmr_res <- mvmr(mvmr, 0, 1)

# MR-PRESSO for mvMR #
# mr_presso(BetaOutcome = "Y_effect", BetaExposure = c("E1_effect", "E2_effect"), SdOutcome = "Y_se", SdExposure = c("E1_se", "E2_se"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)

mvmr_pleiotropy <- as.data.frame(pleiotropy_mvmr(mvmr))
mvmr_strength <- as.data.frame(strength_mvmr(mvmr))

mvmr_results_chr_bio_F <- cbind.data.frame(mvmr_pleiotropy, mvmr_strength)

colnames (mvmr_results_chr_bio_F)[3] <- "F-stat Chronotype" 
colnames (mvmr_results_chr_bio_F)[4] <- "F-stat Testosterone"
rownames (mvmr_results_chr_bio_F)[1] <- "Chrono + Bioavailable Testosterone (F)"

# Save all sensitivity outputs #

mvmr_sensitivity <- rbind.data.frame(mvmr_results_chr_tot_F, mvmr_results_chr_bio_F)
#write.csv(mvmr_sensitivity, "mvMR_sensitivity_BrCa.csv", row.names=T, quote=F)

# Save all MVMR outputs # 

mvmr_results <- rbind.data.frame(or_mv_chr_tot, or_mv_chr_bio)
#write.csv(mvmr_results, "mvMR_results_BrCa.csv", row.names=T, quote=F)
