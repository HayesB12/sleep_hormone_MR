#install.packages("devtools")
library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
#devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
#devtools::install_github("MendelianRandomization")
library(MendelianRandomization)
#devtools::install_github("MRPRESSO")
library(MRPRESSO)

##########################
##########################
### FEMALE
##########################
##########################

# Chronotype <> Total Testosterone
# Chronotype <> Bioavailable Testosterone

# Update script throughout depending on data being used #


# CHRONOTYPE > TOTAL TESTOSTERONE #

####
# Read in exposure instrument & outcome data 
####

exposure_dat <- read_exposure_data("chronotype_female.txt", sep = "\t",snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval")

outcome_dat <- read_outcome_data("outcome_chrono_totaltest_F.txt", sep = "\t", snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")

####
# Harmonise data & run MR
####

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dat <- clump_data(dat)

res_test <- mr(dat)
or_chr_tottest_F <- generate_odds_ratios(res_test)

#write.csv(or_chr_tottest_F, "or_chrono_totaltest_female_bdMR.csv", row.names=F, quote=F)

# SENSITIVITY TESTS (RAW) #

chr_tot_het_raw <- mr_heterogeneity(dat)
chr_tot_het_raw$I2 <- ((chr_tot_het_raw$Q-chr_tot_het_raw$Q_df)/chr_tot_het_raw$Q)*100
chr_tot_het_raw$outcome="Total Testosterone" # Update for different datasets

chr_tot_pleio_raw <- mr_pleiotropy_test(dat)

dat$samplesize.exposure <- 244207 # update samplesize depending on dataset used
dat$samplesize.outcome <- 199956  # update samplesize depending on dataset used

dat_steiger <- steiger_filtering(dat)
total_r2 <- sum(dat_steiger$rsq.exposure) 
F_stat <- (244207-209-1)/209 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

chr_tot_het_raw$r2 <- total_r2
chr_tot_het_raw$F_stat <- F_stat

total_r2_out <- sum(dat$rsq.outcome) 
steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)
steiger_sens$sensitivity_ratio
table(dat$steiger_dir)

#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  SignifThreshold = 0.05)

# Determine SNPs removed by steiger filtering #

dat_removed <- subset(dat, steiger_dir==FALSE)
SNPs <- as.data.frame(dat_removed$SNP)

# Harmonised data with steiger filtered SNPs removed #

dat_filtered <- subset(dat, steiger_dir==TRUE)

# SENSITIVITY TESTS (STEIGER-FILTERED) #  

chr_tot_het_filt <- mr_heterogeneity(dat_filtered)
chr_tot_pleio_filtered <- mr_pleiotropy_test(dat_filtered)

total_r2 <- sum(dat_filtered$rsq.exposure) 
F_stat <- (244207-204-1)/204 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

total_r2_out <- sum(dat_filtered$rsq.outcome) 
steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)
steiger_sens$sensitivity_ratio

chr_tot_het_filt$total_r2<- total_r2
chr_tot_het_filt$F_stat<- F_stat

res_test <- mr(dat_filtered)
or_chr_tottest_F_filt <- generate_odds_ratios(res_test)

#write.csv(chr_tot_het_raw, "sensitivity_chrono_totaltest_F_raw_supplementary.csv", row.names = FALSE, quote = FALSE)
write.csv(chr_tot_het_filt, "sensitivity_chrono_totaltest_F_filt_supplementary.csv", row.names = FALSE, quote = FALSE)

#write.csv(or_chr_tottest_F, "or_chrono_totaltest_F_raw_bdMR.csv", row.names = FALSE, quote = FALSE)
write.csv(or_chr_tottest_F, "or_chrono_totaltest_F_filt_bdMR.csv", row.names = FALSE, quote = FALSE)

# RADIAL MR #

res_single1 <- mr_singlesnp(dat_filtered)
dat1R <- dat_filtered[dat_filtered$SNP%in%res_single1$SNP,]

raddat1 <- format_radial(dat1R$beta.exposure, dat1R$beta.outcome, dat1R$se.exposure, dat1R$se.outcome, dat1R$SNP)
ivwrad1 <- ivw_radial(raddat1, alpha=0.05/204, weights=3)
dim(ivwrad1$outliers)[1] 

eggrad1 <- egger_radial(raddat1, alpha=0.05/204, weights=3)
eggrad1$coef 
dim(eggrad1$outliers)[1] 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
#plot_radial(c(ivwrad1,eggrad1), TRUE, FALSE, TRUE)

#ivwrad1$qstatistic 
#ivwrad1$sortoutliers <- ivwrad1$outliers[order(ivwrad1$outliers$p.value),]
#ivwrad1$sortoutliers$Qsum <- cumsum(ivwrad1$sortoutliers$Q_statistic)
#ivwrad1$sortoutliers$Qdif <- ivwrad1$sortoutliers$Qsum - ivwrad1$qstatistic
#write.csv(ivwrad$sortoutliers, "chrono_outliers_001.csv", row.names=F, quote=F)

#eggrad1$qstatistic 
#eggrad1$sortoutliers <- eggrad1$outliers[order(eggrad1$outliers$p.value),]
#eggrad1$sortoutliers$Qsum <- cumsum(eggrad1$sortoutliers$Q_statistic)
#eggrad1$sortoutliers$Qdif <- eggrad1$sortoutliers$Qsum - eggrad1$qstatistic
#write.csv(eggrad$sortoutliers, "chrono_outliers_egger_01.csv", row.names=F, quote=F)

# REMOVE TOP OUTLIERS #

dat1_outliers <- dat1R[!dat1R$SNP %in% ivwrad1$outliers$SNP,]
mr_dat1_outliers <- mr(dat1_outliers)
or_dat1_outliers <- generate_odds_ratios(mr_dat1_outliers)

#write.csv(dat1_outliers, "chrono_totaltest_F_filtered_outliers_removed.csv", row.names=F, quote=F)
write.csv(or_dat1_outliers, "or_chrono_totaltest_F_filtered_outliers_removed.csv", row.names=F, quote=F)

# # SENSITIVITY TESTS (OUTLIER REMOVED) # #

chr_tot_het_outliers <- mr_heterogeneity(dat1_outliers)
chr_tot_het_outliers$I2 <- ((chr_tot_het_outliers$Q-chr_tot_het_outliers$Q_df)/chr_tot_het_outliers$Q)*100
chr_tot_het_outliers$outcome="Total Testosterone" # Update for different datasets

chr_tot_pleio_outliers <- mr_pleiotropy_test(dat1_outliers)

total_r2 <- sum(dat1_outliers$rsq.exposure) 
F_stat <- (244207-192-1)/192 * total_r2 / (1 - total_r2)

chr_tot_het_F_outliers$total_r2 <- total_r2
chr_tot_het_F_outliers$F_stat<- F_stat

write.csv(chr_tot_het_F_outliers, "sensitivity_chrono_totaltest_F_supplementary_outliers_removed.csv", row.names = FALSE, quote = FALSE)

# TESTOSTERONE > CHRONOTYPE #

####
# Read in exposure instrument & outcome data 
####

exposure_dat <- read_exposure_data("totaltest_female_inv_fin.txt", sep = "\t",
                                   snp_col = "SNP", beta_col = "BETA",
                                   se_col = "SE", effect_allele_col = "ALLELE1",
                                   other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")

outcome_dat <- read_outcome_data("outcome_totaltest_chrono_female.txt", sep = "\t",
                                 snp_col = "SNP", beta_col = "BETA",
                                 se_col = "SE", effect_allele_col = "ALLELE1",
                                 other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")

####
# Harmonise data & run MR
####

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dat <- clump_data(dat)

res_test <- mr(dat)
or_tottest_chr_F <- generate_odds_ratios(res_test)

#write.csv(or_tottest_chr_F, "or_totaltest_chrono_female_bdMR.csv", row.names=F, quote=F)

# # SENSITIVITY TESTS (RAW) # #

tot_chr_het_raw <- mr_heterogeneity(dat)
tot_chr _het_raw$I2 <- ((tot_chr_het_raw$Q-tot_chr_het_raw$Q_df)/tot_chr_het_raw$Q)*100
tot_chr _het_raw$outcome= "Chronotype" # Update for different datasets

tot_chr_pleio_raw <- mr_pleiotropy_test(dat)

dat$samplesize.exposure <- 199956 # update samplesize depending on dataset used
dat$samplesize.outcome <- 244207  # update samplesize depending on dataset used

dat_steiger <- steiger_filtering(dat)
total_r2 <- sum(dat_steiger$rsq.exposure) 
F_stat <- (199956-140-1)/140 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

tot_chr_het_raw$r2 <- total_r2
tot_chr_het_raw$F_stat <- F_stat

total_r2_out <- sum(dat$rsq.outcome) 
steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)
steiger_sens$sensitivity_ratio
table(dat$steiger_dir)

#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  SignifThreshold = 0.05)

# Determine SNPs removed by steiger filtering #

dat_removed <- subset(dat, steiger_dir==FALSE)
SNPs <- as.data.frame(dat_removed$SNP)

# Harmonised data with steiger filtered SNPs removed #

dat_filtered <- subset(dat, steiger_dir==TRUE)

# SENSITIVITY TESTS (STEIGER-FILTERED) #  

tot_chr_het_filt <- mr_heterogeneity(dat_filtered)
tot_chr_pleio_filtered <- mr_pleiotropy_test(dat_filtered)

total_r2 <- sum(dat_filtered$rsq.exposure) 
F_stat <- (199956-139-1)/139 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

total_r2_out <- sum(dat_filtered$rsq.outcome) 
steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)
steiger_sens$sensitivity_ratio

tot_chr_het_filt$total_r2<- total_r2
tot_chr_het_filt$F_stat<- F_stat

res_test <- mr(dat_filtered)
or_tot_chr_F_filt <- generate_odds_ratios(res_test)

#write.csv(tot_chr_het_raw, "sensitivity_totaltest_chrono_F_raw_supplementary.csv", row.names = FALSE, quote = FALSE)
write.csv(tot_chr_het_filt, "sensitivity_totaltest_chrono_F_filt_supplementary.csv", row.names = FALSE, quote = FALSE)

#write.csv(or_tot_chrtest_F, "or_totaltest_chrono_F_raw_bdMR.csv", row.names = FALSE, quote = FALSE)
write.csv(or_tot_chrtest_F, "or_totaltest_chrono_F_filt_bdMR.csv", row.names = FALSE, quote = FALSE)

# RADIAL MR #

res_single1 <- mr_singlesnp(dat_filtered)
dat1R <- dat_filtered[dat_filtered$SNP%in%res_single1$SNP,]

raddat1 <- format_radial(dat1R$beta.exposure, dat1R$beta.outcome, dat1R$se.exposure, dat1R$se.outcome, dat1R$SNP)
ivwrad1 <- ivw_radial(raddat1, alpha=0.05/139, weights=3)
dim(ivwrad1$outliers)[1] 

eggrad1 <- egger_radial(raddat1, alpha=0.05/139, weights=3)
eggrad1$coef 
dim(eggrad1$outliers)[1] 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
#plot_radial(c(ivwrad1,eggrad1), TRUE, FALSE, TRUE)

#ivwrad1$qstatistic 
#ivwrad1$sortoutliers <- ivwrad1$outliers[order(ivwrad1$outliers$p.value),]
#ivwrad1$sortoutliers$Qsum <- cumsum(ivwrad1$sortoutliers$Q_statistic)
#ivwrad1$sortoutliers$Qdif <- ivwrad1$sortoutliers$Qsum - ivwrad1$qstatistic
#write.csv(ivwrad$sortoutliers, "chrono_outliers_001.csv", row.names=F, quote=F)

#eggrad1$qstatistic 
#eggrad1$sortoutliers <- eggrad1$outliers[order(eggrad1$outliers$p.value),]
#eggrad1$sortoutliers$Qsum <- cumsum(eggrad1$sortoutliers$Q_statistic)
#eggrad1$sortoutliers$Qdif <- eggrad1$sortoutliers$Qsum - eggrad1$qstatistic
#write.csv(eggrad$sortoutliers, "chrono_outliers_egger_01.csv", row.names=F, quote=F)

# REMOVE TOP OUTLIERS #

dat1_outliers <- dat1R[!dat1R$SNP %in% ivwrad1$outliers$SNP,]
mr_dat1_outliers <- mr(dat1_outliers)
or_dat1_outliers <- generate_odds_ratios(mr_dat1_outliers)

#write.csv(dat1_outliers, "totaltest_chrono_F_filtered_outliers_removed.csv", row.names=F, quote=F)
write.csv(or_dat1_outliers, "or_totaltest_chrono_F_filtered_outliers_removed.csv", row.names=F, quote=F)

# SENSITIVITY TESTS (OUTLIERS REMOVED) #

tot_chr_het_outliers <- mr_heterogeneity(dat1_outliers)
tot_chr_het_outliers$I2 <- ((tot_chr_het_outliers$Q-tot_chr_het_outliers$Q_df)/tot_chr_het_outliers$Q)*100
tot_chr_het_outliers$outcome="Chronotype" # Update for different datasets

tot_chr_pleio_outliers <- mr_pleiotropy_test(dat1_outliers)

total_r2 <- sum(dat1_outliers$rsq.exposure) 
F_stat <- (199956-130-1)/130 * total_r2 / (1 - total_r2)

tot_chr_het_F_outliers$total_r2 <- total_r2
tot_chr_het_F_outliers$F_stat<- F_stat

write.csv(tot_chr_het_F_outliers, "sensitivity_totaltest_chrono_F_supplementary_outliers_removed.csv", row.names = FALSE, quote = FALSE)

##########################
##########################
### MALE
##########################
##########################

# Chronotype <> Total Testosterone
# Chronotype <> Bioavailable Testosterone

# Update script throughout depending on data being used #


# CHRONOTYPE > BIOAVAILABLE TESTOSTERONE #

####
# Read in exposure instrument & outcome data 
####

exposure_dat <- read_exposure_data("chronotype_male.txt", sep = "\t",snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")

outcome_dat <- read_outcome_data("outcome_chrono_biotest_M.txt", sep = "\t", snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")

####
# Harmonise data & run MR
####

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dat <- clump_data(dat)

res_test <- mr(dat)
or_chr_biotest_F <- generate_odds_ratios(res_test)

#write.csv(or_chr_biotest_M, "or_chrono_biotest_male_bdMR.csv", row.names=F, quote=F)

# SENSITIVITY TESTS (RAW) #

chr_bio_het_raw <- mr_heterogeneity(dat)
chr_bio_het_raw$I2 <- ((chr_bio_het_raw$Q-chr_bio_het_raw$Q_df)/chr_bio_het_raw$Q)*100
chr_bio_het_raw$outcome="Bioavailable Testosterone" # Update for different datasets

chr_bio_pleio_raw <- mr_pleiotropy_test(dat)

dat$samplesize.exposure <- 205527 # update samplesize depending on dataset used
dat$samplesize.outcome <- 184204  # update samplesize depending on dataset used

dat_steiger <- steiger_filtering(dat)
total_r2 <- sum(dat_steiger$rsq.exposure) 
F_stat <- (205527-212-1)/212 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

chr_bio_het_raw$r2 <- total_r2
chr_bio_het_raw$F_stat <- F_stat

total_r2_out <- sum(dat$rsq.outcome) 
steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)
steiger_sens$sensitivity_ratio
table(dat$steiger_dir)

#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  SignifThreshold = 0.05)

# Determine SNPs removed by steiger filtering #

dat_removed <- subset(dat, steiger_dir==FALSE)
SNPs <- as.data.frame(dat_removed$SNP)

# Harmonised data with steiger filtered SNPs removed #

dat_filtered <- subset(dat, steiger_dir==TRUE)

# SENSITIVITY TESTS (STEIGER-FILTERED) #  

chr_bio_het_filt <- mr_heterogeneity(dat_filtered)
chr_bio_pleio_filtered <- mr_pleiotropy_test(dat_filtered)

total_r2 <- sum(dat_filtered$rsq.exposure) 
F_stat <- (205527-202-1)/202 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

total_r2_out <- sum(dat_filtered$rsq.outcome) 
steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)
steiger_sens$sensitivity_ratio

chr_bio_het_filt$total_r2<- total_r2
chr_bio_het_filt$F_stat<- F_stat

res_test <- mr(dat_filtered)
or_chr_biotest_F_filt <- generate_odds_ratios(res_test)

#write.csv(chr_bio_het_raw, "sensitivity_chrono_biotest_M_raw_supplementary.csv", row.names = FALSE, quote = FALSE)
write.csv(chr_bio_het_filt, "sensitivity_chrono_biotest_M_filt_supplementary.csv", row.names = FALSE, quote = FALSE)

#write.csv(or_chr_biotest_F, "or_chrono_biotest_M_raw_bdMR.csv", row.names = FALSE, quote = FALSE)
write.csv(or_chr_biotest_F, "or_chrono_biotest_M_filt_bdMR.csv", row.names = FALSE, quote = FALSE)

# RADIAL MR #

res_single1 <- mr_singlesnp(dat_filtered)
dat1R <- dat_filtered[dat_filtered$SNP%in%res_single1$SNP,]

raddat1 <- format_radial(dat1R$beta.exposure, dat1R$beta.outcome, dat1R$se.exposure, dat1R$se.outcome, dat1R$SNP)
ivwrad1 <- ivw_radial(raddat1, alpha=0.05/195, weights=3)
dim(ivwrad1$outliers)[1] 

eggrad1 <- egger_radial(raddat1, alpha=0.05/195, weights=3)
eggrad1$coef 
dim(eggrad1$outliers)[1] 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
#plot_radial(c(ivwrad1,eggrad1), TRUE, FALSE, TRUE)

#ivwrad1$qstatistic 
#ivwrad1$sortoutliers <- ivwrad1$outliers[order(ivwrad1$outliers$p.value),]
#ivwrad1$sortoutliers$Qsum <- cumsum(ivwrad1$sortoutliers$Q_statistic)
#ivwrad1$sortoutliers$Qdif <- ivwrad1$sortoutliers$Qsum - ivwrad1$qstatistic
#write.csv(ivwrad$sortoutliers, "chrono_outliers_001.csv", row.names=F, quote=F)

#eggrad1$qstatistic 
#eggrad1$sortoutliers <- eggrad1$outliers[order(eggrad1$outliers$p.value),]
#eggrad1$sortoutliers$Qsum <- cumsum(eggrad1$sortoutliers$Q_statistic)
#eggrad1$sortoutliers$Qdif <- eggrad1$sortoutliers$Qsum - eggrad1$qstatistic
#write.csv(eggrad$sortoutliers, "chrono_outliers_egger_01.csv", row.names=F, quote=F)

# REMOVE TOP OUTLIERS #

dat1_outliers <- dat1R[!dat1R$SNP %in% ivwrad1$outliers$SNP,]
mr_dat1_outliers <- mr(dat1_outliers)
or_dat1_outliers <- generate_odds_ratios(mr_dat1_outliers)

#write.csv(dat1_outliers, "chrono_biotest_M_filtered_outliers_removed.csv", row.names=F, quote=F)
write.csv(or_dat1_outliers, "or_chrono_biotest_M_filtered_outliers_removed.csv", row.names=F, quote=F)

# SENSITIVITY TESTS (OUTLIERS REMOVED) #

chr_bio_het_outliers <- mr_heterogeneity(dat1_outliers)
chr_bio_het_outliers$I2 <- ((chr_bio_het_outliers$Q-chr_bio_het_outliers$Q_df)/chr_bio_het_outliers$Q)*100
chr_bio_het_outliers$outcome="Bioavailable Testosterone" # Update for different datasets

chr_bio_pleio_outliers <- mr_pleiotropy_test(dat1_outliers)

total_r2 <- sum(dat1_outliers$rsq.exposure) 
F_stat <- (205527-193-1)/193 * total_r2 / (1 - total_r2)

chr_bio_het_M_outliers$total_r2 <- total_r2
chr_bio_het_M_outliers$F_stat<- F_stat

write.csv(chr_bio_het_M_outliers, "sensitivity_chrono_biotest_M_supplementary_outliers_removed.csv", row.names = FALSE, quote = FALSE)

# TESTOSTERONE > CHRONOTYPE #

####
# Read in exposure instrument & outcome data 
####

exposure_dat <- read_exposure_data("biotest_female.txt", sep = "\t",
                                   snp_col = "SNP", beta_col = "BETA",
                                   se_col = "SE", effect_allele_col = "ALLELE1",
                                   other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")

outcome_dat <- read_outcome_data("outcome_biotest_chrono_male.txt", sep = "\t",
                                 snp_col = "SNP", beta_col = "BETA",
                                 se_col = "SE", effect_allele_col = "ALLELE1",
                                 other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")

####
# Harmonise data & run MR
####

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dat <- clump_data(dat)

res_test <- mr(dat)
or_biotest_chr_M <- generate_odds_ratios(res_test)

#write.csv(or_biotest_chr_M, "or_biotest_chrono_male_bdMR.csv", row.names=F, quote=F)

# SENSITIVITY TESTS (RAW) #

bio_chr_het_raw <- mr_heterogeneity(dat)
bio_chr _het_raw$I2 <- ((bio_chr_het_raw$Q-bio_chr_het_raw$Q_df)/bio_chr_het_raw$Q)*100
bio_chr _het_raw$outcome= "Chronotype" # Update for different datasets

bio_chr_pleio_raw <- mr_pleiotropy_test(dat)

dat$samplesize.exposure <- 184204 # update samplesize depending on dataset used
dat$samplesize.outcome <- 205527  # update samplesize depending on dataset used

dat_steiger <- steiger_filtering(dat)
total_r2 <- sum(dat_steiger$rsq.exposure) 
F_stat <- (184204-70-1)/70 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

bio_chr_het_raw$r2 <- total_r2
bio_chr_het_raw$F_stat <- F_stat

total_r2_out <- sum(dat$rsq.outcome) 
steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)
steiger_sens$sensitivity_ratio
table(dat$steiger_dir)

#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  SignifThreshold = 0.05)

# Determine SNPs removed by steiger filtering #

dat_removed <- subset(dat, steiger_dir==FALSE)
SNPs <- as.data.frame(dat_removed$SNP)

# Harmonised data with steiger filtered SNPs removed #

dat_filtered <- subset(dat, steiger_dir==TRUE)

# SENSITIVITY TESTS (STEIGER-FILTERED) #  

bio_chr_het_filt <- mr_heterogeneity(dat_filtered)
bio_chr_pleio_filtered <- mr_pleiotropy_test(dat_filtered)

total_r2 <- sum(dat_filtered$rsq.exposure) 
F_stat <- (184204-69-1)/69 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

total_r2_out <- sum(dat_filtered$rsq.outcome) 
steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)
steiger_sens$sensitivity_ratio

bio_chr_het_filt$total_r2<- total_r2
bio_chr_het_filt$F_stat<- F_stat

res_test <- mr(dat_filtered)
or_bio_chr_M_filt <- generate_odds_ratios(res_test)

#write.csv(bio_chr_het_raw, "sensitivity_biotest_chrono_M_raw_supplementary.csv", row.names = FALSE, quote = FALSE)
write.csv(bio_chr_het_filt, "sensitivity_biotest_chrono_M_filt_supplementary.csv", row.names = FALSE, quote = FALSE)

#write.csv(or_bio_chr_M, "or_biotest_chrono_M_raw_bdMR.csv", row.names = FALSE, quote = FALSE)
write.csv(or_bio_chr_M, "or_biotest_chrono_M_filt_bdMR.csv", row.names = FALSE, quote = FALSE)

# RADIAL MR #

res_single1 <- mr_singlesnp(dat_filtered)
dat1R <- dat_filtered[dat_filtered$SNP%in%res_single1$SNP,]

raddat1 <- format_radial(dat1R$beta.exposure, dat1R$beta.outcome, dat1R$se.exposure, dat1R$se.outcome, dat1R$SNP)
ivwrad1 <- ivw_radial(raddat1, alpha=0.05/66, weights=3)
dim(ivwrad1$outliers)[1] 

eggrad1 <- egger_radial(raddat1, alpha=0.05/66, weights=3)
eggrad1$coef 
dim(eggrad1$outliers)[1] 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
#plot_radial(c(ivwrad1,eggrad1), TRUE, FALSE, TRUE)

#ivwrad1$qstatistic 
#ivwrad1$sortoutliers <- ivwrad1$outliers[order(ivwrad1$outliers$p.value),]
#ivwrad1$sortoutliers$Qsum <- cumsum(ivwrad1$sortoutliers$Q_statistic)
#ivwrad1$sortoutliers$Qdif <- ivwrad1$sortoutliers$Qsum - ivwrad1$qstatistic
#write.csv(ivwrad$sortoutliers, "chrono_outliers_001.csv", row.names=F, quote=F)

#eggrad1$qstatistic 
#eggrad1$sortoutliers <- eggrad1$outliers[order(eggrad1$outliers$p.value),]
#eggrad1$sortoutliers$Qsum <- cumsum(eggrad1$sortoutliers$Q_statistic)
#eggrad1$sortoutliers$Qdif <- eggrad1$sortoutliers$Qsum - eggrad1$qstatistic
#write.csv(eggrad$sortoutliers, "chrono_outliers_egger_01.csv", row.names=F, quote=F)

# REMOVE TOP OUTLIERS #

dat1_outliers <- dat1R[!dat1R$SNP %in% ivwrad1$outliers$SNP,]
mr_dat1_outliers <- mr(dat1_outliers)
or_dat1_outliers <- generate_odds_ratios(mr_dat1_outliers)

#write.csv(dat1_outliers, "biotest_chrono_M_filtered_outliers_removed.csv", row.names=F, quote=F)
write.csv(or_dat1_outliers, "or_biotest_chrono_M_filtered_outliers_removed.csv", row.names=F, quote=F)

# SENSITIVITY TESTS (OUTLIERS REMOVED) #

bio_chr_het_outliers <- mr_heterogeneity(dat1_outliers)
bio_chr_het_outliers$I2 <- ((bio_chr_het_outliers$Q-bio_chr_het_outliers$Q_df)/bio_chr_het_outliers$Q)*100
bio_chr_het_outliers$outcome= "Chronotype" # Update for different datasets

bio_chr_pleio_outliers <- mr_pleiotropy_test(dat1_outliers)

total_r2 <- sum(dat1_outliers$rsq.exposure) 
F_stat <- (184204-66-1)/66 * total_r2 / (1 - total_r2)

bio_chr_het_M_outliers$total_r2 <- total_r2
bio_chr_het_M_outliers$F_stat<- F_stat

write.csv(bio_chr_het_M_outliers, "sensitivity_biotest_chrono_M_supplementary_outliers_removed.csv", row.names = FALSE, quote = FALSE)
