
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

# Chronotype > Breast Cancer (overall)
# Total Testosterone > Breast Cancer (overall)
# Bioavailable Testosterone > Breast Cancer (overall)
# SHBG > Breast Cancer (overall)
# Oestradiol > Breast Cancer (overall)

# Update script throughout depending on data being used #


####
# Read in exposure instrument & outcome data 
####

exposure_dat <- read_exposure_data("chronotype_female.txt", sep = "\t", ,snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col = "eaf", pval_col = "pval")

outcome_dat <- read_outcome_data("outcome_chrono_BrCa_overall.txt", sep = "\t",
                                 snp_col = "SNP", beta_col = "Beta.meta",
                                 se_col = "sdE.meta", effect_allele_col = "Effect.Onco",
                                 other_allele_col = "Baseline.Onco", eaf_col = "EAFcontrols.Onco", pval_col = "p.meta")

####
# Harmonise data & run MR
####

dat <- harmonise_data(exposure_dat, outcome_dat, action= 2)
dat <- clump_data(dat)

res_chronotype <- mr(dat)
or_chronotype <- generate_odds_ratios(res_chronotype)

#write.csv(or_chronotype, "or_chronotype_BrCa_overall.csv", row.names=F, quote=F)

# SENSITIVITY TESTS #

het_overall <- mr_heterogeneity(dat)
het_overall$I2 <- ((het_overall$Q-het_overall$Q_df)/het_overall$Q)*100
het_overall$outcome="Overall Breast Cancer" # Update for different datasets

pleio_overall <- mr_pleiotropy_test(dat)

dat$samplesize.exposure <- 244207 # update samplesize depending on dataset used
dat$samplesize.outcome <- 247173  # update samplesize depending on dataset used

dat_steiger <- steiger_filtering(dat)
total_r2 <- sum(dat_steiger$rsq.exposure) 
F_stat <- (244207-198-1)/198 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

##total_r2_out <- sum(dat$rsq.outcome) 
##steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)

het_overall$r2 <- total_r2
het_overall$F_stat <- F_stat

#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  SignifThreshold = 0.05)

####
# RADIAL MR 
####

res_single1 <- mr_singlesnp(dat)
dat1R <- dat[dat$SNP%in%res_single1$SNP,]

raddat1 <- format_radial(dat1R$beta.exposure, dat1R$beta.outcome, dat1R$se.exposure, dat1R$se.outcome, dat1R$SNP)
ivwrad1 <- ivw_radial(raddat1, alpha=0.05/190, weights=3) # Where alpha=0.05/nSNP
dim(ivwrad1$outliers)[1] 

eggrad1 <- egger_radial(raddat1, alpha=0.05/190, weights=3) # Where alpha=0.05/nSNP
eggrad1$coef 
dim(eggrad1$outliers)[1] 

# If outliers detected # 

SNPs <- as.data.frame(ivwrad1$outliers$SNP)

#plot_radial(ivwrad1, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad1,eggrad1), TRUE, FALSE, TRUE)

#ivwrad1$qstatistic 
#ivwrad1$sortoutliers <- ivwrad1$outliers[order(ivwrad1$outliers$p.value),]
#ivwrad1$sortoutliers$Qsum <- cumsum(ivwrad1$sortoutliers$Q_statistic)
#ivwrad1$sortoutliers$Qdif <- ivwrad1$sortoutliers$Qsum - ivwrad1$qstatistic
#write.csv(ivwrad$sortoutliers, "chronotype_outliers_001.csv", row.names=F, quote=F)

#eggrad1$qstatistic 
#eggrad1$sortoutliers <- eggrad1$outliers[order(eggrad1$outliers$p.value),]
#eggrad1$sortoutliers$Qsum <- cumsum(eggrad1$sortoutliers$Q_statistic)
#eggrad1$sortoutliers$Qdif <- eggrad1$sortoutliers$Qsum - eggrad1$qstatistic
#write.csv(eggrad$sortoutliers, "chronotype_outliers_egger_01.csv", row.names=F, quote=F)

# REMOVE TOP OUTLIERS #

dat1_outliers <- dat1R[!dat1R$SNP %in% ivwrad1$outliers$SNP,]
mr_dat1_outliers <- mr(dat1_outliers)
or_dat1_outliers <- generate_odds_ratios(mr_dat1_outliers)

#write.csv(or_dat1_outliers, "or_chronotype_BrCa_overall_outliers_removed.csv", row.names=F, quote=F)
#write.csv(ivwrad1$outliers$SNP, "chronotype_BrCa_outlier_snps_removed.csv", row.names=F, quote=F)

##########################
##########################
### MALE
##########################
##########################

# Chronotype > Prostate Cancer
# Total Testosterone > Prostate Cancer
# Bioavailable Testosterone > Prostate Cancer
# SHBG > Prostate Cancer
# Oestradiol > Prostate Cancer

# Update script throughout depending on data being used #


####
# Read in exposure instrument & outcome data 
####

exposure_dat <- read_exposure_data("chronotype_male.txt", sep = "\t",snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM")

outcome_dat <- extract_outcome_data(exposure_dat$SNP, 'ieu-a-1174') # SNPs extracted frm prostate cancer GWAS in MR-BASE

dat <- harmonise_data(exposure_dat, outcome_dat, action=2)
dat <- clump_data(dat)

res_chronotype <- mr(dat)
or_chronotype <- generate_odds_ratios(res_chronotype)

write.csv(or_chronotype, "or_chronotype_PrCa_uvMR.csv", row.names=F, quote=F)

# SENSITIVITY TESTS #

het_overall <- mr_heterogeneity(dat)
het_overall$I2 <- ((het_overall$Q-het_overall$Q_df)/het_overall$Q)*100
het_overall$outcome="Overall Prostate Cancer" # Update for different datasets

pleio_overall <- mr_pleiotropy_test(dat)

dat$samplesize.exposure <- 180386 # update samplesize depending on dataset used
dat$samplesize.outcome <- 140254  # update samplesize depending on dataset used

dat_steiger <- steiger_filtering(dat)
total_r2 <- sum(dat_steiger$rsq.exposure) 
F_stat <- (180386-203-1)/203 * total_r2 / (1 - total_r2) # update samplesize and nSNP depending on data used

##total_r2_out <- sum(dat$rsq.outcome) 
##steiger_sens <- steiger_sensitivity(total_r2, total_r2_out)

het_overall$r2 <- total_r2
het_overall$F_stat <- F_stat

#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  SignifThreshold = 0.05)

####
# RADIAL MR 
####

res_single1 <- mr_singlesnp(dat)
dat1R <- dat[dat$SNP%in%res_single1$SNP,]

raddat1 <- format_radial(dat1R$beta.exposure, dat1R$beta.outcome, dat1R$se.exposure, dat1R$se.outcome, dat1R$SNP)
ivwrad1 <- ivw_radial(raddat1, alpha=0.05/195, weights=3) # Where alpha=0.05/nSNP
dim(ivwrad1$outliers)[1] 

eggrad1 <- egger_radial(raddat1, alpha=0.05/195, weights=3) # Where alpha=0.05/nSNP
eggrad1$coef 
dim(eggrad1$outliers)[1] 

# If outliers detected # 

SNPs <- as.data.frame(ivwrad1$outliers$SNP)

#plot_radial(ivwrad1, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad1,eggrad1), TRUE, FALSE, TRUE)

#ivwrad1$qstatistic 
#ivwrad1$sortoutliers <- ivwrad1$outliers[order(ivwrad1$outliers$p.value),]
#ivwrad1$sortoutliers$Qsum <- cumsum(ivwrad1$sortoutliers$Q_statistic)
#ivwrad1$sortoutliers$Qdif <- ivwrad1$sortoutliers$Qsum - ivwrad1$qstatistic
#write.csv(ivwrad$sortoutliers, "chronotype_outliers_001.csv", row.names=F, quote=F)

#eggrad1$qstatistic 
#eggrad1$sortoutliers <- eggrad1$outliers[order(eggrad1$outliers$p.value),]
#eggrad1$sortoutliers$Qsum <- cumsum(eggrad1$sortoutliers$Q_statistic)
#eggrad1$sortoutliers$Qdif <- eggrad1$sortoutliers$Qsum - eggrad1$qstatistic
#write.csv(eggrad$sortoutliers, "chronotype_outliers_egger_01.csv", row.names=F, quote=F)

# REMOVE TOP OUTLIERS #

dat1_outliers <- dat1R[!dat1R$SNP %in% ivwrad1$outliers$SNP,]
mr_dat1_outliers <- mr(dat1_outliers)
or_dat1_outliers <- generate_odds_ratios(mr_dat1_outliers)

#write.csv(or_dat1_outliers, "or_chronotype_PrCa_outliers_removed.csv", row.names=F, quote=F)
#write.csv(ivwrad1$outliers$SNP, "chronotype_PrCa_outlier_snps_removed.csv", row.names=F, quote=F)