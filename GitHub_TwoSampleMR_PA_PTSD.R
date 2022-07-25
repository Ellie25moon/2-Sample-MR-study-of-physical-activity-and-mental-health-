####################################################################################
##################### Physical Activity -> PTSD ####################################
####################################################################################

### set working directory
setwd("~/Documents/Study1_GWAS")

### load packages 
library(TwoSampleMR)
library(LDlinkR)
library (dplyr)
library(mr.raps)
library(MRPRESSO)
library(ggpubr)
library(writexl)
library(readr)

### Upload GWAS summary stats for the outcome and save as csv file  

# note: (i) edit name of the summary stats file
# (ii) you might need to use a different function to read the data depending on the format of the file
ptsd_full<- read_tsv("", col_names = T) 

head(ptsd_full)
colnames(ptsd_full)

## select and rename columns needed for the analysis 
# note: (i) columns needed to run MR: SNP ID, Effect Allele, Other Allele, BETA, SE, P-value, Allele Frequency
# (i) edit column names of 'ptsd_full' as needed
ptsd <- data.frame(
  SNP = ptsd_full$RSID,
  beta = ptsd_full$BETA,
  se = ptsd_full$SE,
  effect_allele = ptsd_full$REF,
  other_allele = ptsd_full$ALT,
  eaf = ptsd_full$AF, 
  pval = ptsd_full$PVALUE
)

write.csv(ptsd, file= "ptsd.csv")


##### (a) Self-reported moderate-vigorous activity 

### 1. Prepare exposure data 

# first genetic instrument (G1)
exp_modvigpaself1 <- extract_instruments("ebi-a-GCST006097", p1=5e-8)

# second genetic instrument (G2)
exp_modvigpaself2 <- extract_instruments("ebi-a-GCST006097", p1=10e-6)

# clumping 
exp_modvigpaself1_clu<- clump_data(exp_modvigpaself1)
exp_modvigpaself2_clu<- clump_data(exp_modvigpaself2)


### 2. Prepare outcome data 

## G1

out_modvigpaself1_ptsd1 <- read_outcome_data(
  snps = exp_modvigpaself1_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_modvigpaself1_ptsd1)
nrow(exp_modvigpaself1_clu)
# if nrow out = exp, run code on line 80 and then move to G2 (line 139)
# if nrow out != exp -> run code from line 83

out_modvigpaself1_ptsd1_proxy <- out_modvigpaself1_ptsd1


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g1<- exp_modvigpaself1_clu[!exp_modvigpaself1_clu$SNP %in% out_modvigpaself1_ptsd1$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g1$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy1<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy1<- ptsd_proxy1 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_modvigpaself1_ptsd1[1,11] ,
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy1 <- ptsd_proxy1[!ptsd_proxy1$SNP %in% out_modvigpaself1_ptsd1$SNP, ]

# combine with outcome data
out_modvigpaself1_ptsd1_proxy<- rbind(out_modvigpaself1_ptsd1,ptsd_proxy1)


### G2 

out_modvigpaself1_ptsd2 <- read_outcome_data(
  snps = exp_modvigpaself2_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_modvigpaself1_ptsd2)
nrow(exp_modvigpaself2_clu)
# if nrow out = exp, run code on line 159 and then move to part 3 (data harmonization, line 218)
# if nrow out != exp -> run code from line 162

out_modvigpaself1_ptsd2_proxy <- out_modvigpaself1_ptsd2


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g2<- exp_modvigpaself2_clu[!exp_modvigpaself2_clu$SNP %in% out_modvigpaself1_ptsd2$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g2$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy2<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy2<- ptsd_proxy2 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_modvigpaself1_ptsd2[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy2 <- ptsd_proxy2[!ptsd_proxy2$SNP %in% out_modvigpaself1_ptsd2$SNP, ]

# combine with outcome data
out_modvigpaself1_ptsd2_proxy<- rbind(out_modvigpaself1_ptsd2,ptsd_proxy2)


### 3. Harmonise data 

dat_modvigpaself_ptsd1 <- harmonise_data(
  exposure_dat = exp_modvigpaself1_clu, 
  outcome_dat = out_modvigpaself1_ptsd1_proxy, 
  action =2)

dat_modvigpaself_ptsd2 <- harmonise_data(
  exposure_dat = exp_modvigpaself2_clu, 
  outcome_dat = out_modvigpaself1_ptsd2_proxy, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets

## Add column for total sample size of outcome GWAS  
# Note: edit sample size accordingly 
dat_modvigpaself_ptsd1$samplesize.outcome<- XXX
dat_modvigpaself_ptsd2$samplesize.outcome<- XXX

## Pruning 
dat_modvigpaself_ptsd1<-power_prune(dat_modvigpaself_ptsd1,method=1,dist.outcome="binary")
dat_modvigpaself_ptsd2<-power_prune(dat_modvigpaself_ptsd2,method=1,dist.outcome="binary")

nrow(dat_modvigpaself_ptsd1)
nrow(dat_modvigpaself_ptsd2)

## Edit names of outcome and exposure columns 
dat_modvigpaself_ptsd1$exposure[dat_modvigpaself_ptsd1$exposure == "Moderate to vigorous physical activity levels || id:ebi-a-GCST006097"] <- "Self-reported moderate-vigorous PA (G1)"
dat_modvigpaself_ptsd1$outcome[dat_modvigpaself_ptsd1$outcome == "outcome"] <- "PTSD"

dat_modvigpaself_ptsd2$exposure[dat_modvigpaself_ptsd2$exposure ==" || id:ebi-a-GCST006097"] <- "Self-reported moderate-vigorous PA (G2)"
dat_modvigpaself_ptsd2$outcome[dat_modvigpaself_ptsd2$outcome == "outcome"] <- "PTSD"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_modvigpaself_ptsd1<-mr(dat_modvigpaself_ptsd1, method_list=c("mr_ivw_mre",
                                                                    "mr_egger_regression",
                                                                    "mr_weighted_median", 
                                                                    "mr_weighted_mode"))
results_modvigpaself_ptsd1

results_modvigpaself_ptsd2<-mr(dat_modvigpaself_ptsd2,  method_list=c("mr_ivw_mre", "mr_egger_regression",
                                                                     "mr_weighted_median", 
                                                                     "mr_weighted_mode"))
results_modvigpaself_ptsd2

## raps 

# g1
raps<- mr.raps(dat_modvigpaself_ptsd1$beta.exposure, dat_modvigpaself_ptsd1$beta.outcome, dat_modvigpaself_ptsd1$se.exposure, 
               dat_modvigpaself_ptsd1$se.outcome)

raps_df<- data.frame(id.exposure = results_modvigpaself_ptsd1[1,1],
                     id.outcome = results_modvigpaself_ptsd1[1,2], 
                     outcome = results_modvigpaself_ptsd1[1,3],
                     exposure = results_modvigpaself_ptsd1[1,4], 
                     method = "raps", 
                     nsnp = results_modvigpaself_ptsd1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_modvigpaself_ptsd1<-rbind(results_modvigpaself_ptsd1,raps_df )
results_modvigpaself_ptsd1$analysis = "modvigpaself_ptsd1"
results_modvigpaself_ptsd1 <- generate_odds_ratios(results_modvigpaself_ptsd1)

#g2
raps2<- mr.raps(dat_modvigpaself_ptsd2$beta.exposure, dat_modvigpaself_ptsd2$beta.outcome, dat_modvigpaself_ptsd2$se.exposure, 
                dat_modvigpaself_ptsd2$se.outcome)

raps_df2<- data.frame(id.exposure = results_modvigpaself_ptsd2[1,1],
                      id.outcome = results_modvigpaself_ptsd2[1,2], 
                      outcome = results_modvigpaself_ptsd2[1,3],
                      exposure = results_modvigpaself_ptsd2[1,4], 
                      method = "raps", 
                      nsnp = results_modvigpaself_ptsd2[1,6],
                      b = raps[["beta.hat"]], 
                      se = raps[["beta.se"]], 
                      pval = raps[["beta.p.value"]])

results_modvigpaself_ptsd2<-rbind(results_modvigpaself_ptsd2,raps_df2 )
results_modvigpaself_ptsd2$analysis = "modvigpaself_ptsd2"
results_modvigpaself_ptsd2 <- generate_odds_ratios(results_modvigpaself_ptsd2)

# heterogeneity stats 

hetero_modvigpaself_ptsd1<- mr_heterogeneity(dat_modvigpaself_ptsd1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_modvigpaself_ptsd1$analysis<- "modvigpaself_ptsd1"

hetero_modvigpaself_ptsd2<- mr_heterogeneity(dat_modvigpaself_ptsd2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_modvigpaself_ptsd2$analysis<- "modvigpaself_ptsd2"

# check intercept of Egger regression 

intercept_modvigpaself_ptsd1<- mr_pleiotropy_test(dat_modvigpaself_ptsd1)
intercept_modvigpaself_ptsd1$analysis<- "modvigpaself_ptsd1"

intercept_modvigpaself_ptsd2<- mr_pleiotropy_test(dat_modvigpaself_ptsd2)
intercept_modvigpaself_ptsd2$analysis<- "modvigpaself_ptsd2"


# steiger directionality test 

steiger_modvigpaself_ptsd1 <- directionality_test(dat_modvigpaself_ptsd1)
steiger_modvigpaself_ptsd1$analysis <- "modvigpaself_ptsd1"

steiger_modvigpaself_ptsd2 <- directionality_test(dat_modvigpaself_ptsd2)
steiger_modvigpaself_ptsd2$analysis <- "modvigpaself_ptsd2"

# MR presso 

#g1
presso<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                  SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_modvigpaself_ptsd1, 
                  NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso[["Main MR results"]]
presso$`MR-PRESSO results`$`Global Test`$RSSobs
presso$`MR-PRESSO results`$`Global Test`$Pvalue

presso_modvigpaself_ptsd1<- presso[["Main MR results"]]
presso_modvigpaself_ptsd1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                       presso$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_modvigpaself_ptsd1$analysis<- "modvigpaself_ptsd1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_modvigpaself_ptsd2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_modvigpaself_ptsd2<- presso2[["Main MR results"]]
presso_modvigpaself_ptsd2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                       presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_modvigpaself_ptsd2$analysis<- "modvigpaself_ptsd2"

### 6.Graphs 

#Scatterplot

#g1
res <- mr(dat_modvigpaself_ptsd1, method_list=c("mr_ivw_mre",
                                               "mr_egger_regression",
                                               "mr_weighted_median", 
                                               "mr_weighted_mode"))
p1g1<-mr_scatter_plot(res, dat_modvigpaself_ptsd1)

#g2
res2 <- mr(dat_modvigpaself_ptsd2, method_list=c("mr_ivw_mre",
                                                "mr_egger_regression",
                                                "mr_weighted_median", 
                                                "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_modvigpaself_ptsd2)


# combine 

plot_modvigpaself_ptsd12 <- ggarrange(p1g1[[paste(exp_modvigpaself1[1,7],".",out_modvigpaself1_ptsd1[1,11], sep="")]],
                                     p1g2[[paste(exp_modvigpaself2[1,7],".",out_modvigpaself1_ptsd2[1,11], sep="")]],
                                     labels = c("(a) Scatter plots", ""),
                                     hjust = -0.05,
                                     ncol = 2, nrow = 1,
                                     common.legend =F)
plot_modvigpaself_ptsd12
ggsave("plot_modvigpaself_ptsd12.jpeg", device = "jpeg",dpi = 300, width =15, height =7, limitsize = F)


##### (b) Accelerometer-based physical activity (mean acceleration) 

### 1.Prepare exposure data 

# first genetic instrument (G1)
exp_accpamean1 <- extract_instruments("ebi-a-GCST006099", p1=5e-8)

# second genetic instrument (G2)
exp_accpamean2 <- extract_instruments("ebi-a-GCST006099", p1=10e-6)

# clumping 
exp_accpamean1_clu<- clump_data(exp_accpamean1)
exp_accpamean2_clu<- clump_data(exp_accpamean2)


### 2. Prepare outcome data 

## G1

out_accpamean1_ptsd1 <- read_outcome_data(
  snps = exp_accpamean1_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_accpamean1_ptsd1)
nrow(exp_accpamean1_clu)
# if nrow out = exp, run code on line 436 and then move to G2 (line 495)
# if nrow out != exp -> run code from line 439

out_accpamean1_ptsd1_proxy <- out_accpamean1_ptsd1


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g1<- exp_accpamean1_clu[!exp_accpamean1_clu$SNP %in% out_accpamean1_ptsd1$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g1$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy1<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy1<- ptsd_proxy1 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_accpamean1_ptsd1[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy1 <- ptsd_proxy1[!ptsd_proxy1$SNP %in% out_accpamean1_ptsd1$SNP, ]

# combine with outcome data
out_accpamean1_ptsd1_proxy<- rbind(out_accpamean1_ptsd1,ptsd_proxy1)


### G2 

out_accpamean1_ptsd2 <- read_outcome_data(
  snps = exp_accpamean2_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_accpamean1_ptsd2)
nrow(exp_accpamean2_clu)
# if nrow out = exp, run code on line 515 and then move to part 3 (data harmonization, line 574)
# if nrow out != exp -> run code from line 518

out_accpamean1_ptsd2_proxy <- out_accpamean1_ptsd2


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g2<- exp_accpamean2_clu[!exp_accpamean2_clu$SNP %in% out_accpamean1_ptsd2$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g2$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy2<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy2<- ptsd_proxy2 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_accpamean1_ptsd2[1,11] ,  
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy2 <- ptsd_proxy2[!ptsd_proxy2$SNP %in% out_accpamean1_ptsd2$SNP, ]

# combine with outcome data
out_accpamean1_ptsd2_proxy<- rbind(out_accpamean1_ptsd2,ptsd_proxy2)


### 3. Harmonise data 

dat_accpamean_ptsd1 <- harmonise_data(
  exposure_dat = exp_accpamean1_clu, 
  outcome_dat = out_accpamean1_ptsd1_proxy, 
  action =2)

dat_accpamean_ptsd2 <- harmonise_data(
  exposure_dat = exp_accpamean2_clu, 
  outcome_dat = out_accpamean1_ptsd2_proxy, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets

## Add column for total sample size of outcome GWAS  
# Note: edit sample size accordingly 
dat_accpamean_ptsd1$samplesize.outcome<- XXX
dat_accpamean_ptsd2$samplesize.outcome<- XXX

## Pruning 
dat_accpamean_ptsd1<-power_prune(dat_accpamean_ptsd1,method=1,dist.outcome="binary")
dat_accpamean_ptsd2<-power_prune(dat_accpamean_ptsd2,method=1,dist.outcome="binary")

nrow(dat_accpamean_ptsd1)
nrow(dat_accpamean_ptsd2)

## Edit names of outcome and exposure columns 
dat_accpamean_ptsd1$exposure[dat_accpamean_ptsd1$exposure == "Accelerometer-based physical activity measurement (average acceleration) || id:ebi-a-GCST006099"] <- "Accelerometer-based PA, mean acceleration (G1)"
dat_accpamean_ptsd1$outcome[dat_accpamean_ptsd1$outcome == "outcome"] <- "PTSD"

dat_accpamean_ptsd2$exposure[dat_accpamean_ptsd2$exposure ==" || id:ebi-a-GCST006099"] <- "Accelerometer-based PA, mean acceleration (G2)"
dat_accpamean_ptsd2$outcome[dat_accpamean_ptsd2$outcome == "outcome"] <- "PTSD"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_accpamean_ptsd1<-mr(dat_accpamean_ptsd1, method_list=c("mr_ivw_mre",
                                                                     "mr_egger_regression",
                                                                     "mr_weighted_median", 
                                                                     "mr_weighted_mode"))
results_accpamean_ptsd1

results_accpamean_ptsd2<-mr(dat_accpamean_ptsd2,  method_list=c("mr_ivw_mre", "mr_egger_regression",
                                                                      "mr_weighted_median", 
                                                                      "mr_weighted_mode"))
results_accpamean_ptsd2

## raps 

# g1
raps<- mr.raps(dat_accpamean_ptsd1$beta.exposure, dat_accpamean_ptsd1$beta.outcome, dat_accpamean_ptsd1$se.exposure, 
               dat_accpamean_ptsd1$se.outcome)

raps_df<- data.frame(id.exposure = results_accpamean_ptsd1[1,1],
                     id.outcome = results_accpamean_ptsd1[1,2], 
                     outcome = results_accpamean_ptsd1[1,3],
                     exposure = results_accpamean_ptsd1[1,4], 
                     method = "raps", 
                     nsnp = results_accpamean_ptsd1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_accpamean_ptsd1<-rbind(results_accpamean_ptsd1,raps_df )
results_accpamean_ptsd1$analysis = "accpamean_ptsd1"
results_accpamean_ptsd1 <- generate_odds_ratios(results_accpamean_ptsd1)

#g2
raps2<- mr.raps(dat_accpamean_ptsd2$beta.exposure, dat_accpamean_ptsd2$beta.outcome, dat_accpamean_ptsd2$se.exposure, 
                dat_accpamean_ptsd2$se.outcome)

raps_df2<- data.frame(id.exposure = results_accpamean_ptsd2[1,1],
                      id.outcome = results_accpamean_ptsd2[1,2], 
                      outcome = results_accpamean_ptsd2[1,3],
                      exposure = results_accpamean_ptsd2[1,4], 
                      method = "raps", 
                      nsnp = results_accpamean_ptsd2[1,6],
                      b = raps[["beta.hat"]], 
                      se = raps[["beta.se"]], 
                      pval = raps[["beta.p.value"]])

results_accpamean_ptsd2<-rbind(results_accpamean_ptsd2,raps_df2 )
results_accpamean_ptsd2$analysis = "accpamean_ptsd2"
results_accpamean_ptsd2 <- generate_odds_ratios(results_accpamean_ptsd2)

# heterogeneity stats 

hetero_accpamean_ptsd1<- mr_heterogeneity(dat_accpamean_ptsd1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_accpamean_ptsd1$analysis<- "accpamean_ptsd1"

hetero_accpamean_ptsd2<- mr_heterogeneity(dat_accpamean_ptsd2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_accpamean_ptsd2$analysis<- "accpamean_ptsd2"

# check intercept of Egger regression 

intercept_accpamean_ptsd1<- mr_pleiotropy_test(dat_accpamean_ptsd1)
intercept_accpamean_ptsd1$analysis<- "accpamean_ptsd1"

intercept_accpamean_ptsd2<- mr_pleiotropy_test(dat_accpamean_ptsd2)
intercept_accpamean_ptsd2$analysis<- "accpamean_ptsd2"


# steiger directionality test 

steiger_accpamean_ptsd1 <- directionality_test(dat_accpamean_ptsd1)
steiger_accpamean_ptsd1$analysis <- "accpamean_ptsd1"

steiger_accpamean_ptsd2 <- directionality_test(dat_accpamean_ptsd2)
steiger_accpamean_ptsd2$analysis <- "accpamean_ptsd2"

# MR presso 

#g1
presso<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                  SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_accpamean_ptsd1, 
                  NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso[["Main MR results"]]
presso$`MR-PRESSO results`$`Global Test`$RSSobs
presso$`MR-PRESSO results`$`Global Test`$Pvalue

presso_accpamean_ptsd1<- presso[["Main MR results"]]
presso_accpamean_ptsd1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                        presso$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_accpamean_ptsd1$analysis<- "accpamean_ptsd1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_accpamean_ptsd2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_accpamean_ptsd2<- presso2[["Main MR results"]]
presso_accpamean_ptsd2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                        presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_accpamean_ptsd2$analysis<- "accpamean_ptsd2"


### 6.Graphs 

#Scatterplot

#g1
res <- mr(dat_accpamean_ptsd1, method_list=c("mr_ivw_mre",
                                                "mr_egger_regression",
                                                "mr_weighted_median", 
                                                "mr_weighted_mode"))
p1g1<-mr_scatter_plot(res, dat_accpamean_ptsd1)

#g2
res2 <- mr(dat_accpamean_ptsd2, method_list=c("mr_ivw_mre",
                                                 "mr_egger_regression",
                                                 "mr_weighted_median", 
                                                 "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_accpamean_ptsd2)


# combine 

plot_accpamean_ptsd12 <- ggarrange(p1g1[[paste(exp_accpamean1[1,7],".",out_accpamean1_ptsd1[1,11], sep="")]],
                                      p1g2[[paste(exp_accpamean2[1,7],".",out_accpamean1_ptsd2[1,11], sep="")]],
                                      labels = c("(a) Scatter plots", ""),
                                      hjust = -0.05,
                                      ncol = 2, nrow = 1,
                                      common.legend =F)
plot_accpamean_ptsd12
ggsave("plot_accpamean_ptsd12.jpeg", device = "jpeg",dpi = 300, width =15, height =7, limitsize = F)


##### (c) Accelerometer-based moderate physical activity 

### 1. Prepare exposure data 

# load excel file 
accmodpa<-read.csv("moderatepa.csv")
names(accmodpa)

# create df and change column names 
accmodpa_df <- data.frame(
  SNP = accmodpa$SNP,
  beta = accmodpa$BETA,
  se = accmodpa$SE,
  effect_allele = accmodpa$ALLELE1,
  other_allele = accmodpa$ALLELE0,
  eaf = accmodpa$A1FREQ, 
  pval = accmodpa$P_BOLT_LMM_INF)

# first exposure instrument (G1)
exp_accmodpa1 <- format_data(accmodpa_df, type="exposure")
exp_accmodpa1<- exp_accmodpa1[exp_accmodpa1$pval.exposure <=5e-8, ]
summary(exp_accmodpa1$pval.exposure)

# second genetic instrument (G2)
exp_accmodpa2 <- format_data(accmodpa_df, type="exposure")
exp_accmodpa2<- exp_accmodpa2[exp_accmodpa2$pval.exposure <=10e-6, ]
summary(exp_accmodpa2$pval.exposure)

# clumping 
exp_accmodpa1_clu<- clump_data(exp_accmodpa1)
exp_accmodpa2_clu<- clump_data(exp_accmodpa2)

# add column for sample size 
exp_accmodpa1_clu$samplesize.exposure<- 91105
exp_accmodpa2_clu$samplesize.exposure<- 91105

### 2. Prepare outcome data 

## G1

out_accmodpa1_ptsd1 <- read_outcome_data(
  snps = exp_accmodpa1_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_accmodpa1_ptsd1)
nrow(exp_accmodpa1_clu)
# if nrow out = exp -> run code on line 814 and then move to G2 (line 873)
# if nrow out != exp -> run code from line 817

out_accmodpa1_ptsd1_proxy <- out_accmodpa1_ptsd1


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g1<- exp_accmodpa1_clu[!exp_accmodpa1_clu$SNP %in% out_accmodpa1_ptsd1$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g1$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy1<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy1<- ptsd_proxy1 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_accmodpa1_ptsd1[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy1 <- ptsd_proxy1[!ptsd_proxy1$SNP %in% out_accmodpa1_ptsd1$SNP, ]

# combine with outcome data
out_accmodpa1_ptsd1_proxy<- rbind(out_accmodpa1_ptsd1,ptsd_proxy1)


### G2 

out_accmodpa1_ptsd2 <- read_outcome_data(
  snps = exp_accmodpa2_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_accmodpa1_ptsd2)
nrow(exp_accmodpa2_clu)
# if nrow out = exp, run code on line 894 and then move to part 3 (data harmonization, line 952)
# if nrow out != exp -> run code from line 896

out_accmodpa1_ptsd2_proxy <- out_accmodpa1_ptsd2


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g2<- exp_accmodpa2_clu[!exp_accmodpa2_clu$SNP %in% out_accmodpa1_ptsd2$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g2$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy2<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy2<- ptsd_proxy2 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_accmodpa1_ptsd2[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy2 <- ptsd_proxy2[!ptsd_proxy2$SNP %in% out_accmodpa1_ptsd2$SNP, ]

# combine with outcome data
out_accmodpa1_ptsd2_proxy<- rbind(out_accmodpa1_ptsd2,ptsd_proxy2)


### 3. Harmonise data 

dat_accmodpa_ptsd1 <- harmonise_data(
  exposure_dat = exp_accmodpa1_clu, 
  outcome_dat = out_accmodpa1_ptsd1_proxy, 
  action =2)

dat_accmodpa_ptsd2 <- harmonise_data(
  exposure_dat = exp_accmodpa2_clu, 
  outcome_dat = out_accmodpa1_ptsd2_proxy, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets

## Add column for total sample size of outcome GWAS  
# Note: edit sample size accordingly 
dat_accmodpa_ptsd1$samplesize.outcome<- XXX
dat_accmodpa_ptsd2$samplesize.outcome<- XXX

## Pruning 
dat_accmodpa_ptsd1<-power_prune(dat_accmodpa_ptsd1,method=1,dist.outcome="binary")
dat_accmodpa_ptsd2<-power_prune(dat_accmodpa_ptsd2,method=1,dist.outcome="binary")

nrow(dat_accmodpa_ptsd1)
nrow(dat_accmodpa_ptsd2)

## Edit names of outcome and exposure columns 
dat_accmodpa_ptsd1$exposure[dat_accmodpa_ptsd1$exposure == "exposure"] <- "Accelerometer-based moderate PA (G1)"
dat_accmodpa_ptsd1$outcome[dat_accmodpa_ptsd1$outcome == "outcome"] <- "PTSD"

dat_accmodpa_ptsd2$exposure[dat_accmodpa_ptsd2$exposure =="exposure"] <- "Accelerometer-based moderate PA (G2)"
dat_accmodpa_ptsd2$outcome[dat_accmodpa_ptsd2$outcome == "outcome"] <- "PTSD"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_accmodpa_ptsd1<-mr(dat_accmodpa_ptsd1, #method_list=c("mr_ivw_mre","mr_egger_regression",
                                                                     #"mr_weighted_median", 
                                                                     #"mr_weighted_mode"
                           ) # note: only one method available for G1 due to low number of SNPs
results_accmodpa_ptsd1

results_accmodpa_ptsd2<-mr(dat_accmodpa_ptsd2,  method_list=c("mr_ivw_mre", "mr_egger_regression",
                                                                      "mr_weighted_median", 
                                                                      "mr_weighted_mode"))
results_accmodpa_ptsd2

## raps 

# g1
raps<- mr.raps(dat_accmodpa_ptsd1$beta.exposure, dat_accmodpa_ptsd1$beta.outcome, dat_accmodpa_ptsd1$se.exposure, 
               dat_accmodpa_ptsd1$se.outcome)

raps_df<- data.frame(id.exposure = results_accmodpa_ptsd1[1,1],
                     id.outcome = results_accmodpa_ptsd1[1,2], 
                     outcome = results_accmodpa_ptsd1[1,3],
                     exposure = results_accmodpa_ptsd1[1,4], 
                     method = "raps", 
                     nsnp = results_accmodpa_ptsd1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_accmodpa_ptsd1<-rbind(results_accmodpa_ptsd1,raps_df )
results_accmodpa_ptsd1$analysis = "accmodpa_ptsd1"
results_accmodpa_ptsd1 <- generate_odds_ratios(results_accmodpa_ptsd1)

#g2
raps2<- mr.raps(dat_accmodpa_ptsd2$beta.exposure, dat_accmodpa_ptsd2$beta.outcome, dat_accmodpa_ptsd2$se.exposure, 
                dat_accmodpa_ptsd2$se.outcome)

raps_df2<- data.frame(id.exposure = results_accmodpa_ptsd2[1,1],
                      id.outcome = results_accmodpa_ptsd2[1,2], 
                      outcome = results_accmodpa_ptsd2[1,3],
                      exposure = results_accmodpa_ptsd2[1,4], 
                      method = "raps", 
                      nsnp = results_accmodpa_ptsd2[1,6],
                      b = raps[["beta.hat"]], 
                      se = raps[["beta.se"]], 
                      pval = raps[["beta.p.value"]])

results_accmodpa_ptsd2<-rbind(results_accmodpa_ptsd2,raps_df2 )
results_accmodpa_ptsd2$analysis = "accmodpa_ptsd2"
results_accmodpa_ptsd2 <- generate_odds_ratios(results_accmodpa_ptsd2)

# heterogeneity stats (# note: method not available for G1 due to low number of SNPs)

#hetero_accmodpa_ptsd1<- mr_heterogeneity(dat_accmodpa_ptsd1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
#hetero_accmodpa_ptsd1$analysis<- "accmodpa_ptsd1"

hetero_accmodpa_ptsd2<- mr_heterogeneity(dat_accmodpa_ptsd2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_accmodpa_ptsd2$analysis<- "accmodpa_ptsd2"

# check intercept of Egger regression (# note: method not available for G1 due to low number of SNPs)

#intercept_accmodpa_ptsd1<- mr_pleiotropy_test(dat_accmodpa_ptsd1)
#intercept_accmodpa_ptsd1$analysis<- "accmodpa_ptsd1"

intercept_accmodpa_ptsd2<- mr_pleiotropy_test(dat_accmodpa_ptsd2)
intercept_accmodpa_ptsd2$analysis<- "accmodpa_ptsd2"


# steiger directionality test 

steiger_accmodpa_ptsd1 <- directionality_test(dat_accmodpa_ptsd1)
steiger_accmodpa_ptsd1$analysis <- "accmodpa_ptsd1"

steiger_accmodpa_ptsd2 <- directionality_test(dat_accmodpa_ptsd2)
steiger_accmodpa_ptsd2$analysis <- "accmodpa_ptsd2"

# MR presso 

#g1  (# note: method not available for G1 due to low number of SNPs)
#presso<-mr_presso(BetaOutcome = "beta.outcome", 
#                  BetaExposure = "beta.exposure", 
#                  SdOutcome = "se.outcome", 
#                  SdExposure = "se.exposure", 
#                  OUTLIERtest = TRUE,
#                  DISTORTIONtest = TRUE, 
#                  data = dat_accmodpa_ptsd1, 
#                  NbDistribution = 1000,  
#                  SignifThreshold = 0.05)
#
#presso[["Main MR results"]]
#presso$`MR-PRESSO results`$`Global Test`$RSSobs
#presso$`MR-PRESSO results`$`Global Test`$Pvalue

#presso_accmodpa_ptsd1<- presso[["Main MR results"]]
#presso_accmodpa_ptsd1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
#                                        presso$`MR-PRESSO results`$`Global Test`$Pvalue)
#presso_accmodpa_ptsd1$analysis<- "accmodpa_ptsd1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_accmodpa_ptsd2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_accmodpa_ptsd2<- presso2[["Main MR results"]]
presso_accmodpa_ptsd2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                        presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_accmodpa_ptsd2$analysis<- "accmodpa_ptsd2"

### 6.Graphs 

#Scatterplot

#g1 (# note: method not available for G1 due to low number of SNPs)
#res <- mr(dat_accmodpa_ptsd1#, method_list=c("mr_ivw_mre", 
#                                                "mr_egger_regression",
#                                                "mr_weighted_median", 
#                                                "mr_weighted_mode"))
#p1g1<-mr_scatter_plot(res, dat_accmodpa_ptsd1)

#g2
res2 <- mr(dat_accmodpa_ptsd2, method_list=c("mr_ivw_mre",
                                                 "mr_egger_regression",
                                                 "mr_weighted_median", 
                                                 "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_accmodpa_ptsd2)


# combine 

plot_accmodpa_ptsd2 <- ggarrange( p1g2[[paste(exp_accmodpa2[1,11],".",out_accmodpa1_ptsd2[1,11], sep="")]],
                                  labels = c("(a) Scatter plot"),
                                  hjust = 0,
                                  ncol = 1, nrow = 1,
                                  common.legend =F)
plot_accmodpa_ptsd2
ggsave("plot_accmodpa_ptsd2.jpeg", device = "jpeg",dpi = 300, width =7, height =7, limitsize = F)


##### (d) Accelerometer-based walking

### 1. Prepare exposure data 

# load excel file 
accwalk<-read.csv("walking.csv")
names(accwalk)
head(accwalk)

# create df and change column names 
accwalk_df <- data.frame(
  SNP = accwalk$SNP,
  beta = accwalk$BETA,
  se = accwalk$SE,
  effect_allele = accwalk$ALLELE1,
  other_allele = accwalk$ALLELE0,
  eaf = accwalk$A1FREQ, 
  pval = accwalk$P_BOLT_LMM_INF)

# first exposure instrument (G1)
exp_accwalk1 <- format_data(accwalk_df, type="exposure")
exp_accwalk1<- exp_accwalk1[exp_accwalk1$pval.exposure <=5e-8, ]
summary(exp_accwalk1$pval.exposure)

# second genetic instrument (G2)
exp_accwalk2 <- format_data(accwalk_df, type="exposure")
exp_accwalk2<- exp_accwalk2[exp_accwalk2$pval.exposure <=10e-6, ]
summary(exp_accwalk2$pval.exposure)

# clumping 
exp_accwalk1_clu<- clump_data(exp_accwalk1)
exp_accwalk2_clu<- clump_data(exp_accwalk2)

# add column for sample size 
exp_accwalk1_clu$samplesize.exposure<- 91105
exp_accwalk2_clu$samplesize.exposure<- 91105

### 2. Prepare outcome data 

## G1

out_accwalk1_ptsd1 <- read_outcome_data(
  snps = exp_accwalk1_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_accwalk1_ptsd1)
nrow(exp_accwalk1_clu)
# if nrow out = exp, run code on line 1191 and then move to G2 (line 1250)
# if nrow out != exp -> run code from line 1194

out_accwalk1_ptsd1_proxy <- out_accwalk1_ptsd1


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g1<- exp_accwalk1_clu[!exp_accwalk1_clu$SNP %in% out_accwalk1_ptsd1$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g1$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy1<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy1<- ptsd_proxy1 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_accwalk1_ptsd1[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy1 <- ptsd_proxy1[!ptsd_proxy1$SNP %in% out_accwalk1_ptsd1$SNP, ]

# combine with outcome data
out_accwalk1_ptsd1_proxy<- rbind(out_accwalk1_ptsd1,ptsd_proxy1)


### G2 

out_accwalk1_ptsd2 <- read_outcome_data(
  snps = exp_accwalk2_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_accwalk1_ptsd2)
nrow(exp_accwalk2_clu)
# if nrow out = exp, run code on line 1270 and then move to part 3 (data harmonization, line 1329)
# if nrow out != exp -> run code from line 1275

out_accwalk1_ptsd2_proxy <- out_accwalk1_ptsd2


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g2<- exp_accwalk2_clu[!exp_accwalk2_clu$SNP %in% out_accwalk1_ptsd2$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g2$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy2<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy2<- ptsd_proxy2 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_accwalk1_ptsd2[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy2 <- ptsd_proxy2[!ptsd_proxy2$SNP %in% out_accwalk1_ptsd2$SNP, ]

# combine with outcome data
out_accwalk1_ptsd2_proxy<- rbind(out_accwalk1_ptsd2,ptsd_proxy2)


### 3. Harmonise data 

dat_accwalk_ptsd1 <- harmonise_data(
  exposure_dat = exp_accwalk1_clu, 
  outcome_dat = out_accwalk1_ptsd1_proxy, 
  action =2)

dat_accwalk_ptsd2 <- harmonise_data(
  exposure_dat = exp_accwalk2_clu, 
  outcome_dat = out_accwalk1_ptsd2_proxy, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets

## Add column for total sample size of outcome GWAS  
# Note: edit sample size accordingly 
dat_accwalk_ptsd1$samplesize.outcome<- XXX
dat_accwalk_ptsd2$samplesize.outcome<- XXX

## Pruning 
dat_accwalk_ptsd1<-power_prune(dat_accwalk_ptsd1,method=1,dist.outcome="binary")
dat_accwalk_ptsd2<-power_prune(dat_accwalk_ptsd2,method=1,dist.outcome="binary")

nrow(dat_accwalk_ptsd1)
nrow(dat_accwalk_ptsd2)

## Edit names of outcome and exposure columns 
dat_accwalk_ptsd1$exposure[dat_accwalk_ptsd1$exposure == "exposure"] <- "Accelerometer-based walking (G1)"
dat_accwalk_ptsd1$outcome[dat_accwalk_ptsd1$outcome == "outcome"] <- "PTSD"

dat_accwalk_ptsd2$exposure[dat_accwalk_ptsd2$exposure =="exposure"] <- "Accelerometer-based walking (G2)"
dat_accwalk_ptsd2$outcome[dat_accwalk_ptsd2$outcome == "outcome"] <- "PTSD"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_accwalk_ptsd1<-mr(dat_accwalk_ptsd1 #, method_list=c("mr_ivw_mre","mr_egger_regression",
                                                                    # "mr_weighted_median", 
                                                                     #"mr_weighted_mode")
)
results_accwalk_ptsd1

results_accwalk_ptsd2<-mr(dat_accwalk_ptsd2,  method_list=c("mr_ivw_mre", "mr_egger_regression",
                                                                      "mr_weighted_median", 
                                                                      "mr_weighted_mode"))
results_accwalk_ptsd2

## raps 

# g1
raps<- mr.raps(dat_accwalk_ptsd1$beta.exposure, dat_accwalk_ptsd1$beta.outcome, dat_accwalk_ptsd1$se.exposure, 
               dat_accwalk_ptsd1$se.outcome)

raps_df<- data.frame(id.exposure = results_accwalk_ptsd1[1,1],
                     id.outcome = results_accwalk_ptsd1[1,2], 
                     outcome = results_accwalk_ptsd1[1,3],
                     exposure = results_accwalk_ptsd1[1,4], 
                     method = "raps", 
                     nsnp = results_accwalk_ptsd1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_accwalk_ptsd1<-rbind(results_accwalk_ptsd1,raps_df )
results_accwalk_ptsd1$analysis = "accwalk_ptsd1"
results_accwalk_ptsd1 <- generate_odds_ratios(results_accwalk_ptsd1)

#g2
raps2<- mr.raps(dat_accwalk_ptsd2$beta.exposure, dat_accwalk_ptsd2$beta.outcome, dat_accwalk_ptsd2$se.exposure, 
                dat_accwalk_ptsd2$se.outcome)

raps_df2<- data.frame(id.exposure = results_accwalk_ptsd2[1,1],
                      id.outcome = results_accwalk_ptsd2[1,2], 
                      outcome = results_accwalk_ptsd2[1,3],
                      exposure = results_accwalk_ptsd2[1,4], 
                      method = "raps", 
                      nsnp = results_accwalk_ptsd2[1,6],
                      b = raps[["beta.hat"]], 
                      se = raps[["beta.se"]], 
                      pval = raps[["beta.p.value"]])

results_accwalk_ptsd2<-rbind(results_accwalk_ptsd2,raps_df2 )
results_accwalk_ptsd2$analysis = "accwalk_ptsd2"
results_accwalk_ptsd2 <- generate_odds_ratios(results_accwalk_ptsd2)

# heterogeneity stats (note: method not available for G1 due to low number of SNPs)

#hetero_accwalk_ptsd1<- mr_heterogeneity(dat_accwalk_ptsd1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
#hetero_accwalk_ptsd1$analysis<- "accwalk_ptsd1"

hetero_accwalk_ptsd2<- mr_heterogeneity(dat_accwalk_ptsd2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_accwalk_ptsd2$analysis<- "accwalk_ptsd2"

# check intercept of Egger regression (note: method not available for G1 due to low number of SNPs)

#intercept_accwalk_ptsd1<- mr_pleiotropy_test(dat_accwalk_ptsd1)
#intercept_accwalk_ptsd1$analysis<- "accwalk_ptsd1"

intercept_accwalk_ptsd2<- mr_pleiotropy_test(dat_accwalk_ptsd2)
intercept_accwalk_ptsd2$analysis<- "accwalk_ptsd2"


# steiger directionality test 

steiger_accwalk_ptsd1 <- directionality_test(dat_accwalk_ptsd1)
steiger_accwalk_ptsd1$analysis <- "accwalk_ptsd1"

steiger_accwalk_ptsd2 <- directionality_test(dat_accwalk_ptsd2)
steiger_accwalk_ptsd2$analysis <- "accwalk_ptsd2"

# MR presso 

#g1 ( note: method not available for G1 due to low number of SNPs)
#presso<-mr_presso(BetaOutcome = "beta.outcome", 
#                  BetaExposure = "beta.exposure", 
#                  SdOutcome = "se.outcome", 
#                  SdExposure = "se.exposure", 
#                  OUTLIERtest = TRUE,
#                  DISTORTIONtest = TRUE, 
#                  data = dat_accwalk_ptsd1, 
#                  NbDistribution = 1000,  
#                  SignifThreshold = 0.05)
#
#presso[["Main MR results"]]
#presso$`MR-PRESSO results`$`Global Test`$RSSobs
#presso$`MR-PRESSO results`$`Global Test`$Pvalue

#presso_accwalk_ptsd1<- presso[["Main MR results"]]
#presso_accwalk_ptsd1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
#                                        presso$`MR-PRESSO results`$`Global Test`$Pvalue)
#presso_accwalk_ptsd1$analysis<- "accwalk_ptsd1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_accwalk_ptsd2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_accwalk_ptsd2<- presso2[["Main MR results"]]
presso_accwalk_ptsd2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                        presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_accwalk_ptsd2$analysis<- "accwalk_ptsd2"

### 6.Graphs 

#Scatterplot

#g1
#res <- mr(dat_accwalk_ptsd1, method_list=c("mr_ivw_mre",
#                                                "mr_egger_regression",
#                                                "mr_weighted_median", 
#                                                "mr_weighted_mode"))
#p1g1<-mr_scatter_plot(res, dat_accwalk_ptsd1)

#g2
res2 <- mr(dat_accwalk_ptsd2, method_list=c("mr_ivw_mre",
                                                 "mr_egger_regression",
                                                 "mr_weighted_median", 
                                                 "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_accwalk_ptsd2)


# combine 

plot_accwalk_ptsd2 <- ggarrange( p1g2[[paste(exp_accwalk2[1,11],".",out_accwalk1_ptsd2[1,11], sep="")]],
                                  labels = c("(a) Scatter plot"),
                                  hjust = 0,
                                  ncol = 1, nrow = 1,
                                  common.legend =F)
plot_accwalk_ptsd2
ggsave("plot_accwalk_ptsd2.jpeg", device = "jpeg",dpi = 300, width =7, height =7, limitsize = F)

##### (e) Accelerometer-based sedentary behavior 

### 1. Prepare exposure data 

# load excel file 
accsed<-read.csv("sedentarybehaviour.csv")
names(accsed)
head(accsed)

# create df and change column names 
accsed_df <- data.frame(
  SNP = accsed$SNP,
  beta = accsed$BETA,
  se = accsed$SE,
  effect_allele = accsed$ALLELE1,
  other_allele = accsed$ALLELE0,
  eaf = accsed$A1FREQ, 
  pval = accsed$P_BOLT_LMM_INF)

# first exposure instrument (G1)
exp_accsed1 <- format_data(accsed_df, type="exposure")
exp_accsed1<- exp_accsed1[exp_accsed1$pval.exposure <=5e-8, ]
summary(exp_accsed1$pval.exposure)

# second genetic instrument (G2)
exp_accsed2 <- format_data(accsed_df, type="exposure")
exp_accsed2<- exp_accsed2[exp_accsed2$pval.exposure <=10e-6, ]
summary(exp_accsed2$pval.exposure)

# clumping 
exp_accsed1_clu<- clump_data(exp_accsed1)
exp_accsed2_clu<- clump_data(exp_accsed2)

# add column for sample size 
exp_accsed1_clu$samplesize.exposure<- 91105
exp_accsed2_clu$samplesize.exposure<- 91105

### 2. Prepare outcome data 

## G1

out_accsed1_ptsd1 <- read_outcome_data(
  snps = exp_accsed1_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_accsed1_ptsd1)
nrow(exp_accsed1_clu)
# if nrow out = exp, run code on line 1567 and then move to G2 (line 1626)
# if nrow out != exp -> run code from line 1570

out_accsed1_ptsd1_proxy <- out_accsed1_ptsd1


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g1<- exp_accsed1_clu[!exp_accsed1_clu$SNP %in% out_accsed1_ptsd1$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g1$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 
 
ptsd_proxy1<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy1<- ptsd_proxy1 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_accsed1_ptsd1[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy1 <- ptsd_proxy1[!ptsd_proxy1$SNP %in% out_accsed1_ptsd1$SNP, ]

# combine with outcome data
out_accsed1_ptsd1_proxy<- rbind(out_accsed1_ptsd1,ptsd_proxy1)


### G2 

out_accsed1_ptsd2 <- read_outcome_data(
  snps = exp_accsed2_clu$SNP,
  filename = "ptsd.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(out_accsed1_ptsd2)
nrow(exp_accsed2_clu)
# if nrow out = exp, run code on line 1646 and then move to part 3 (data harmonization, line 1705)
# if nrow out != exp -> run code from line 1652

out_accsed1_ptsd2_proxy <- out_accsed1_ptsd2


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g2<- exp_accsed2_clu[!exp_accsed2_clu$SNP %in% out_accsed1_ptsd2$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_g2$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "abf3700533f3")
  df = rbind(df, trial)  
}

# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

ptsd_proxy2<- filter(ptsd, SNP %in% df_high$RS_Number) 

# rename column names 

ptsd_proxy2<- ptsd_proxy2 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         eaf.outcome = eaf,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_accsed1_ptsd2[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         eaf.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         data_source.outcome)

# remove duplicates    
ptsd_proxy2 <- ptsd_proxy2[!ptsd_proxy2$SNP %in% out_accsed1_ptsd2$SNP, ]

# combine with outcome data
out_accsed1_ptsd2_proxy<- rbind(out_accsed1_ptsd2,ptsd_proxy2)


### 3. Harmonise data 

dat_accsed_ptsd1 <- harmonise_data(
  exposure_dat = exp_accsed1_clu, 
  outcome_dat = out_accsed1_ptsd1_proxy, 
  action =2)

dat_accsed_ptsd2 <- harmonise_data(
  exposure_dat = exp_accsed2_clu, 
  outcome_dat = out_accsed1_ptsd2_proxy, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets

## Add column for total sample size of outcome GWAS  
# Note: edit sample size accordingly 
dat_accsed_ptsd1$samplesize.outcome<- XXX
dat_accsed_ptsd2$samplesize.outcome<- XXX

## Pruning 
dat_accsed_ptsd1<-power_prune(dat_accsed_ptsd1,method=1,dist.outcome="binary")
dat_accsed_ptsd2<-power_prune(dat_accsed_ptsd2,method=1,dist.outcome="binary")

nrow(dat_accsed_ptsd1)
nrow(dat_accsed_ptsd2)

## Edit names of outcome and exposure columns 
dat_accsed_ptsd1$exposure[dat_accsed_ptsd1$exposure == "exposure"] <- "Accelerometer-based sedentary behaviour (G1)"
dat_accsed_ptsd1$outcome[dat_accsed_ptsd1$outcome == "outcome"] <- "PTSD"

dat_accsed_ptsd2$exposure[dat_accsed_ptsd2$exposure == "exposure"] <- "Accelerometer-based sedentary behaviour (G2)"
dat_accsed_ptsd2$outcome[dat_accsed_ptsd2$outcome == "outcome"] <- "PTSD"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_accsed_ptsd1<-mr(dat_accsed_ptsd1, method_list=c("mr_ivw_mre",
                                                                     "mr_egger_regression",
                                                                     "mr_weighted_median", 
                                                                     "mr_weighted_mode"))
results_accsed_ptsd1

results_accsed_ptsd2<-mr(dat_accsed_ptsd2,  method_list=c("mr_ivw_mre", "mr_egger_regression",
                                                                      "mr_weighted_median", 
                                                                      "mr_weighted_mode"))
results_accsed_ptsd2

## raps 

# g1
raps<- mr.raps(dat_accsed_ptsd1$beta.exposure, dat_accsed_ptsd1$beta.outcome, dat_accsed_ptsd1$se.exposure, 
               dat_accsed_ptsd1$se.outcome)

raps_df<- data.frame(id.exposure = results_accsed_ptsd1[1,1],
                     id.outcome = results_accsed_ptsd1[1,2], 
                     outcome = results_accsed_ptsd1[1,3],
                     exposure = results_accsed_ptsd1[1,4], 
                     method = "raps", 
                     nsnp = results_accsed_ptsd1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_accsed_ptsd1<-rbind(results_accsed_ptsd1,raps_df )
results_accsed_ptsd1$analysis = "accsed_ptsd1"
results_accsed_ptsd1 <- generate_odds_ratios(results_accsed_ptsd1)

#g2
raps2<- mr.raps(dat_accsed_ptsd2$beta.exposure, dat_accsed_ptsd2$beta.outcome, dat_accsed_ptsd2$se.exposure, 
                dat_accsed_ptsd2$se.outcome)

raps_df2<- data.frame(id.exposure = results_accsed_ptsd2[1,1],
                      id.outcome = results_accsed_ptsd2[1,2], 
                      outcome = results_accsed_ptsd2[1,3],
                      exposure = results_accsed_ptsd2[1,4], 
                      method = "raps", 
                      nsnp = results_accsed_ptsd2[1,6],
                      b = raps[["beta.hat"]], 
                      se = raps[["beta.se"]], 
                      pval = raps[["beta.p.value"]])

results_accsed_ptsd2<-rbind(results_accsed_ptsd2,raps_df2 )
results_accsed_ptsd2$analysis = "accsed_ptsd2"
results_accsed_ptsd2 <- generate_odds_ratios(results_accsed_ptsd2)

# heterogeneity stats 

hetero_accsed_ptsd1<- mr_heterogeneity(dat_accsed_ptsd1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_accsed_ptsd1$analysis<- "accsed_ptsd1"

hetero_accsed_ptsd2<- mr_heterogeneity(dat_accsed_ptsd2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_accsed_ptsd2$analysis<- "accsed_ptsd2"

# check intercept of Egger regression 

intercept_accsed_ptsd1<- mr_pleiotropy_test(dat_accsed_ptsd1)
intercept_accsed_ptsd1$analysis<- "accsed_ptsd1"

intercept_accsed_ptsd2<- mr_pleiotropy_test(dat_accsed_ptsd2)
intercept_accsed_ptsd2$analysis<- "accsed_ptsd2"


# steiger directionality test 

steiger_accsed_ptsd1 <- directionality_test(dat_accsed_ptsd1)
steiger_accsed_ptsd1$analysis <- "accsed_ptsd1"

steiger_accsed_ptsd2 <- directionality_test(dat_accsed_ptsd2)
steiger_accsed_ptsd2$analysis <- "accsed_ptsd2"

# MR presso 

#g1
presso<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                  SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_accsed_ptsd1, 
                  NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso[["Main MR results"]]
presso$`MR-PRESSO results`$`Global Test`$RSSobs
presso$`MR-PRESSO results`$`Global Test`$Pvalue

presso_accsed_ptsd1<- presso[["Main MR results"]]
presso_accsed_ptsd1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                        presso$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_accsed_ptsd1$analysis<- "accsed_ptsd1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_accsed_ptsd2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_accsed_ptsd2<- presso2[["Main MR results"]]
presso_accsed_ptsd2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                        presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_accsed_ptsd2$analysis<- "accsed_ptsd2"

### 6.Graphs 

#Scatterplot

#g1
res <- mr(dat_accsed_ptsd1, method_list=c("mr_ivw_mre",
                                                "mr_egger_regression",
                                                "mr_weighted_median", 
                                                "mr_weighted_mode"))
p1g1<-mr_scatter_plot(res, dat_accsed_ptsd1)

#g2
res2 <- mr(dat_accsed_ptsd2, method_list=c("mr_ivw_mre",
                                                 "mr_egger_regression",
                                                 "mr_weighted_median", 
                                                 "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_accsed_ptsd2)


# combine 

plot_accsed_ptsd12 <- ggarrange(p1g1[[paste(exp_accsed1[1,11],".",out_accsed1_ptsd1[1,11],sep="")]],
                                p1g2[[paste(exp_accsed2[1,11],".",out_accsed1_ptsd2[1,11], sep="")]],
                                labels = c("(a) Scatter plots", ""),
                                hjust = -0.05,
                                ncol = 2, nrow = 1,
                                common.legend =F)
plot_accsed_ptsd12
ggsave("plot_accsed_ptsd12.jpeg", device = "jpeg",dpi = 300, width =15, height =7, limitsize = F)


### Dataframes with results of MR and sensitivity analyses for all exposures 

results_pa_ptsd<- rbind(results_modvigpaself_ptsd1,
                       results_modvigpaself_ptsd2,
                       results_accpamean_ptsd1,
                       results_accpamean_ptsd2,
                       results_accmodpa_ptsd1, 
                       results_accmodpa_ptsd2,
                       results_accwalk_ptsd1,
                       results_accwalk_ptsd2,
                       results_accsed_ptsd1,
                       results_accsed_ptsd2)

hetero_pa_ptsd<-rbind(hetero_modvigpaself_ptsd1,
                     hetero_modvigpaself_ptsd2,
                     hetero_accpamean_ptsd1,
                     hetero_accpamean_ptsd2,
                     hetero_accmodpa_ptsd2,
                     hetero_accwalk_ptsd2,
                     hetero_accsed_ptsd1,
                     hetero_accsed_ptsd2)

intercept_pa_ptsd<-rbind(intercept_modvigpaself_ptsd1,
                        intercept_modvigpaself_ptsd2,
                        intercept_accpamean_ptsd1,
                        intercept_accpamean_ptsd2,
                        intercept_accmodpa_ptsd2,
                        intercept_accwalk_ptsd2,
                        intercept_accsed_ptsd1,
                        intercept_accsed_ptsd2)

steiger_pa_ptsd<-rbind(steiger_modvigpaself_ptsd1,
                      steiger_modvigpaself_ptsd2,
                      steiger_accpamean_ptsd1,
                      steiger_accpamean_ptsd2,
                      steiger_accmodpa_ptsd1, 
                      steiger_accmodpa_ptsd2,
                      steiger_accwalk_ptsd1,
                      steiger_accwalk_ptsd2,
                      steiger_accsed_ptsd1,
                      steiger_accsed_ptsd2)

presso_pa_ptsd<- rbind(presso_modvigpaself_ptsd1,
                      presso_modvigpaself_ptsd2,
                      presso_accpamean_ptsd1,
                      presso_accpamean_ptsd2,
                      presso_accmodpa_ptsd2,
                      presso_accwalk_ptsd2,
                      presso_accsed_ptsd1,
                      presso_accsed_ptsd2)


write_xlsx(list(results_pa_ptsd = results_pa_ptsd,
                hetero_pa_ptsd = hetero_pa_ptsd,
                intercept_pa_ptsd = intercept_pa_ptsd,
                steiger_pa_ptsd = steiger_pa_ptsd,
                presso_pa_ptsd = presso_pa_ptsd), "results_pa_ptsd_all.xlsx")


save.image("~/Documents/Study1_GWAS/PA_ptsd.RData")

