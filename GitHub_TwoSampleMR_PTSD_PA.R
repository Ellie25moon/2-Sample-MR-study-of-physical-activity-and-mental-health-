####################################################################################
######################### PTSD -> Physical Activity ###############################
####################################################################################

setwd("~/Documents/Study1_GWAS")

# load packages 
library(TwoSampleMR)
library(LDlinkR)
library (dplyr)
library(mr.raps)
library(MRPRESSO)
library(ggpubr)
library(writexl)

##### (a) Self-reported moderate-vigorous activity 

### 1. Prepare exposure data 

# upload excel file 
ptsd_df<-read.csv("ptsd.csv")
names(ptsd_df)

# add column for total sample size of PTSD GWAS
# Note: edit sample size accordingly 
ptsd_df$samplesize.exposure<- XXX

# first genetic instrument (G1)
exp_ptsd1 <- format_data(ptsd_df, type="exposure",samplesize_col = "samplesize.exposure")
exp_ptsd1<- exp_ptsd1[exp_ptsd1$pval.exposure <=5e-8, ]
nrow(exp_ptsd1)

# second genetic instrument (G2)
exp_ptsd2 <- format_data(ptsd_df, type="exposure", samplesize_col = "samplesize.exposure")
exp_ptsd2<- exp_ptsd2[exp_ptsd2$pval.exposure <=10e-6, ]
nrow(exp_ptsd2)

# clumping 
exp_ptsd1_clu<- clump_data(exp_ptsd1)
exp_ptsd2_clu<- clump_data(exp_ptsd2)

### 2. Prepare outcome data 
out_modvigpaself1 <- extract_outcome_data(
  snps = exp_ptsd1_clu$SNP,
  outcomes = 'ebi-a-GCST006097'
)

nrow(out_modvigpaself1)

out_modvigpaself2 <- extract_outcome_data(
  snps = exp_ptsd2_clu$SNP,
  outcomes = 'ebi-a-GCST006097'
)

nrow(out_modvigpaself2)

# Note: proxy SNPs are added by default when using internal outcome data  

### 3. Harmonise data 
dat_ptsd_modvigpaself1 <- harmonise_data(
  exposure_dat = exp_ptsd1_clu, 
  outcome_dat = out_modvigpaself1, 
  action =2)

dat_ptsd_modvigpaself2 <- harmonise_data(
  exposure_dat = exp_ptsd2_clu, 
  outcome_dat = out_modvigpaself2, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets

dat_ptsd_modvigpaself1<-power_prune(dat_ptsd_modvigpaself1,method=2,dist.outcome="continuous")

dat_ptsd_modvigpaself2<-power_prune(dat_ptsd_modvigpaself2,method=2,dist.outcome="continuous")

nrow(dat_ptsd_modvigpaself1)
nrow(dat_ptsd_modvigpaself2)

## Edit names of outcome and exposure
dat_ptsd_modvigpaself1$exposure[dat_ptsd_modvigpaself1$exposure == "exposure"] <- "PTSD (G1)"
dat_ptsd_modvigpaself1$outcome[dat_ptsd_modvigpaself1$outcome == "Moderate to vigorous physical activity levels || id:ebi-a-GCST006097"] <- "Self-reported moderate-vigorous PA"

dat_ptsd_modvigpaself2$exposure[dat_ptsd_modvigpaself2$exposure == "exposure"] <- "PTSD (G2)"
dat_ptsd_modvigpaself2$outcome[dat_ptsd_modvigpaself2$outcome == "Moderate to vigorous physical activity levels || id:ebi-a-GCST006097"] <- "Self-reported moderate-vigorous PA"


##### 5.MR analyses 

# Ivw, egger, weighted_median, weighted mode 

results_ptsd_modvigpaself1<-mr(dat_ptsd_modvigpaself1, method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                                    "mr_weighted_median", 
                                                                    "mr_weighted_mode"))
results_ptsd_modvigpaself1

results_ptsd_modvigpaself2<-mr(dat_ptsd_modvigpaself2,  method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                                     "mr_weighted_median", 
                                                                     "mr_weighted_mode"))
results_ptsd_modvigpaself2

# raps 

# g1
raps<- mr.raps(dat_ptsd_modvigpaself1$beta.exposure, dat_ptsd_modvigpaself1$beta.outcome, dat_ptsd_modvigpaself1$se.exposure, 
               dat_ptsd_modvigpaself1$se.outcome)

raps_df<- data.frame(id.exposure = results_ptsd_modvigpaself1[1,1],
                     id.outcome = results_ptsd_modvigpaself1[1,2], 
                     outcome = results_ptsd_modvigpaself1[1,3],
                     exposure = results_ptsd_modvigpaself1[1,4], 
                     method = "raps", 
                     nsnp = results_ptsd_modvigpaself1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_ptsd_modvigpaself1<-rbind(results_ptsd_modvigpaself1,raps_df )
results_ptsd_modvigpaself1$analysis = "ptsd_modvigpaself1"
results_ptsd_modvigpaself1 <- generate_odds_ratios(results_ptsd_modvigpaself1)

#g2
raps2<- mr.raps(dat_ptsd_modvigpaself2$beta.exposure, dat_ptsd_modvigpaself2$beta.outcome, dat_ptsd_modvigpaself2$se.exposure, 
                dat_ptsd_modvigpaself2$se.outcome)

raps_df2<- data.frame(id.exposure = results_ptsd_modvigpaself2[1,1],
                      id.outcome = results_ptsd_modvigpaself2[1,2], 
                      outcome = results_ptsd_modvigpaself2[1,3],
                      exposure = results_ptsd_modvigpaself2[1,4], 
                      method = "raps", 
                      nsnp = results_ptsd_modvigpaself2[1,6],
                      b = raps2[["beta.hat"]], 
                      se = raps2[["beta.se"]], 
                      pval = raps2[["beta.p.value"]])

results_ptsd_modvigpaself2<-rbind(results_ptsd_modvigpaself2,raps_df2 )
results_ptsd_modvigpaself2$analysis = "ptsd_modvigpaself2"
results_ptsd_modvigpaself2 <- generate_odds_ratios(results_ptsd_modvigpaself2)


# heterogeneity stats 

hetero_ptsd_modvigpaself1<- mr_heterogeneity(dat_ptsd_modvigpaself1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_modvigpaself1$analysis<- "ptsd_modvigpaself1"

hetero_ptsd_modvigpaself2<- mr_heterogeneity(dat_ptsd_modvigpaself2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_modvigpaself2$analysis<- "ptsd_modvigpaself2"

# check intercept of Egger regression 

intercept_ptsd_modvigpaself1<- mr_pleiotropy_test(dat_ptsd_modvigpaself1)
intercept_ptsd_modvigpaself1$analysis<- "ptsd_modvigpaself1"

intercept_ptsd_modvigpaself2<- mr_pleiotropy_test(dat_ptsd_modvigpaself2)
intercept_ptsd_modvigpaself2$analysis<- "ptsd_modvigpaself2"


# steiger directionality test 

steiger_ptsd_modvigpaself1 <- directionality_test(dat_ptsd_modvigpaself1)
steiger_ptsd_modvigpaself1$analysis <- "ptsd_modvigpaself1"

steiger_ptsd_modvigpaself2 <- directionality_test(dat_ptsd_modvigpaself2)
steiger_ptsd_modvigpaself2$analysis <- "ptsd_modvigpaself2"


# MR presso 

#g1

presso<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                 SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_ptsd_modvigpaself1, 
                 NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso[["Main MR results"]]
presso$`MR-PRESSO results`$`Global Test`$RSSobs
presso$`MR-PRESSO results`$`Global Test`$Pvalue

presso_ptsd_modvigpaself1<- presso[["Main MR results"]]
presso_ptsd_modvigpaself1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                       presso$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_modvigpaself1$analysis<- "ptsd_modvigpaself1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_ptsd_modvigpaself2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_ptsd_modvigpaself2<- presso2[["Main MR results"]]
presso_ptsd_modvigpaself2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                       presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_modvigpaself2$analysis<- "ptsd_modvigpaself2"


### 6. Graphs 

#Scatterplot

#g1

res <- mr(dat_ptsd_modvigpaself1,method_list=c("mr_ivw_mre",
                                              "mr_egger_regression",
                                              "mr_weighted_median", 
                                              "mr_weighted_mode"))

p1g1<-mr_scatter_plot(res, dat_ptsd_modvigpaself1)

#g2
res2 <- mr(dat_ptsd_modvigpaself2, method_list=c("mr_ivw_mre",
                                                "mr_egger_regression",
                                                "mr_weighted_median", 
                                                "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_ptsd_modvigpaself2)


# combine 
plot_ptsd_modvigpaself12 <- ggarrange(p1g1[[paste(exp_ptsd1[1,12],".",out_modvigpaself1[1,12], sep="")]],
                                      p1g2[[paste(exp_ptsd2[1,12],".",out_modvigpaself2[1,12], sep="")]],
                                      labels = c("(a) Scatter plots", ""),
                                      hjust = -0.05,
                                      ncol = 2, nrow = 1,
                                      common.legend =F)
plot_ptsd_modvigpaself12
ggsave("plot_ptsd_modvigpaself12.jpeg", device = "jpeg",dpi = 300, width =15, height =7, limitsize = F)


##### (b) Accelerometer-based physical activity (mean acceleration) 

### 1.Prepare exposure data 
head(ptsd_df)

# first exposure instrument (G1)
head(exp_ptsd1)

# second genetic instrument (G2)
head(exp_ptsd2)

# clumping 
head(exp_ptsd1_clu)
head(exp_ptsd2_clu) 

### 2.Prepare outcome data 
out_accpamean1 <- extract_outcome_data(
  snps = exp_ptsd1_clu$SNP,
  outcomes = 'ebi-a-GCST006099'
)
nrow(out_accpamean1)

out_accpamean2 <- extract_outcome_data(
  snps = exp_ptsd2_clu$SNP,
  outcomes = 'ebi-a-GCST006099'
)
nrow(out_accpamean2)
# Note: proxy SNPs are added by default when using internal outcome data  

### 3.Harmonise data 
dat_ptsd_accpamean1 <- harmonise_data(
  exposure_dat = exp_ptsd1_clu, 
  outcome_dat = out_accpamean1, 
  action =2)

dat_ptsd_accpamean2 <- harmonise_data(
  exposure_dat = exp_ptsd2_clu, 
  outcome_dat = out_accpamean2, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets

dat_ptsd_accpamean1<-power_prune(dat_ptsd_accpamean1,method=2,dist.outcome="continuous")

dat_ptsd_accpamean2<-power_prune(dat_ptsd_accpamean2,method=2,dist.outcome="continuous")


## Edit names of outcome and exposure
dat_ptsd_accpamean1$exposure[dat_ptsd_accpamean1$exposure == "exposure"] <- "PTSD (G1)"
dat_ptsd_accpamean1$outcome[dat_ptsd_accpamean1$outcome == "Accelerometer-based physical activity measurement (average acceleration) || id:ebi-a-GCST006099"] <- "Accelerometer-based PA, mean acceleration"

dat_ptsd_accpamean2$exposure[dat_ptsd_accpamean2$exposure == "exposure"] <- "PTSD (G2)"
dat_ptsd_accpamean2$outcome[dat_ptsd_accpamean2$outcome == "Accelerometer-based physical activity measurement (average acceleration) || id:ebi-a-GCST006099"] <- "Accelerometer-based PA, mean acceleration"

nrow(dat_ptsd_accpamean1)
nrow(dat_ptsd_accpamean2)

##### 5.MR analyses 

# Ivw, egger, weighted_median, weighted mode 

results_ptsd_accpamean1<-mr(dat_ptsd_accpamean1,method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                             "mr_weighted_median", 
                                                             "mr_weighted_mode"))
results_ptsd_accpamean1

results_ptsd_accpamean2<-mr(dat_ptsd_accpamean2,  method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                               "mr_weighted_median", 
                                                               "mr_weighted_mode"))
results_ptsd_accpamean2

# raps 

# g1
raps<- mr.raps(dat_ptsd_accpamean1$beta.exposure, dat_ptsd_accpamean1$beta.outcome, dat_ptsd_accpamean1$se.exposure, 
               dat_ptsd_accpamean1$se.outcome)

raps_df<- data.frame(id.exposure = results_ptsd_accpamean1[1,1],
                     id.outcome = results_ptsd_accpamean1[1,2], 
                     outcome = results_ptsd_accpamean1[1,3],
                     exposure = results_ptsd_accpamean1[1,4], 
                     method = "raps", 
                     nsnp = results_ptsd_accpamean1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_ptsd_accpamean1<-rbind(results_ptsd_accpamean1,raps_df )
results_ptsd_accpamean1$analysis = "ptsd_accpamean1"
results_ptsd_accpamean1 <- generate_odds_ratios(results_ptsd_accpamean1)

#g2
raps2<- mr.raps(dat_ptsd_accpamean2$beta.exposure, dat_ptsd_accpamean2$beta.outcome, dat_ptsd_accpamean2$se.exposure, 
                dat_ptsd_accpamean2$se.outcome)

raps_df2<- data.frame(id.exposure = results_ptsd_accpamean2[1,1],
                      id.outcome = results_ptsd_accpamean2[1,2], 
                      outcome = results_ptsd_accpamean2[1,3],
                      exposure = results_ptsd_accpamean2[1,4], 
                      method = "raps", 
                      nsnp = results_ptsd_accpamean2[1,6],
                      b = raps2[["beta.hat"]], 
                      se = raps2[["beta.se"]], 
                      pval = raps2[["beta.p.value"]])

results_ptsd_accpamean2<-rbind(results_ptsd_accpamean2,raps_df2 )
results_ptsd_accpamean2$analysis = "ptsd_accpamean2"
results_ptsd_accpamean2 <- generate_odds_ratios(results_ptsd_accpamean2)


# heterogeneity stats 

hetero_ptsd_accpamean1<- mr_heterogeneity(dat_ptsd_accpamean1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_accpamean1$analysis<- "ptsd_accpamean1"

hetero_ptsd_accpamean2<- mr_heterogeneity(dat_ptsd_accpamean2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_accpamean2$analysis<- "ptsd_accpamean2"

# check intercept of Egger regression 

intercept_ptsd_accpamean1<- mr_pleiotropy_test(dat_ptsd_accpamean1)
intercept_ptsd_accpamean1$analysis<- "ptsd_accpamean1"

intercept_ptsd_accpamean2<- mr_pleiotropy_test(dat_ptsd_accpamean2)
intercept_ptsd_accpamean2$analysis<- "ptsd_accpamean2"


# steiger directionality test 

steiger_ptsd_accpamean1 <- directionality_test(dat_ptsd_accpamean1)
steiger_ptsd_accpamean1$analysis <- "ptsd_accpamean1"

steiger_ptsd_accpamean2 <- directionality_test(dat_ptsd_accpamean2)
steiger_ptsd_accpamean2$analysis <- "ptsd_accpamean2"


# MR presso 

#g1

presso<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                  SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_ptsd_accpamean1, 
                  NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso[["Main MR results"]]
presso$`MR-PRESSO results`$`Global Test`$RSSobs
presso$`MR-PRESSO results`$`Global Test`$Pvalue

presso_ptsd_accpamean1<- presso[["Main MR results"]]
presso_ptsd_accpamean1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                   presso$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_accpamean1$analysis<- "ptsd_accpamean1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_ptsd_accpamean2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_ptsd_accpamean2<- presso2[["Main MR results"]]
presso_ptsd_accpamean2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                    presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_accpamean2$analysis<- "ptsd_accpamean2"

### 6.Graphs 

#Scatterplot

#g1
res <- mr(dat_ptsd_accpamean1, method_list=c("mr_ivw_mre",
                                            "mr_egger_regression",
                                            "mr_weighted_median", 
                                            "mr_weighted_mode"))

p1g1<-mr_scatter_plot(res, dat_ptsd_accpamean1)

#g2
res2 <- mr(dat_ptsd_accpamean2, method_list=c("mr_ivw_mre",
                                             "mr_egger_regression",
                                             "mr_weighted_median", 
                                             "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_ptsd_accpamean2)


# combine 
plot_ptsd_accpamean12 <- ggarrange(p1g1[[paste(exp_ptsd1[1,12],".",out_accpamean1[1,12], sep="")]],
                                   p1g2[[paste(exp_ptsd2[1,12],".",out_accpamean2[1,12], sep="")]],
                                   labels = c("(a) Scatter plots", ""),
                                   hjust = -0.05,
                                   ncol = 2, nrow = 1,
                                   common.legend =F)
plot_ptsd_accpamean12
ggsave("plot_ptsd_accpamean12.jpeg", device = "jpeg",dpi = 300, width =15, height =7, limitsize = F)

##### (c) Accelerometer-based moderate physical activity 

### 1.Prepare exposure data 
head(ptsd_df)

# first exposure instrument (G1)
head(exp_ptsd1)

# second genetic instrument (G2)
head(exp_ptsd2)

# clumping 
head(exp_ptsd1_clu)
head(exp_ptsd2_clu) 

### 2.Prepare outcome data 

### G1
out_accmodpa1 <- read_outcome_data(
  snps = exp_ptsd1_clu$SNP,
  filename = "moderatepa.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF"
  #units_col = "Units",
  #gene_col = "Gene",
)
nrow(exp_ptsd1_clu)
nrow(out_accmodpa1)
# if nrow out = exp -> run code on line 479 and then move to G2 (line 536)
# if nrow out != exp -> run code from line 482

out_accmodpa1_proxy <- out_accmodpa1


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g1<- exp_ptsd1_clu[!exp_ptsd1_clu$SNP %in% out_accmodpa1$SNP, ]

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
accmodpa_proxy1<- filter(accmodpa, SNP %in% df_high$RS_Number) 

# rename column names 

accmodpa_proxy1<- accmodpa_proxy1 %>%
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
         id.outcome = out_accmodpa1[1,11] ,
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
accmodpa_proxy1 <- accmodpa_proxy1[!accmodpa_proxy1$SNP %in% out_accmodpa1$SNP, ]

# combine with outcome data
out_accmodpa1_proxy<- rbind(out_accmodpa1,accmodpa_proxy1)

### G2
out_accmodpa2 <- read_outcome_data(
  snps = exp_ptsd2_clu$SNP,
  filename = "moderatepa.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF"
  #units_col = "Units",
  #gene_col = "Gene",
)

nrow(exp_ptsd2_clu)
nrow(out_accmodpa2)
# if nrow out = exp, run code on line 557 and then move to Part 3 (data harmonisation, line 615)
# if nrow out != exp -> run code from line 560

out_accmodpa2_proxy <- out_accmodpa2


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g2<- exp_ptsd2_clu[!exp_ptsd2_clu$SNP %in% out_accmodpa2$SNP, ]

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
accmodpa_proxy2<- filter(accmodpa, SNP %in% df_high$RS_Number) 

# rename column names 

accmodpa_proxy2<- accmodpa_proxy2 %>%
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
         id.outcome = out_accmodpa2[1,11] ,
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
accmodpa_proxy2 <- accmodpa_proxy2[!accmodpa_proxy2$SNP %in% out_accmodpa2$SNP, ]

# combine with outcome data
out_accmodpa2_proxy<- rbind(out_accmodpa2,accmodpa_proxy2)


### 3. Harmonise data 
dat_ptsd_accmodpa1 <- harmonise_data(
  exposure_dat = exp_ptsd1_clu, 
  outcome_dat = out_accmodpa1_proxy, 
  action =2)

dat_ptsd_accmodpa2 <- harmonise_data(
  exposure_dat = exp_ptsd2_clu, 
  outcome_dat = out_accmodpa2_proxy, 
  action =2)

### Drop duplicate exposure-outcome summary sets

dat_ptsd_accmodpa1<-power_prune(dat_ptsd_accmodpa1,method=2,dist.outcome="continuous")

dat_ptsd_accmodpa2<-power_prune(dat_ptsd_accmodpa2,method=2,dist.outcome="continuous")


## Edit names of outcome and exposure
dat_ptsd_accmodpa1$exposure[dat_ptsd_accmodpa1$exposure == "exposure"] <- "PTSD (G1)"
dat_ptsd_accmodpa1$outcome[dat_ptsd_accmodpa1$outcome == "outcome"] <- "Accelerometer-based moderate PA"

dat_ptsd_accmodpa2$exposure[dat_ptsd_accmodpa2$exposure == "exposure"] <- "PTSD (G2)"
dat_ptsd_accmodpa2$outcome[dat_ptsd_accmodpa2$outcome == "outcome"] <- "Accelerometer-based moderate PA"

dat_ptsd_accmodpa1$samplesize.outcome<- 91105
dat_ptsd_accmodpa2$samplesize.outcome<- 91105

nrow(dat_ptsd_accmodpa1)
nrow(dat_ptsd_accmodpa2)

##### 5.MR analyses 

# Ivw, egger, weighted_median, weighted mode 

results_ptsd_accmodpa1<-mr(dat_ptsd_accmodpa1,method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                           "mr_weighted_median", 
                                                           "mr_weighted_mode"))
results_ptsd_accmodpa1

results_ptsd_accmodpa2<-mr(dat_ptsd_accmodpa2,  method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                             "mr_weighted_median", 
                                                             "mr_weighted_mode"))
results_ptsd_accmodpa2

# raps 

# g1
raps<- mr.raps(dat_ptsd_accmodpa1$beta.exposure, dat_ptsd_accmodpa1$beta.outcome, dat_ptsd_accmodpa1$se.exposure, 
               dat_ptsd_accmodpa1$se.outcome)

raps_df<- data.frame(id.exposure = results_ptsd_accmodpa1[1,1],
                     id.outcome = results_ptsd_accmodpa1[1,2], 
                     outcome = results_ptsd_accmodpa1[1,3],
                     exposure = results_ptsd_accmodpa1[1,4], 
                     method = "raps", 
                     nsnp = results_ptsd_accmodpa1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_ptsd_accmodpa1<-rbind(results_ptsd_accmodpa1,raps_df )
results_ptsd_accmodpa1$analysis = "ptsd_accmodpa1"
results_ptsd_accmodpa1 <- generate_odds_ratios(results_ptsd_accmodpa1)

#g2
raps2<- mr.raps(dat_ptsd_accmodpa2$beta.exposure, dat_ptsd_accmodpa2$beta.outcome, dat_ptsd_accmodpa2$se.exposure, 
                dat_ptsd_accmodpa2$se.outcome)

raps_df2<- data.frame(id.exposure = results_ptsd_accmodpa2[1,1],
                      id.outcome = results_ptsd_accmodpa2[1,2], 
                      outcome = results_ptsd_accmodpa2[1,3],
                      exposure = results_ptsd_accmodpa2[1,4], 
                      method = "raps", 
                      nsnp = results_ptsd_accmodpa2[1,6],
                      b = raps2[["beta.hat"]], 
                      se = raps2[["beta.se"]], 
                      pval = raps2[["beta.p.value"]])

results_ptsd_accmodpa2<-rbind(results_ptsd_accmodpa2,raps_df2 )
results_ptsd_accmodpa2$analysis = "ptsd_accmodpa2"
results_ptsd_accmodpa2 <- generate_odds_ratios(results_ptsd_accmodpa2)


# heterogeneity stats 

hetero_ptsd_accmodpa1<- mr_heterogeneity(dat_ptsd_accmodpa1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_accmodpa1$analysis<- "ptsd_accmodpa1"

hetero_ptsd_accmodpa2<- mr_heterogeneity(dat_ptsd_accmodpa2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_accmodpa2$analysis<- "ptsd_accmodpa2"

# check intercept of Egger regression 

intercept_ptsd_accmodpa1<- mr_pleiotropy_test(dat_ptsd_accmodpa1)
intercept_ptsd_accmodpa1$analysis<- "ptsd_accmodpa1"

intercept_ptsd_accmodpa2<- mr_pleiotropy_test(dat_ptsd_accmodpa2)
intercept_ptsd_accmodpa2$analysis<- "ptsd_accmodpa2"


# steiger directionality test 

steiger_ptsd_accmodpa1 <- directionality_test(dat_ptsd_accmodpa1)
steiger_ptsd_accmodpa1$analysis <- "ptsd_accmodpa1"

steiger_ptsd_accmodpa2 <- directionality_test(dat_ptsd_accmodpa2)
steiger_ptsd_accmodpa2$analysis <- "ptsd_accmodpa2"


# MR presso 

#g1

presso<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                  SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_ptsd_accmodpa1, 
                  NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso[["Main MR results"]]
presso$`MR-PRESSO results`$`Global Test`$RSSobs
presso$`MR-PRESSO results`$`Global Test`$Pvalue

presso_ptsd_accmodpa1<- presso[["Main MR results"]]
presso_ptsd_accmodpa1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                   presso$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_accmodpa1$analysis<- "ptsd_accmodpa1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_ptsd_accmodpa2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_ptsd_accmodpa2<- presso2[["Main MR results"]]
presso_ptsd_accmodpa2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                   presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_accmodpa2$analysis<- "ptsd_accmodpa2"

### 6.Graphs 

#Scatterplot

#g1
res <- mr(dat_ptsd_accmodpa1, method_list=c("mr_ivw_mre",
                                           "mr_egger_regression",
                                           "mr_weighted_median", 
                                           "mr_weighted_mode"))

p1g1<-mr_scatter_plot(res, dat_ptsd_accmodpa1)

#g2
res2 <- mr(dat_ptsd_accmodpa2, method_list=c("mr_ivw_mre",
                                            "mr_egger_regression",
                                            "mr_weighted_median", 
                                            "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_ptsd_accmodpa2)


# combine 
plot_ptsd_accmodpa12 <- ggarrange(p1g1[[paste(exp_ptsd1[1,12],".",out_accmodpa1[1,11], sep="")]],
                                  p1g2[[paste(exp_ptsd2[1,12],".",out_accmodpa2[1,11], sep="")]],
                                  labels = c("(a) Scatter plots", ""),
                                  hjust = -0.05,
                                  ncol = 2, nrow = 1,
                                  common.legend =F)
plot_ptsd_accmodpa12
ggsave("plot_ptsd_accmodpa12.jpeg", device = "jpeg",dpi = 300, width =15, height =7, limitsize = F)


##### (d) Accelerometer-based walking

### 1.Prepare exposure data 
head(ptsd_df)

# first exposure instrument (G1)
head(exp_ptsd1)

# second genetic instrument (G2)
head(exp_ptsd2)

# clumping 
head(exp_ptsd1_clu)
head(exp_ptsd2_clu) 

### 2.Prepare outcome data 

### G1
out_accwalk1 <- read_outcome_data(
  snps = exp_ptsd1_clu$SNP,
  filename = "walking.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF"
  #units_col = "Units",
  #gene_col = "Gene",
)
nrow(exp_ptsd1_clu)
nrow(out_accwalk1)
# if nrow out = exp -> run code on line 832 and then move to G2 (line 890)
# if nrow out != exp -> run code from line 835

out_accwalk1_proxy <- out_accwalk1


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g1<- exp_ptsd1_clu[!exp_ptsd1_clu$SNP %in% out_accwalk1$SNP, ]

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
accwalk_proxy1<- filter(accwalk, SNP %in% df_high$RS_Number) 

# rename column names 

accwalk_proxy1<- accwalk_proxy1 %>%
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
         id.outcome = out_accwalk1[1,11] ,
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
accwalk_proxy1 <- accwalk_proxy1[!accwalk_proxy1$SNP %in% out_accwalk1$SNP, ]

# combine with outcome data
out_accwalk1_proxy<- rbind(out_accwalk1,accwalk_proxy1)


### G2
out_accwalk2 <- read_outcome_data(
  snps = exp_ptsd2_clu$SNP,
  filename = "walking.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF"
  #units_col = "Units",
  #gene_col = "Gene",
)

nrow(exp_ptsd2_clu)
nrow(out_accwalk2)
# if nrow out = exp -> run code on line 911 and then move to Part 3 (data harmonisation, line 969)
# if nrow out != exp -> run code from line 914

out_accwalk2_proxy <- out_accwalk2


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g2<- exp_ptsd2_clu[!exp_ptsd2_clu$SNP %in% out_accwalk2$SNP, ]

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
accwalk_proxy2<- filter(accwalk, SNP %in% df_high$RS_Number) 

# rename column names 

accwalk_proxy2<- accwalk_proxy2 %>%
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
         id.outcome = out_accwalk2[1,11] ,
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
accwalk_proxy2 <- accwalk_proxy2[!accwalk_proxy2$SNP %in% out_accwalk2$SNP, ]

# combine with outcome data
out_accwalk2_proxy<- rbind(out_accwalk2,accwalk_proxy2)


### 3.Harmonise data 
dat_ptsd_accwalk1 <- harmonise_data(
  exposure_dat = exp_ptsd1_clu, 
  outcome_dat = out_accwalk1_proxy, 
  action =2)

dat_ptsd_accwalk2 <- harmonise_data(
  exposure_dat = exp_ptsd2_clu, 
  outcome_dat = out_accwalk2_proxy, 
  action =2)

### 4.Drop duplicate exposure-outcome summary sets

dat_ptsd_accwalk1<-power_prune(dat_ptsd_accwalk1,method=2,dist.outcome="continuous")

dat_ptsd_accwalk2<-power_prune(dat_ptsd_accwalk2,method=2,dist.outcome="continuous")


## Edit names of outcome and exposure
dat_ptsd_accwalk1$exposure[dat_ptsd_accwalk1$exposure == "exposure"] <- "PTSD (G1)"
dat_ptsd_accwalk1$outcome[dat_ptsd_accwalk1$outcome == "outcome"] <- "Accelerometer-based walking"

dat_ptsd_accwalk2$exposure[dat_ptsd_accwalk2$exposure == "exposure"] <- "PTSD (G2)"
dat_ptsd_accwalk2$outcome[dat_ptsd_accwalk2$outcome == "outcome"] <- "Accelerometer-based walking"

dat_ptsd_accwalk1$samplesize.outcome<- 91105
dat_ptsd_accwalk2$samplesize.outcome<- 91105

nrow(dat_ptsd_accwalk1)
nrow(dat_ptsd_accwalk2)

##### 5.MR analyses 

# Ivw, egger, weighted_median, weighted mode 

results_ptsd_accwalk1<-mr(dat_ptsd_accwalk1,method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                         "mr_weighted_median", 
                                                         "mr_weighted_mode"))
results_ptsd_accwalk1

results_ptsd_accwalk2<-mr(dat_ptsd_accwalk2,  method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                           "mr_weighted_median", 
                                                           "mr_weighted_mode"))
results_ptsd_accwalk2

# raps 

# g1
raps<- mr.raps(dat_ptsd_accwalk1$beta.exposure, dat_ptsd_accwalk1$beta.outcome, dat_ptsd_accwalk1$se.exposure, 
               dat_ptsd_accwalk1$se.outcome)

raps_df<- data.frame(id.exposure = results_ptsd_accwalk1[1,1],
                     id.outcome = results_ptsd_accwalk1[1,2], 
                     outcome = results_ptsd_accwalk1[1,3],
                     exposure = results_ptsd_accwalk1[1,4], 
                     method = "raps", 
                     nsnp = results_ptsd_accwalk1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_ptsd_accwalk1<-rbind(results_ptsd_accwalk1,raps_df )
results_ptsd_accwalk1$analysis = "ptsd_accwalk1"
results_ptsd_accwalk1 <- generate_odds_ratios(results_ptsd_accwalk1)

#g2
raps2<- mr.raps(dat_ptsd_accwalk2$beta.exposure, dat_ptsd_accwalk2$beta.outcome, dat_ptsd_accwalk2$se.exposure, 
                dat_ptsd_accwalk2$se.outcome)

raps_df2<- data.frame(id.exposure = results_ptsd_accwalk2[1,1],
                      id.outcome = results_ptsd_accwalk2[1,2], 
                      outcome = results_ptsd_accwalk2[1,3],
                      exposure = results_ptsd_accwalk2[1,4], 
                      method = "raps", 
                      nsnp = results_ptsd_accwalk2[1,6],
                      b = raps2[["beta.hat"]], 
                      se = raps2[["beta.se"]], 
                      pval = raps2[["beta.p.value"]])

results_ptsd_accwalk2<-rbind(results_ptsd_accwalk2,raps_df2 )
results_ptsd_accwalk2$analysis = "ptsd_accwalk2"
results_ptsd_accwalk2 <- generate_odds_ratios(results_ptsd_accwalk2)


# heterogeneity stats 
hetero_ptsd_accwalk1<- mr_heterogeneity(dat_ptsd_accwalk1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_accwalk1$analysis<- "ptsd_accwalk1"

hetero_ptsd_accwalk2<- mr_heterogeneity(dat_ptsd_accwalk2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_accwalk2$analysis<- "ptsd_accwalk2"

# check intercept of Egger regression 
intercept_ptsd_accwalk1<- mr_pleiotropy_test(dat_ptsd_accwalk1)
intercept_ptsd_accwalk1$analysis<- "ptsd_accwalk1"

intercept_ptsd_accwalk2<- mr_pleiotropy_test(dat_ptsd_accwalk2)
intercept_ptsd_accwalk2$analysis<- "ptsd_accwalk2"


# steiger directionality test 

steiger_ptsd_accwalk1 <- directionality_test(dat_ptsd_accwalk1)
steiger_ptsd_accwalk1$analysis <- "ptsd_accwalk1"

steiger_ptsd_accwalk2 <- directionality_test(dat_ptsd_accwalk2)
steiger_ptsd_accwalk2$analysis <- "ptsd_accwalk2"


# MR presso 

#g1

presso<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                  SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_ptsd_accwalk1, 
                  NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso[["Main MR results"]]
presso$`MR-PRESSO results`$`Global Test`$RSSobs
presso$`MR-PRESSO results`$`Global Test`$Pvalue

presso_ptsd_accwalk1<- presso[["Main MR results"]]
presso_ptsd_accwalk1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                  presso$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_accwalk1$analysis<- "ptsd_accwalk1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_ptsd_accwalk2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_ptsd_accwalk2<- presso2[["Main MR results"]]
presso_ptsd_accwalk2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                  presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_accwalk2$analysis<- "ptsd_accwalk2"

### 6.Graphs 

#Scatterplot

#g1
res <- mr(dat_ptsd_accwalk1, method_list=c("mr_ivw_mre",
                                          "mr_egger_regression",
                                          "mr_weighted_median", 
                                          "mr_weighted_mode"))

p1g1<-mr_scatter_plot(res, dat_ptsd_accwalk1)

#g2
res2 <- mr(dat_ptsd_accwalk2, method_list=c("mr_ivw_mre",
                                           "mr_egger_regression",
                                           "mr_weighted_median", 
                                           "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_ptsd_accwalk2)


# combine 
plot_ptsd_accwalk12 <-  ggarrange(p1g1[[paste(exp_ptsd1[1,12],".",out_accwalk1[1,11], sep="")]],
                                  p1g2[[paste(exp_ptsd2[1,12],".",out_accwalk2[1,11], sep="")]],
                                  labels = c("(a) Scatter plots", ""),
                                  hjust = -0.05,
                                  ncol = 2, nrow = 1,
                                  common.legend =F)
plot_ptsd_accwalk12
ggsave("plot_ptsd_accwalk12.jpeg", device = "jpeg",dpi = 300, width =15, height =7, limitsize = F)

##### (e) Accelerometer-based sedentary behavior 

### 1.Prepare exposure data 
head(ptsd_df)

# first exposure instrument (G1)
head(exp_ptsd1)

# second genetic instrument (G2)
head(exp_ptsd2)

# clumping 
head(exp_ptsd1_clu)
head(exp_ptsd2_clu) 

### 2.Prepare outcome data 

### G1
out_accsed1 <- read_outcome_data(
  snps = exp_ptsd1_clu$SNP,
  filename = "sedentarybehaviour.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF"
  #units_col = "Units",
  #gene_col = "Gene",
)
nrow(exp_ptsd1_clu)
nrow(out_accsed1)
# if nrow out = exp -> run code on line 1183 and then move to G2 (line 1241)
# if nrow out != exp -> run code from line 1186

out_accsed1_proxy <- out_accsed1


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g1<- exp_ptsd1_clu[!exp_ptsd1_clu$SNP %in% out_accsed1$SNP, ]

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
accsed_proxy1<- filter(accsed, SNP %in% df_high$RS_Number) 

# rename column names 

accsed_proxy1<- accsed_proxy1 %>%
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
         id.outcome = out_accsed1[1,11] ,
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
accsed_proxy1 <- accsed_proxy1[!accsed_proxy1$SNP %in% out_accsed1$SNP, ]

# combine with outcome data
out_accsed1_proxy<- rbind(out_accsed1,accsed_proxy1)


### G2
out_accsed2 <- read_outcome_data(
  snps = exp_ptsd2_clu$SNP,
  filename = "sedentarybehaviour.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF"
  #units_col = "Units",
  #gene_col = "Gene",
)

nrow(exp_ptsd2_clu)
nrow(out_accsed2)
# if nrow out = exp -> run code on line 1262 and then move to Part 3 (data harmonisation, line 1320)
# if nrow out != exp -> run code from line 1265

out_accsed2_proxy <- out_accsed2


## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_g2<- exp_ptsd2_clu[!exp_ptsd2_clu$SNP %in% out_accsed2$SNP, ]

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
accsed_proxy2<- filter(accsed, SNP %in% df_high$RS_Number) 

# rename column names 

accsed_proxy2<- accsed_proxy2 %>%
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
         id.outcome = out_accsed2[1,11] ,
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
accsed_proxy2 <- accsed_proxy2[!accsed_proxy2$SNP %in% out_accsed2$SNP, ]

# combine with outcome data
out_accsed2_proxy<- rbind(out_accsed2,accsed_proxy2)


### 3.Harmonise data 
dat_ptsd_accsed1 <- harmonise_data(
  exposure_dat = exp_ptsd1_clu, 
  outcome_dat = out_accsed1_proxy, 
  action =2)

dat_ptsd_accsed2 <- harmonise_data(
  exposure_dat = exp_ptsd2_clu, 
  outcome_dat = out_accsed2_proxy, 
  action =2)

### 4.Drop duplicate exposure-outcome summary sets

dat_ptsd_accsed1<-power_prune(dat_ptsd_accsed1,method=2,dist.outcome="continuous")

dat_ptsd_accsed2<-power_prune(dat_ptsd_accsed2,method=2,dist.outcome="continuous")


## Edit names of outcome and exposure
dat_ptsd_accsed1$exposure[dat_ptsd_accsed1$exposure == "exposure"] <- "PTSD (G1)"
dat_ptsd_accsed1$outcome[dat_ptsd_accsed1$outcome == "outcome"] <- "Accelerometer-based sedentary behaviour"

dat_ptsd_accsed2$exposure[dat_ptsd_accsed2$exposure == "exposure"] <- "PTSD (G2)"
dat_ptsd_accsed2$outcome[dat_ptsd_accsed2$outcome == "outcome"] <- "Accelerometer-based sedentary behaviour"

dat_ptsd_accsed1$samplesize.outcome<- 91105
dat_ptsd_accsed2$samplesize.outcome<- 91105

nrow(dat_ptsd_accsed1)
nrow(dat_ptsd_accsed2)

##### 5.MR analyses 

# Ivw, egger, weighted_median, weighted mode 

results_ptsd_accsed1<-mr(dat_ptsd_accsed1,method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                       "mr_weighted_median", 
                                                       "mr_weighted_mode"))
results_ptsd_accsed1

results_ptsd_accsed2<-mr(dat_ptsd_accsed2,  method_list=c( "mr_ivw_mre", "mr_egger_regression",
                                                         "mr_weighted_median", 
                                                         "mr_weighted_mode"))
results_ptsd_accsed2

# raps 

# g1
raps<- mr.raps(dat_ptsd_accsed1$beta.exposure, dat_ptsd_accsed1$beta.outcome, dat_ptsd_accsed1$se.exposure, 
               dat_ptsd_accsed1$se.outcome)

raps_df<- data.frame(id.exposure = results_ptsd_accsed1[1,1],
                     id.outcome = results_ptsd_accsed1[1,2], 
                     outcome = results_ptsd_accsed1[1,3],
                     exposure = results_ptsd_accsed1[1,4], 
                     method = "raps", 
                     nsnp = results_ptsd_accsed1[1,6],
                     b = raps[["beta.hat"]], 
                     se = raps[["beta.se"]], 
                     pval = raps[["beta.p.value"]])

results_ptsd_accsed1<-rbind(results_ptsd_accsed1,raps_df )
results_ptsd_accsed1$analysis = "ptsd_accsed1"
results_ptsd_accsed1 <- generate_odds_ratios(results_ptsd_accsed1)

#g2
raps2<- mr.raps(dat_ptsd_accsed2$beta.exposure, dat_ptsd_accsed2$beta.outcome, dat_ptsd_accsed2$se.exposure, 
                dat_ptsd_accsed2$se.outcome)

raps_df2<- data.frame(id.exposure = results_ptsd_accsed2[1,1],
                      id.outcome = results_ptsd_accsed2[1,2], 
                      outcome = results_ptsd_accsed2[1,3],
                      exposure = results_ptsd_accsed2[1,4], 
                      method = "raps", 
                      nsnp = results_ptsd_accsed2[1,6],
                      b = raps2[["beta.hat"]], 
                      se = raps2[["beta.se"]], 
                      pval = raps2[["beta.p.value"]])

results_ptsd_accsed2<-rbind(results_ptsd_accsed2,raps_df2 )
results_ptsd_accsed2$analysis = "ptsd_accsed2"
results_ptsd_accsed2 <- generate_odds_ratios(results_ptsd_accsed2)


# heterogeneity stats 

hetero_ptsd_accsed1<- mr_heterogeneity(dat_ptsd_accsed1, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_accsed1$analysis<- "ptsd_accsed1"

hetero_ptsd_accsed2<- mr_heterogeneity(dat_ptsd_accsed2, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_ptsd_accsed2$analysis<- "ptsd_accsed2"

# check intercept of Egger regression 

intercept_ptsd_accsed1<- mr_pleiotropy_test(dat_ptsd_accsed1)
intercept_ptsd_accsed1$analysis<- "ptsd_accsed1"

intercept_ptsd_accsed2<- mr_pleiotropy_test(dat_ptsd_accsed2)
intercept_ptsd_accsed2$analysis<- "ptsd_accsed2"


# steiger directionality test 

steiger_ptsd_accsed1 <- directionality_test(dat_ptsd_accsed1)
steiger_ptsd_accsed1$analysis <- "ptsd_accsed1"

steiger_ptsd_accsed2 <- directionality_test(dat_ptsd_accsed2)
steiger_ptsd_accsed2$analysis <- "ptsd_accsed2"


# MR presso 

#g1

presso<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                 SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_ptsd_accsed1, 
                  NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso[["Main MR results"]]
presso$`MR-PRESSO results`$`Global Test`$RSSobs
presso$`MR-PRESSO results`$`Global Test`$Pvalue

presso_ptsd_accsed1<- presso[["Main MR results"]]
presso_ptsd_accsed1$globaltest<-c(presso$`MR-PRESSO results`$`Global Test`$RSSobs,
                                 presso$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_accsed1$analysis<- "ptsd_accsed1"

#g2
presso2<-mr_presso(BetaOutcome = "beta.outcome", 
                   BetaExposure = "beta.exposure", 
                   SdOutcome = "se.outcome", 
                   SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE, 
                   data = dat_ptsd_accsed2, 
                   NbDistribution = 1000,  
                   SignifThreshold = 0.05)

presso_ptsd_accsed2<- presso2[["Main MR results"]]
presso_ptsd_accsed2$globaltest<-c(presso2$`MR-PRESSO results`$`Global Test`$RSSobs,
                                 presso2$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_ptsd_accsed2$analysis<- "ptsd_accsed2"

### 6.Graphs 

#Scatterplot

#g1
res <- mr(dat_ptsd_accsed1, method_list=c("mr_ivw_mre",
                                         "mr_egger_regression",
                                         "mr_weighted_median", 
                                         "mr_weighted_mode"))

p1g1<-mr_scatter_plot(res, dat_ptsd_accsed1)

#g2
res2 <- mr(dat_ptsd_accsed2, method_list=c("mr_ivw_mre",
                                          "mr_egger_regression",
                                          "mr_weighted_median", 
                                          "mr_weighted_mode"))
p1g2<-mr_scatter_plot(res2, dat_ptsd_accsed2)


# combine 
plot_ptsd_accsed12 <-  ggarrange(p1g1[[paste(exp_ptsd1[1,12],".",out_accsed1[1,11], sep="")]],
                                 p1g2[[paste(exp_ptsd2[1,12],".",out_accsed2[1,11], sep="")]],
                                 labels = c("(a) Scatter plots", ""),
                                 hjust = -0.05,
                                 ncol = 2, nrow = 1,
                                 common.legend =F)
plot_ptsd_accsed12
ggsave("plot_ptsd_accsed12.jpeg", device = "jpeg",dpi = 300, width =15, height =7, limitsize = F)

### Dataframes with results of MR and sensitivity analyses for all exposures 

results_ptsd_pa<- rbind(results_ptsd_modvigpaself1, 
                       results_ptsd_modvigpaself2, 
                       results_ptsd_accpamean1, 
                       results_ptsd_accpamean2, 
                       results_ptsd_accmodpa1, 
                       results_ptsd_accmodpa2, 
                       results_ptsd_accwalk1,
                       results_ptsd_accwalk2,
                       results_ptsd_accsed1,
                       results_ptsd_accsed2)

hetero_ptsd_pa<-rbind(hetero_ptsd_modvigpaself1,
                     hetero_ptsd_modvigpaself2,
                     hetero_ptsd_accpamean1,
                     hetero_ptsd_accpamean2,
                     hetero_ptsd_accmodpa1,
                     hetero_ptsd_accmodpa2,
                     hetero_ptsd_accwalk1,
                     hetero_ptsd_accwalk2,
                     hetero_ptsd_accsed1, 
                     hetero_ptsd_accsed2)

intercept_ptsd_pa<-rbind(intercept_ptsd_modvigpaself1,
                        intercept_ptsd_modvigpaself2,
                        intercept_ptsd_accpamean1,
                        intercept_ptsd_accpamean2,
                        intercept_ptsd_accmodpa1,
                        intercept_ptsd_accmodpa2,
                        intercept_ptsd_accwalk1,
                        intercept_ptsd_accwalk2,
                        intercept_ptsd_accsed1,
                        intercept_ptsd_accsed2)

steiger_ptsd_pa<-rbind(steiger_ptsd_modvigpaself1, 
                      steiger_ptsd_modvigpaself2, 
                      steiger_ptsd_accpamean1, 
                      steiger_ptsd_accpamean2, 
                      steiger_ptsd_accmodpa1, 
                      steiger_ptsd_accmodpa2, 
                      steiger_ptsd_accwalk1,
                      steiger_ptsd_accwalk2,
                      steiger_ptsd_accsed1,
                      steiger_ptsd_accsed2)

presso_ptsd_pa<- rbind(presso_ptsd_modvigpaself1,
                      presso_ptsd_modvigpaself2,
                      presso_ptsd_accpamean1,
                      presso_ptsd_accpamean2,
                      presso_ptsd_accmodpa1,
                      presso_ptsd_accmodpa2,
                      presso_ptsd_accwalk1,
                      presso_ptsd_accwalk2,
                      presso_ptsd_accsed1, 
                      presso_ptsd_accsed2)


write_xlsx(list(results_ptsd_pa = results_ptsd_pa,
                hetero_ptsd_pa = hetero_ptsd_pa,
                intercept_ptsd_pa = intercept_ptsd_pa,
                steiger_ptsd_pa = steiger_ptsd_pa,
                presso_ptsd_pa = presso_ptsd_pa), "results_ptsd_pa_all.xlsx")


save.image("~/Documents/Study1_GWAS/ptsd_PA.RData")

