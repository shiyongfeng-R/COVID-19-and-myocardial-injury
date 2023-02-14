library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(forestplot)
library(forestploter)
library(haven)
library(autoReg)
# First step:import the COVID19 datasets
Covid_population <- read.delim("~/Desktop/NALFD AND heart failure/COVID19数据源7/Covid vs. populationCOVID19_HGI_C2_ALL_leave_23andme_20220403_GRCh37.tsv")
View(Covid_population)
hVSnothosp <- read.delim("~/Desktop/NALFD AND heart failure/COVID19数据源7/Hospitalized covid vs. not hospitalized covidCOVID19_HGI_B1_ALL_leave_23andme_20220403_GRCh37.tsv")
View(hVSnothosp)
hvspop <- read.delim("~/Desktop/NALFD AND heart failure/COVID19数据源7/Hospitalized covid vs. populationCOVID19_HGI_B2_ALL_leave_23andme_20220403_GRCh37.tsv")
View(hvspop)
Verysere <- read.delim("~/Desktop/NALFD AND heart failure/COVID19数据源7/Very severe respiratory confirmed covid vs. populationCOVID19_HGI_A2_ALL_leave_23andme_20220403_GRCh37.tsv")
View(Verysere)

# Prepare the standard data
a <- Verysere
a <- hVSnothosp
a <- Covid_population
a <- hvspop


# COVID-19 version7 transformation
names(a)
names(a)[5] <- "snp1"
names(a)[c(18,4,3,7,8,17,9)] <- c("SNP","effect_allele","other_allele","beta","se","eaf","pval")
names(a)

# downloaded datasets were changed into standard exposure data

b <- subset(a,pval<5e-08) # normal value of p-value

b <- subset(a,pval<1e-09) # tightened threshold value of P-value 

c <- format_data(b,type = "exposure")


exp_dat_clumped <- clump_data(c,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 1,clump_p2 = 1,pop = "EUR")


# download outcome dataset

a <- th

names(a)
names(a)[c(5,6,7,11,12,14)] <- c("effect_allele","other_allele","eaf","beta","se","pval")
bmi_exp_dat_clumped <- exp_dat

d <- merge(bmi_exp_dat_clumped,a,by.x = "SNP",by.y ="SNP" )
d <- merge(bmi_exp_dat_clumped,t,by.x = "SNP",by.y ="SNP" )
write.csv(d,file = "outcome.csv")

outcome_dat <- read_outcome_data(snps = bmi_exp_dat_clumped$SNP,filename = "outcome.csv",sep = ",",snp_col = "SNP",
                                 beta_col = "beta",se_col = "se",effect_allele_col = "effect_allele",other_allele_col = "other_allele",
                                 eaf_col = "eaf",pval_col = "pval",samplesize_col = "samplesize")

## prepare outcome datasets
exp_dat <- exp_dat_clumped

out_dat <- extract_outcome_data(snps = exp_dat$SNP,outcomes="prot-a-3066") # TNI


dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat = out_dat)


dat <- subset(dat,pval.outcome>5e-08)

res1<- generate_odds_ratios(mr_res = mr(dat))
res1

res1<- generate_odds_ratios(mr_res = mr(dat2,method_list = "mr_ivw"))
res1<- generate_odds_ratios(mr_res = mr(dat2,method_list =  c("mr_ivw","mr_egger_regression","mr_weighted_median")))
res1
write.csv(res1,file = "result2.csv")

# mv MR

id_exposure <- c("ebi-a-GCST011494","ebi-a-GCST006414","ukb-b-12493","ebi-a-GCST006867","finn-b-I9_NONISCHCARDMYOP") 
#("ieu-b-38","ukb-b-12493","ebi-a-GCST011365")
id_exposure <- c("ebi-a-GCST006414","ukb-b-12493")
HF_riskfactor <- c("ebi-a-GCST011365","ieu-b-40","ebi-a-GCST006867","ukb-b-12493")
id_exposure <- c("ebi-a-GCST006414","ebi-a-GCST011494","ebi-a-GCST011365","ebi-a-GCST011365")
id_outcome <- "ebi-a-GCST009541"

exposure_dat <- mv_extract_exposures(id_exposure)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
res <- mv_multiple(mvdat)
res
res <- data.frame(res)
names(res) <- c("id.exposure","exposure","id.outcome","outcome","nsnp","b","se","pval")
dat1 <- generate_odds_ratios(res)
mv_multiple(mvdat,plots = T)

# 计算中介占比
a <- dat$b
a <- dat$lo_ci
a <- dat$up_ci
b <- a[1]*a[2]/(a[3]+a[1]*a[2])
b

dat <- harmonise_data(exposure_dat = bmi_exp_dat_clumped,outcome_dat = outcome_dat)
dat2 <- subset(dat,pval.outcome>5e-08)
res1<- generate_odds_ratios(mr_res = mr(dat2))
res1

write.csv(res1,file = "result2.csv")
# 绘图
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat2))

mr_scatter_plot(mr_results = mr(dat2,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat2)
mr_heterogeneity(dat2)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat2))
mr_pleiotropy_test(dat2)
mr_forest_plot(mr_singlesnp(dat2))
write.csv(m1,file = "m1.csv")
write.csv(m11,file = "m11.csv")

run_mr_presso(dat2,NbDistribution = 10000)
# 计算 F
n=2942817
n=1163698
n=337159
n=462933
n=33394
n=977323
n=1030836
n=395795
n=681275
n=655666
n=463010
k=21
k=18
k=446
k=6
k=64
k=96
dat2 <- exp_dat
# mean F
dat2$F <- 2*(1-dat2$eaf.exposure)*dat2$eaf.exposure*(dat2$beta.exposure/(dat2$se.exposure*dat2$samplesize.exposure^0.5))*(dat2$beta.exposure/(dat2$se.exposure*dat2$samplesize.exposure^0.5))
dat2$F1 <- (dat2$F/(1-dat2$F))*(dat2$samplesize.exposure-2)
mean(dat2$F1)

# TOTAL F
dat <- dat2[dat2$mr_keep,]
dat$F <- 2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure/(dat$se.exposure*(dat$samplesize.exposure)^0.5))^2
R2 <- sum(dat$F)
value <- R2*(n-k-1)/((1-R2)*k)
value


#绘制森林图数据准备
dataset <- dis
dataset <- res1
dataset$pval <- sprintf("%0.4f",dataset$pval)

dataset <- dataset[dataset$method=="Inverse variance weighted",]
dataset$Discovery <- paste(rep(" ",40),collapse = "")
dataset$Replication <- paste(rep(" ",40),collapse = "")

dataset$'' <- paste(rep(" ",40),collapse = "")


dataset$or <- sprintf("%0.2f",dataset$or)
dataset$or_lci95 <- sprintf("%0.2f",dataset$or_lci95)
dataset$or_uci95 <- sprintf("%0.2f",dataset$or_uci95)

dataset$'OR(95%CI)' <- paste(dataset$or,"(",dataset$or_lci95,"-",dataset$or_uci95,")")
names(dataset)
names(dataset)[c(3,4,5,6,9)] <- c("Outcome","Exposure","Method","nSNP","Pval")
dataset$Method[grepl("Inverse",dataset$Method)] <- "IVW"
#绘制森林图perfect
tm <- forest_theme(base_size = 10,
                   refline_col = "green",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")


plot <- forest(dataset[,c(4,3,5,6,15,16,9)],est=dataset$or,lower=dataset$or_lci95,
               upper=dataset$or_uci95,
               ci_column = 5,
               ref_line = 1,
               vert_line =NULL,
               xlim = c(0.5,1.5),
               ticks_at = c(0.5,1,1.5),
               theme = tm, arrow_lab=NULL,xlab = "OR(95%CI)")
add_underline(plot, row = NULL, col = NULL, part = "header")


# 修改exposure and outcome
dataset$Exposure[grepl("Atrial",dataset$Exposure)] <- "AF"
dataset$Exposure[grepl("Heart",dataset$Exposure)] <- "HF"
dataset$Outcome[grepl("Heart",dataset$Outcome)] <- "HF"
dataset$Exposure[grepl("hypertension",dataset$Exposure)] <- "Hypertension"
dataset$Outcome[grepl("Atrial",dataset$Outcome)] <- "AF"
dataset$Outcome[grepl("diabetes",dataset$Outcome)] <- "Diabetes"
dataset$Outcome[grepl("body",dataset$Outcome)] <- "BMI"
dataset$Outcome[grepl("Myocardial",dataset$Outcome)] <- "MI"
dataset$Outcome[grepl("hypertension",dataset$Outcome)] <- "Hypertension"
dataset$Exposure[grepl("Heart",dataset$Exposure)] <- "HF"
#绘制森林图数据准备
dataset <- dis3
dataset <- dis 
dataset <- res1
dataset$outcome <- c("Waist circumference","Apo-A1","Apo-B","HDL-C","LDL-C","TG","SBP","BMI")
dataset <- dataset[dataset$method=="Inverse variance weighted",]
dataset <- dataset[-c(1,2,3),]
dataset$pval <- round(dataset$pval,3)

dataset$Replication <- paste(rep(" ",40),collapse = "")
dataset$Discovery <- paste(rep(" ",40),collapse = "")
dataset$b <- round(dataset$b,2)
dataset$lo_ci <- round(dataset$lo_ci,2)
dataset$up_ci <- round(dataset$up_ci,2)
dataset$'beta(95%CI)' <- paste(dataset$b,"(",dataset$lo_ci,"-",dataset$up_ci,")")
names(dataset)
names(dataset)[c(3,6,9)] <- c("Outcome","nSNP","Pval")
#绘制森林图perfect
tm <- forest_theme(base_size = 10,
                   refline_col = "green",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")


plot <- forest(dataset[,c(3,15,6,16,9)],est=dataset$b,lower=dataset$lo_ci,
               upper=dataset$up_ci,
               ci_column = 2,
               ref_line = 0,
               vert_line =NULL,
               xlim = c(-2,2),
               ticks_at = c(-2,-1,0,1,2),
               theme = tm, arrow_lab=NULL,xlab = "beta(95%CI)")
add_underline(plot, row = NULL, col = NULL, part = "header")

library(forestplot)
install.packages("haven")
library(haven)
library(dplyr)
library(autoReg)

ft <- gaze(method~.,data=dis) %>% myft()
ft
