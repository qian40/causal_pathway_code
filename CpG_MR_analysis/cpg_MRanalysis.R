library(TwoSampleMR)
library(dplyr)
library(ieugwasr)
library(plinkbinr)

# load pruning result
load("meqtl_clumpresult.Rdata")

load("meqtl.Rdata")
colnames(result) = c("CpG","SNP","p-value")
all_mrdata = left_join(result,meqtl)
rm(meqtl)

# filter CpG with over 3 meQTLs
count_result = dplyr::count(result,CpG)
final_result = subset(count_result,count_result$n>=3)
use_cpg = unique(final_result$CpG)
mrdata = subset(all_mrdata,all_mrdata$CpG %in% use_cpg)
use_snp = unique(mrdata$SNP)

freqdata = read.table("freqdata.frq",head=TRUE)
colnames(freqdata) = c("chr","snpid","A1","A2","MAF","NCHROBS")
snpmap = read.table("snploc.map")
colnames(snpmap) = c("chr","SNP","unknown","pos")
freqdata$SNP = snpmap$SNP
use_freq = subset(freqdata,freqdata$SNP %in% use_snp)
rm(freqdata)
rm(snpmap)

# MR analysis for ADHD
adhd = read.table("GWAS/ADHD.meta",header=TRUE)
use_adhd = subset(adhd,adhd$SNP %in% use_snp)
rm(adhd)
adhd_mrdata = subset(mrdata,mrdata$SNP %in% use_adhd$SNP)
adhd_freq = subset(use_freq,use_freq$SNP %in% use_adhd$SNP)
count_result = dplyr::count(adhd_mrdata,CpG)
final_result = subset(count_result,count_result$n>=3)
adhd_cpg = final_result$CpG
adhd_mrdata = subset(adhd_mrdata,adhd_mrdata$CpG %in% adhd_cpg)
adhd_snp = unique(adhd_mrdata$SNP)
use_adhd = subset(use_adhd,use_adhd$SNP %in% adhd_snp)
adhd_freq = subset(adhd_freq,adhd_freq$SNP %in% adhd_snp)
adhd_freq = select(adhd_freq,c("SNP","A1","A2","MAF"))
adhd_mrdata = left_join(adhd_mrdata,adhd_freq)
use_adhd$beta = log(use_adhd$OR)
use_adhd$se = abs(use_adhd$beta/qnorm(use_adhd$P/2,lower.tail=F))
use_adhd$phenotype = "ADHD"
# MR analysis
exposure_dat = format_data(
  adhd_mrdata,
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
outcome_dat = format_data(
  use_adhd,
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "FRQ_U_34194",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/7568)
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% positive_cpg)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# MR analysis for autism
autism = read.table("GWAS/autism",header=TRUE,sep="\t",quote="")
use_autism = subset(autism,autism$SNP %in% use_snp)
rm(autism)
autism_mrdata = subset(mrdata,mrdata$SNP %in% use_autism$SNP)
autism_freq = subset(use_freq,use_freq$SNP %in% use_autism$SNP)
count_result = dplyr::count(autism_mrdata,CpG)
final_result = subset(count_result,count_result$n>=3)
autism_cpg = final_result$CpG
autism_mrdata = subset(autism_mrdata,autism_mrdata$CpG %in% autism_cpg)
autism_snp = unique(autism_mrdata$SNP)
use_autism = subset(use_autism,use_autism$SNP %in% autism_snp)
autism_freq = subset(autism_freq,autism_freq$SNP %in% autism_snp)
autism_freq = select(autism_freq,c("SNP","A1","A2","MAF"))
autism_mrdata = left_join(autism_mrdata,autism_freq)
use_autism$beta = log(use_autism$OR)
use_autism$se = abs(use_autism$beta/qnorm(use_autism$P/2,lower.tail=F))
use_autism$phenotype = "autism"
# MR analysis
exposure_dat = format_data(
  autism_mrdata,
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
outcome_dat = format_data(
  use_autism,
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/7644)
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% positive_cpg)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# MR analysis for bipolar
bipolar = read.table("GWAS/bipolar.tsv",sep="\t",quote="")
colnames(bipolar)[3:10] = c("SNP","A1","A2","beta","se","P","unknown","MAF")
use_bipolar = subset(bipolar,bipolar$SNP %in% use_snp)
rm(bipolar)
bipolar_mrdata = subset(mrdata,mrdata$SNP %in% use_bipolar$SNP)
bipolar_freq = subset(use_freq,use_freq$SNP %in% use_bipolar$SNP)
count_result = dplyr::count(bipolar_mrdata,CpG)
final_result = subset(count_result,count_result$n>=3)
bipolar_cpg = final_result$CpG
bipolar_mrdata = subset(bipolar_mrdata,bipolar_mrdata$CpG %in% bipolar_cpg)
bipolar_snp = unique(bipolar_mrdata$SNP)
use_bipolar = subset(use_bipolar,use_bipolar$SNP %in% bipolar_snp)
bipolar_freq = subset(bipolar_freq,bipolar_freq$SNP %in% bipolar_snp)
bipolar_freq = select(bipolar_freq,c("SNP","A1","A2","MAF"))
bipolar_mrdata = left_join(bipolar_mrdata,bipolar_freq)
use_bipolar$phenotype = "bipolar"
# MR analysis
exposure_dat = format_data(
  bipolar_mrdata,
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
outcome_dat = format_data(
  use_bipolar,
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/7466)
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% positive_cpg)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# MR analysis for depression
depression = read.table("GWAS/depression2.txt",header=TRUE,quote="")
colnames(depression)[1:6] = c("SNP","A1","A2","MAF","beta","se")
use_depression = subset(depression,depression$SNP %in% use_snp)
rm(depression)
use_depression$A1 = toupper(use_depression$A1)
use_depression$A2 = toupper(use_depression$A2)
depression_mrdata = subset(mrdata,mrdata$SNP %in% use_depression$SNP)
depression_freq = subset(use_freq,use_freq$SNP %in% use_depression$SNP)
count_result = dplyr::count(depression_mrdata,CpG)
final_result = subset(count_result,count_result$n>=3)
depression_cpg = final_result$CpG
depression_mrdata = subset(depression_mrdata,depression_mrdata$CpG %in% depression_cpg)
depression_snp = unique(depression_mrdata$SNP)
use_depression = subset(use_depression,use_depression$SNP %in% depression_snp)
depression_freq = subset(depression_freq,depression_freq$SNP %in% depression_snp)
depression_freq = select(depression_freq,c("SNP","A1","A2","MAF"))
depression_mrdata = left_join(depression_mrdata,depression_freq)
use_depression$phenotype = "depression"
# MR analysis
exposure_dat = format_data(
  depression_mrdata,
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
outcome_dat = format_data(
  use_depression,
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/6790)
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% positive_cpg)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# MR analysis for insomnia
insomnia = read.table("GWAS/insomnia",header=TRUE,quote="")
use_insomnia = subset(insomnia,insomnia$SNP %in% use_snp)
rm(insomnia)
insomnia_mrdata = subset(mrdata,mrdata$SNP %in% use_insomnia$SNP)
insomnia_freq = subset(use_freq,use_freq$SNP %in% use_insomnia$SNP)
count_result = dplyr::count(insomnia_mrdata,CpG)
final_result = subset(count_result,count_result$n>=3)
insomnia_cpg = final_result$CpG
insomnia_mrdata = subset(insomnia_mrdata,insomnia_mrdata$CpG %in% insomnia_cpg)
insomnia_snp = unique(insomnia_mrdata$SNP)
use_insomnia = subset(use_insomnia,use_insomnia$SNP %in% insomnia_snp)
insomnia_freq = subset(insomnia_freq,insomnia_freq$SNP %in% insomnia_snp)
insomnia_freq = select(insomnia_freq,c("SNP","A1","A2","MAF"))
insomnia_mrdata = left_join(insomnia_mrdata,insomnia_freq)
use_insomnia$beta = log(use_insomnia$OR)
use_insomnia$se = abs(use_insomnia$beta/qnorm(use_insomnia$P/2,lower.tail=F))
use_insomnia$phenotype = "insomnia"
# MR analysis
exposure_dat = format_data(
  insomnia_mrdata,
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
outcome_dat = format_data(
  use_insomnia,
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/7391)
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% positive_cpg)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# MR analysis for schizo
schizo = read.table("GWAS/schizophrenia.tsv",header=TRUE,sep="\t",quote="")
colnames(schizo)[2] = "SNP"
colnames(schizo)[6] = "MAF"
colnames(schizo)[9:11] = c("beta","se","P")
use_schizo = subset(schizo,schizo$SNP %in% use_snp)
rm(schizo)
schizo_mrdata = subset(mrdata,mrdata$SNP %in% use_schizo$SNP)
schizo_freq = subset(use_freq,use_freq$SNP %in% use_schizo$SNP)
count_result = dplyr::count(schizo_mrdata,CpG)
final_result = subset(count_result,count_result$n>=3)
schizo_cpg = final_result$CpG
schizo_mrdata = subset(schizo_mrdata,schizo_mrdata$CpG %in% schizo_cpg)
schizo_snp = unique(schizo_mrdata$SNP)
use_schizo = subset(use_schizo,use_schizo$SNP %in% schizo_snp)
schizo_freq = subset(schizo_freq,schizo_freq$SNP %in% schizo_snp)
schizo_freq = select(schizo_freq,c("SNP","A1","A2","MAF"))
schizo_mrdata = left_join(schizo_mrdata,schizo_freq)
use_schizo$phenotype = "schizo"
# MR analysis
exposure_dat = format_data(
  schizo_mrdata,
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
outcome_dat = format_data(
  use_schizo,
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/7466)
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% positive_cpg)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)
