library(TwoSampleMR)
library(dplyr)
library(ieugwasr)
library(plinkbinr)

causal_cpg = read.csv("causal_cpg_result.CSV")
causal_cpg_adhd = subset(causal_cpg,causal_cpg$Phenotype=="ADHD")$CpG
causal_cpg_autism = subset(causal_cpg,causal_cpg$Phenotype=="Autism")$CpG
causal_cpg_bipolar = subset(causal_cpg,causal_cpg$Phenotype=="Bipolar")$CpG
causal_cpg_depression = subset(causal_cpg,causal_cpg$Phenotype=="Depression")$CpG
causal_cpg_insomnia = subset(causal_cpg,causal_cpg$Phenotype=="Insomnia")$CpG
causal_cpg_schizo = subset(causal_cpg,causal_cpg$Phenotype=="Schizophrenia")$CpG

# filter meqtl
cis = subset(cis,cis$discovery.population != "SA")
meqtl = subset(cis,cis$cpg %in% causal_cpg$CpG)
meqtl = data.frame(rsid=meqtl$snp,A1=meqtl$A1,A2=meqtl$A2,MAF=meqtl$eaf.eu,cpg=meqtl$cpg,beta=meqtl$beta.eur,se=meqtl$se.eur,`p-value`=meqtl$p.eur)
meqtl$pair = paste(meqtl$cpg,meqtl$rsid,sep="_")

# LDpruning
clump_data = data.frame(rsid=meqtl$rsid,id=meqtl$cpg,pval=meqtl$p.value,pair=meqtl$pair)
cpg_list = unique(clump_data$id)
result = data.frame()
error_list = c()
for (i in 1:52){
  delayedAssign("do.next", {next})
  tmp_cpg = cpg_list[i]
  tmp_clump_data = subset(clump_data,clump_data$id==tmp_cpg)
  tmp_result = tryCatch({
    ld_clump(dat=tmp_clump_data,clump_kb=2000,clump_r2=0.01,clump_p=1,plink_bin=get_plink_exe(),bfile="LDpruning/EUR")
  }, error = function(i){
    error_list = c(error_list,i)
    force(do.next)
  })
  result = rbind(result,tmp_result)
}
clump_result = result
colnames(clump_result)[2] = "cpg"

# preprocess MRdata
mrdata = left_join(clump_result,meqtl)
mrdata = mrdata[-10]

# MR analysis for ADHD
adhd = read.table("GWAS/ADHD.meta",header=TRUE)
use_snp = unique(clump_result$rsid)
use_adhd = subset(adhd,adhd$SNP %in% use_snp)
rm(adhd)
adhd_mrdata = subset(mrdata,mrdata$rsid %in% use_adhd$SNP)
adhd_snp = unique(adhd_mrdata$rsid)
use_adhd = subset(use_adhd,use_adhd$SNP %in% adhd_snp)
use_adhd$beta = log(use_adhd$OR)
use_adhd$se = abs(use_adhd$beta/qnorm(use_adhd$P/2,lower.tail=F))
use_adhd$phenotype = "ADHD"
over_3 = subset(count(adhd_mrdata,cpg),n>=1)
over3_adhd_mrdata = subset(adhd_mrdata,adhd_mrdata$cpg %in% over_3$cpg)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_adhd_mrdata,
  type = "exposure",
  phenotype_col = "cpg",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
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
mr_model = mr(subset(harmonise_dat,harmonise_dat$exposure %in% causal_cpg_adhd))

# MR analysis for bipolar
bipolar = read.table("GWAS/bipolar.tsv",sep="\t",quote="")
colnames(bipolar)[3:10] = c("SNP","A1","A2","beta","se","P","unknown","MAF")
use_snp = unique(clump_result$rsid)
use_bipolar = subset(bipolar,bipolar$SNP %in% use_snp)
rm(bipolar)
bipolar_mrdata = subset(mrdata,mrdata$rsid %in% use_bipolar$SNP)
bipolar_snp = unique(bipolar_mrdata$rsid)
use_bipolar = subset(use_bipolar,use_bipolar$SNP %in% bipolar_snp)
use_bipolar$phenotype = "bipolar"
over_3 = subset(count(bipolar_mrdata,cpg),n>=1)
over3_bipolar_mrdata = subset(bipolar_mrdata,bipolar_mrdata$cpg %in% over_3$cpg)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_bipolar_mrdata,
  type = "exposure",
  phenotype_col = "cpg",
  snp_col = "rsid",
  beta_col ="beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
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
mr_model = mr(subset(harmonise_dat,harmonise_dat$exposure %in% causal_cpg_bipolar))

# MR analysis for schizo
schizo = read.table("GWAS/schizophrenia.tsv",header=TRUE,sep="\t",quote="")
colnames(schizo)[2] = "SNP"
colnames(schizo)[6] = "MAF"
colnames(schizo)[9:11] = c("beta","se","P")
use_snp = unique(clump_result$rsid)
use_schizo = subset(schizo,schizo$SNP %in% use_snp)
rm(schizo)
schizo_mrdata = subset(mrdata,mrdata$rsid %in% use_schizo$SNP)
schizo_snp = unique(schizo_mrdata$rsid)
use_schizo = subset(use_schizo,use_schizo$SNP %in% schizo_snp)
use_schizo$phenotype = "schizo"
over_3 = subset(count(schizo_mrdata,cpg),n>=1)
over3_schizo_mrdata = subset(schizo_mrdata,schizo_mrdata$cpg %in% over_3$cpg)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_schizo_mrdata,
  type = "exposure",
  phenotype_col = "cpg",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
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
mr_model = mr(subset(harmonise_dat,harmonise_dat$exposure %in% causal_cpg_schizo))

# MR analysis for depression
depression = read.table("GWAS/depression2.txt",header=TRUE,quote="")
colnames(depression)[1:6] = c("SNP","A1","A2","MAF","beta","se")
use_snp = unique(clump_result$rsid)
use_depression = subset(depression,depression$SNP %in% use_snp)
rm(depression)
use_depression$A1 = toupper(use_depression$A1)
use_depression$A2 = toupper(use_depression$A2)
depression_mrdata = subset(mrdata,mrdata$rsid %in% use_depression$SNP)
depression_snp = unique(depression_mrdata$rsid)
use_depression = subset(use_depression,use_depression$SNP %in% depression_snp)
use_depression$phenotype = "depression"
over_3 = subset(count(depression_mrdata,cpg),n>=1)
over3_depression_mrdata = subset(depression_mrdata,depression_mrdata$cpg %in% over_3$cpg)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_depression_mrdata,
  type = "exposure",
  phenotype_col = "cpg",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
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
mr_model = mr(subset(harmonise_dat,harmonise_dat$exposure %in% causal_cpg_depression))

# MR analysis for insomnia
insomnia = read.table("GWAS/insomnia",header=TRUE,quote="")
use_insomnia = subset(insomnia,insomnia$SNP %in% use_snp)
rm(insomnia)
insomnia_mrdata = subset(mrdata,mrdata$rsid %in% use_insomnia$SNP)
insomnia_snp = unique(insomnia_mrdata$rsid)
use_insomnia = subset(use_insomnia,use_insomnia$SNP %in% insomnia_snp)
use_insomnia$beta = log(use_insomnia$OR)
use_insomnia$se = abs(use_insomnia$beta/qnorm(use_insomnia$P/2,lower.tail=F))
use_insomnia$phenotype = "insomnia"
over_3 = subset(count(insomnia_mrdata,cpg),n>=1)
over3_insomnia_mrdata = subset(insomnia_mrdata,insomnia_mrdata$cpg %in% over_3$cpg)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_insomnia_mrdata,
  type = "exposure",
  phenotype_col = "cpg",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
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
mr_model = mr(subset(harmonise_dat,harmonise_dat$exposure %in% causal_cpg_insomnia))

# MR analysis for autism
autism = read.table("GWAS/autism",header=TRUE,sep="\t",quote="")
use_autism = subset(autism,autism$SNP %in% use_snp)
rm(autism)
autism_mrdata = subset(mrdata,mrdata$rsid %in% use_autism$SNP)
autism_snp = unique(autism_mrdata$rsid)
use_autism = subset(use_autism,use_autism$SNP %in% autism_snp)
use_autism$beta = log(use_autism$OR)
use_autism$se = abs(use_autism$beta/qnorm(use_autism$P/2,lower.tail=F))
use_autism$phenotype = "autism"
over_3 = subset(count(autism_mrdata,cpg),n>=1)
over3_autism_mrdata = subset(autism_mrdata,autism_mrdata$cpg %in% over_3$cpg)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_autism_mrdata,
  type = "exposure",
  phenotype_col = "cpg",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
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
mr_model = mr(subset(harmonise_dat,harmonise_dat$exposure %in% causal_cpg_autism))
