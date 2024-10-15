# eQTL result
eqtl = read.table("gtex_brain/Brain_Spinal_cord_cervical_c-1.v8.EUR.signif_pairs.txt",sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)

library(biomaRt)
# transform ENSG... to gene symbol
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensembl_attrs = c("ensembl_gene_id", "hgnc_symbol")
ensembl_idtype = "ensembl_gene_id"
ensembl_ids = unique(eqtl$phenotype_id)
ensembl_ids_cleaned = gsub("\\.\\d+", "", ensembl_ids)
gene_info = getBM(attributes=ensembl_attrs, filters=ensembl_idtype, values=ensembl_ids_cleaned, mart=ensembl)
colnames(gene_info) = c("ensembl_id", "gene_name")
eqtl$ensembl_ids = gsub("\\.\\d+", "", eqtl$phenotype_id)
colnames(eqtl)[13] = "ensembl_id"
library(dplyr)
eqtl = left_join(eqtl,gene_info)

# transform chr1_777550_T_C_b38 to rs id
snp_list = unique(eqtl$variant_id)
snploc = read.table("snploc.map",sep="\t",quote="",header=FALSE,stringsAsFactors=FALSE)
snp_info <- data.frame(do.call(rbind, strsplit(snp_list, "_"))) %>%
  rename(chr = X1, pos = X2, ref = X3, alt = X4)
snp_info$variant_id = snp_list
snp_A1_A2 = data.frame(variant_id = snp_info$variant_id, A1 = snp_info$ref, A2 = snp_info$alt)
snp_info$chr = gsub("chr","",snp_info$chr)
snp_info$chr = as.integer(snp_info$chr)
snp_info$pos = as.integer(snp_info$pos)
snp_info = na.omit(snp_info)
colnames(snploc) = c("chr","rsid","unknown","pos")
snp_to_rsid = left_join(snp_info,snploc)
rm(snploc)
library(tidyr)
snp_to_rsid = snp_to_rsid %>% unite(variant_id, chr:X5, sep = "_")
snp_to_rsid = snp_to_rsid[,-3]
snp_to_rsid$variant_id = paste("chr", snp_to_rsid$variant_id, sep = "")
eqtl = left_join(eqtl,snp_to_rsid)
eqtl$gene_name[eqtl$gene_name==""] = NA
eqtl = na.omit(eqtl)
colnames(eqtl)[14] = "gene_name"
eqtl = left_join(eqtl, snp_A1_A2)

causal_gene = unique(read.csv("causal_gene_result.CSV")$Gene)
eqtl_filter = subset(eqtl, eqtl$gene_name %in% causal_gene)

# LD clump
library(ieugwasr)
library(plinkbinr)
clump_data = data.frame(rsid=eqtl_filter$rsid,id=eqtl_filter$gene_name,pval=eqtl_filter$pval_nominal)
gene_list = unique(clump_data$id)
result = data.frame()
error_list = c()
for (i in 1:length(gene_list)){
  delayedAssign("do.next", {next})
  tmp_gene = gene_list[i]
  tmp_clump_data = subset(clump_data,clump_data$id==tmp_gene)
  tmp_result = tryCatch({
    ld_clump(dat=tmp_clump_data,clump_kb=2000,clump_r2=0.01,clump_p=1,plink_bin=get_plink_exe(),bfile="LDpruning/EUR")
  }, error = function(i){
    error_list = c(error_list,i)
    force(do.next)
  })
  result = rbind(result,tmp_result)
}

clump_result = result
colnames(clump_result)[2] = "gene_name"

# preprocess MR data
mrdata = left_join(clump_result,eqtl)

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
count_result = dplyr::count(adhd_mrdata,gene_name)
over_3 = subset(count_result,count_result$n>=1)
over3_adhd_mrdata = subset(adhd_mrdata,adhd_mrdata$gene_name %in% gene_adhd)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_adhd_mrdata,
  type = "exposure",
  phenotype_col = "gene_name",
  snp_col = "rsid",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval_nominal"
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
significant = subset(mr_model,mr_model$pval<0.05/length(over_3$gene_name))
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
postive_gene = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% postive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

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
# 筛选大于等于3
count_result = dplyr::count(bipolar_mrdata,gene_name)
over_3 = subset(count_result,count_result$n>=1)
over3_bipolar_mrdata = subset(bipolar_mrdata,bipolar_mrdata$gene_name %in% gene_bipolar)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_bipolar_mrdata,
  type = "exposure",
  phenotype_col = "gene_name",
  snp_col = "rsid",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval_nominal"
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
significant = subset(mr_model,mr_model$pval<0.05/length(over_3$gene_name))
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
postive_gene = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% postive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

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
count_result = dplyr::count(schizo_mrdata,gene_name)
over_3 = subset(count_result,count_result$n>=1)
over3_schizo_mrdata = subset(schizo_mrdata,schizo_mrdata$gene_name %in% gene_schizo)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_schizo_mrdata,
  type = "exposure",
  phenotype_col = "gene_name",
  snp_col = "rsid",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval_nominal"
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
significant = subset(mr_model,mr_model$pval<0.05)
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
postive_gene = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% postive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

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
count_result = dplyr::count(depression_mrdata,gene_name)
over_3 = subset(count_result,count_result$n>=1)
over3_depression_mrdata = subset(depression_mrdata,depression_mrdata$gene_name %in% gene_depression)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_depression_mrdata,
  type = "exposure",
  phenotype_col = "gene_name",
  snp_col = "rsid",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval_nominal"
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
significant = subset(mr_model,mr_model$pval<0.05/length(over_3$gene_name))
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
postive_gene = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% postive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

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
count_result = dplyr::count(insomnia_mrdata,gene_name)
over_3 = subset(count_result,count_result$n>=1)
over3_insomnia_mrdata = subset(insomnia_mrdata,insomnia_mrdata$gene_name %in% gene_insomnia)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_insomnia_mrdata,
  type = "exposure",
  phenotype_col = "gene_name",
  snp_col = "rsid",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval_nominal"
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
significant = subset(mr_model,mr_model$pval<0.05/length(unique(mr_model$exposure)))
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
postive_gene = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% postive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

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
count_result = dplyr::count(autism_mrdata,gene_name)
over_3 = subset(count_result,count_result$n>=1)
over3_autism_mrdata = subset(autism_mrdata,autism_mrdata$gene_name %in% gene_autism)
# MR analysis
library(TwoSampleMR)
exposure_dat = format_data(
  over3_autism_mrdata,
  type = "exposure",
  phenotype_col = "gene_name",
  snp_col = "rsid",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval_nominal"
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
significant = subset(mr_model,mr_model$pval<0.05/length(over_3$gene_name))
ivw = subset(significant,significant$method=="Inverse variance weighted")
wm = subset(significant,significant$method=="Weighted median")
# heterogeneity test and pleiotropy test
postive_gene = unique(significant$exposure)
positive_harmonise_dat = subset(harmonise_dat,harmonise_dat$exposure %in% postive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)
