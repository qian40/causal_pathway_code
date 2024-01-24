causal_cpg = read.csv("causal_cpg_result.CSV")
causal_gene = read.csv("causal_gene_result.CSV")

# causal CpG
cpg_adhd = subset(causal_cpg,causal_cpg$Phenotype=="ADHD")$CpG
cpg_autism = subset(causal_cpg,causal_cpg$Phenotype=="Autism")$CpG
cpg_bipolar = subset(causal_cpg,causal_cpg$Phenotype=="Bipolar")$CpG
cpg_depression = subset(causal_cpg,causal_cpg$Phenotype=="Depression")$CpG
cpg_insomnia = subset(causal_cpg,causal_cpg$Phenotype=="Insomnia")$CpG
cpg_schizo = subset(causal_cpg,causal_cpg$Phenotype=="Schizophrenia")$CpG

# causal Gene
gene_adhd = subset(causal_gene,causal_gene$Phenotype=="ADHD")$Gene
gene_autism = subset(causal_gene,causal_gene$Phenotype=="Autism")$Gene
gene_bipolar = subset(causal_gene,causal_gene$Phenotype=="Bipolar")$Gene
gene_depression = subset(causal_gene,causal_gene$Phenotype=="Depression")$Gene
gene_insomnia = subset(causal_gene,causal_gene$Phenotype=="Insomnia")$Gene
gene_schizo = subset(causal_gene,causal_gene$Phenotype=="Schizophrenia")$Gene

# meQTL of causal CpGs as IVs
use_result = subset(result,result$id %in% causal_cpg$CpG)
colnames(use_result) = c("CpG","SNP","pval")
use_meqtl = left_join(use_result,my_meqtl)
use_meqtl = use_meqtl[,-7]

use_snp = unique(use_meqtl$SNP)
# chr_id = 1,5,6,7,8,13,22（chr id which has causal CpG）
snpchr22 = as.data.frame(snpchr22)
snpchr22 = snpchr22[ids]
tmp_snp_data = snpchr22[rownames(snpchr22) %in% use_snp, ]
snp_data = rbind(snp_data,tmp_snp_data)
rm(snpchr22)

# causal gene as outcome
use_gene = unique(causal_gene$Gene)
gene_pos = subset(gene_pos,gene_pos$gene %in% use_gene)
gene_data = gene_df[, colnames(gene_df) %in% c("Subject",gene_pos$probe_id)]

# linear regression
result = data.frame()
for (probe in colnames(gene_data)[-1]){
  for (i in 1:184){
    lmmodel = lm(gene_data[probe][,1]~t(snp_data[i,])+., data=gene_cov)
    result = rbind(result, c(probe, rownames(snp_data[i,]), summary(lmmodel)$coefficients[2,]), stringsAsFactors = F)
  }
}
colnames(result) = c("probe_id","SNP","beta","se","t-value","p-value")

result = left_join(result,gene_pos[,1:2])

snp_info = data.frame(SNP=use_meqtl$SNP,A1=use_meqtl$A1,A2=use_meqtl$A2,MAF=use_meqtl$MAF)
snp_info = unique(snp_info)
result = left_join(result,snp_info)

# normalize beta and se
snp_sd = apply(snp_data,1,sd,na.rm=TRUE)
snp_sd = data.frame(SNP = names(snp_sd), snp_sd = snp_sd)
result = left_join(result,snp_sd)
gene_sd = apply(gene_data[,-1],2,sd,na.rm=TRUE)
gene_sd = data.frame(gene = names(gene_sd), gene_sd = gene_sd)
colnames(gene_sd)[1] = "probe_id"
result = left_join(result,gene_sd)
result$adjust_beta = as.numeric(result$beta)*result$snp_sd/result$gene_sd
result$adjust_se = as.numeric(result$se)*result$snp_sd/result$gene_sd

# filter duplicated
gene_pos_unique = gene_pos[!duplicated(gene_pos$gene),]
result_unique = subset(result,result$probe_id %in% gene_pos_unique$probe_id)

# MR analysis
library(TwoSampleMR)
# ADHD
exposure_dat = format_data(
  subset(use_meqtl,use_meqtl$CpG %in% cpg_adhd),
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$gene %in% gene_adhd),
  type = "outcome",
  phenotype_col = "gene",
  snp_col = "SNP",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_adhd)*length(gene_adhd)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# autism
exposure_dat = format_data(
  subset(use_meqtl,use_meqtl$CpG %in% cpg_autism),
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$gene %in% gene_autism),
  type = "outcome",
  phenotype_col = "gene",
  snp_col = "SNP",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_autism)*length(gene_autism)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# bipolar
exposure_dat = format_data(
  subset(use_meqtl,use_meqtl$CpG %in% cpg_bipolar),
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$gene %in% gene_bipolar),
  type = "outcome",
  phenotype_col = "gene",
  snp_col = "SNP",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_bipolar)*length(gene_bipolar)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# schizo
exposure_dat = format_data(
  subset(use_meqtl,use_meqtl$CpG %in% cpg_schizo),
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$gene %in% gene_schizo),
  type = "outcome",
  phenotype_col = "gene",
  snp_col = "SNP",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_schizo)*length(gene_schizo)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# insomnia
exposure_dat = format_data(
  subset(use_meqtl,use_meqtl$CpG %in% cpg_insomnia),
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$gene %in% gene_insomnia),
  type = "outcome",
  phenotype_col = "gene",
  snp_col = "SNP",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_insomnia)*length(gene_insomnia)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# depression
exposure_dat = format_data(
  subset(use_meqtl,use_meqtl$CpG %in% cpg_depression),
  type = "exposure",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "std",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$gene %in% gene_depression),
  type = "outcome",
  phenotype_col = "gene",
  snp_col = "SNP",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_depression)*length(gene_depression)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# eQTL of causal genes as IVs
use_eqtl = subset(clump_result_all,clump_result_all$gene %in% causal_gene$Gene)

use_snp = unique(use_eqtl$rsid)

ids = as.character(ordermethy$ID)

# chr_id = 1,2,3,5,6,7,8,10,11,12,14,15,16,17,19,20,22（chr id which has causal CpG）
snpchr22 = as.data.frame(snpchr22)
snpchr22 = snpchr22[ids]
tmp_snp_data = snpchr22[rownames(snpchr22) %in% use_snp, ]
snp_data = rbind(snp_data,tmp_snp_data)
rm(snpchr22)

# causal cpg as outcome
use_cpg = unique(causal_cpg$CpG)
cpg_data = ordermethy[, colnames(ordermethy) %in% c("ID",use_cpg)]

# linear regression
result = data.frame()
for (cpg in colnames(cpg_data)[-1]){
  for (i in 1:132){
    lmmodel = lm(cpg_data[cpg][,1]~t(snp_data[i,])+., data=ordercvrt_noid)
    result = rbind(result, c(cpg, rownames(snp_data[i,]), summary(lmmodel)$coefficients[2,]), stringsAsFactors = F)
  }
}
colnames(result) = c("CpG","SNP","beta","se","t-value","p-value")

snp_info = data.frame(SNP=use_eqtl$rsid,A1=use_eqtl$A1,A2=use_eqtl$A2,MAF=use_eqtl$MAF)
snp_info = unique(snp_info)
result = left_join(result,snp_info)

# MR analysis
library(TwoSampleMR)
# ADHD
exposure_dat = format_data(
  subset(use_eqtl,use_eqtl$gene %in% gene_adhd),
  type = "exposure",
  phenotype_col = "gene",
  snp_col = "rsid",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$CpG %in% cpg_adhd),
  type = "outcome",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)

# autism
exposure_dat = format_data(
  subset(use_eqtl,use_eqtl$gene %in% gene_autism),
  type = "exposure",
  phenotype_col = "gene",
  snp_col = "rsid",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$CpG %in% cpg_autism),
  type = "outcome",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)

# bipolar
exposure_dat = format_data(
  subset(use_eqtl,use_eqtl$gene %in% gene_bipolar),
  type = "exposure",
  phenotype_col = "gene",
  snp_col = "rsid",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$CpG %in% cpg_bipolar),
  type = "outcome",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_bipolar)*length(gene_bipolar)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# schizo
exposure_dat = format_data(
  subset(use_eqtl,use_eqtl$gene %in% gene_schizo),
  type = "exposure",
  phenotype_col = "gene",
  snp_col = "rsid",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$CpG %in% cpg_schizo),
  type = "outcome",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_schizo)*length(gene_schizo)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# insomnia
exposure_dat = format_data(
  subset(use_eqtl,use_eqtl$gene %in% gene_insomnia),
  type = "exposure",
  phenotype_col = "gene",
  snp_col = "rsid",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$CpG %in% cpg_insomnia),
  type = "outcome",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_insomnia)*length(gene_insomnia)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# depression
exposure_dat = format_data(
  subset(use_eqtl,use_eqtl$gene %in% gene_depression),
  type = "exposure",
  phenotype_col = "gene",
  snp_col = "rsid",
  beta_col = "adjust_beta",
  se_col = "adjust_se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "pval"
)
outcome_dat = format_data(
  subset(result2,result2$CpG %in% cpg_depression),
  type = "outcome",
  phenotype_col = "CpG",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p-value"
)
harmonise_dat = harmonise_data(exposure_dat,outcome_dat,action=2)
mr_model = mr(harmonise_dat)
significant = subset(mr_model,mr_model$pval<0.05/(length(cpg_depression)*length(gene_depression)))
# heterogeneity test and pleiotropy test
positive_cpg = unique(significant$exposure)
positive_gene = unique(significant$outcome)
positive_harmonise_dat = subset(harmonise_dat,exposure %in% positive_cpg & outcome %in% positive_gene)
het_test = mr_heterogeneity(positive_harmonise_dat)
ple_test = mr_pleiotropy_test(positive_harmonise_dat)

# calculate mediation proportion
causal_mediation = read.csv("causal_mediation.CSV")
cpg_effect = data.frame(CpG=causal_cpg$CpG,beta0=causal_cpg$beta,Phenotype=causal_cpg$Phenotype)
gene_effect = data.frame(Gene=causal_gene$Gene,beta2=causal_gene$beta,Phenotype=causal_gene$Phenotype)
cpgtogene_effect = data.frame(CpG=causal_mediation$CpG,Gene=causal_mediation$Gene,Phenotype=causal_mediation$Phenotype,beta1=causal_mediation$beta)
mediation_df = left_join(left_join(cpgtogene_effect,cpg_effect),gene_effect)
mediation_df$mediation_effect = mediation_df$beta1*mediation_df$beta2
mediation_df$direct_effect = mediation_df$beta0-mediation_df$mediation_effect
mediation_df$mediation_proportion = mediation_df$mediation_effect/mediation_df$beta0
