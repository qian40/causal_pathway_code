causal_gene = unique(read.csv("causal_gene_result.CSV")$Gene)
causal_cpg = unique(read.csv("causal_cpg_result.CSV")$CpG)
mediation = read.csv("causal_mediation.CSV")

# colocalization
library(coloc)
use_eqtl = subset(eqtl,eqtl$gene %in% causal_gene)
use_meqtl = subset(my_meqtl,my_meqtl$CpG %in% causal_cpg)

# BTN3A2
BTN3A2_eqtl = subset(use_eqtl,use_eqtl$gene=="BTN3A2")
use_BTN3A2_eqtl = data.frame(gene_name=BTN3A2_eqtl$gene,SNP=BTN3A2_eqtl$rsid,beta=BTN3A2_eqtl$adjust_beta,se=BTN3A2_eqtl$adjust_se,pval=BTN3A2_eqtl$`p-value`)
rm(BTN3A2_eqtl)
BTN3A2_cpg = mediation$CpG[which(mediation$Gene=="BTN3A2")]
coloc_snp1 = data.frame()
for (i in 1:27){
  tmp_cpg = BTN3A2_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_BTN3A2_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=591, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp1 = rbind(coloc_snp1,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(dim(coloc_snp1))
}

# IER3
IER3_eqtl = subset(use_eqtl,use_eqtl$gene=="IER3")
use_IER3_eqtl = data.frame(gene_name=IER3_eqtl$gene,SNP=IER3_eqtl$rsid,beta=IER3_eqtl$adjust_beta,se=IER3_eqtl$adjust_se,pval=IER3_eqtl$`p-value`)
rm(IER3_eqtl)
IER3_cpg = mediation$CpG[which(mediation$Gene=="IER3")]
coloc_snp2 = data.frame()
for (i in 1:18){
  tmp_cpg = IER3_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_IER3_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=591, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp2 = rbind(coloc_snp2,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(tmp_cpg)
  print(dim(coloc_snp2))
}

# HLA-DPB1
HLA_DPB1_eqtl = subset(use_eqtl,use_eqtl$gene=="HLA-DPB1")
use_HLA_DPB1_eqtl = data.frame(gene_name=HLA_DPB1_eqtl$gene,SNP=HLA_DPB1_eqtl$rsid,beta=HLA_DPB1_eqtl$adjust_beta,se=HLA_DPB1_eqtl$adjust_se,pval=HLA_DPB1_eqtl$`p-value`)
rm(HLA_DPB1_eqtl)
HLA_DPB1_cpg = mediation$CpG[which(mediation$Gene=="HLA-DPB1")]
coloc_snp3 = data.frame()
for (i in 1:15){
  tmp_cpg = HLA_DPB1_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_HLA_DPB1_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=591, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp3 = rbind(coloc_snp3,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(tmp_cpg)
  print(dim(coloc_snp3))
}

# HCG9
HCG9_eqtl = subset(use_eqtl,use_eqtl$gene=="HCG9")
use_HCG9_eqtl = data.frame(gene_name=HCG9_eqtl$gene,SNP=HCG9_eqtl$rsid,beta=HCG9_eqtl$adjust_beta,se=HCG9_eqtl$adjust_se,pval=HCG9_eqtl$`p-value`)
rm(HCG9_eqtl)
HCG9_cpg = mediation$CpG[which(mediation$Gene=="HCG9")]
coloc_snp4 = data.frame()
for (i in 1:3){
  tmp_cpg = HCG9_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_HCG9_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=591, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp4 = rbind(coloc_snp4,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(tmp_cpg)
  print(dim(coloc_snp4))
}

# CCS
CCS_eqtl = subset(use_eqtl,use_eqtl$gene=="CCS")
use_CCS_eqtl = data.frame(gene_name=CCS_eqtl$gene,SNP=CCS_eqtl$rsid,beta=CCS_eqtl$adjust_beta,se=CCS_eqtl$adjust_se,pval=CCS_eqtl$`p-value`)
rm(CCS_eqtl)
CCS_cpg = mediation$CpG[which(mediation$Gene=="CCS")]
coloc_snp5 = data.frame()
for (i in 1:1){
  tmp_cpg = CCS_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_CCS_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=591, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp5 = rbind(coloc_snp5,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(tmp_cpg)
  print(dim(coloc_snp5))
}

# MRPL2
MRPL2_eqtl = subset(use_eqtl,use_eqtl$gene=="MRPL2")
use_MRPL2_eqtl = data.frame(gene_name=MRPL2_eqtl$gene,SNP=MRPL2_eqtl$rsid,beta=MRPL2_eqtl$adjust_beta,se=MRPL2_eqtl$adjust_se,pval=MRPL2_eqtl$`p-value`)
rm(MRPL2_eqtl)
MRPL2_cpg = mediation$CpG[which(mediation$Gene=="MRPL2")]
coloc_snp6 = data.frame()
for (i in 1:1){
  tmp_cpg = MRPL2_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_MRPL2_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=591, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp6 = rbind(coloc_snp6,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(tmp_cpg)
  print(dim(coloc_snp6))
}

# TRIM10
TRIM10_eqtl = subset(use_eqtl,use_eqtl$gene=="TRIM10")
use_TRIM10_eqtl = data.frame(gene_name=TRIM10_eqtl$gene,SNP=TRIM10_eqtl$rsid,beta=TRIM10_eqtl$adjust_beta,se=TRIM10_eqtl$adjust_se,pval=TRIM10_eqtl$`p-value`)
rm(TRIM10_eqtl)
TRIM10_cpg = mediation$CpG[which(mediation$Gene=="TRIM10")]
coloc_snp7 = data.frame()
for (i in 1:1){
  tmp_cpg = TRIM10_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_TRIM10_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=591, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp7 = rbind(coloc_snp7,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(tmp_cpg)
  print(dim(coloc_snp7))
}

# HLA_DQA1
HLA_DQA1_eqtl = subset(use_eqtl,use_eqtl$gene=="HLA-DQA1")
use_HLA_DQA1_eqtl = data.frame(gene_name=HLA_DQA1_eqtl$gene,SNP=HLA_DQA1_eqtl$rsid,beta=HLA_DQA1_eqtl$adjust_beta,se=HLA_DQA1_eqtl$adjust_se,pval=HLA_DQA1_eqtl$`p-value`)
rm(HLA_DQA1_eqtl)
HLA_DQA1_cpg = mediation$CpG[which(mediation$Gene=="HLA-DQA1")]
coloc_snp8 = data.frame()
for (i in 1:21){
  tmp_cpg = HLA_DQA1_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_HLA_DQA1_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=591, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp8 = rbind(coloc_snp8,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(tmp_cpg)
  print(dim(coloc_snp8))
}
