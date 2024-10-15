causal_gene = unique(read.csv("causal_gene_result.CSV")$Gene)
causal_cpg = unique(read.csv("causal_cpg_result.CSV")$CpG)
mediation = read.csv("causal_mediation.CSV")
load("meqtl.RData")

# colocalization
library(coloc)
use_eqtl = subset(eqtl,eqtl$gene %in% causal_gene)
use_meqtl = subset(my_meqtl,my_meqtl$CpG %in% causal_cpg)

# ATXN1
ATXN1_eqtl = subset(use_eqtl,use_eqtl$gene=="ATXN1")
use_ATXN1_eqtl = data.frame(gene_name=ATXN1_eqtl$gene,SNP=ATXN1_eqtl$rsid,beta=ATXN1_eqtl$beta,se=ATXN1_eqtl$std,pval=ATXN1_eqtl$`p-value`)
rm(ATXN1_eqtl)
ATXN1_cpg = mediation$CpG[which(mediation$Gene=="ATXN1")]
coloc_snp1 = data.frame()
for (i in 1:1){
  tmp_cpg = ATXN1_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_ATXN1_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=545, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp1 = rbind(coloc_snp1,tmp_result$results %>% filter(SNP.PP.H4>0.95))
  print(i)
  print(dim(coloc_snp1))
}

# PGBD1
PGBD1_eqtl = subset(use_eqtl,use_eqtl$gene=="PGBD1")
use_PGBD1_eqtl = data.frame(gene_name=PGBD1_eqtl$gene,SNP=PGBD1_eqtl$rsid,beta=PGBD1_eqtl$beta,se=PGBD1_eqtl$std,pval=PGBD1_eqtl$`p-value`)
rm(PGBD1_eqtl)
PGBD1_cpg = mediation$CpG[which(mediation$Gene=="PGBD1")]
coloc_snp1 = data.frame()
for (i in 1:6){
  tmp_cpg = PGBD1_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_PGBD1_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
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

# MRPL2
MRPL2_eqtl = subset(use_eqtl,use_eqtl$gene=="MRPL2")
use_MRPL2_eqtl = data.frame(gene_name=MRPL2_eqtl$gene,SNP=MRPL2_eqtl$rsid,beta=MRPL2_eqtl$beta,se=MRPL2_eqtl$std,pval=MRPL2_eqtl$`p-value`)
rm(MRPL2_eqtl)
MRPL2_cpg = mediation$CpG[which(mediation$Gene=="MRPL2")]
coloc_snp1 = data.frame()
for (i in 1:2){
  tmp_cpg = MRPL2_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_MRPL2_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=545, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp1 = rbind(coloc_snp1,tmp_result$results %>% filter(SNP.PP.H4>0.90))
  print(i)
  print(dim(coloc_snp1))
}

# MAD1L1
MAD1L1_eqtl = subset(use_eqtl,use_eqtl$gene=="MAD1L1")
use_MAD1L1_eqtl = data.frame(gene_name=MAD1L1_eqtl$gene,SNP=MAD1L1_eqtl$rsid,beta=MAD1L1_eqtl$beta,se=MAD1L1_eqtl$std,pval=MAD1L1_eqtl$`p-value`)
rm(MAD1L1_eqtl)
MAD1L1_cpg = mediation$CpG[which(mediation$Gene=="MAD1L1")]
coloc_snp1 = data.frame()
for (i in 1:1){
  tmp_cpg = MAD1L1_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_MAD1L1_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=545, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp1 = rbind(coloc_snp1,tmp_result$results %>% filter(SNP.PP.H4>0.90))
  print(i)
  print(dim(coloc_snp1))
}

# NT5DC1
NT5DC1_eqtl = subset(use_eqtl,use_eqtl$gene=="NT5DC1")
use_NT5DC1_eqtl = data.frame(gene_name=NT5DC1_eqtl$gene,SNP=NT5DC1_eqtl$rsid,beta=NT5DC1_eqtl$beta,se=NT5DC1_eqtl$std,pval=NT5DC1_eqtl$`p-value`)
rm(NT5DC1_eqtl)
NT5DC1_cpg = mediation$CpG[which(mediation$Gene=="NT5DC1")]
coloc_snp1 = data.frame()
for (i in 1:1){
  tmp_cpg = NT5DC1_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_NT5DC1_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=545, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp1 = rbind(coloc_snp1,tmp_result$results %>% filter(SNP.PP.H4>0.90))
  print(i)
  print(dim(coloc_snp1))
}

# HLA-DRB1
HLA_DRB1_eqtl = subset(use_eqtl,use_eqtl$gene=="HLA-DRB1")
use_HLA_DRB1_eqtl = data.frame(gene_name=HLA_DRB1_eqtl$gene,SNP=HLA_DRB1_eqtl$rsid,beta=HLA_DRB1_eqtl$beta,se=HLA_DRB1_eqtl$std,pval=HLA_DRB1_eqtl$`p-value`)
rm(HLA_DRB1_eqtl)
HLA_DRB1_cpg = mediation$CpG[which(mediation$Gene=="HLA-DRB1")]
coloc_snp1 = data.frame()
for (i in 1:19){
  tmp_cpg = HLA_DRB1_cpg[i]
  tmp_meqtl = subset(use_meqtl,use_meqtl$CpG==tmp_cpg)
  input = merge(tmp_meqtl,use_HLA_DRB1_eqtl,by="SNP",all=FALSE,suffixes=c("_meqtl","_eqtl"))
  if (dim(input)[1]==0){
    next
  }
  tmp_result = coloc.abf(dataset1=list(pvalues=as.numeric(input$`p-value`), type="quant", beta=as.numeric(input$beta_meqtl), varbeta=as.numeric(input$std)^2, N=1274, snp=input$SNP), 
                         dataset2=list(pvalues=as.numeric(input$pval), type="quant", beta=as.numeric(input$beta_eqtl), varbeta=as.numeric(input$se)^2, N=545, snp=input$SNP), 
                         MAF=input$MAF, p1=2e-11, p2=1e-7, p12=2e-12)
  coloc_snp1 = rbind(coloc_snp1,tmp_result$results %>% filter(SNP.PP.H4>0.90))
  print(i)
  print(dim(coloc_snp1))
}
