# gene position annotation
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_names = unique(gene_map$SYMBOL)
gene_info = getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position", "end_position"),
                   filters = "external_gene_name",
                   values = gene_names,
                   mart = ensembl)

# filter
gene_location = subset(gene_info,gene_info$chromosome_name %in% 1:22)
gene_location_tmp = subset(gene_info,!gene_info$chromosome_name %in% c(1:22,"X","Y"))
gene_location_tmp$chromosome_name = gsub("CHR_HSCHR(\\d+).*", "\\1", gene_location_tmp$chromosome_name)
gene_location = rbind(gene_location, subset(gene_location_tmp,gene_location_tmp$chromosome_name %in% 1:22))
gene_location_tmp2 = data.frame("gene"=gene_location$hgnc_symbol,"chr"=gene_location$chromosome_name,"pos"=gene_location$start_position)
gene_location_tmp3 = subset(gene_location_tmp2,gene_location_tmp2$gene!="")
gene_location_tmp4 = gene_location_tmp3[!duplicated(gene_location_tmp3$gene),]

# combine
gene_map_tmp = data.frame("probe_id"=gene_map$..PROBE_ID,"gene"=gene_map$SYMBOL)
gene_map_tmp2 = subset(gene_map_tmp,gene_map_tmp$gene!="")
gene_pos = left_join(gene_map_tmp2,gene_location_tmp4)
gene_pos = na.omit(gene_pos)
save(gene_pos,file="gene_pos.Rdata")

# filter residual
gene_res = gene_res_data.filter[c("Subject",gene_pos$probe_id)]
ids = as.character(gene_res$Subject)

snp_pos$chr = as.character(gsub("chr","",snp_pos$chr))

# eqtl of chr 1
cis_result = data.frame()
load("./snpchr/snpchr1.Rdata")
load("./prefilter.Rdata")
snpchr1 = as.data.frame(snpchr1)
use_snp = snpchr1[ids]
rm(snpchr1)
gene_chr_i = subset(gene_pos,gene_pos$chr==1)
use_gene = gene_res[gene_chr_i$probe_id]
rownames(use_gene) = gene_res$Subject
snp_chr_i = subset(snp_pos,snp_pos$chr=="1")
for (i in 1:dim(gene_chr_i)[1]){
  tmp_probe = gene_chr_i$probe_id[i]
  tmp_gene_pos = gene_chr_i$pos[i]
  tmp_cis_snp = subset(snp_chr_i,abs(snp_chr_i$pos-tmp_gene_pos)<1e6)
  tmp_cis_snp_id = tmp_cis_snp$snpid
  # tmp_cisdata = use_snp[tmp_cis_snp_id,]
  tmp_cisdata = subset(use_snp,rownames(use_snp) %in% tmp_cis_snp_id)
  tmp_cisdata = as.data.frame(t(tmp_cisdata))
  print(i)
  if (dim(tmp_cisdata)[2] == 0){
    print("next")
    next
  }
  for (j in 1:dim(tmp_cisdata)[2]){
    lmmodel = lm(use_gene[,i]~tmp_cisdata[,j])
    if (summary(lmmodel)$coefficients[2,4] > 1e-3){
      next
    }
    cis_result = rbind(cis_result, c(tmp_probe, colnames(tmp_cisdata)[j], summary(lmmodel)$coefficients[2,]), stringsAsFactors = F)
  }
}
save(cis_result,file="cis_result1.Rdata")

# combine
load("cis_result1.Rdata")
eqtl_my = cis_result
colnames(eqtl_my) = c("probe","rsid","beta","se","t-value","p-value")
for (i in 2:22){
  load(paste0("cis_result",i,".Rdata"))
  print(dim(cis_result))
  colnames(cis_result) = c("probe","rsid","beta","se","t-value","p-value")
  eqtl_my = rbind(eqtl_my,cis_result)
}
dim(subset(eqtl_my,eqtl_my$`p-value`<1e-4))
probe_gene = data.frame("probe"=gene_pos$probe_id,"gene"=gene_pos$gene)
eqtl_my = left_join(eqtl_my,probe_gene)
eqtl_my$pair = paste(eqtl_my$gene,eqtl_my$rsid,sep="_")
eqtl$pair = paste(eqtl$gene_name,eqtl$rsid,sep="_")
eqtl_dupli = subset(eqtl_my,duplicated(eqtl_my$pair))

# preprocess
gene_cov = subset(gene_cov,gene_cov$Subject %in% ids)
rownames(gene_cov) = gene_cov$Subject
gene_cov = gene_cov[-1]
gene_df = subset(gene_df,gene_df$Subject %in% ids)

# recalculate
load("eqtl_mapping/cis_result1.Rdata")
load("snpchr/snpchr1.Rdata")
load("recalculate.Rdata")

cis_final = data.frame()
chr = "1"
snpchr1 = as.data.frame(snpchr1)
snpchr1 = snpchr1[ids]
for (i in 1:dim(cis_result)[1]){
  if (i %% 100 == 0){
    print(paste0(i,"/",dim(cis_result)[1]))
  }
  probe = cis_result[i,1]
  snp = cis_result[i,2]
  lmmodel = lm(gene_df[probe][,1]~t(subset(snpchr1, rownames(snpchr1)==snp))+., data=gene_cov)
  cis_final = rbind(cis_final, c(probe, snp, chr, summary(lmmodel)$coefficients[2,]), stringsAsFactors = F)
}
colnames(cis_final) = c("probe","SNP","chr","beta","std","t-value","p-value")
save(cis_final, file="eqtl_mapping/cis_final1.Rdata")

# combine
load("eqtl_mapping/cis_final1.Rdata")
eqtl_my = cis_final
for (i in 2:22){
  load(paste0("eqtl_mapping/cis_final",i,".Rdata"))
  print(dim(cis_final))
  eqtl_my = rbind(eqtl_my,cis_final)
}
eqtl_my$beta = as.numeric(eqtl_my$beta)
eqtl_my$std = as.numeric(eqtl_my$std)
eqtl_my$`t-value` = as.numeric(eqtl_my$`t-value`)
eqtl_my$`p-value` = as.numeric(eqtl_my$`p-value`)
eqtl_my = subset(eqtl_my,eqtl_my$`p-value`<1e-3)
eqtl_my = left_join(eqtl_my,probe_gene)
eqtl_my$pair = paste(eqtl_my$gene,eqtl_my$SNP,sep="_")
eqtl_valid = data.frame(eqtl$gene_name,eqtl$rsid,eqtl$pair,eqtl$slope,eqtl$slope_se,eqtl$pval_nominal,eqtl$variant_id)
colnames(eqtl_valid) = c("gene","rsid","pair","beta","se","p-value","variant_id")
colnames(eqtl_my)[2] = "rsid"
colnames(eqtl_my)[5] = "se"
