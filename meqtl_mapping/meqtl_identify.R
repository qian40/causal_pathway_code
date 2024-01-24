# load data
load("ordercvrt.Rdata")
load("ordercvrt_noid.Rdata")
load("ordermethy.Rdata")
load("pos.Rdata")
load("cisresult.Rdata")
for (i in 1:22){
  load(paste0("snpchr", i, ".Rdata"))
}
snps_list = list(snpchr1,snpchr2,snpchr3,snpchr4,snpchr5,snpchr6,snpchr7,snpchr8,snpchr9,snpchr10,snpchr11,snpchr12,snpchr13,snpchr14,snpchr15,snpchr16,snpchr17,snpchr18,snpchr19,snpchr20,snpchr21,snpchr22)
names(snps_list) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
rm(snpchr1,snpchr2,snpchr3,snpchr4,snpchr5,snpchr6,snpchr7,snpchr8,snpchr9,snpchr10,snpchr11,snpchr12,snpchr13,snpchr14,snpchr15,snpchr16,snpchr17,snpchr18,snpchr19,snpchr20,snpchr21,snpchr22)

all_data = merge(merge(cvrt, methy, by="ID"), snp1, by="ID")
snps = snps$ColumnSubsample(colnames(snps) %in% filter_id)

# save residual
res_data = all_data

# compute lm residual: methy~cvrt
for (i in 23:372603){
  lmmodel = lm(all_data[,i]~all_data[,2]+all_data[,3]+all_data[,4]+
                 all_data[,5]+all_data[,6]+all_data[,7]+all_data[,8]+
                 all_data[,9]+all_data[,10]+all_data[,11]+all_data[,12]+
                 all_data[,13]+all_data[,14]+all_data[,15]+all_data[,16]+
                 all_data[,17]+all_data[,18]+all_data[,19]+all_data[,20]+
                 all_data[,21]+all_data[,22])
  res = residuals(lmmodel)
  res_data[,i] = res
  if (i%%10000 == 0){
    print(i)
  }
}

# find cis, threshold=1e6
findcis = function(cpg_name){
  tmp_cpg = cpg_pos[cpg_pos$geneid==cpg_name,]
  chr_num = tmp_cpg$chr
  tmp_cpg_pos = tmp_cpg$left
  tmp_snp = snp_pos[snp_pos$chr==chr_num & abs(snp_pos$pos-tmp_cpg_pos)<1e6,]
  cis_snpid = tmp_snp$snpid
  cisdata = snps_list[[chr_num]][rownames(snps_list[[chr_num]]) %in% cis_snpid, ]
  if (is.null(dim(cisdata))){
    cisdata = as.data.frame(t(cisdata))
    rownames(cisdata) = rownames(snps_list[[chr_num]])[rownames(snps_list[[chr_num]]) %in% cis_snpid]
  }
  cisdata = as.data.frame(t(cisdata))
  return(cisdata)
}

# linear regression to find cis-meQTL, results in cis_result
cismeqtl = function(){
  cis_result = data.frame()
  threshold = 1e-6
  for (i in 1:372603){
    cpg_name = colnames(res_data)[i]
    tmp_cis = findcis(cpg_name)
    if (dim(tmp_cis)[2] == 0){
      next
    }
    for (j in 1:dim(tmp_cis)[2]){
      lmmodel = lm(res_data[,i]~tmp_cis[,j])
      if (summary(lmmodel)$coefficients[2,4] > threshold){
        next
      }
      cis = colnames(tmp_cis)[j]      
      cis_result = rbind(cis_result, c(cpg_name, cis, summary(lmmodel)$coefficients[2,4]), stringsAsFactors = F)
    }
    if (i %% 100 == 0){
      print(i)
    }
    if (i %% 10000 == 0){
      save(cis_result, i, file="cis_result2.Rdata")
    }
  }
  colnames(cis_result) = c("CpG","SNP", "p-value")
}

# recalculate
cisfinal = data.frame()
threshold = 3.1e-11
for (i in 1:15689116){
  cpg = cisresult[i,1]
  snp = cisresult[i,2]
  chr_num = cpg_pos[cpg_pos$geneid==cpg,]$chr
  lmmodel = lm(ordermethy[cpg][,1]~t(subset(snps_list[[chr_num]], rownames(snps_list[[chr_num]])==snp))+., data=ordercvrt_noid)
  if (summary(lmmodel)$coefficients[2,4] > threshold){
    next
  }
  cisfinal = rbind(cisfinal, c(cpg, snp, chr_num, summary(lmmodel)$coefficients[2,4]), stringsAsFactors = F)
}
