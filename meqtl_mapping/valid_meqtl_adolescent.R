# load data
load("./data/meqtl.Rdata")
valid_meqtl = read.table("./data/15up.ALL.M.tab", header=TRUE)
snp_info = read.table("./data/ariesmqtlsnps.bim", header=FALSE)
colnames(snp_info) = c("chr","SNP","none","mapinfo","A1","A2")
valid_meqtl$pair = paste(valid_meqtl$gene,valid_meqtl$SNP,sep="_")
valid_meqtl7 = valid_meqtl
valid_meqtl = subset(valid_meqtl,`p.value`<1e-14)
valid_meqtl = left_join(valid_meqtl,select(snp_info,SNP,A1,A2))

# filter replicated data
repli_pair = intersect(valid_meqtl$pair,my_meqtl$pair)
repli_data = subset(my_meqtl,pair %in% repli_pair)
repli_valid_data = subset(valid_meqtl,pair %in% repli_pair)

# compare beta direction
use_repli_data = select(repli_data,c("pair","beta","A1","A2","tvalue"))
use_repli_valid_data = select(repli_valid_data,c("pair","beta","A1","A2","t.stat"))
colnames(use_repli_valid_data)[2:5] = c("beta_valid","valid_A1","valid_A2","tvalue_valid")
all_repli_data = left_join(use_repli_data,use_repli_valid_data)
reverse_index = (all_repli_data$A1==all_repli_data$valid_A2 & all_repli_data$A2==all_repli_data$valid_A1)
reverse_index[is.na(reverse_index)] = FALSE
obverse_index = (all_repli_data$A1==all_repli_data$valid_A1 & all_repli_data$A2==all_repli_data$valid_A2)
obverse_index[is.na(obverse_index)] = FALSE
other_index = !(reverse_index|obverse_index)
all_repli_data$beta_adjust = all_repli_data$beta_valid
all_repli_data[reverse_index,]$beta_adjust = -all_repli_data[reverse_index,]$beta_adjust
all_repli_data$tvalue_adjust = all_repli_data$tvalue_valid
all_repli_data[reverse_index,]$tvalue_adjust = -all_repli_data[reverse_index,]$tvalue_adjust
all_repli_data$tmp = as.numeric(all_repli_data$beta)*as.numeric(all_repli_data$beta_adjust)
all_repli_data$tmp2 = as.numeric(all_repli_data$tvalue)*as.numeric(all_repli_data$tvalue_adjust)
final_repli_data = subset(all_repli_data,!other_index)
sum(final_repli_data$tmp>0)
sum(final_repli_data$tmp2>0)

# scatter plot
set.seed(1)
plot_index = sample(1:1382244,100000)
# plot beta
plot(final_repli_data$beta[plot_index],final_repli_data$beta_adjust[plot_index],pch=16,cex=0.5,xlab="our_beta",ylab="valid_beta")
abline(h=0)
abline(v=0)
abline(a=0,b=1,col="red")
cor(as.numeric(final_repli_data$beta),as.numeric(final_repli_data$beta_adjust))
# plot t-value
plot(final_repli_data$tvalue[plot_index],final_repli_data$tvalue_adjust[plot_index],pch=16,cex=0.5,xlab="our_tval",ylab="valid_tval")
abline(h=0)
abline(v=0)
abline(a=0,b=1,col="red")
cor(as.numeric(final_repli_data$tvalue),as.numeric(final_repli_data$tvalue_adjust))
