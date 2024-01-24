# load data
load("eqtl_mapping/eqtl_valid.Rdata")
load("eqtl_mapping/eqtl_my_FDRcorrect.Rdata")
# replication
replicate_pair = intersect(eqtl_my$pair,eqtl_valid$pair)
replicate_eqtl_my = subset(eqtl_my,eqtl_my$pair %in% replicate_pair)
replicate_eqtl_my$my_beta = replicate_eqtl_my$adjust_beta
replicate_eqtl_valid = subset(eqtl_valid,eqtl_valid$pair %in% replicate_pair)
# remove duplicated data
replicate_eqtl_valid = subset(replicate_eqtl_valid,!(replicate_eqtl_valid$pair=="RSPH3_rs2143653"&replicate_eqtl_valid$A2=="G"))

# calculate tvalue
replicate_eqtl_my$tvalue = as.numeric(replicate_eqtl_my$adjust_beta)/replicate_eqtl_my$adjust_se
replicate_eqtl_valid$tvalue_valid = as.numeric(replicate_eqtl_valid$beta)/as.numeric(replicate_eqtl_valid$se)
repli_data = replicate_eqtl_my
repli_valid_data = replicate_eqtl_valid
# compare beta direction
use_repli_data = select(repli_data,c("pair","adjust_beta","A1","A2","tvalue"))
use_repli_valid_data = select(repli_valid_data,c("pair","beta","A1","A2","tvalue_valid"))
colnames(use_repli_data)[2] = "beta"
colnames(use_repli_valid_data)[2:5] = c("beta_valid","valid_A1","valid_A2","tvalue_valid")
all_repli_data = left_join(use_repli_data,use_repli_valid_data)
reverse_index = (all_repli_data$A1==all_repli_data$valid_A2 & all_repli_data$A2==all_repli_data$valid_A1)
reverse_index[is.na(reverse_index)] = FALSE
obverse_index = (all_repli_data$A1==all_repli_data$valid_A1 & all_repli_data$A2==all_repli_data$valid_A2)
obverse_index[is.na(obverse_index)] = FALSE
other_index = !(reverse_index|obverse_index)
all_repli_data$beta_adjust = all_repli_data$beta_valid
all_repli_data[obverse_index,]$beta_adjust = -all_repli_data[obverse_index,]$beta_adjust
all_repli_data$tvalue_adjust = all_repli_data$tvalue_valid
all_repli_data[obverse_index,]$tvalue_adjust = -all_repli_data[obverse_index,]$tvalue_adjust
all_repli_data$tmp = as.numeric(all_repli_data$beta)*as.numeric(all_repli_data$beta_adjust)
all_repli_data$tmp2 = as.numeric(all_repli_data$tvalue)*as.numeric(all_repli_data$tvalue_adjust)
final_repli_data = subset(all_repli_data,!other_index)
sum(final_repli_data$tmp>0)
sum(final_repli_data$tmp2>0)

# scatter plot
set.seed(1)
plot_index = sample(1:263403,100000)
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
