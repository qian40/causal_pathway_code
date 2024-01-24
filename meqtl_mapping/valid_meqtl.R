# load data
load("./data/cis_cosmo.RData")
load("./data/meqtl.Rdata")

# divide ancestry
cis_EUR = subset(cis,cis$discovery.population=="Eur"|cis$discovery.population=="Both")
cis_SA = subset(cis,cis$discovery.population=="SA"|cis$discovery.population=="Both")
valid_meqtl = select(cis_EUR,c("snp","snp.chr","A1","A2","eaf.eu","cpg","beta.eur","se.eur","p.eur"))
# valid_meqtl = select(cis_SA,c("snp","snp.chr","A1","A2","eaf.sa","cpg","beta.sa","se.sa","p.sa"))
colnames(valid_meqtl) = c("snp","chr","A1","A2","eaf","cpg","beta","se","p")
valid_meqtl$pair = paste(valid_meqtl$cpg,valid_meqtl$snp,sep="_")

# calculate t-value
my_meqtl$tvalue = as.numeric(my_meqtl$beta)/as.numeric(my_meqtl$std)
valid_meqtl$tvalue_valid = as.numeric(valid_meqtl$beta)/as.numeric(valid_meqtl$se)

# filter replicated data
index = my_meqtl$pair %in% valid_meqtl$pair
repli_pair = my_meqtl$pair[index]
repli_data = subset(my_meqtl,pair %in% repli_pair)
repli_valid_data = subset(valid_meqtl,pair %in% repli_pair)

# compare beta direction
use_repli_data = select(repli_data,c("pair","beta","A1","A2","MAF","tvalue"))
use_repli_valid_data = select(repli_valid_data,c("pair","beta","A1","A2","eaf","tvalue_valid"))
colnames(use_repli_valid_data)[2:4] = c("beta_valid","valid_A1","valid_A2")
all_repli_data = left_join(use_repli_data,use_repli_valid_data)
reverse_index = (all_repli_data$A1==all_repli_data$valid_A2 & all_repli_data$A2==all_repli_data$valid_A1)
obverse_index = (all_repli_data$A1==all_repli_data$valid_A1 & all_repli_data$A2==all_repli_data$valid_A2)
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
set.seed(2)
plot_index = sample(1:2957864,100000)
# plot_index = sample(1:3621033,100000)
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