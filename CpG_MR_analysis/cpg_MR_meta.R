meta_df = read.csv("causal_cpg_result.CSV")
sum(is.na(meta_df$beta_validation))
meta_df$tmp = meta_df$beta*meta_df$beta_validation
sum(meta_df$tmp>0,na.rm = TRUE)

sig_gene = meta_df
colnames(sig_gene)[6:8] = c("imagen_beta","imagen_se","imagen_p")
colnames(sig_gene)[12:14] = c("gtex_beta","gtex_se","gtex_p")
sig_gene$imagen_z = sig_gene$imagen_beta/sig_gene$imagen_se
sig_gene$gtex_z = sig_gene$gtex_beta/sig_gene$gtex_se
sig_gene$imagen_nobs = 1274
sig_gene$gtex_nobs = 1731

sig_gene=sig_gene%>%
  mutate(beta_meta=(gtex_beta/(gtex_se**2)+imagen_beta/(imagen_se**2))/(1/gtex_se**2+1/imagen_se**2),
         z_meta=(gtex_z*gtex_nobs+imagen_z*imagen_nobs)/sqrt(gtex_nobs**2+imagen_nobs**2),
         n_meta=gtex_nobs+imagen_nobs)%>%
  mutate(p_meta=case_when(z_meta>0~1-pnorm(z_meta),
                          z_meta<0~pnorm(z_meta)))

sum(sig_gene$p_meta<sig_gene$imagen_p,na.rm = TRUE)
