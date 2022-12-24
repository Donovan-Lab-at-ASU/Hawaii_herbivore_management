# -----------------------------------------------------------------------
# Spatial summaries of herbivore biomass in Hawai‘i
# -----------------------------------------------------------------------
# Model and code by Mary Donovan - marydonovan@asu.edu
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

# pulls in posteriors from full models and calculates summaries for various spatial groups (moku, mma, etc)

library(rjags)
library(dplyr)
library(raster)

# predictor matrix - rows corresponds to the 100x100m pixels within the study domain
predictor_full_wMMA_noNA <- readRDS(file="Data/full_predictor_matrix_wMMA_noNA_sand.Rdata")

# model files
mod_out <- readRDS("model_out/TotalHerbivore_gamma_model_out.Rdata")

# moku means
grepgo <- grep('b0.moku',colnames(mod_out[[1]]))
b0moku_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])

# betas
grepgo <- grep('beta',colnames(mod_out[[1]]))
beta_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
beta_matrix <- as.matrix(beta_post)

# 3 moku with no data, but we still want to calculate a posterior, so need to update the index and add columns of zeros to the moku intercept matrix (and thus giving those intercepts a value of zero)
b0moku_post <- cbind(b0moku_post,rep(0,15000),rep(0,15000),rep(0,15000))
predictor_full_wMMA_noNA$moku_ind <- ifelse(predictor_full_wMMA_noNA$moku_ind == -36, 41, predictor_full_wMMA_noNA$moku_ind)
predictor_full_wMMA_noNA$moku_ind <- ifelse(predictor_full_wMMA_noNA$moku_ind == -35, 42, predictor_full_wMMA_noNA$moku_ind)
predictor_full_wMMA_noNA$moku_ind <- ifelse(predictor_full_wMMA_noNA$moku_ind == -13, 43, predictor_full_wMMA_noNA$moku_ind)

# Remove a few rows for pixels determined to be not fit in our domain (estuaries, on land, etc) [done by hand by J Lecky Oct 2022]
predictor_full_wMMA_noNA <- predictor_full_wMMA_noNA %>% filter(is.na(bad_pt))

# by moku -----------------------------------------------------------------
moku_summary <- data.frame(moku_ind=unique(predictor_full_wMMA_noNA$moku_ind),
                           mean=NA,up80=NA,down80=NA,up95=NA,down95=NA,up50=NA,down50=NA,
                           mean_noSand=NA,up80_noSand=NA,down80_noSand=NA,up95_noSand=NA,down95_noSand=NA,up50_noSand=NA,down50_noSand=NA)

for(k in 1:nrow(moku_summary)){
  # subset predictor matrix to moku
  predictor_full_sub <- predictor_full_wMMA_noNA %>% filter(moku_ind==k)
  pred_post_full <- matrix(NA,15000,nrow(predictor_full_sub))
  
  # calculate the posterior
  for(i in 1:nrow(predictor_full_sub)){
    pred_post_full[,i] <- b0moku_post[,predictor_full_sub$moku_ind[i]] + beta_matrix %*% t(as.matrix(predictor_full_sub[i,1:28]))
  }
  
  # summarise posterior
  moku_summary$mean[moku_summary$moku_ind==k] <- exp(mean(pred_post_full))
  moku_summary$down80[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full,0.10))
  moku_summary$up80[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full,0.90))
  moku_summary$down95[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full,0.025))
  moku_summary$up95[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full,0.975))
  moku_summary$down50[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full,0.25))
  moku_summary$up50[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full,0.75))
  
  pred_post_full_noSand <- pred_post_full
  pred_post_full_noSand[,which(predictor_full_sub$HIMARC_hard_soft==1)] <- NA
  
  moku_summary$mean_noSand[moku_summary$moku_ind==k] <- exp(mean(pred_post_full_noSand,na.rm=T))
  moku_summary$down80_noSand[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full_noSand,0.10,na.rm=T))
  moku_summary$up80_noSand[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full_noSand,0.90,na.rm=T))
  moku_summary$down95_noSand[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full_noSand,0.025,na.rm=T))
  moku_summary$up95_noSand[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full_noSand,0.975,na.rm=T))
  moku_summary$down50_noSand[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full_noSand,0.25,na.rm=T))
  moku_summary$up50_noSand[moku_summary$moku_ind==k] <- exp(quantile(pred_post_full_noSand,0.75,na.rm=T))
  
}

write.csv(moku_summary, 'outputs/posterior_summary_moku_TotalHerbivores.csv',row.names=F)


# mma_comparison ----------------------------------------------------------

mlcd_allMHI <- data.frame(mean=NA,up80=NA,down80=NA,up95=NA,down95=NA,up50=NA,down50=NA,
                          mean_noSand=NA,up80_noSand=NA,down80_noSand=NA,up95_noSand=NA,down95_noSand=NA,up50_noSand=NA,down50_noSand=NA)

# subset predictor matrix to moku
predictor_full_sub <- predictor_full_wMMA_noNA %>% filter(Full_NoTak=='Yes')
pred_post_full <- matrix(NA,15000,nrow(predictor_full_sub))

# calculate the posterior
for(i in 1:nrow(predictor_full_sub)){
  pred_post_full[,i] <- b0moku_post[,predictor_full_sub$moku_ind[i]] + beta_matrix %*% t(as.matrix(predictor_full_sub[i,1:28]))
}

# summarise posterior
mlcd_allMHI$mean <- exp(mean(pred_post_full))
mlcd_allMHI$down80 <- exp(quantile(pred_post_full,0.10))
mlcd_allMHI$up80 <- exp(quantile(pred_post_full,0.90))
mlcd_allMHI$down95 <- exp(quantile(pred_post_full,0.025))
mlcd_allMHI$up95 <- exp(quantile(pred_post_full,0.975))
mlcd_allMHI$down50 <- exp(quantile(pred_post_full,0.25))
mlcd_allMHI$up50 <- exp(quantile(pred_post_full,0.75))

pred_post_full_noSand <- pred_post_full
pred_post_full_noSand[,which(predictor_full_sub$HIMARC_hard_soft==1)] <- NA

mlcd_allMHI$mean_noSand <- exp(mean(pred_post_full_noSand,na.rm=T))
mlcd_allMHI$down80_noSand <- exp(quantile(pred_post_full_noSand,0.10,na.rm=T))
mlcd_allMHI$up80_noSand <- exp(quantile(pred_post_full_noSand,0.90,na.rm=T))
mlcd_allMHI$down95_noSand <- exp(quantile(pred_post_full_noSand,0.025,na.rm=T))
mlcd_allMHI$up95_noSand <- exp(quantile(pred_post_full_noSand,0.975,na.rm=T))
mlcd_allMHI$down50_noSand <- exp(quantile(pred_post_full_noSand,0.25,na.rm=T))
mlcd_allMHI$up50_noSand <- exp(quantile(pred_post_full_noSand,0.75,na.rm=T))


write.csv(mlcd_allMHI, 'outputs/posterior_summary_mlcdMHIwide_TotalHerbivores.csv',row.names=F)

# barplot -----------------------------------------------------------------
moku_label <- read.csv('data/moku_label-4.csv', encoding = "UTF-8")

moku_order_use <- read.csv('model_out/moku_order_TotalHerbivore.csv')
colnames(moku_order_use)[2] <- 'moku_ind'
moku_order_use <- rbind(moku_order_use,data.frame(moku=c("KAUPO","KAHIKINUI","MANA"),moku_ind=c(41,42,43))) # add rows fo 3 moku with no data

mlcd <- mlcd_allMHI

moku_sum <- left_join(moku_summary, moku_order_use, by='moku_ind')
moku_sum <- left_join(moku_sum, moku_label, by='moku')
moku_sum <- moku_sum %>% arrange(desc(mean_noSand))

moku_sum$colgo <- rgb(227, 221, 209,max=255)#rgb(213, 229, 232,max=255)#'#DAA14C'
moku_sum$colgo[moku_sum$mean_noSand >= quantile(moku_sum$mean_noSand,0.75)] <- rgb(92, 148, 160, max=255)#'#5C94A0'
moku_sum$colgo[moku_sum$mean_noSand <= quantile(moku_sum$mean_noSand,0.25)] <- '#964336'

indlab <- expression("Herbivore Biomass"~~bgroup("(",'g '*m^{-2},")"))

png(file=paste0('outputs/Figure3.png'),height=2300,width=1700,res=200)
par(mfrow=c(1,1),mar=c(5,9,1,1),oma=c(0,0,0,0))
bp <- barplot((rev(moku_sum$mean_noSand)),border=F,names.arg=rep("",nrow(moku_sum))
              ,col='white',space=.3,xlim=c(0,max(moku_sum$up50_noSand)),
              ylab="",main="",cex.axis=1.2, horiz=T,xaxt='n')
rect(xleft=mlcd$down50_noSand,ybottom=bp[43],xright=mlcd$up50_noSand,ytop=bp[1],col='grey90',border=F)
barplot((rev(moku_sum$mean_noSand)),border=T,names.arg=rep("",nrow(moku_sum))
        ,col=(rev(moku_sum$colgo)),space=.3,xlim=c(0,max(moku_sum$up80_noSand)),
        ylab="",main="",cex.axis=1.2, horiz=T,xaxt='n',add=T)
axis(1,las=1,cex.axis=1.2,mgp=c(1,0.7,0))
plotrix::plotCI(rev(moku_sum$mean_noSand),bp, ui=rev(moku_sum$up50_noSand), li=rev(moku_sum$down50_noSand),add=T,err="x",pch=NA,sfrac=0.003)
text(par("usr")[3] +1.5,bp, labels = rev(moku_sum$Moku_olelo), pos = 2, xpd = TRUE)
mtext(indlab,side=1,outer=F,cex=1.3,line=3)
xb <- ((mlcd$up50_noSand - mlcd$down50_noSand)/2) + mlcd$down50_noSand
text(xb,bp[2]+0.5,'Predicted biomass in MLCDs',col='grey10',pos=1,font=3,cex=0.9)
dev.off()


# predicted raster --------------------------------------------------------

# moku mean 
grepgo <- grep('b0', colnames(mod_out[[1]]))
b0_summary <- summary(mcmc.list(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo]))
b0_moku_out <- data.frame(b0_summary$quantiles[grep('b0.m',rownames(b0_summary$quantiles)),])
moku_mean <- b0_moku_out$X50.

mas1 = raster("Data/moku_rasters/Moku_1.tif")
mas2 = raster("Data/moku_rasters/Moku_2.tif")
mas3 = raster("Data/moku_rasters/Moku_3.tif")
mas4 = raster("Data/moku_rasters/Moku_4.tif")
mas5 = raster("Data/moku_rasters/Moku_5.tif")
mas6 = raster("Data/moku_rasters/Moku_6.tif")
mas7 = raster("Data/moku_rasters/Moku_7.tif")
mas8 = raster("Data/moku_rasters/Moku_8.tif")
mas9 = raster("Data/moku_rasters/Moku_9.tif")
mas10 = raster("Data/moku_rasters/Moku_10.tif")
mas11 = raster("Data/moku_rasters/Moku_11.tif")
mas12 = raster("Data/moku_rasters/Moku_12.tif")
mas13 = raster("Data/moku_rasters/Moku_13.tif")
mas14 = raster("Data/moku_rasters/Moku_14.tif")
mas15 = raster("Data/moku_rasters/Moku_15.tif")
mas16 = raster("Data/moku_rasters/Moku_16.tif")
mas17 = raster("Data/moku_rasters/Moku_17.tif")
mas18 = raster("Data/moku_rasters/Moku_18.tif")
mas19 = raster("Data/moku_rasters/Moku_19.tif")
mas20 = raster("Data/moku_rasters/Moku_20.tif")
mas21 = raster("Data/moku_rasters/Moku_21.tif")
mas22 = raster("Data/moku_rasters/Moku_22.tif")
mas23 = raster("Data/moku_rasters/Moku_23.tif")
mas24 = raster("Data/moku_rasters/Moku_24.tif")
mas25 = raster("Data/moku_rasters/Moku_25.tif")
mas26 = raster("Data/moku_rasters/Moku_26.tif")
mas27 = raster("Data/moku_rasters/Moku_27.tif")
mas28 = raster("Data/moku_rasters/Moku_28.tif")
mas29 = raster("Data/moku_rasters/Moku_29.tif")
mas30 = raster("Data/moku_rasters/Moku_30.tif")
mas31 = raster("Data/moku_rasters/Moku_31.tif")
mas32 = raster("Data/moku_rasters/Moku_32.tif")
mas33 = raster("Data/moku_rasters/Moku_33.tif")
mas34 = raster("Data/moku_rasters/Moku_34.tif")
mas35 = raster("Data/moku_rasters/Moku_35.tif")
mas36 = raster("Data/moku_rasters/Moku_36.tif")
mas37 = raster("Data/moku_rasters/Moku_37.tif")
mas38 = raster("Data/moku_rasters/Moku_38.tif")
mas39 = raster("Data/moku_rasters/Moku_39.tif")
mas40 = raster("Data/moku_rasters/Moku_40.tif")

outras0 = ((mas1 / 1) * moku_mean[1]) + ((mas2 / 2) * moku_mean[2]) + ((mas3 / 3) * moku_mean[3]) + ((mas4 / 4) * moku_mean[4]) + ((mas5 / 5) * moku_mean[5]) + ((mas6 / 6) * moku_mean[6]) + ((mas7 / 7) * moku_mean[7]) + ((mas8 / 8) * moku_mean[8]) + ((mas9 / 9) * moku_mean[9]) + ((mas10 / 10) * moku_mean[10])

outras1 = outras0 + ((mas11 / 11) * moku_mean[11]) + ((mas12 / 12) * moku_mean[12]) + ((mas13 / 13) * moku_mean[13]) + ((mas14 / 14) * moku_mean[14]) + ((mas15 / 15) * moku_mean[15]) + ((mas16 / 16) * moku_mean[16]) + ((mas17 / 17) * moku_mean[17]) + ((mas18 / 18) * moku_mean[18]) + ((mas19 / 19) * moku_mean[19]) + ((mas20 / 20) * moku_mean[20])

outras2 = outras1 + ((mas21 / 21) * moku_mean[21]) + ((mas22 / 22) * moku_mean[22]) + ((mas23 / 23) * moku_mean[23]) + ((mas24 / 24) * moku_mean[24]) + ((mas25 / 25) * moku_mean[25]) + ((mas26 / 26) * moku_mean[26]) + ((mas27 / 27) * moku_mean[27]) + ((mas28 / 28) * moku_mean[28]) + ((mas29 / 29) * moku_mean[29]) + ((mas30 / 30) * moku_mean[30])

outras3 = outras2 + ((mas31 / 31) * moku_mean[31]) + ((mas32 / 32) * moku_mean[32]) + ((mas33 / 33) * moku_mean[33]) + ((mas34 / 34) * moku_mean[34]) + ((mas35 / 35) * moku_mean[35]) + ((mas36 / 36) * moku_mean[36]) + ((mas37 / 37) * moku_mean[37]) + ((mas38 / 38) * moku_mean[38]) + ((mas39 / 39) * moku_mean[39]) + ((mas40 / 40) * moku_mean[40])

writeRaster(outras3, paste0('outputs/TotalHerbivore_moku_mean_v2.tif'),
            format="GTiff",
            overwrite=TRUE,
            NAflag=-9999)

moku_mean <- outras3


# read in predictor rasters 

ras1 = raster("Data/all_predictors_CommonExtent_100m/01_Analysis_Mask_100m.tif")
ras2 = raster("Data/all_predictors_CommonExtent_100m/02_hab_100m_boulder_ready.tif")
ras3 = raster("Data/all_predictors_CommonExtent_100m/03_hab_100m_coral_ready.tif")
ras4 = raster("Data/all_predictors_CommonExtent_100m/04_hab_100m_pavement_ready.tif")
ras5 = raster("Data/all_predictors_CommonExtent_100m/05_CHL_ANOM_F_ready.tif")
ras6 = raster("Data/all_predictors_CommonExtent_100m/06_CHL_LTM_ready.tif")
ras7 = raster("Data/all_predictors_CommonExtent_100m/07_CHL_ANOM_M_ready.tif")
ras8 = raster("Data/all_predictors_CommonExtent_100m/08_WAV_ANOM_F_ready.tif")
ras9 = raster("Data/all_predictors_CommonExtent_100m/09_WAV_ANOM_M_ready.tif")
ras10 = raster("Data/all_predictors_CommonExtent_100m/10_PAR_LTM_ready.tif")
ras11 = raster("Data/all_predictors_CommonExtent_100m/11_SST_LTM_ready.tif")
ras12 = raster("Data/all_predictors_CommonExtent_100m/12_SST_STD_ready.tif")
ras13 = raster("Data/all_predictors_CommonExtent_100m/13_OTP_MHI_Fishing_Commercial_Line_ready.tif")
ras14 = raster("Data/all_predictors_CommonExtent_100m/14_OTP_MHI_Fishing_Commercial_Spear_ready.tif")
ras15 = raster("Data/all_predictors_CommonExtent_100m/OTP_MHI_Fishing_combined_net.tif")
ras16 = raster("Data/all_predictors_CommonExtent_100m/15_OTP_MHI_Fishing_NonCommercial_BoatBased_Spear_ready.tif")
ras17 = raster("Data/all_predictors_CommonExtent_100m/17_OTP_MHI_Fishing_NonCommercial_ShoreBased_Line_ready.tif")
ras18 = raster("Data/all_predictors_CommonExtent_100m/18_OTP_MHI_Fishing_NonCommercial_ShoreBased_Spear_ready.tif")
ras19 = raster("Data/all_predictors_CommonExtent_100m/OTP_MHI_Fishing_NonCommercial_ShoreBased_Net_ready.tif")
ras20 = raster("Data/all_predictors_CommonExtent_100m/19_Fishing_Aquarium_AnnAvg_NumHa_ready.tif")
ras21 = raster("Data/all_predictors_CommonExtent_100m/20_OTP_MHI_HabitatModification_ready.tif")
ras22 = raster("Data/all_predictors_CommonExtent_100m/21_OTP_MHI_OSDS_TotalEffluent_ready.tif")
ras23 = raster("Data/all_predictors_CommonExtent_100m/22_OTP_MHI_Sedimentation_ready.tif")
ras24 = raster("Data/all_predictors_CommonExtent_100m/23_Ag_Mean_Decay_01_ready.tif")
ras25 = raster("Data/all_predictors_CommonExtent_100m/24_Golf_Mean_Decay_01_ready.tif")
ras26 = raster("Data/all_predictors_CommonExtent_100m/25_LBSP_Urban_runoff_01_ready.tif")
ras27 = raster("Data/all_predictors_CommonExtent_100m/26_Depth100m_FocMEAN60m_Mosaic_rank_ready.tif") * -1
ras28 = raster("Data/all_predictors_CommonExtent_100m/27_Slope2_100m_FocMEAN60m_combined_transform_ready.tif")

# center and multiply by beta 
scale_ab <- read.csv('model_out/predictor_scale_TotalHerbivore.csv') # for transforming to predictor scale used in model

grepgo <- grep('beta',colnames(mod_out[[1]]))
beta_out <- summary(mcmc.list(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo]))
betas <- beta_out$quantiles[,'50%']

ras1_p = ras1 * betas[1]
ras2_p = ras2 * betas[2]
ras3_p = ras3 * betas[3]
ras4_p = ras4 * betas[4]
ras5_p = (scale_ab[5,2] + scale_ab[5,3]*ras5) * betas[5]
ras6_p = (scale_ab[6,2] + scale_ab[6,3]*ras6) * betas[6]
ras7_p = (scale_ab[7,2] + scale_ab[7,3]*ras7) * betas[7]
ras8_p = (scale_ab[8,2] + scale_ab[8,3]*ras8) * betas[8]
ras9_p = (scale_ab[9,2] + scale_ab[9,3]*ras9) * betas[9]
ras10_p = (scale_ab[10,2] + scale_ab[10,3]*ras10) * betas[10]
ras11_p = (scale_ab[11,2] + scale_ab[11,3]*ras11) * betas[11]
ras12_p = (scale_ab[12,2] + scale_ab[12,3]*ras12) * betas[12]
ras13_p = (scale_ab[13,2] + scale_ab[13,3]*ras13) * betas[13]
ras14_p = (scale_ab[14,2] + scale_ab[14,3]*ras14) * betas[14]
ras15_p = (scale_ab[15,2] + scale_ab[15,3]*ras15) * betas[15]
ras16_p = (scale_ab[16,2] + scale_ab[16,3]*ras16) * betas[16]
ras17_p = (scale_ab[17,2] + scale_ab[17,3]*ras17) * betas[17]
ras18_p = (scale_ab[18,2] + scale_ab[18,3]*ras18) * betas[18]
ras19_p = (scale_ab[19,2] + scale_ab[19,3]*ras19) * betas[19]
ras20_p = (scale_ab[20,2] + scale_ab[20,3]*ras20) * betas[20]
ras21_p = (scale_ab[21,2] + scale_ab[21,3]*ras21) * betas[21]
ras22_p = (scale_ab[22,2] + scale_ab[22,3]*ras22) * betas[22]
ras23_p = (scale_ab[23,2] + scale_ab[23,3]*ras23) * betas[23]
ras24_p = (scale_ab[24,2] + scale_ab[24,3]*ras24) * betas[24]
ras25_p = (scale_ab[25,2] + scale_ab[25,3]*ras25) * betas[25]
ras26_p = (scale_ab[26,2] + scale_ab[26,3]*ras26) * betas[26]
ras27_p = (scale_ab[27,2] + scale_ab[27,3]*ras27) * betas[27]
ras28_p = (scale_ab[28,2] + scale_ab[28,3]*ras28) * betas[28]


# calc prediction 

out_ras = raster::overlay(moku_mean,ras1_p,ras2_p,ras3_p,ras4_p,ras5_p,ras6_p,ras7_p,ras8_p,ras9_p,ras10_p,ras11_p,ras12_p,ras13_p,ras14_p,ras15_p,ras16_p,ras17_p,ras18_p,ras19_p,ras20_p,ras21_p,ras22_p,ras23_p,ras24_p,ras25_p,ras26_p,ras27_p,ras28_p,fun=sum)

out_ras_exp <- exp(out_ras)

writeRaster(out_ras_exp, paste0('outputs/TotalHerbivore_pred100m.tif'),
            format="GTiff",
            overwrite=TRUE,
            NAflag=-9999)


