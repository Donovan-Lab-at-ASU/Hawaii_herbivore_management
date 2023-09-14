# -----------------------------------------------------------------------
# Estimating fishing effects on herbivores and translation to benthos
# -----------------------------------------------------------------------
# Model and code by Mary Donovan - marydonovan@asu.edu
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

library(dplyr) # v1.0.7
library(rjags) # v4.10 linked to JAGS 4.3.0
library(ggplot2) # v3.4.0
library(vioplot) # v0.4.0


# inputs ------------------------------------------------------------------

# fish models
hgd_mod <- readRDS("model_out/Hgd_gamma.Rdata")
hbrow_mod <- readRDS("model_out/Hbrow_gamma.Rdata")
hscex_mod <- readRDS("model_out/Hscex_gamma.Rdata")

# predictor matrix - rows corresponds to the 100x100m pixels within the study domain
predictor_full <- readRDS(file="Data/full_predictor_matrix.Rdata")

# 3 moku with no data, but we still want to calculate a posterior, so need to update the index and add columns of zeros to the moku intercept matrix (and thus giving those intercepts a value of zero)
predictor_full$moku_ind <- ifelse(predictor_full$moku_ind == -36, 41, predictor_full$moku_ind)
predictor_full$moku_ind <- ifelse(predictor_full$moku_ind == -35, 42, predictor_full$moku_ind)
predictor_full$moku_ind <- ifelse(predictor_full$moku_ind == -13, 43, predictor_full$moku_ind)

predictor_full$moku <- ifelse(predictor_full$moku_ind == 41,"KAUPO", predictor_full$moku)
predictor_full$moku <- ifelse(predictor_full$moku_ind == 42,"KAHIKINUI", predictor_full$moku)
predictor_full$moku <- ifelse(predictor_full$moku_ind == 43,"MANA", predictor_full$moku)

# Remove a few rows for pixels determined to be not fit in our domain (estuaries, on land, etc) 
predictor_full <- predictor_full %>% filter(is.na(bad_pt))

pred_full <- predictor_full

# predict fish data -------------------------------------------------------

############## HSCEX
mod_out <- hscex_mod

# betas
grepgo <- grep('beta',colnames(mod_out[[1]]))
beta_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
beta_matrix <- as.matrix(beta_post)

# moku means
grepgo <- grep('b0.moku',colnames(mod_out[[1]]))
b0moku_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
b0moku_post <- cbind(b0moku_post,rep(0,15000),rep(0,15000),rep(0,15000)) # for 3 missing moku

# create new levels of human impacts
predictor_full_minF <- pred_full
for(k in 13:20){predictor_full_minF[,k] <- min(predictor_full_minF[,k])}

pixel_out <- data.frame(rowid=seq(1:nrow(pred_full)),
                        pred=NA,pred_up=NA,pred_down=NA,
                        minHF=NA,minF_up=NA,minF_down=NA
)

pred_post_full <- matrix(NA,15000,nrow(pred_full))
for(i in 1:nrow(pred_full)){
  pred_post_full[,i] <- b0moku_post[,pred_full$moku_ind[i]] + beta_matrix %*% t(as.matrix(pred_full[i,1:28]))
}

for(i in 1:ncol(pred_post_full)){
  pixel_out$pred[i] <- mean(pred_post_full[,i])
  pixel_out$pred_up[i] <- quantile(pred_post_full[,i],0.75)
  pixel_out$pred_down[i] <- quantile(pred_post_full[,i],0.25)
}

#min
pred_post_full <- matrix(NA,15000,nrow(predictor_full_minF))
for(i in 1:nrow(predictor_full_minF)){
  pred_post_full[,i] <- b0moku_post[,predictor_full_minF$moku_ind[i]] + beta_matrix %*% t(as.matrix(predictor_full_minF[i,1:28]))
}

for(i in 1:ncol(pred_post_full)){
  pixel_out$minHF[i] <- mean(pred_post_full[,i])
  pixel_out$minF_up[i] <- quantile(pred_post_full[,i],0.75)
  pixel_out$minF_down[i] <- quantile(pred_post_full[,i],0.25)
}


pixel_out_hscex <- pixel_out %>% rename_at(vars(colnames(pixel_out)), ~paste0("hscex_",colnames(pixel_out)))
saveRDS(pixel_out_hscex, file="outputs/pixel_out_hscex.Rdata")
pixel_out_hscex <- readRDS("outputs/pixel_out_hscex.Rdata")

############## HGD
mod_out <- hgd_mod

# betas
grepgo <- grep('beta',colnames(mod_out[[1]]))
beta_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
beta_matrix <- as.matrix(beta_post)

# moku means
grepgo <- grep('b0.moku',colnames(mod_out[[1]]))
b0moku_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
b0moku_post <- cbind(b0moku_post,rep(0,15000),rep(0,15000),rep(0,15000)) # for 3 missing moku

# create new levels of human impacts
predictor_full_minF <- pred_full
for(k in 13:20){predictor_full_minF[,k] <- min(predictor_full_minF[,k])}

pixel_out <- data.frame(rowid=seq(1:nrow(pred_full)),
                        pred=NA,pred_up=NA,pred_down=NA,
                        minHF=NA,minF_up=NA,minF_down=NA
)

pred_post_full <- matrix(NA,15000,nrow(pred_full))
for(i in 1:nrow(pred_full)){
  pred_post_full[,i] <- b0moku_post[,pred_full$moku_ind[i]] + beta_matrix %*% t(as.matrix(pred_full[i,1:28]))
}

for(i in 1:ncol(pred_post_full)){
  pixel_out$pred[i] <- mean(pred_post_full[,i])
  pixel_out$pred_up[i] <- quantile(pred_post_full[,i],0.75)
  pixel_out$pred_down[i] <- quantile(pred_post_full[,i],0.25)
}

#min
pred_post_full <- matrix(NA,15000,nrow(predictor_full_minF))
for(i in 1:nrow(predictor_full_minF)){
  pred_post_full[,i] <- b0moku_post[,predictor_full_minF$moku_ind[i]] + beta_matrix %*% t(as.matrix(predictor_full_minF[i,1:28]))
}

for(i in 1:ncol(pred_post_full)){
  pixel_out$minHF[i] <- mean(pred_post_full[,i])
  pixel_out$minF_up[i] <- quantile(pred_post_full[,i],0.75)
  pixel_out$minF_down[i] <- quantile(pred_post_full[,i],0.25)
}

pixel_out_hgd <- pixel_out %>% rename_at(vars(colnames(pixel_out)), ~paste0("hgd_",colnames(pixel_out)))
saveRDS(pixel_out_hgd, file="outputs/pixel_out_hgd.Rdata")
pixel_out_hgd <- readRDS("outputs/pixel_out_hgd.Rdata")



############## HBROW
mod_out <- hbrow_mod

# betas
grepgo <- grep('beta',colnames(mod_out[[1]]))
beta_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
beta_matrix <- as.matrix(beta_post)

# moku means
grepgo <- grep('b0.moku',colnames(mod_out[[1]]))
b0moku_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
b0moku_post <- cbind(b0moku_post,rep(0,15000),rep(0,15000),rep(0,15000)) # for 3 missing moku

# create new levels of human impacts
predictor_full_minF <- pred_full
for(k in 13:20){predictor_full_minF[,k] <- min(predictor_full_minF[,k])}

pixel_out <- data.frame(rowid=seq(1:nrow(pred_full)),
                        pred=NA,pred_up=NA,pred_down=NA,
                        minHF=NA,minF_up=NA,minF_down=NA
                        # ,halfF=NA,halfF_up=NA,halfF_down=NA
)

pred_post_full <- matrix(NA,15000,nrow(pred_full))
for(i in 1:nrow(pred_full)){
  pred_post_full[,i] <- b0moku_post[,pred_full$moku_ind[i]] + beta_matrix %*% t(as.matrix(pred_full[i,1:28]))
}

for(i in 1:ncol(pred_post_full)){
  pixel_out$pred[i] <- mean(pred_post_full[,i])
  pixel_out$pred_up[i] <- quantile(pred_post_full[,i],0.75)
  pixel_out$pred_down[i] <- quantile(pred_post_full[,i],0.25)
}

#min
pred_post_full <- matrix(NA,15000,nrow(predictor_full_minF))
for(i in 1:nrow(predictor_full_minF)){
  pred_post_full[,i] <- b0moku_post[,predictor_full_minF$moku_ind[i]] + beta_matrix %*% t(as.matrix(predictor_full_minF[i,1:28]))
}

for(i in 1:ncol(pred_post_full)){
  pixel_out$minHF[i] <- mean(pred_post_full[,i])
  pixel_out$minF_up[i] <- quantile(pred_post_full[,i],0.75)
  pixel_out$minF_down[i] <- quantile(pred_post_full[,i],0.25)
}


pixel_out_hbrow <- pixel_out %>% rename_at(vars(colnames(pixel_out)), ~paste0("hbrow_",colnames(pixel_out)))
saveRDS(pixel_out_hbrow, file="outputs/pixel_out_hbrow.Rdata")
pixel_out_hbrow <- readRDS("outputs/pixel_out_hbrow.Rdata")



# predict benthic response ------------------------------------------------

# set up driver data 
beta_names <- read.csv('model_out/beta_names_TotalHerbivore.csv')
colnames(pred_full)[1:26] <- beta_names$x

pred_full.go <- pred_full %>% select(OTP_CHL_ANOM_F:OTP_SST_STD
                                     ,OTP_HabitatModification:LBSP2_Urban_runoff,
                                     beta_depth,beta_rugosity
) 

herb_full <- data.frame(Scraper_log=pixel_out_hscex$hscex_pred,Grazer_log=pixel_out_hgd$hgd_pred,Browser_log=pixel_out_hbrow$hbrow_pred)
herb_full$Scraper_x_Grazer <- herb_full$Scraper_log*herb_full$Grazer_log
herb_full$Scraper_x_Browser <- herb_full$Scraper_log*herb_full$Browser_log
herb_full$Grazer_x_Browser <- herb_full$Grazer_log*herb_full$Browser_log
herb_full$Scaper_x_Grazer_x_Browser <- herb_full$Scraper_log*herb_full$Grazer_log*herb_full$Browser_log

pred_full.go <- cbind(herb_full,pred_full.go)

pred_ab <- read.csv('benthic_all_drivers_out/model_out/predictor_scale_ratio.CM.csv')

pred_full_scaled <- pred_full.go
for(i in 1:7){pred_full_scaled[,i] <- pred_ab$a[i] + pred_ab$b[i]*pred_full.go[,i]}


herbs_all_min <- data.frame(Scraper_log=pixel_out_hscex$hscex_minHF,Grazer_log=pixel_out_hgd$hgd_minHF,Browser_log=pixel_out_hbrow$hbrow_minHF)
herbs_all_min$Scraper_x_Grazer <- herbs_all_min$Scraper_log*herbs_all_min$Grazer_log
herbs_all_min$Scraper_x_Browser <- herbs_all_min$Scraper_log*herbs_all_min$Browser_log
herbs_all_min$Grazer_x_Browser <- herbs_all_min$Grazer_log*herbs_all_min$Browser_log
herbs_all_min$Scaper_x_Grazer_x_Browser <- herbs_all_min$Scraper_log*herbs_all_min$Grazer_log*herbs_all_min$Browser_log
pred_full.go <- pred_full %>% select(OTP_CHL_ANOM_F:OTP_SST_STD
                                     ,OTP_HabitatModification:LBSP2_Urban_runoff,
                                     beta_depth,beta_rugosity) 
pred_full.go <- cbind(herbs_all_min,pred_full.go)
pred_full_scaled_min <- pred_full.go
for(i in 1:7){pred_full_scaled_min[,i] <- pred_ab$a[i] + pred_ab$b[i]*pred_full.go[,i]}


# predict benthic

pred <- data.frame(pred_all=rep(NA,nrow(herbs_pred)),pred_all_down=rep(NA,nrow(herbs_pred)),
                   pred_all_up=rep(NA,nrow(herbs_pred)),pred_min_all=rep(NA,nrow(herbs_pred)),
                   pred_min_all_down=rep(NA,nrow(herbs_pred)),pred_min_all_up=rep(NA,nrow(herbs_pred))
                  )

# pred_all
X <- model.matrix(~1 + as.matrix(pred_full_scaled))
X <- X[,2:ncol(X)] # no need for intercept

zmPrp <- readRDS(file='benthic_all_drivers_out/model_out/benthic_model_cm_alldrivers_out.Rdata')
grepgo <- c(grep('beta',colnames(zmPrp[[1]]))) 
beta_d <- rbind(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo])
beta_matrix <- as.matrix(beta_d)

# predict from model output ignoring random intercepts
for(i in 1:nrow(pred_full_scaled)){
  mcmc_go <- beta_matrix %*% X[i,]
  pred$pred_all[i] <- median(mcmc_go)
  pred$pred_all_down[i] <- quantile(mcmc_go,0.25)
  pred$pred_all_up[i] <- quantile(mcmc_go,0.75)
  
}

# pred_min_all
X <- model.matrix(~1 + as.matrix(pred_full_scaled_min))
X <- X[,2:ncol(X)] # no need for intercept

# predict from model output ignoring random intercepts
for(i in 1:nrow(pred_full_scaled_min)){
  mcmc_go <- beta_matrix %*% X[i,]
  pred$pred_min_all[i] <- median(mcmc_go)
  pred$pred_min_all_down[i] <- quantile(mcmc_go,0.25)
  pred$pred_min_alll_up[i] <- quantile(mcmc_go,0.75)
  
}

pred_copy <- pred

pred$hard_soft <- pred_full$HIMARC_hard_soft
pred$x <- pred_full$x; pred$y <- pred_full$y
pred$moku <- pred_full$moku
pred$moku[pred_full$moku_ind==41] <- "KAUPO"
pred$moku[pred_full$moku_ind==42] <- "KAHIKINUI"
pred$moku[pred_full$moku_ind==43] <- "MANA"

saveRDS(pred, "benthic_all_drivers_out/predbenthic_bigmatrix.Rdata")


# total herbivores % potential biomass ------------------------------------

h_mod <- readRDS("model_out/H_gamma.Rdata") # total herbivores

mod_out <- h_mod

# betas
grepgo <- grep('beta',colnames(mod_out[[1]]))
beta_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
beta_matrix <- as.matrix(beta_post)

# moku means
grepgo <- grep('b0.moku',colnames(mod_out[[1]]))
b0moku_post <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])
b0moku_post <- cbind(b0moku_post,rep(0,15000),rep(0,15000),rep(0,15000)) # for 3 missing moku

# create new levels of human impacts
predictor_full_minF <- pred_full
for(k in 13:20){predictor_full_minF[,k] <- min(predictor_full_minF[,k])}

pixel_out <- data.frame(rowid=seq(1:nrow(pred_full)),
                        pred=NA,pred_up=NA,pred_down=NA,
                        minHF=NA,minF_up=NA,minF_down=NA
)

pred_post_full <- matrix(NA,15000,nrow(pred_full))
pred_post_full_min <- matrix(NA,15000,nrow(pred_full))
pred_post_full_ratio <- matrix(NA,15000,nrow(pred_full))
for(i in 1:nrow(pred_full)){
  pred_post_full[,i] <- b0moku_post[,pred_full$moku_ind[i]] + beta_matrix %*% t(as.matrix(pred_full[i,1:28]))
  pred_post_full_min[,i] <- b0moku_post[,pred_full$moku_ind[i]] + beta_matrix %*% t(as.matrix(predictor_full_minF[i,1:28]))
  pred_post_full_ratio[,i] <- pred_post_full[,i]/pred_post_full_min[,i]
}

ratio_post <- pred_post_full_ratio
# saveRDS(ratio_post, file="outputs/ratio_post.Rdata")

moku_ratio <- data.frame(moku=unique(pred_full$moku),mean=NA,up=NA,down=NA,median=NA)
for(k in 1:nrow(moku_ratio)){
  moku_go <- which(pred_full$moku==moku_ratio$moku[k])
  ratio_post_sel <- ratio_post[,moku_go]
  moku_ratio$mean[k] <- mean(ratio_post_sel)
  moku_ratio$up[k] <- quantile(ratio_post_sel,0.75)
  moku_ratio$down[k] <- quantile(ratio_post_sel,0.25)
  moku_ratio$median[k] <- median(ratio_post_sel)
}

ratio.mean <- rep(NA,nrow(pred_full))
for(i in 1:ncol(pred_post_full_ratio)){
  ratio.mean[i] <- mean(pred_post_full_ratio[,i])
}
saveRDS(ratio.mean, file="outputs/ratio.mean.Rdata")
ratio.mean <- readRDS("outputs/ratio.mean.Rdata")

# benthic logistic regression ---------------------------------------------

temp <- data.frame(ratio.mean,ben_sig=ifelse(pred$ben_rp_sig=='black',1,0),hard_soft=pred$hard_soft,moku=pred_full$moku)

data_go <- temp %>% filter(ratio.mean < 5) %>% mutate(ratio.mean=ifelse(ratio.mean < 0,0,ratio.mean)) %>% filter(hard_soft==0) # removes one outlier pixel and soft bottom pixels
mod <- glm(ben_sig ~ ratio.mean, data=data_go,family='binomial')
p <- 0.99
(log(p/(1-p)) - coef(mod)[1])/coef(mod)[2]

# set up plot
data_go$y1 <- data_go$ben_sig
data_go$y2 <- ifelse(data_go$y1==0,-0.01,1.01)
data_go$y <- as.logical(data_go$y1)
data_go$x1 <- data_go$ratio.mean
data_go$x1 <- ifelse(data_go$x1 > 1, 1, data_go$x1)

a <- density(data_go$ratio.mean[data_go$ben_sig==1])
b <- data.frame(a$x, a$y)

a2 <- density(data_go$ratio.mean[data_go$ben_sig==0])
b2 <- data.frame(a2$x, a2$y)

##### summarise across moku
moku_under <- data_go %>% mutate(ratio_over = ifelse(ratio.mean > 0.80, 'over','under')) %>% group_by(moku,ratio_over) %>% tally() %>% tidyr::pivot_wider(names_from=ratio_over,values_from=n,values_fill=0) %>% mutate(sum = rowSums(across(where(is.numeric)))) %>% mutate(prop=under/sum) %>% arrange(prop)
moku_under %>% print(n=43)
sum(moku_under$under)/sum(moku_under$sum)

# Figure 3 ----------------------------------------------------------------
moku_label <- read.csv('data/moku_label-4.csv', encoding = "UTF-8")

moku_order_use <- read.csv('model_out/moku_order_TotalHerbivore.csv')
colnames(moku_order_use)[2] <- 'moku_ind'
moku_order_use <- rbind(moku_order_use,data.frame(moku=c("KAUPO","KAHIKINUI","MANA"),moku_ind=c(41,42,43))) # add rows fo 3 moku with no data

mlcd <- read.csv('outputs/posterior_summary_mlcdMHIwide_TotalHerbivores.csv')
moku_summary <- read.csv('outputs/posterior_summary_moku_TotalHerbivores.csv')

moku_sum <- left_join(moku_summary, moku_order_use, by='moku_ind')
moku_sum <- left_join(moku_sum, moku_label, by='moku')
moku_sum <- moku_sum %>% arrange(desc(mean_noSand))

moku_sum$colgo <- rgb(227, 221, 209,max=255)#rgb(213, 229, 232,max=255)#'#DAA14C'
moku_sum$colgo[moku_sum$mean_noSand >= quantile(moku_sum$mean_noSand,0.75)] <- rgb(92, 148, 160, max=255)#'#5C94A0'
moku_sum$colgo[moku_sum$mean_noSand <= quantile(moku_sum$mean_noSand,0.25)] <- '#964336'

indlab <- expression("Herbivore Biomass"~~bgroup("(",'g '*m^{-2},")"))
moku_sum$moku_lab <- paste0(moku_sum$Moku_olelo," (",seq(1:43),")")


pdf(file='outputs/Figure3_new.pdf',height=10.2,width=16,encoding = "CP1257.enc")
# png(file=paste0('outputs/Figure3_new.png'),height=2100,width=3500,res=200)

# Create the layout
nf <- layout(matrix(c(1,1,2,2), nrow=2, byrow=TRUE), heights=c(2.25,1.75) )


par(mar=c(0,5,0.25,1),oma=c(0,0,0,0),bg=NA)
bp <- barplot(((moku_sum$mean_noSand)),border=F,names.arg=rep("",nrow(moku_sum))
              ,col='white',space=.3,ylab="",main="",cex.axis=1.3, horiz=F,xaxt='n',yaxt="n"
              ,ylim=c(0,max(moku_sum$up50_noSand))
              ,xlim=c(0,56.2))
# rect(xleft=bp[1],xright=bp[43],ybottom=mlcd$down50_noSand,ytop=mlcd$up50_noSand,col='grey90',border=F)
barplot(((moku_sum$mean_noSand)),border=T,names.arg=rep("",nrow(moku_sum))
        ,col=((moku_sum$colgo)),space=.3,ylim=c(0,max(moku_sum$up50_noSand)),
        ylab="",main="",cex.axis=1.3, horiz=F,xaxt='n',add=T)
plotrix::plotCI(bp,(moku_sum$mean_noSand), ui=(moku_sum$up50_noSand), li=(moku_sum$down50_noSand),add=T,err="y",pch=NA,sfrac=0.003)
# text(bp+0.25, par("usr")[3] - 4, labels = moku_sum$moku_lab, srt = 55, pos = 2, xpd = TRUE)
mtext(indlab,side=2,outer=F,cex=1.2,line=2.25)
xb <- ((mlcd$up50_noSand - mlcd$down50_noSand)/2) + mlcd$down50_noSand
moku_sum$num_col <- c(rep("white",11),rep("black",21),rep("white",11))
text(bp,3.5,seq(1:43),col=moku_sum$num_col,font=2)
# text(bp[2]+0.5,xb,'Predicted biomass in MLCDs',col='grey10',pos=1,font=3,cex=0.9)

# datplot <- datplot %>% left_join(data.frame(moku=moku_sum$moku,bp=bp[,1]),by='moku')
par(mar=c(10,5,0.25,1))
vioplot(datplot$ratio.mean~datplot$bp,pchMed=NA,drawRect=F,col='white',xaxt="n",frame.plot=F,
        ylab="",cex.axis=1.3,cex.lab=3,xlab="",lwd=1.2)
# abline(h=0.8,col="red")
segments(0,0.8,43.5,0.8,col='red',lwd=1.5)
mtext("% Potential Biomass",side=2,outer=F,cex=1.2,line=3,at=0.5)
# text(bp*1, par("usr")[3] - 0, labels = moku_sum$moku_lab, srt = 55, pos = 2, xpd = TRUE)
text(bp*0.78, par("usr")[3] - 0, labels = moku_sum$moku_lab, srt = 55, pos = 2, xpd = T,cex=1.3)
dev.off()


# Figure 4 ----------------------------------------------------------------

library(ggbreak)
# 
p <- ggplot(data_go, aes(x=x1,y=y1)) +
  # geom_point(aes(x=x1,y=y2),cex = 1,color="black",alpha=0.5,pch=1) +
  stat_smooth(method = "glm", se = T, method.args = list(family = binomial)) +
  # geom_density(
  #   data = data_go[!data_go$y, ], aes(x = x1, y = after_stat(- density - 2)),
  #   alpha = 0.2, fill = "white"
  # ) +
  # geom_polygon(
  #   data = data_go[data_go$y, ], aes(x = x1, y = after_stat(density)),
  #   alpha = 0.2, fill = "blue", stat = "density"
  # ) +
  # geom_path(
  #   data = data_go[data_go$y, ], aes(x = x1, y = after_stat(2 + density)),
  #   stat = "density", color = "black"
# ) +
# geom_hline(yintercept=2, colour="white", size=0.75) +
# geom_hline(yintercept=-2, colour="white", size=0.75) +
geom_line(data=b,aes(x=a.x,y=a.y + 2),color='black',lwd=1) +
  geom_line(data=b2,aes(x=a2.x,y=-1*a2.y - 2),color='black',lwd=1) +
  theme_bw() + ylim(-10,10) + xlim(0,1.02) +
  scale_y_break(c(-1.93,-0.02),scale=4,ticklabels=c(0,0.5,1),space=0,expand=F) + scale_y_break(c(1.02,1.93),scale=1,space=0,expand=F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Probability Benthos Unaffected") + xlab("% Potential Biomass") +
  # geom_vline(aes(xintercept=(log(p/(1-p)) - coef(mod)[1])/coef(mod)[2])) +
  # geom_segment(aes(x = 0.8, y = 0.99, xend = 0.8, yend = 0.0),arrow=arrow(length = unit(0.3,"cm")),color="red") +
  # geom_segment(aes(x = 0, y = 0.99, xend = 0.8, yend = 0.99),arrow=arrow(length = unit(0.3,"cm")),color="red",lineend = "round") +
  # geom_segment(aes(x = 0.8, y = -8, xend = 0.8, yend = -10),arrow=arrow(length = unit(0.3,"cm")),color="red") +
  theme(axis.text = element_text(size = 15,colour="black"),axis.title = element_text(size = 15)) +
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black")) +
  theme(axis.ticks.length=unit(.25, "cm"))
p

