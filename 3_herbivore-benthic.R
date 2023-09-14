# -----------------------------------------------------------------------
# Models of herbivore-benthic relationships in Hawai‘i
# -----------------------------------------------------------------------
# Model and code by Mary Donovan - marydonovan@asu.edu
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

library(plotrix) # v3.7.8
library(dplyr) # v1.0.7
library(rjags) # v4.10 linked to JAGS 4.3.0
library(parallel)  # v4.0.3
library(bayesplot) # v1.8.0

# input data --------------------------------------------------------------
data <- read.csv("data/donovan_etal_HIMARCdata_for_herbivore_publication.csv")
data <- data[complete.cases(data),] # subset for only those observations with both fish and benthic data
data$Scraper_log <- log(data$Scraper+1)
data$Grazer_log <- log(data$Grazer+1)
data$Browser_log <- log(data$Browser+1)

# herbivore-driver model --------------------------------------------------
data$ratio.CM <- log((data$Coral + data$Calg + 1)/(data$Macro + 1))

varuse <- 'ratio.CM' # response variable

cat("model{
    # Likelihood
    for(i in 1:N){
      y[i] ~ dnorm(mu[i],tau_y)
      mu[i] <- b0.y[year[i]] + b0.d[dataset[i]] + b0.moku[moku[i]] + inprod(beta,X[i,])
      y.new[i] ~ dnorm(mu[i],tau_y)
    }
    tau_y ~ dgamma(0.1,0.1)
    sigma_y <- 1/sqrt(tau_y)
    for(b in 1:B){ beta[b] ~ dnorm(0,0.001) }

    ## year
    for(z in 1:(Y-1)){
      b0.y[z] ~ dnorm(0,tau_year) # zero centered prior
    }
    tau_year ~ dgamma(0.1,0.1)
    sigma_year <- 1/sqrt(tau_year)
    b0.y[Y] <- -sum(b0.y[1:(Y-1)]) # sum to zero constraint
    
    ## dataset
    for(d in 1:(D-1)){
      b0.d[d] ~ dnorm(0,tau_dataset)
    }
    tau_dataset ~ dgamma(0.1,0.1)
    sigma_dataset <- 1/sqrt(tau_dataset)
    b0.d[D] <- -sum(b0.d[1:(D-1)])
    
    ## moku
    for(m in 1:(M-1)){
      b0.moku[m] ~ dnorm(0,tau_moku)
    }
    tau_moku ~ dgamma(0.1,0.1)
    sigma_moku <- 1/sqrt(tau_moku)
    b0.moku[M] <- -sum(b0.moku[1:(M-1)])
    
    # Posterior predictive checks
    y.mean <- mean(y)
    ynew.mean <- mean(y.new)
    pval.mean <- step(y.mean-ynew.mean)

    sd.y <- sd(y)
    sd.y.new <- sd(y.new)
    pval.sd <- step(sd.y.new-sd.y)

    # R-squared
    varF <- sd(y)^2
    varE <- sd(y - y.new)^2
    R2 <- varF/(varF+varE)
    
    }",file='benthic_model.jags')

# response variable
response <- data[,varuse]
# response <- log(response)

# predictors
predict.go <- data %>% select(OTP_CHL_ANOM_F:OTP_SST_STD
                              ,OTP_HabitatModification:LBSP2_Urban_runoff
)
predict.go.save <- predict.go
for(i in 1:ncol(predict.go)) predict.go[,i] <- scale(predict.go[,i])[,1]
predict.go$depth <- scale(data$depth_predicted)[,1]
predict.go.save$depth <- data$depth_predicted
predict.go$rugosity <- scale(data$rugosity_predicted)[,1]
predict.go.save$rugosity <- data$rugosity_predicted
X <- model.matrix(~1+as.matrix(predict.go))
X <- X[,2:ncol(X)] # no need for intercept

scaled_pred <- data[,c('Scraper_log','Grazer_log','Browser_log')]
scaled_pred$Scraper_x_Grazer <- scaled_pred$Scraper_log*scaled_pred$Grazer_log
scaled_pred$Scraper_x_Browser <- scaled_pred$Scraper_log*scaled_pred$Browser_log
scaled_pred$Grazer_x_Browser <- scaled_pred$Grazer_log*scaled_pred$Browser_log
scaled_pred$Scaper_x_Grazer_x_Browser <- scaled_pred$Scraper_log*scaled_pred$Grazer_log*scaled_pred$Browser_log
predict.go.save <- cbind(scaled_pred,predict.go.save)
for(i in 1:ncol(scaled_pred)) scaled_pred[,i] <- scale(scaled_pred[,i])[,1]

X <- cbind(model.matrix(~1 + scaled_pred$Scraper_log + scaled_pred$Grazer_log + scaled_pred$Browser_log
                        + scaled_pred$Scraper_x_Grazer 
                        + scaled_pred$Scraper_x_Browser
                        + scaled_pred$Grazer_x_Browser 
                        + scaled_pred$Scaper_x_Grazer_x_Browser
),X)
X <- X[,2:ncol(X)] # no need for intercept
write.csv(colnames(X),paste0('model_out/beta_names_benthic.csv'),row.names = F)

# export conversion from scaled predictors
predictor_scale_ab <- data.frame(pred_ord=seq(1:23),pred=colnames(predict.go.save),a=NA,b=NA)
for(i in 1:23){
  pred_lm <- lm(X[,i]~predict.go.save[,i])
  predictor_scale_ab$a[i] <- summary(pred_lm)$coefficients[1,1]
  predictor_scale_ab$b[i] <- summary(pred_lm)$coefficients[2,1]
}

write.csv(predictor_scale_ab, paste0('model_out/predictor_scale_',varuse,'.csv'),row.names=F)


# moku levels
data$moku <- as.factor(as.character(data$moku))
mok_levs <- levels(data$moku)
write.csv(mok_levs,paste0('model_out/mok_levs_benthic.csv'),row.names = F)
nmoku <- length(unique(data$moku))

# create index
data$year.f <- as.numeric(as.factor(as.character(data$Year)))
data$moku.f <- as.numeric(as.factor(as.character(data$moku)))
data$dataset.f <- as.numeric(as.factor(as.character(data$Dataset)))

# save for export
year <- data %>% distinct(year.f)
dataset <- data %>% distinct(dataset.f)
moku <- data %>% distinct(moku.f) 

# model inputs
nX <- ncol(X)
N=length(response)

initFunc <- function(){return(list(
  beta=rnorm(nX,0,1)
  ,tau_year=runif(1,0.05,0.5),
  tau_dataset=runif(1,0.05,0.5),
  tau_moku=runif(1,0.05,0.5)
))} # initial conditions

n.adapt <- 500; n.update <- 5000; n.iter <- 2000

mdat <- list(N=length(response),
             y=response,
             X=X,
             B=ncol(X),
             year=data$year.f,
             dataset=data$dataset.f,
             moku=data$moku.f,
             Y=length(unique(data$year.f)),
             D=length(unique(data$dataset.f)),
             M=length(unique(data$moku.f))
)

write.csv(data.frame(
  response=response,
  year=data$Year,
  year.f=data$year.f,
  dataset=data$Dataset,
  dataset.f=data$dataset.f,
  moku=data$moku,
  moku.f=data$moku.f,
  replicate_id=data$replicate_id,
  Scraper_log=data$Scraper_log,
  Grazer_log=data$Grazer_log,
  Browser_log=data$Browser_log,
  X
), paste0('model_out/data_in_benthic.csv'), row.names=F)

# run
mod <- jags.model("benthic_model.jags", data= mdat, n.chains=3, n.adapt=n.adapt,inits = initFunc())
update(mod, n.iter = n.update)
zmCore <- coda.samples(mod, c("pval.mean","pval.sd","R2","sigma_y",
                              "beta","b0.moku","b0.y","b0.d","y.new"
),n.iter=n.iter, n.thin=1)
zmPrp <- zmCore

# export 
saveRDS(zmPrp, file='model_out/benthic_model_cm_alldrivers_out.Rdata')

grepgo <- grep('beta',colnames(zmPrp[[1]]))
beta_out <- summary(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))
beta_out
write.csv(
  data.frame(beta_out$quantiles,c(colnames(X)))
  , 'model_out/benthic_betas.csv', row.names = T)

grepgo <- grep(c('R2'),colnames(zmPrp[[1]]))
R2_sum <- summary(mcmc.list(zmPrp[[1]][,grepgo],
                            zmPrp[[2]][,grepgo],
                            zmPrp[[3]][,grepgo]))
R2_sum

# posterior check
grepgo <- grep('y.new', colnames(zmPrp[[1]]))
yrep <- rbind(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo])
old <- bayesplot_theme_set(ggplot2::theme_minimal())

moku_label <- read.csv('data/moku_label-4.csv', encoding = "UTF-8")
data <- dplyr::left_join(data,moku_label,by='moku') # modify moku label

png(file="outputs/FigureS6_revised.png",height=2500,width=3500,res=300)
bayesplot::ppc_stat_grouped(response, yrep, group = data$Moku_olelo, stat = "mean")
dev.off()


# gelman-rubin convergence
grepgo <- grep('beta', colnames(zmPrp[[1]]))
gel_check <- gelman.diag(mcmc.list(
  zmPrp[[1]][, grepgo],zmPrp[[2]][, grepgo],zmPrp[[3]][, grepgo]),multivariate = F)[[1]][,1]
gel_check

grepgo <- grep('b0', colnames(zmPrp[[1]]))
gel_check <- c(gel_check,gelman.diag(mcmc.list(
  zmPrp[[1]][, grepgo],zmPrp[[2]][, grepgo],zmPrp[[3]][, grepgo]),multivariate = F)[[1]][,1])
gel_check
max(gel_check)

temp <- data.frame(gel_check,var=names(gel_check))
temp



# coef plot ---------------------------------------------------------------

labs <- c(
  'Scrapers',
  'Grazers',
  'Browsers',
  'Scrapers:Grazers',
  'Scrapers:Browsers',
  'Grazers:Browsers',
  'Scrapers:Grazers:Browsers',
  expression('Chl-'*italic(a)~'Anomaly Freq' ), 
  expression('Chl-'*italic(a)~'Mean' ), 
  expression('Chl-'*italic(a)~'Anomaly Max' ), 
  'Wave Anomaly Freq',
  'Wave Anomaly Max',
  'Irradiance Mean',
  'Temperature Mean',
  'Temperature SD',
  'Habitat Modification',
  'OSDS Effluent',
  'Sedimentation',
  'Agriculture Runoff',
  'Golf Course Runoff',
  'Urban Runoff',
  'Depth',
  'Rugosity')


png(file=paste0('outputs/Figure4a.png'),height=2800,width=4000,res=300)
par(mfrow=c(1,1),mar=c(5,13,2,2),mgp=c(3.2,1,0),oma=c(1,7,0,0))
plot(beta_out$quantiles[,'50%'],
     rev(seq(from=1,to=23)),
     xlim=c(-1.1,1.16),
     xlab="Effect on log-Calcified:Macroalgae",ylab='',yaxt='n',cex.axis=1.9,cex.lab=1.8,bty='l')

axis(2,at=rev(seq(from=1,to=23)),labels=labs,las=2,cex.axis=1.6,cex.lab=1.8)

colgo <- rev(c(rep(rgb(160,32,240,255,max=255),7),rep('#6ba9ed',8),rep('#e5a913',6),rep('#303740',2)))
rect(-1.5, -1, 1.16, 2.5, border=NA, col=rgb(48,55,64,30,max=255))
rect(-1.5, 2.5, 1.16, 8.5, border=NA, col=rgb(229,169,19,30,max=255))
rect(-1.5, 8.5, 1.16, 16.5, border=NA, col=rgb(107,169,217,30,max=255))
rect(-1.5, 16.5, 1.16, 27.5, border=NA, col=rgb(160,32,240,30,max=255))

abline(v=0,lty=2,lwd=1.5)

plotrix::plotCI(beta_out$quantiles[,'50%'],
                rev(seq(1:23)),
                ui=beta_out$quantiles[,'97.5%'],
                li=beta_out$quantiles[,'2.5%'],
                err='x',add=T,pch=NA,cex=1.5,lwd=3,sfrac=0,
                col=rev(colgo))

plotrix::plotCI(beta_out$quantiles[,'50%'],
                rev(seq(1:23)),
                ui=beta_out$quantiles[,'75%'],
                li=beta_out$quantiles[,'25%'],
                err='x',add=T,pch=NA,cex=1.5,lwd=7,sfrac=0,
                col=rev(colgo))

points(beta_out$quantiles[,'50%'],
       rev(seq(1:23)),
       pch=21,bg=rev(colgo),cex=2.5,lwd=2)

text(-1.15,23.25,'Herbivores',srt=0,cex=1.75,font=2,col=rgb(160,32,240,255,max=255),pos=4)
text(-1.15,15.75,'Oceanography',srt=0,cex=1.75,font=2,col='#6ba9ed',pos=4)
text(-1.15,7.75,'Pollution',srt=0,cex=1.75,font=2,col='#e5a913',pos=4)
text(-1.15,1.75,'Habitat',srt=0,cex=1.75,font=2,col='#303740',pos=4)

dev.off()

# # herbivore interactions --------------------------------------------------

# beta mcmc matrix
grepgo <- c(grep('beta',colnames(zmPrp[[1]])))
beta_d <- rbind(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo])
beta_matrix <- as.matrix(beta_d)

png(file='outputs/FigureS2_revised.png',height=2000,width=2600,res=300)
par(mfrow=c(3,3),mgp=c(2,0.8,0),mar=c(4,3,2,1),oma=c(0,2,0,15),xpd=F)

# predict
n_pred <- 100

x_go <- data.frame(pred_var = c('Scraper_log','Grazer_log','Browser_log'))

X_new <- data.frame(
  Scraper_log = seq(from=min(scaled_pred$Scraper_log),to=max(scaled_pred$Scraper_log),length=n_pred),
  Grazer_log = rep(mean(scaled_pred$Grazer_log), n_pred),
  Browser_log = rep(mean(scaled_pred$Browser_log), n_pred))


X_new_list <- list(
  'no_brow_low_Grazer_log' = X_new,
  'no_brow_high_Grazer_log' = X_new,
  'no_brow_no_Grazer_log' = X_new,
  'low_brow_no_Grazer_log' = X_new,
  'high_brow_no_Grazer_log' = X_new
)
X_new_list[[1]]$Browser_log <- 0; X_new_list[[1]]$Grazer_log <- rep(quantile(scaled_pred$Grazer_log,0.10), n_pred)
X_new_list[[2]]$Browser_log <- 0; X_new_list[[2]]$Grazer_log <- rep(quantile(scaled_pred$Grazer_log,0.90), n_pred)
X_new_list[[3]]$Browser_log <- 0; X_new_list[[3]]$Grazer_log <- 0
X_new_list[[4]]$Browser_log <- rep(quantile(scaled_pred$Browser_log,0.10), n_pred); X_new_list[[4]]$Grazer_log <- 0
X_new_list[[4]]$Browser_log <- rep(quantile(scaled_pred$Browser_log,0.90), n_pred); X_new_list[[4]]$Grazer_log <- 0

pred_out <- list(5)
for(k in 1:5){
  
  X_new <- X_new_list[[k]]
  
  # convert to model matrix
  # X_new_mat <- as.matrix(model.matrix( ~ X_new$Scraper_log + X_new$Grazer_log + X_new$Browser_log
  #                                      + X_new$Scraper_log*X_new$Grazer_log + X_new$Scraper_log*X_new$Browser_log
  #                                      + X_new$Grazer_log*X_new$Browser_log))
  X_new_mat <- as.matrix(model.matrix( ~ X_new$Scraper_log + X_new$Grazer_log + X_new$Browser_log
                                       + X_new$Scraper_log*X_new$Grazer_log + X_new$Scraper_log*X_new$Browser_log
                                       + X_new$Grazer_log*X_new$Browser_log + X_new$Scraper_log*X_new$Grazer_log*X_new$Browser_log))
  X_new_mat <- as.matrix(X_new_mat[,2:ncol(X_new_mat)])
  X_new_mat <- cbind(X_new_mat,matrix(0,nrow=n_pred,ncol=16))
  
  # empty vectors for outputs
  pred <- data.frame(x = X_new$Scraper_log,
                     pred = rep(NA, n_pred), pred_up = rep(NA, n_pred), pred_down = rep(NA, n_pred))
  
  # predict from model output based on quantiles ignoring random intercepts
  for(i in 1:n_pred){
    mcmc_go <- beta_matrix %*% X_new_mat[i,]
    pred$pred[i] <- median(mcmc_go)
    pred$pred_up[i] <- quantile(mcmc_go,0.975)
    pred$pred_down[i] <- quantile(mcmc_go,0.025)
  }
  pred_out[[k]] <- pred
  
}

plot((pred_out[[1]]$x),pred_out[[1]]$pred,ylim=c(-3,4),type='n',
     xlab='Scrapers',ylab='',cex.axis=1.3,cex.lab=1.3,main="Browsers = 0")
abline(h=0,lty=4)
points((pred_out[[1]]$x),pred_out[[1]]$pred,type='l',lwd=3,col='lightblue')
points((pred_out[[2]]$x),pred_out[[2]]$pred,type='l',lwd=3,col='blue')
# legend('topleft',legend=c('low Grazers','high Grazers'),col=c('lightblue','blue'),lwd=1,bty='n',cex=0.7)

plot((pred_out[[3]]$x),pred_out[[3]]$pred,ylim=c(-3,4),type='n',
     xlab='Scrapers',ylab='',cex.axis=1.3,cex.lab=1.3,main="Browsers = 0, Grazers = 0")
abline(h=0,lty=4)
points((pred_out[[3]]$x),pred_out[[3]]$pred,type='l',lwd=3,col='red')

plot((pred_out[[4]]$x),pred_out[[4]]$pred,ylim=c(-3,4),type='n',
     xlab='Scrapers',ylab='',cex.axis=1.3,cex.lab=1.3,main="Grazers = 0")
abline(h=0,lty=4)
points((pred_out[[4]]$x),pred_out[[4]]$pred,type='l',lwd=3,col='lightgreen')
points((pred_out[[5]]$x),pred_out[[5]]$pred,type='l',lwd=3,col='darkgreen')
# legend('topleft',legend=c('low Browsers','high Browsers'),col=c('lightgreen','darkgreen'),lwd=1,bty='n',cex=0.7)


X_new <- data.frame(
  Scraper_log = rep(mean(scaled_pred$Scraper_log), n_pred),
  Grazer_log = seq(from=min(scaled_pred$Grazer_log),to=max(scaled_pred$Grazer_log),length=n_pred),
  Browser_log = rep(mean(scaled_pred$Browser_log), n_pred)
)

X_new_list <- list(
  'no_brow_low_Scraper_log' = X_new,
  'no_brow_high_Scraper_log' = X_new,
  'no_brow_no_Scraper_log' = X_new,
  'low_brow_no_Scraper_log' = X_new,
  'high_brow_no_Scraper_log' = X_new
)
X_new_list[[1]]$Browser_log <- 0; X_new_list[[1]]$Scraper_log <- rep(quantile(scaled_pred$Scraper_log,0.10), n_pred)
X_new_list[[2]]$Browser_log <- 0; X_new_list[[2]]$Scraper_log <- rep(quantile(scaled_pred$Scraper_log,0.90), n_pred)
X_new_list[[3]]$Browser_log <- 0; X_new_list[[3]]$Scraper_log <- 0
X_new_list[[4]]$Browser_log <- rep(quantile(scaled_pred$Browser_log,0.10), n_pred); X_new_list[[4]]$Scraper_log <- 0
X_new_list[[4]]$Browser_log <- rep(quantile(scaled_pred$Browser_log,0.90), n_pred); X_new_list[[4]]$Scraper_log <- 0

pred_out <- list(5)
for(k in 1:5){
  
  X_new <- X_new_list[[k]]
  
  # convert to model matrix
  # X_new_mat <- as.matrix(model.matrix( ~ X_new$Scraper_log + X_new$Grazer_log + X_new$Browser_log
  #                                      + X_new$Scraper_log*X_new$Grazer_log + X_new$Scraper_log*X_new$Browser_log
  #                                      + X_new$Grazer_log*X_new$Browser_log))
  X_new_mat <- as.matrix(model.matrix( ~ X_new$Scraper_log + X_new$Grazer_log + X_new$Browser_log
                                       + X_new$Scraper_log*X_new$Grazer_log + X_new$Scraper_log*X_new$Browser_log
                                       + X_new$Grazer_log*X_new$Browser_log + X_new$Scraper_log*X_new$Grazer_log*X_new$Browser_log))
  X_new_mat <- as.matrix(X_new_mat[,2:ncol(X_new_mat)])
  X_new_mat <- cbind(X_new_mat,matrix(0,nrow=n_pred,ncol=16))
  
  # empty vectors for outputs
  pred <- data.frame(x = X_new$Grazer_log,
                     pred = rep(NA, n_pred), pred_up = rep(NA, n_pred), pred_down = rep(NA, n_pred))
  
  # predict from model output based on quantiles ignoring random intercepts
  for(i in 1:n_pred){
    mcmc_go <- beta_matrix %*% X_new_mat[i,]
    pred$pred[i] <- median(mcmc_go)
    pred$pred_up[i] <- quantile(mcmc_go,0.975)
    pred$pred_down[i] <- quantile(mcmc_go,0.025)
  }
  pred_out[[k]] <- pred
  
}


plot((pred_out[[1]]$x),pred_out[[1]]$pred,ylim=c(-3,4),type='n',
     xlab='Grazers',ylab='',cex.axis=1.3,cex.lab=1.3,main="Browsers = 0")
abline(h=0,lty=4)
points((pred_out[[1]]$x),pred_out[[1]]$pred,type='l',lwd=3,col='pink')
points((pred_out[[2]]$x),pred_out[[2]]$pred,type='l',lwd=3,col='red')
# legend('topleft',legend=c('low Scrapers','high Scrapers'),col=c('pink','red'),lwd=1,bty='n',cex=0.7)

plot((pred_out[[3]]$x),pred_out[[3]]$pred,ylim=c(-3,4),type='n',
     xlab='Grazers',ylab='',cex.axis=1.3,cex.lab=1.3,main="Browsers = 0, Scrapers = 0")
abline(h=0,lty=4)
points((pred_out[[3]]$x),pred_out[[3]]$pred,type='l',lwd=3,col='blue')

plot((pred_out[[4]]$x),pred_out[[4]]$pred,ylim=c(-3,4),type='n',
     xlab='Grazers',ylab='',cex.axis=1.3,cex.lab=1.3,main="Scrapers = 0")
abline(h=0,lty=4)
points((pred_out[[4]]$x),pred_out[[4]]$pred,type='l',lwd=3,col='lightgreen')
points((pred_out[[5]]$x),pred_out[[5]]$pred,type='l',lwd=3,col='darkgreen')
# legend('topleft',legend=c('low Browsers','high Browsers'),col=c('lightgreen','darkgreen'),lwd=1,bty='n',cex=0.7)




X_new <- data.frame(
  Scraper_log = rep(mean(scaled_pred$Scraper_log), n_pred),
  Grazer_log = rep(mean(scaled_pred$Grazer_log), n_pred),
  Browser_log = seq(from=min(scaled_pred$Browser_log),to=max(scaled_pred$Browser_log),length=n_pred)
)

X_new_list <- list(
  'no_Grazer_log_low_Scraper_log' = X_new,
  'no_Grazer_log_high_Scraper_log' = X_new,
  'no_Grazer_log_no_Scraper_log' = X_new,
  'low_Grazer_log_no_Scraper_log' = X_new,
  'high_Grazer_log_no_Scraper_log' = X_new
)
X_new_list[[1]]$Grazer_log <- 0; X_new_list[[1]]$Scraper_log <- rep(quantile(scaled_pred$Scraper_log,0.10), n_pred)
X_new_list[[2]]$Grazer_log <- 0; X_new_list[[2]]$Scraper_log <- rep(quantile(scaled_pred$Scraper_log,0.90), n_pred)
X_new_list[[3]]$Grazer_log <- 0; X_new_list[[3]]$Scraper_log <- 0
X_new_list[[4]]$Grazer_log <- rep(quantile(scaled_pred$Grazer_log,0.10), n_pred); X_new_list[[4]]$Scraper_log <- 0
X_new_list[[4]]$Grazer_log <- rep(quantile(scaled_pred$Grazer_log,0.90), n_pred); X_new_list[[4]]$Scraper_log <- 0

pred_out <- list(5)
for(k in 1:5){
  
  X_new <- X_new_list[[k]]
  
  # convert to model matrix
  # X_new_mat <- as.matrix(model.matrix( ~ X_new$Scraper_log + X_new$Grazer_log + X_new$Browser_log
  #                                      + X_new$Scraper_log*X_new$Grazer_log + X_new$Scraper_log*X_new$Browser_log
  #                                      + X_new$Grazer_log*X_new$Browser_log))
  X_new_mat <- as.matrix(model.matrix( ~ X_new$Scraper_log + X_new$Grazer_log + X_new$Browser_log
                                       + X_new$Scraper_log*X_new$Grazer_log + X_new$Scraper_log*X_new$Browser_log
                                       + X_new$Grazer_log*X_new$Browser_log + X_new$Scraper_log*X_new$Grazer_log*X_new$Browser_log))
  X_new_mat <- as.matrix(X_new_mat[,2:ncol(X_new_mat)])
  X_new_mat <- cbind(X_new_mat,matrix(0,nrow=n_pred,ncol=16))
  
  # empty vectors for outputs
  pred <- data.frame(x = X_new$Browser_log,
                     pred = rep(NA, n_pred), pred_up = rep(NA, n_pred), pred_down = rep(NA, n_pred))
  
  # predict from model output based on quantiles ignoring random intercepts
  for(i in 1:n_pred){
    mcmc_go <- beta_matrix %*% X_new_mat[i,]
    pred$pred[i] <- median(mcmc_go)
    pred$pred_up[i] <- quantile(mcmc_go,0.975)
    pred$pred_down[i] <- quantile(mcmc_go,0.025)
  }
  pred_out[[k]] <- pred
  
}


plot((pred_out[[1]]$x),pred_out[[1]]$pred,ylim=c(-3,4),type='n',
     xlab='Browsers',ylab='',cex.axis=1.3,cex.lab=1.3,main='Grazers = 0')
abline(h=0,lty=4)
# segments(0,0,max((pred_out[[1]]$x)),lty=2)
points((pred_out[[1]]$x),pred_out[[1]]$pred,type='l',lwd=3,col='pink')
points((pred_out[[2]]$x),pred_out[[2]]$pred,type='l',lwd=3,col='red')
# legend('topleft',legend=c('low Scrapers','high Scrapers'),col=c('pink','red'),lwd=1,bty='n',cex=0.7)

plot((pred_out[[3]]$x),pred_out[[3]]$pred,ylim=c(-3,4),type='n',
     xlab='Browsers',ylab='',cex.axis=1.3,cex.lab=1.3,main='Grazers = 0, Scrapers = 0')
abline(h=0,lty=4)
# segments(0,0,max((pred_out[[3]]$x)),lty=2)
points((pred_out[[3]]$x),pred_out[[3]]$pred,type='l',lwd=3,col='darkgreen')

plot((pred_out[[4]]$x),pred_out[[4]]$pred,ylim=c(-3,4),type='n',
     xlab='Browsers',ylab='',cex.axis=1.3,cex.lab=1.3,main='Scrapers = 0')
abline(h=0,lty=4)
# segments(0,0,max((pred_out[[4]]$x)),lty=2)
points((pred_out[[4]]$x),pred_out[[4]]$pred,type='l',lwd=3,col='lightblue')
points((pred_out[[5]]$x),pred_out[[5]]$pred,type='l',lwd=3,col='blue')
# legend('topleft',legend=c('low Grazers','high Grazers'),col=c('lightblue','blue'),lwd=1,bty='n',cex=0.7)

mtext('log Calcified:Macroalgae',side=2,outer=T)

# legend(5,1,legend=c('low Grazers','high Grazers','low Scrapers','high Scrapers','low Browsers','high Browsers'),col=c('lightblue','blue','pink','red','lightgreen','darkgreen'),lwd=1,bty='n',cex=1, xpd=T)
dev.off()



png(file="outputs/FigS2_legend.png",height=2000,width=2600,res=300)
plot(seq(1:3),seq(1:3),pch=NA,xaxt="n",yaxt="n",ylab="",xlab="",bty="n")
legend("topright",legend=c('low Grazers','high Grazers','low Scrapers','high Scrapers','low Browsers','high Browsers'),col=c('lightblue','blue','pink','red','lightgreen','darkgreen'),lwd=5,bty='n',cex=2, xpd=T)
dev.off()