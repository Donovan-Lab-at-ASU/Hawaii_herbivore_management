# -----------------------------------------------------------------------
# Hierarchical Bayesian model of drivers of herbivore biomass in Hawai‘i
# -----------------------------------------------------------------------
# Model and code by Mary Donovan - marydonovan@asu.edu
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

var_in <- commandArgs(trailingOnly = TRUE) # this code is written to run as a loop on a high performance computer that utilizes multiple cores
var_opts <- c('TotalHerbivore','Scraper','Grazer','Browser')
varuse <- var_opts[as.numeric(var_in)]
varuse

# Initialization ----------------------------------------------------------
library(rjags) # v4.10 linked to JAGS 4.3.0
library(parallel) # v4.0.3
library(dplyr) # v1.0.7
library(tidyr) # v1.1.2

# data --------------------------------------------------------------------
data <- read.csv("data/donovan_etal_HIMARCdata_for_herbivore_publication.csv")

# model -------------------------------------------------------------------
cat("model{
    # Likelihood
    for(i in 1:N){
      y[i] ~ dgamma(shape, shape / exp(mu[i]))
      mu[i] <- b0.y[year[i]] + b0.d[dataset[i]] + b0.moku[moku[i]] + inprod(beta,X[i,]) + beta_depth*depth[i] + beta_rugosity*rugosity[i]
      y.new[i] ~ dgamma(shape, shape / exp(mu[i]))
    }
    shape ~ dunif(0,100)
    for(b in 1:B){ beta[b] ~ dnorm(0,0.001) }
    beta_depth ~ dnorm(0,0.001)
    beta_rugosity ~ dnorm(0,0.001)

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
    
    }",file='herbivore_driver_model.jags') # write out model to working directory


# set up ------------------------------------------------------------------
# set up final inputs

response <- data[,varuse]
response[response==0] <- 0.001 # gamma isn't defined at zero

# moku factor levels
data$moku <- as.factor(as.character(data$moku))
mok_levs <- levels(data$moku)
nmoku <- length(unique(data$moku))

# order habitat factor
levels(as.factor(data$habitat))
data$habitat2 <- factor(data$habitat,levels(as.factor(data$habitat))[c(4,1,2,3)])
levels(as.factor(data$habitat2))

# create predictor matrix
predict.go <- data %>% select(OTP_CHL_ANOM_F:LBSP2_Urban_runoff)
X <- cbind(model.matrix(~1+ data$habitat2), as.matrix(data.frame(predict.go)))
Xsave <- X
for(i in 5:ncol(X)) X[,i] <- scale(X[,i])[,1]
write.csv(colnames(X),paste0('model_out/beta_names_',varuse,'.csv'),row.names = F)

# depth and rugosity inputs
depth.go <- scale(data$depth_predicted)[,1]
# summary(depth.go)

rugosity.go <- scale(data$rugosity_predicted)[,1]
# summary(rugosity.go)

# export conversion from scaled predictors
predictor_scale_ab <- data.frame(pred_ord=seq(1:28),a=c(rep(1,4),rep(NA,24)),b=c(rep(1,4),rep(NA,24)))
for(i in 5:26){
  pred_lm <- lm(X[,i]~Xsave[,i])
  predictor_scale_ab$a[i] <- summary(pred_lm)$coefficients[1,1]
  predictor_scale_ab$b[i] <- summary(pred_lm)$coefficients[2,1]
}
pred_lm <- lm(depth.go~data$depth_predicted)
predictor_scale_ab$a[27] <- summary(pred_lm)$coefficients[1,1]
predictor_scale_ab$b[27] <- summary(pred_lm)$coefficients[2,1]
pred_lm <- lm(rugosity.go~data$rugosity_predicted)
predictor_scale_ab$a[28] <- summary(pred_lm)$coefficients[1,1]
predictor_scale_ab$b[28] <- summary(pred_lm)$coefficients[2,1]

write.csv(predictor_scale_ab, paste0('model_out/predictor_scale_',varuse,'.csv'),row.names=F)

# final factor levels
data$year.f <- as.numeric(as.factor(as.character(data$Year)))
data$moku.f <- as.numeric(as.factor(as.character(data$moku)))
data$dataset.f <- as.numeric(as.factor(as.character(data$Dataset)))

year <- data %>% distinct(year.f)
dataset <- data %>% distinct(dataset.f)
moku <- data %>% distinct(moku.f) 

# export moku index
write.csv(data %>% distinct(moku.f,moku), paste0('model_out/moku_order_',varuse,'.csv'),row.names=F)


# run ---------------------------------------------------------------------
nX <- ncol(X)
N=length(response)

# model initial conditions
initFunc <- function(){return(list(
  shape=0.65,
  beta=rnorm(nX,0,1),
  beta_depth=0.5
))}

# model inputs
n.adapt <- 500; n.update <- 2000; n.iter <- 5000

mdat <- list(N=length(response),
             y=response,
             X=X,
             B=ncol(X),
             year=data$year.f,
             dataset=data$dataset.f,
             moku=data$moku.f,
             Y=length(unique(data$year.f)),
             D=length(unique(data$dataset.f)),
             M=length(unique(data$moku.f)),
             depth=depth.go,
             rugosity=rugosity.go
)

# save model data
write.csv(data.frame(
  response=response,
  year=data$Year,
  year.f=data$year.f,
  dataset=data$Dataset,
  dataset.f=data$dataset.f,
  moku=data$moku,
  moku.f=data$moku.f,
  X,
  depth=depth.go,
  rugosity.go
), paste0('model_out/data_in','_',varuse,'.csv'), row.names=F)

# run chains in parallel
cl <- makeCluster(3)
clusterExport(cl, c("mdat","initFunc","n.adapt","n.update","n.iter","nX","N"))

out <- clusterEvalQ(cl,{
  library(rjags)
  mod <- jags.model("herbivore_driver_model.jags", data= mdat, n.chains=1, n.adapt=n.adapt,inits = initFunc())
  update(mod, n.iter = n.update)
  zmCore <- coda.samples(mod, c("pval.mean","pval.sd","R2","shape",
                                "beta","beta_depth","beta_rugosity",'sigma_year','sigma_moku','mu_year','mu_moku',
                                "b0.y","b0.d","b0.moku","b0","y.new"
  ),n.iter=n.iter, n.thin=1)
  return(as.mcmc(zmCore))
})
zmPrp <- mcmc.list(out)
stopCluster(cl)

# export ------------------------------------------------------------------
saveRDS(zmPrp, file=paste0('model_out/',varuse,"_gamma_model_out.Rdata"))


grepgo <- grep('y.new', colnames(zmPrp[[1]]))
ynew_summary <- summary(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))

ynew_out <- ynew_summary$quantiles
ynew_out <- data.frame(ynew_out,response)
write.csv(ynew_out,paste0('model_out/',varuse,"_gamma_","ynew.csv"), row.names = T)

grepgo <- grep('beta',colnames(zmPrp[[1]]))
beta_out <- summary(mcmc.list(zmPrp[[1]][,grepgo],zmPrp[[2]][,grepgo],zmPrp[[3]][,grepgo]))
write.csv(
  cbind(beta_out$quantiles,c(colnames(X),'depth','rugosity'))
  , paste0('model_out/',varuse,"_gamma_","beta_quantiles.csv"), row.names = T)

grepgo <- grep(c('pval.mean'),colnames(zmPrp[[1]]))
pval_mean_sum <- summary(mcmc.list(zmPrp[[1]][,grepgo],
                                   zmPrp[[2]][,grepgo],
                                   zmPrp[[3]][,grepgo]))
grepgo <- grep(c('pval.sd'),colnames(zmPrp[[1]]))
pval_sd_sum <- summary(mcmc.list(zmPrp[[1]][,grepgo],
                                 zmPrp[[2]][,grepgo],
                                 zmPrp[[3]][,grepgo]))
grepgo <- grep(c('R2'),colnames(zmPrp[[1]]))
R2_sum <- summary(mcmc.list(zmPrp[[1]][,grepgo],
                            zmPrp[[2]][,grepgo],
                            zmPrp[[3]][,grepgo]))

pos_out <- data.frame(
  value = c('pval_mean','pval_sd','R2'),
  mean = c(pval_mean_sum$statistics['Mean'],pval_sd_sum$statistics['Mean'],R2_sum$statistics['Mean'])
)
pos_out
write.csv(pos_out, paste0('model_out/',varuse,"_gamma_",'model_checks','.csv'), row.names = F)

grepgo <- grep('beta', colnames(zmPrp[[1]]))
gel_check <- gelman.diag(mcmc.list(
  zmPrp[[1]][, grepgo],zmPrp[[2]][, grepgo],zmPrp[[3]][, grepgo]),multivariate = F)[[1]][,1]
# gel_check

grepgo <- grep('b0', colnames(zmPrp[[1]]))
gel_check <- c(gel_check,gelman.diag(mcmc.list(
  zmPrp[[1]][, grepgo],zmPrp[[2]][, grepgo],zmPrp[[3]][, grepgo]),multivariate = F)[[1]][,1])
# gel_check
# max(gel_check)

temp <- data.frame(gel_check,var=names(gel_check))
write.csv(temp, paste0('model_out/',varuse,"_gamma_",'gelcheck.csv'), row.names = T)
