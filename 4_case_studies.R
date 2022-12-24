# -----------------------------------------------------------------------
# Case studies of herbivore management in Hawai‘i
# -----------------------------------------------------------------------
# Model and code by Mary Donovan - marydonovan@asu.edu
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

library(plotrix) # v3.7.8
library(dplyr) # v1.0.7
library(boot) # v1.3.25
library(broom) # v0.7.3


# data --------------------------------------------------------------------

# Kahekili
kahe <- read.csv("data/kahekili.csv")
kahe$not_zero <- ifelse(kahe$Parrotfish==0,0,1)
kahe$mgmt <- ifelse(kahe$Year < 2010, '1_before', '2_after'); kahe$mgmt <- as.factor(kahe$mgmt)
kahe <- kahe %>% filter(!Year==2010)

# West Hawaii
whap_parrot <- read.csv("data/westhawaii.csv")
whap_parrot$Date <- as.Date(whap_parrot$Date, format = "%Y-%m-%d")
whap_parrot$RoundID_fac <- as.numeric(as.factor(whap_parrot$RoundID))
whap_parrot$mgmt_area_fac <- as.factor(ifelse(whap_parrot$mgmt_area=='Open','1_Open','2_Closed'))
whap_parrot$mmgmt_gear_fac <- as.factor(ifelse(whap_parrot$mgmt_gear=='Before','1_Before','2_After'))
whap_parrot$Site <- as.factor(whap_parrot$Site)
whap_parrot$bio_use <- ifelse(whap_parrot$biomass_g_m2<1,1,whap_parrot$biomass_g_m2)
whap_parrot$RoundID_fac_s <- scale(whap_parrot$RoundID_fac)[,1]

# Maui
maui_fahu <- read.csv('Data/fahu_parrot.csv') #created in github/HIMARC/indicators/fahu_parrotfish.R
maui_fahu_op <- maui_fahu
maui_fahu_op$zero_one <- ifelse(maui_fahu_op$parrot==0,0,1)
maui_fahu_op$mgmt <- as.factor(maui_fahu_op$open_protected)
maui_fahu_op <- subset(maui_fahu_op, Year==2018 | Year==2019)
maui_fahu_op$year_num <- as.numeric(as.factor(maui_fahu_op$Year))
maui_fahu_op$Site_code <- as.factor(maui_fahu_op$Site_code)


# hurdle functions --------------------------------------------------------

hurdle_fun <- function(data, i){
  dat_boot <- data[i, ]
  m1 <- glm(zero_one ~ mgmt, data = dat_boot, family = binomial(link = logit))
  m2 <- glm(parrot ~ mgmt, data = subset(dat_boot, zero_one == 1), family = Gamma(link = log))
  open_out <- exp(log(plogis(coef(m1)[1])) + log(exp(coef(m2)[1])))
  closed_out <- exp(log(plogis(coef(m1)[1]+coef(m1)[2])) + log(exp(coef(m2)[1]+coef(m1)[2])))
  return(c(open_out,closed_out))
}

hurdle_fun1 <- function(data, i){
  dat_boot <- data[i, ]
  m1 <- glm(zero_one ~ mgmt, data = dat_boot, family = binomial(link = logit))
  m2 <- glm(parrot ~ mgmt, data = subset(dat_boot, zero_one == 1), family = Gamma(link = log))
  open_out <- exp(log(plogis(coef(m1)[1])) + log(exp(coef(m2)[1])))
  closed_out <- exp(log(plogis(coef(m1)[1]+coef(m1)[2])) + log(exp(coef(m2)[1]+coef(m1)[2])))
  return(open_out)
}

nothurdle_fun <- function(data, i){
  dat_boot <- data[i, ]
  # m1 <- glm(zero_one ~ mgmt, data = dat_boot, family = binomial(link = logit))
  m2 <- glm(parrot ~ mgmt, data = subset(dat_boot, zero_one == 1), family = Gamma(link = log))
  open_out <- exp(coef(m2)[1])
  closed_out <- exp(coef(m2)[1]+coef(m2)[2])
  return(c(open_out,closed_out))
}

# analysis ----------------------------------------------------------------

temp <- maui_fahu_op
temp$mgmt <- ifelse(maui_fahu_op$mgmt=='OPEN',"2_open","1_protected")
maui_boot <- boot(temp, hurdle_fun, R = 1000)
maui_ci <- tidy(maui_boot,conf.int=TRUE,conf.method="perc",conf=0.95)

temp <- kahe
temp$mgmt <- ifelse(temp$mgmt=='1_before',"2_before","1_after")
temp$zero_one <- ifelse(temp$Parrotfish==0,0,1)
temp$parrot <- temp$Parrotfish
kahe_boot <- boot(temp, hurdle_fun, R = 1000)
kahe_ci <- tidy(kahe_boot,conf.int=TRUE,conf.method="perc",conf=0.95)

temp <- subset(whap_parrot, mgmt_area=="Open" & Year < 2011 | mgmt_area=="Open" & Year > 2015)
# temp <- subset(temp, biomass_g_m2 < 200)
temp$mgmt <- as.factor(ifelse(temp$Year < 2011, '2_before','1_after'))
# temp$mgmt <- as.factor(ifelse(temp$mgmt_area == 'MPA',"2_closed",'1_open'))
temp$zero_one <- ifelse(temp$biomass_g_m2==0,0,1)
temp$parrot <- temp$biomass_g_m2
wh_boot <- boot(temp, hurdle_fun, R = 1000)
wh_ci <- tidy(wh_boot,conf.int=TRUE,conf.method="perc",conf=0.95)
wh_ci

temp <- subset(whap_parrot, mgmt_area=="MPA"| mgmt_area=="FRA")
temp <- subset(temp, !Site_Name=="Kealakekua Bay MLCD")
temp <- subset(temp, !Site_Name=="Old Kona Airport MLCD")
temp <- subset(temp, Year < 2011 | Year > 2015)
# temp <- subset(temp, biomass_g_m2 < 200)
temp$mgmt <- as.factor(ifelse(temp$Year < 2011, '1_before','2_after'))
# temp$mgmt <- as.factor(ifelse(temp$mgmt_area == 'MPA',"2_closed",'1_open'))
temp$zero_one <- ifelse(temp$biomass_g_m2==0,0,1)
temp$parrot <- temp$biomass_g_m2
wh_boot1 <- boot(temp, nothurdle_fun, R = 1000)
wh_ci1 <- tidy(wh_boot1,conf.int=TRUE,conf.method="perc",conf=0.95)
wh_ci1

temp <- subset(whap_parrot, Site_Name=="Kealakekua Bay MLCD" | Site_Name=="Old Kona Airport MLCD")
temp <- subset(temp, Year < 2011 | Year > 2015)
# temp <- subset(temp, biomass_g_m2 < 200)
temp$mgmt <- as.factor(ifelse(temp$Year < 2011, '1_before','2_after'))
# temp$mgmt <- as.factor(ifelse(temp$mgmt_area == 'MPA',"2_closed",'1_open'))
temp$zero_one <- ifelse(temp$biomass_g_m2==0,0,1)
temp$parrot <- temp$biomass_g_m2
wh_boot1 <- boot(temp, nothurdle_fun, R = 1000)
wh_ci2 <- tidy(wh_boot1,conf.int=TRUE,conf.method="perc",conf=0.95)
wh_ci2


# figure ------------------------------------------------------------------

png(file=paste0('outputs/Figure5.png'),height=1400,width=2800,res=300)
par(mfrow=c(1,1),mar=c(4,5,2,2),mgp=c(2.2,0.7,0),oma=c(1,0,0,0))
plot(c(0.3,0.5,1,1.5,2.0),seq(1:5),type='n',ylim=c(0,100),bty="o",xaxt="n",ylab=expression("Parrotfish Biomass"~~bgroup("(",'g '*m^{-2},")")),xlab="")
plotrix::plotCI(c(0.38,0.52),maui_ci$statistic[c(2,1)],ui=maui_ci$conf.high[c(2,1)],li=maui_ci$conf.low[c(2,1)],add=T,pch=c(21,23),col='black',pt.bg="#8BBAAA",lwd=1.5,sfrac=0)
plotrix::plotCI(c(0.77,.93),kahe_ci$statistic[c(2,1)],ui=kahe_ci$conf.high[c(2,1)],li=kahe_ci$conf.low[c(2,1)],add=T,pch=c(21,23),col='black',pt.bg="#8BBAAA",lwd=1.5,sfrac=0)
plotrix::plotCI(c(1.18,1.34),wh_ci$statistic[c(2,1)],ui=wh_ci$conf.high[c(2,1)],li=wh_ci$conf.low[c(2,1)],add=T,pch=21,col='black',pt.bg="#8BBAAA",lwd=1.5,sfrac=0)
plotrix::plotCI(c(1.5,1.66),wh_ci1$statistic[c(1,2)],ui=wh_ci1$conf.high[c(1,2)],li=wh_ci1$conf.low[c(1,2)],add=T,pch=23,col='black',pt.bg="#8BBAAA",lwd=1.5,sfrac=0)
plotrix::plotCI(c(1.82,1.98),wh_ci2$statistic[c(1,2)],ui=wh_ci2$conf.high[c(1,2)],li=wh_ci2$conf.low[c(1,2)],add=T,pch=23,col='black',pt.bg="#8BBAAA",lwd=1.5,sfrac=0)

axis(1,at=c(0.36,0.52),label=c('Out','In'))
axis(1,at=c(0.44),label=c('Maui'),line=1,tick=F)
axis(1,at=c(0.36,0.52),label=c('',''),line=1.7,tick=T,tck=0)

axis(1,at=c(0.77,.93),label=c('Before','After'))
axis(1,at=c(0.85),label=c('Kahekili'),line=1,tick=F)
axis(1,at=c(0.77,.93),label=c('',''),line=1.7,tick=T,tck=0)

axis(1,at=c(1.18,1.34,1.5,1.66,1.82,1.98),label=c('Before','After','Before','After','Before','After'))
axis(1,at=c(1.26,1.58,1.9),label=c('Out','FMA','MLCD'),line=1,tick=F)
axis(1,at=c(1.58),label=c('West Hawai‘i'),line=2,tick=F)
axis(1,at=c(1.18,1.34),label=c('',''),line=1.7,tick=T,tck=0)
axis(1,at=c(1.5,1.66),label=c('',''),line=1.7,tick=T,tck=0)
axis(1,at=c(1.82,1.98),label=c('',''),line=1.7,tick=T,tck=0)
axis(1,at=c(1.26,1.9),label=c('',''),line=2.7,tick=T,tck=0)

abline(v=0.645)
abline(v=1.055)
text(0.24-.02,99,"A) Marine Reserves",pos=4)
text(0.645-.0105,99,"B) Herbivore Reserve",pos=4)
text(1.055-.01,99,"C) Spearfishing Ban",pos=4)

dev.off()

# supplemental table ------------------------------------------------------

temp <- maui_fahu_op
temp$mgmt <- ifelse(maui_fahu_op$mgmt=='OPEN',"1_open","2_protected")
m1 <- glm(zero_one ~ mgmt, data = temp, family = binomial(link = logit))
summary(m1)
m2 <- glm(parrot ~ mgmt, data = subset(temp, zero_one == 1), family = Gamma(link = log))
summary(m2)

temp <- kahe
temp$mgmt <- ifelse(temp$mgmt=='1_before',"1_before","2_after")
temp$zero_one <- ifelse(temp$Parrotfish==0,0,1)
temp$parrot <- temp$Parrotfish
m1 <- glm(zero_one ~ mgmt, data = temp, family = binomial(link = logit))
summary(m1)
m2 <- glm(parrot ~ mgmt, data = subset(temp, zero_one == 1), family = Gamma(link = log))
summary(m2)


temp <- subset(whap_parrot, mgmt_area=="Open" & Year < 2011 | mgmt_area=="Open" & Year > 2015)
temp$mgmt <- as.factor(ifelse(temp$Year < 2011, '1_before','2_after'))
temp$zero_one <- ifelse(temp$biomass_g_m2==0,0,1)
temp$parrot <- temp$biomass_g_m2
m1 <- glm(zero_one ~ mgmt, data = temp, family = binomial(link = logit))
summary(m1)
m2 <- glm(parrot ~ mgmt, data = subset(temp, zero_one == 1), family = Gamma(link = log))
summary(m2)

temp <- subset(whap_parrot, mgmt_area=="MPA"| mgmt_area=="FRA")
temp <- subset(temp, !Site_Name=="Kealakekua Bay MLCD")
temp <- subset(temp, !Site_Name=="Old Kona Airport MLCD")
temp <- subset(temp, Year < 2011 | Year > 2015)
# temp <- subset(temp, biomass_g_m2 < 200)
temp$mgmt <- as.factor(ifelse(temp$Year < 2011, '1_before','2_after'))
# temp$mgmt <- as.factor(ifelse(temp$mgmt_area == 'MPA',"2_closed",'1_open'))
temp$zero_one <- ifelse(temp$biomass_g_m2==0,0,1)
temp$parrot <- temp$biomass_g_m2
m1 <- glm(zero_one ~ mgmt, data = temp, family = binomial(link = logit))
summary(m1)
m2 <- glm(parrot ~ mgmt, data = subset(temp, zero_one == 1), family = Gamma(link = log))
summary(m2)

p <- predict(m2, newdata = data.frame(mgmt=c('1_before','2_after')))
p[2]/p[1]


temp <- subset(whap_parrot, Site_Name=="Kealakekua Bay MLCD" | Site_Name=="Old Kona Airport MLCD")
temp <- subset(temp, Year < 2011 | Year > 2015)
# temp <- subset(temp, biomass_g_m2 < 200)
temp$mgmt <- as.factor(ifelse(temp$Year < 2011, '1_before','2_after'))
# temp$mgmt <- as.factor(ifelse(temp$mgmt_area == 'MPA',"2_closed",'1_open'))
temp$zero_one <- ifelse(temp$biomass_g_m2==0,0,1)
temp$parrot <- temp$biomass_g_m2
m2 <- glm(parrot ~ mgmt, data = subset(temp, zero_one == 1), family = Gamma(link = log))
summary(m2)
p <- predict(m2, newdata = data.frame(mgmt=c('1_before','2_after')))
p[2]/p[1]
