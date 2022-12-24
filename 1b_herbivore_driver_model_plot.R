# -----------------------------------------------------------------------
# Outputs of model on drivers of herbivore biomass in Hawai‘i
# -----------------------------------------------------------------------
# Model and code by Mary Donovan - marydonovan@asu.edu
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

library(plotrix) # v3.7.8
library(bayesplot) # v1.8.0
library(ggplot2) # v3.3.5


# Figure 1 ----------------------------------------------------------------

mod_dest <- 'model_out'
varuse <- 'TotalHerbivore'

beta_out <- read.csv(paste0(mod_dest,'/',varuse,'_gamma_beta_quantiles.csv'))

labs <- c(
  'Intercept',
  'Boulder',
  'Reef',
  'Pavement',
  expression('Chl-'*italic(a)~'Anomaly Freq' ), 
  expression('Chl-'*italic(a)~'Mean' ), 
  expression('Chl-'*italic(a)~'Anomaly Max' ), 
  'Wave Anomaly Freq',
  'Wave Anomaly Max',
  'Irradiance Mean',
  'Temperature Mean',
  'Temperature SD',
  'Commercial Line Fishing',
  'Commercial Spear Fishing',
  'Combined Boat Net Fishing',
  'Noncomm Boat Spear Fishing',
  'Noncomm Shore Line Fishing',
  'Noncomm Shore Spear Fishing',
  'Noncomm Shore Net Fishing',
  'Aquarium Fishing',
  'Habitat Modification',
  'OSDS Effluent',
  'Sedimentation',
  'Agriculture Runoff',
  'Golf Course Runoff',
  'Urban Runoff',
  'Depth',
  'Rugosity')

png(file=paste0('outputs/Figure1.png'),height=3500,width=4000,res=300)
par(mfrow=c(1,1),mar=c(4,13,2,2),mgp=c(2.7,1,0),oma=c(1,7,0,0))

plot(beta_out$X50.[c(28,27,2:12,21:26,20,13:19)],
     seq(from=1,to=27),
     xlim=c(-0.53,max(beta_out$X97.5.[2:nrow(beta_out)])),
     xlab=expression(beta),ylab='',yaxt='n',cex.axis=1.5,cex.lab=2.1,bty='l')

axis(2,at=seq(from=1,to=27),labels=labs[c(28,27,2:12,21:26,20,13:19)],las=2,cex.axis=1.3,cex.lab=1.8)

abline(v=0,lty=2,lwd=1.5)
segments(0,27.5,0,27.5)

colgo <- rev(c(rep('#d84d32',8),rep('#e5a913',6),rep('#6ba9ed',8),rep('#303740',5)))
rect(-0.6, -1, 0.37, 5.5, border=NA, col=rgb(48,55,64,30,max=255))
rect(-0.6, 5.5, 0.37, 13.5, border=NA, col=rgb(107,169,217,30,max=255))
rect(-0.6, 13.5, 0.37, 19.5, border=NA, col=rgb(229,169,19,30,max=255))
rect(-0.6, 19.5, 0.37, 27.5, border=NA, col=rgb(216,77,50,30,max=255))

plotrix::plotCI(beta_out$X50.[c(28,27,2:12,21:26,20,13:19)],
                seq(1:27),
                ui=beta_out$X97.5.[c(28,27,2:12,21:26,20,13:19)],
                li=beta_out$X2.5.[c(28,27,2:12,21:26,20,13:19)],
                err='x',add=T,pch=NA,col=colgo,cex=1.5,lwd=3,sfrac=0)

plotrix::plotCI(beta_out$X50.[c(28,27,2:12,21:26,20,13:19)],
                seq(1:27),
                ui=beta_out$X75.[c(28,27,2:12,21:26,20,13:19)],
                li=beta_out$X25.[c(28,27,2:12,21:26,20,13:19)],
                err='x',add=T,pch=NA,col=colgo,cex=1.5,lwd=7,sfrac=0)

points(beta_out$X50.[c(28,27,2:12,21:26,20,13:19)],
       seq(1:27),
       pch=21,bg=colgo,cex=2.5,lwd=2)

text(-0.53,23.5,'Fishing',srt=90,cex=1.75,font=2,col='#d84d32')
text(-0.53,16.5,'Pollution',srt=90,cex=1.75,font=2,col='#e5a913')
text(-0.53,9.5,'Oceanography',srt=90,cex=1.75,font=2,col='#6ba9ed')
text(-0.531,3,'Habitat',srt=90,cex=1.75,font=2,col='#303740')

dev.off()


# posterior checks --------------------------------------------------------

data <- read.csv("model_out/data_in_TotalHerbivore.csv")
y <- data$response

mod_out <- readRDS("model_out/TotalHerbivore_gamma_model_out.Rdata")

grepgo <- grep('y.new', colnames(mod_out[[1]]))
yrep <- rbind(mod_out[[1]][,grepgo],mod_out[[2]][,grepgo],mod_out[[3]][,grepgo])

# change the bayesplot theme to theme_minimal and save the old theme
old <- bayesplot_theme_set(theme_minimal())
ppc_stat_grouped(y, yrep, group = data$moku, stat = "mean")

# modify moku label
moku_label <- read.csv('data/moku_label-4.csv', encoding = "UTF-8")
data <- dplyr::left_join(data,moku_label,by='moku')
png(file="outputs/FigureS5.png",height=2500,width=3500,res=300)
ppc_stat_grouped(y, yrep, group = data$Moku_olelo, stat = "mean")
dev.off()

# model checks ------------------------------------------------------------

tot_gel <- read.csv('model_out/TotalHerbivore_gamma_gelcheck.csv')

browser_gel <- read.csv('model_out/Browser_gamma_gelcheck.csv')

grazer_gel <- read.csv('model_out/Grazer_gamma_gelcheck.csv')

scraper_gel <- read.csv('model_out/Scraper_gamma_gelcheck.csv')

gel_out <- data.frame(var=tot_gel$var, total_herbivore=round(tot_gel$gel_check,2), browser=round(browser_gel$gel_check,2), grazer=round(grazer_gel$gel_check,2), scraper=round(scraper_gel$gel_check,2))
head(gel_out)
write.csv(gel_out,"outputs/TableS6.csv",row.names = F)
