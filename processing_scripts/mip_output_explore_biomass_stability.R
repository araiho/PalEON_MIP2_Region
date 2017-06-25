library(maps)
library(corrplot)
library(RCurl)
library(sp)
library(udunits2)
library(fields)

calc.second.deriv <- function(biomassCI,h){
  T <- length(biomassCI) - h
  t <- h + 1
  second.deriv <- sum(((biomassCI[(t:T)+h]-2*biomassCI[(t:T)]+biomassCI[(t:T)-h])/((h*100)^2))^2)
  return(second.deriv)
}

calc.second.deriv.1yr <- function(biomassCI,h){
  T <- length(biomassCI) - h
  t <- h + 1
  second.deriv <- sum(((biomassCI[(t:T)+h]-2*biomassCI[(t:T)]+biomassCI[(t:T)-h])/((h)^2))^2)
  return(second.deriv)
}


map.colors <- function(x,breaks,TITLE){
  breaks <-  c(min(x),quantile(x,c(.025,.25,.5,.75,.975)),max(x))
  colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))
  
  data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
  points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
  legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
  title(TITLE)
}

#### ALL BIOMASS FILES IN kgC/m^2
ed.agb <- readRDS("/Users/paleolab/Downloads/ED2.AGB.rds")
lpj.agb <- readRDS("/Users/paleolab/Downloads/LPJ-GUESS.AGB.rds")
lpj.wsl.agb <- readRDS("/Users/paleolab/Downloads/LPJ-WSL.v1.TotLivBiom.rds")
triffid.agb <- readRDS("/Users/paleolab/Downloads/TRIFFID.TotLivBio_PFT.rds")
load('/Users/paleolab/refab/second.deriv.refab.Rdata')
sd.refab <- second.deriv.unlist
load("/Users/paleolab/Downloads/PalEON_siteInfo_all.RData")
paleon <- read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/EcosystemStability/master/data/paleon_site_info.csv'))
linkages.agb <- readRDS("/Users/paleolab/Downloads/PalEON_regional_LINKAGES.AGB.rds")

linkages.agb.pft <- readRDS("/Users/paleolab/Downloads/PalEON_regional_LINKAGES.AGB.pft.rds")

spp <- c('fir','red maple','sugar maple','yellow birch','hickory',
         'chestnut','beech','white ash','tamarack','aspen','spruce',
         'white pine','white oak','hemlock','elm')


#####
##### REFAB OUTPUT #####
#####
load('~/ReFAB/allPredData.Rda')
load('~/babySTEPPS/biomass.CI13.Rdata')
names(biomassCI) <- unique(x.meta$site.id)[1:182]

refab.lat.lon <- matrix(NA,nrow=length(second.deriv.unlist),ncol=2)

for(i in 1:length(second.deriv.unlist)){
    refab.lat.lon[i,1:2] <- as.numeric(x.meta[x.meta[,1]==names(second.deriv.unlist)[i],c('long','lat')][1,])
}

#####
##### SETTLEMENT BIOMASS ####
#####

#Settlement biomass data product
biomass_dat_est <- read.csv(paste0("~/Downloads/","biomass_prediction_v0.9-10_bam.csv"))

albers = data.frame(biomass_dat_est$x,biomass_dat_est$y)
colnames(albers) = c('x', 'y')

coordinates(albers) <- ~ x + y
proj4string(albers) <-  CRS('+init=epsg:3175')

latlon <- spTransform(albers,CRS('+proj=longlat +ellps=WGS84'))
latlon <- as.matrix(data.frame(latlon))

centers_biomass = latlon
idx_cores = vector(length=nrow(paleon))

for(i in 1:nrow(paleon)){   
  core_site = paleon[i,c('lon','lat')]
  d = rdist(matrix(core_site, ncol=2), as.matrix(centers_biomass))
  if(min(d)<.5){
    idx_cores[i] = which.min(d) 
  }else{
    idx_cores[i] = NA 
  }
}

umw.keep <- which(!is.na(idx_cores))

map('state', xlim=range(paleon$lon[umw.keep])+c(-2, 2), ylim=range(paleon$lat[umw.keep])+c(-1, 1))
points(paleon$lon[umw.keep], paleon$lat[umw.keep], col='black', pch=19)
points(centers_biomass[idx_cores,1], centers_biomass[idx_cores,2], col='blue', pch=8)

tot_biomass <- ud.convert(biomass_dat_est$Total,'Mg/ha','kg/m^2')/2
dat_use <- tot_biomass[idx_cores]

map.colors <- function(x,TITLE){
  breaks <-  c(seq(0,12,1),9999)
  colors <- rev(terrain.colors(length(breaks)-1))
  data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
  points(paleon$lon, paleon$lat, pch=21, cex=1,bg=colors[data_binned],col='lightgray')
  points(paleon$lon[is.na(data_binned)], paleon$lat[is.na(data_binned)], pch=21, cex=1,bg='grey',col='darkgray')
  legend('bottomright',legend=c(signif(breaks[1:(length(breaks)-1)],digits=2),'NA'),pch=rep(21,length(breaks)),pt.bg =c(colors,'grey'))
  title(TITLE)
}

pdf('regional.biomass.map.output.pdf')
par(mfrow=c(3,2),mar=rep(0,4))
map.colors(x = dat_use,TITLE = c('DATA'))
map.colors(x = colMeans(ed.agb[12000:12012,]),TITLE = 'ED')
map.colors(x = lpj.agb[1000,,13],TITLE = 'LPJ-GUESS')
map.colors(x = lpj.wsl.agb[1000,],TITLE = 'LPJ-WSL')
map.colors(x = rowSums(triffid.agb[12012,,1:5]),TITLE = 'JULES')
map.colors(x = linkages.agb[1000,],TITLE = 'LINKAGES')
dev.off()

map.diff.colors <- function(x,TITLE){
  x1 <- dat_use[umw.keep]
  
  breaks <-  c(-999,-5,-1,1,5,10,999)
  colors <- colorRampPalette(c("blue",'yellow','red'))(length(breaks)-1)
  data_binned <-  cut(x - x1, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=range(paleon$lon[umw.keep])+c(-2, 2), ylim=range(paleon$lat[umw.keep])+c(-1, 1))
  points(paleon$lon[umw.keep], paleon$lat[umw.keep], pch=21, cex=2,bg=colors[data_binned])
  points(paleon$lon[umw.keep][is.na(data_binned)], paleon$lat[umw.keep][is.na(data_binned)], pch=21, cex=2,bg='grey')
  legend('topright',cex=.8,title = '= Model - Data (kgC/m^2)',legend=c('NA',signif(breaks[1:(length(breaks)-1)],digits=2)),pch=rep(19,length(breaks)),col=c('grey',colors))
  title(TITLE)
}

pdf('regional.diff.maps.biomass.pdf')
map.diff.colors(x = colMeans(ed.agb[12000:12012,umw.keep]),TITLE = 'ED')
map.diff.colors(x = lpj.agb[1000,umw.keep,13],TITLE = 'LPJ-GUESS')
map.diff.colors(x = lpj.wsl.agb[1000,umw.keep],TITLE = 'LPJ-WSL')
map.diff.colors(x = rowSums(triffid.agb[12012,umw.keep,1:5]),TITLE = 'JULES')
map.diff.colors(x = linkages.agb[1000,umw.keep],TITLE = 'LINKAGES')
dev.off()

data_binned <-  cut(ed.agb[i,], c(breaks), include.lowest = FALSE, labels = FALSE)

map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])

triff.box <- rowSums(triffid.agb[12012,umw.keep,1:5])
triff.box[triff.box<0] <- NA

box.stop <- cbind(dat_use[umw.keep],colMeans(ed.agb[12000:12012,umw.keep]),
                  lpj.agb[1000,umw.keep,13],lpj.wsl.agb[1000,umw.keep],triff.box,
                  linkages.agb[1000,umw.keep])
colnames(box.stop) <- c('Data','ED','LPJ-GUESS','LPJ-WSL','JULES','LINKAGES')

pdf('Biomass.Boxplot.pdf')
boxplot(box.stop,ylab = 'Biomass (kg/m^2)')
dev.off()

par(mfrow=c(4,4))
for(i in seq(1,nrow(paleon),length.out=20)){
data_binned <-  cut(rowSums(lpj.agb[i,,1:13]), c(breaks), include.lowest = FALSE, labels = FALSE)

map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
	
}

par(mfrow=c(4,4))
for(i in seq(1,nrow(paleon),length.out=20)){
  data_binned <-  cut(rowSums(lpj.agb[i,,1:13]), c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
  points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
  
}

#####
##### LINKAGES #####
#####

linkages.agb <- ud.convert(linkages.agb,'Mg/ha','kg/m^2')
linkages.agb[linkages.agb<0]<- NA

time.box <-     seq(1,1000,100)
link.mean.mat <- matrix(NA,254,9)
for(i in 1:254){
  for(t in 1:9){
    link.mean.mat[i,t] <- mean(linkages.agb[time.box[t]:(time.box[t+1]-1),i])
  }
}

second.deriv.keep <- numeric(254)
for(i in 1:254){
  second.deriv.keep[i] <- calc.second.deriv(biomassCI = link.mean.mat[i,],h=1)
}
sd.linkages <- second.deriv.keep

breaks <-  c(min(second.deriv.keep),quantile(second.deriv.keep,c(.025,.25,.5,.75,.975),na.rm = TRUE),max(second.deriv.keep))
#seq(range(second.deriv.keep)[1],range(second.deriv.keep)[2],length.out=10)
colors <- rev(rainbow(length(breaks)-1))

data_binned <-  cut(second.deriv.keep, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('linkages.2nd.deriv.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('square sum of 2nd deriv map -- LINKAGES')
dev.off()

cor.keep <- cor(linkages.agb)
cor.keep[is.na(cor.keep)]<-0

#corrplot(cor.keep,order='FPC') #Don't do all of the sites!

cor.sums <- rowSums(cor.keep)

breaks <-  seq(range(cor.sums,na.rm=TRUE)[1],range(cor.sums,na.rm=TRUE)[2],length.out=10)
colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))

data_binned <-  cut(cor.sums, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('linkages.corr.map.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('correlation map -- linkages')
dev.off()

pdf('linkages.agb.pft.maps.pdf')
for(i in 1:15){
  x = linkages.agb.pft[1000,,i]
  breaks <-  c(seq(0,20,5),9999)
  colors <- rev(terrain.colors(length(breaks)-1))
  data_binned <-  cut(x, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
  points(paleon$lon, paleon$lat, pch=21, cex=1.5,bg=colors[data_binned])
  points(paleon$lon[is.na(data_binned)], paleon$lat[is.na(data_binned)], pch=21, cex=1.5,bg='grey')
  legend('bottomright',title='AGB (kg/m^2)',legend=c(signif(breaks[1:(length(breaks)-1)],digits=2),'NA'),pch=rep(21,length(breaks)),pt.bg =c(colors,'grey'))
  title(spp[i])
}
dev.off()
#####
##### LPJ-GUESS #####
#####
time.box <-     seq(1,1000,100)
lpj.mean.mat <- matrix(NA,254,9)
for(i in 1:254){
  for(t in 1:9){
    lpj.mean.mat[i,t] <- mean(lpj.agb[time.box[t]:(time.box[t+1]-1),i,13])
  }
}

second.deriv.keep <- numeric(254)
for(i in 1:254){
	second.deriv.keep[i] <- calc.second.deriv(biomassCI = lpj.mean.mat[i,],h=1)
}
sd.lpj.guess <- second.deriv.keep

breaks <-  c(min(second.deriv.keep),quantile(second.deriv.keep,c(.025,.25,.5,.75,.975),na.rm = TRUE),max(second.deriv.keep))
#seq(range(second.deriv.keep)[1],range(second.deriv.keep)[2],length.out=10)
colors <- rev(rainbow(length(breaks)-1))

data_binned <-  cut(second.deriv.keep, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('lpj.guess.2nd.deriv.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('square sum of 2nd deriv map -- LPJ')
dev.off()

cor.keep <- cor(lpj.agb[,,13])
cor.keep[is.na(cor.keep)]<-0

#corrplot(cor.keep,order='FPC') #Don't do all of the sites!

cor.sums <- rowSums(cor.keep)

breaks <-  seq(range(cor.sums,na.rm=TRUE)[1],range(cor.sums,na.rm=TRUE)[2],length.out=10)
colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))

data_binned <-  cut(cor.sums, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('lpj.guess.corr.map.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('correlation map -- LPJ')
dev.off()

#####
##### LPJ-WSL #####
#####

time.box <-     seq(1,1000,100)
lpj.wsl.mean.mat <- matrix(NA,254,9)
for(i in 1:254){
  for(t in 1:9){
    lpj.wsl.mean.mat[i,t] <- mean(lpj.wsl.agb[time.box[t]:(time.box[t+1]-1),i])
  }
}

second.deriv.keep <- numeric(254)
for(i in 1:254){
	second.deriv.keep[i] <- calc.second.deriv(biomassCI = lpj.wsl.mean.mat[i,],h=1)
}
sd.lpj.wsl <- second.deriv.keep

breaks <-  c(min(second.deriv.keep),quantile(second.deriv.keep,c(.025,.25,.5,.75,.975)),max(second.deriv.keep))
#seq(range(second.deriv.keep)[1],range(second.deriv.keep)[2],length.out=10)
colors <- rev(rainbow(length(breaks)-1))

data_binned <-  cut(second.deriv.keep, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('lpj.wsl.2nd.deriv.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('square sum of 2nd deriv map -- LPJ-WSL')
dev.off()

cor.keep <- cor(lpj.wsl.agb)
cor.keep[is.na(cor.keep)]<-0

#corrplot(cor.keep,order='FPC') #Don't do all of the sites!

cor.sums <- rowSums(cor.keep)

breaks <-  seq(range(cor.sums,na.rm=TRUE)[1],range(cor.sums,na.rm=TRUE)[2],length.out=10)
colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))

data_binned <-  cut(cor.sums, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('lpj.wsl.corr.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('correlation map -- LPJ-WSL')
dev.off()

#####
##### JULES TRIFFID #####
#####

triffid.agb <- triffid.agb[1:13000,,]

time.box <-     seq(1,13000,1000)
triffid.mean.mat <- matrix(NA,254,11)
for(i in 1:254){
  for(t in 1:11){
    triffid.mean.mat[i,t] <- mean(rowSums(triffid.agb[time.box[t]:(time.box[t+1]-1),i,1:5]))
  }
}

breaks <-  seq(0,42,2)
colors <- rev(terrain.colors(length(breaks)-1))

quartz()
par(mfrow=c(3,3))
for(i in seq(1,nrow(paleon),length.out=20)){
  data_binned <-  cut(rowSums(triffid.agb[i,,1:5]), c(breaks), include.lowest = FALSE, 
                      labels = FALSE)
  
  map('state', xlim=range(paleon$lon)+c(-2, 2),
      ylim=range(paleon$lat)+c(-1, 1))
  points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
  
}

plot(rowSums(triffid.agb[,1,1:5]),ylim=c(0,25),typ='l')
for(i in sample(size = 10,x = 2:254)){
   lines(rowSums(triffid.agb[,i,1:5]),typ='l')
}

second.deriv.keep <- numeric(254)
for(i in 1:254){
	second.deriv.keep[i] <- calc.second.deriv(biomassCI =  triffid.mean.mat[i,] ,h=1)
}
sd.triffid <- second.deriv.keep

breaks <-  seq(range(second.deriv.keep)[1],range(second.deriv.keep)[2],length.out=7)
colors <- rev(rainbow(length(breaks)-1))

data_binned <-  cut(second.deriv.keep, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('jules.2nd.deriv.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('square sum of 2nd deriv map -- JULES TRIFFID')
dev.off()

cor.keep <- cor(triffid.agb[,,1])
cor.keep[is.na(cor.keep)]<-0

#corrplot(cor.keep,order='FPC') #Don't do all of the sites!

cor.sums <- rowSums(cor.keep)

breaks <-  seq(range(cor.sums,na.rm=TRUE)[1],range(cor.sums,na.rm=TRUE)[2],length.out=10)
colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))

data_binned <-  cut(cor.sums, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('jules.corr.map.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('correlation map -- JULES TRIFFID')
dev.off()

#####
##### ED #####
#####

time.box <-     seq(1,12012,1200)
ed.mean.mat <- matrix(NA,254,10)
for(i in 1:254){
  for(t in 1:10){
    ed.mean.mat[i,t] <- mean(ed.agb[time.box[t]:(time.box[t+1]-1),i])
  }
}


second.deriv.keep <- numeric(254)
for(i in 1:254){
	second.deriv.keep[i] <- calc.second.deriv(biomassCI = ed.mean.mat[i,],h=1)
}

sd.ed <- second.deriv.keep

breaks <-  c(min(second.deriv.keep,na.rm=TRUE),quantile(second.deriv.keep,c(.025,.25,.5,.75,.975),na.rm=TRUE),max(second.deriv.keep,na.rm=TRUE))
#seq(range(second.deriv.keep)[1],range(second.deriv.keep)[2],length.out=10)colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))
colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))

data_binned <-  cut(second.deriv.keep, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('ed.2nd.deriv.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('square sum of 2nd deriv map -- ED')
dev.off()

cor.keep <- cor(ed.agb[,1:245])
cor.keep[is.na(cor.keep)]<-0

#corrplot(cor.keep,order='FPC') #Don't do all of the sites!

cor.sums <- rowSums(cor.keep)

breaks <-  seq(range(cor.sums,na.rm=TRUE)[1],range(cor.sums,na.rm=TRUE)[2],length.out=10)
colors <- rev(rainbow(length(breaks)-1,start=0,end=.7))

data_binned <-  cut(cor.sums, c(breaks), include.lowest = FALSE, labels = FALSE)

pdf('ed.corr.map.pdf')
map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1))
points(paleon$lon, paleon$lat, pch=19, cex=1,col=colors[data_binned])
legend('bottomright',legend=signif(breaks[1:(length(breaks)-1)],digits=2),pch=rep(19,length(breaks)-1),col=colors)
title('correlation map -- ED')
dev.off()


#####
##### HISTS #####
#####

# WHOLE DOMAIN
pdf('hist.stability.biomass.overall.pdf')
hist(sd.refab,freq = FALSE,xlim=c(0,.0000015),ylim=c(0,10000000),col='lightgray',main='STABILITY of Whole Domain',xlab='Stability Index')
lines(density(na.omit(sd.ed)),col='red',lwd=2)
lines(density(na.omit(sd.triffid)),col='blue',lwd=2)
lines(density(na.omit(sd.lpj.wsl)),col='green',lwd=2)
lines(density(na.omit(sd.lpj.guess)),col='black',lwd=2)
lines(density(na.omit(sd.linkages)),col='purple',lwd=2)
legend('topright',legend = c('ED','TRIFFID','LPJ-WSL','LPJ-GUESS','LINKAGES','ReFAB Biom.'),pch=rep(19,6),col=c('red','blue','green','black','purple','grey'))
dev.off()

pdf('hist.stability.biomass.umw.pdf')
umw <- which(paleon$umw=='y')
hist(sd.refab,freq = FALSE,xlim=c(0,.0000015),ylim=c(0,15000000),
     col='lightgray',main='STABILITY of UMW',xlab='Stability Index')
lines(density(sd.ed[umw],na.rm=TRUE),col='red',lwd=2)
lines(density(sd.triffid[umw],na.rm=TRUE),col='blue',lwd=2)
lines(density(sd.lpj.wsl[umw],na.rm=TRUE),col='green',lwd=2)
lines(density(sd.lpj.guess[umw],na.rm=TRUE),col='black',lwd=2)
lines(density(sd.linkages[umw],na.rm=TRUE),col='purple',lwd=2)
legend('topright',legend = c('ED','TRIFFID','LPJ-WSL','LPJ-GUESS','LINKAGES','ReFAB Biom.'),pch=rep(19,6),col=c('red','blue','green','black','purple','grey'))
dev.off()

second.deriv.keep.1yr.ed <- numeric(254)
second.deriv.keep.1yr.jules <- numeric(254)
second.deriv.keep.1yr.linkages <- numeric(254)
second.deriv.keep.1yr.lpj.guess <- numeric(254)
second.deriv.keep.1yr.lpj.wsl <- numeric(254)
for(i in 1:254){
  second.deriv.keep.1yr.ed[i] <- calc.second.deriv(biomassCI = ed.agb[seq(1,12012,12),i],h=1)
  second.deriv.keep.1yr.jules[i] <- calc.second.deriv(biomassCI = rowSums(triffid.agb[seq(1,12012,12),i,1:5]),h=1)
  second.deriv.keep.1yr.linkages[i] <- calc.second.deriv(biomassCI = linkages.agb[1:1000,i],h=1)
  second.deriv.keep.1yr.lpj.guess[i] <- calc.second.deriv(biomassCI = lpj.agb[1:1000,i,13],h=1)
  second.deriv.keep.1yr.lpj.wsl[i] <- calc.second.deriv(biomassCI = lpj.wsl.agb[1:1000,i],h=1)
}

pdf('annual.stability.hists.pdf')
hist(1,freq = FALSE,col='white',main='Annual Stability of Whole Domain',xlab='Stability Index',
     xlim=c(0,.000005),ylim=c(0,5000000))
lines(density(na.omit(second.deriv.keep.1yr.ed)),col='red',lwd=3)
lines(density(na.omit(second.deriv.keep.1yr.jules)),col='blue',lwd=3)
lines(density(na.omit(second.deriv.keep.1yr.lpj.wsl)),col='green',lwd=3)
lines(density(na.omit(second.deriv.keep.1yr.lpj.guess)),col='black',lwd=3)
lines(density(na.omit(second.deriv.keep.1yr.linkages)),col='purple',lwd=3)
legend('topright',legend = c('ED','TRIFFID','LPJ-WSL','LPJ-GUESS','LINKAGES'),pch=rep(19,5),col=c('red','blue','green','black','purple'))

hist(1,freq = FALSE,col='white',main='Annual Stability of UMW',xlab='Stability Index',
     xlim=c(0,.000005),ylim=c(0,5000000))
lines(density(na.omit(second.deriv.keep.1yr.ed[umw])),col='red',lwd=3)
lines(density(na.omit(second.deriv.keep.1yr.jules[umw])),col='blue',lwd=3)
lines(density(na.omit(second.deriv.keep.1yr.lpj.wsl[umw])),col='green',lwd=3)
lines(density(na.omit(second.deriv.keep.1yr.lpj.guess[umw])),col='black',lwd=3)
lines(density(na.omit(second.deriv.keep.1yr.linkages[umw])),col='purple',lwd=3)
legend('topright',legend = c('ED','TRIFFID','LPJ-WSL','LPJ-GUESS','LINKAGES'),pch=rep(19,5),col=c('red','blue','green','black','purple'))
dev.off()
all.sd.1yr <- c(second.deriv.keep.1yr.ed, second.deriv.keep.1yr.jules, second.deriv.keep.1yr.lpj.wsl,
                second.deriv.keep.1yr.lpj.guess,second.deriv.keep.1yr.linkages)
map.colors.sd.1yr <- function(x,all.sd.1yr,TITLE){
  breaks <-  c(quantile(all.sd.1yr,c(.025,.25,.5,.75,.975),na.rm = TRUE),999)
  colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
  data_binned <-  cut(x, c(breaks), include.lowest = TRUE, labels = FALSE)
  
  map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1), col='gray')
  points(paleon$lon, paleon$lat, pch=21, cex=1,bg=colors[data_binned],col='gray')
  points(paleon$lon[is.na(data_binned)], paleon$lat[is.na(data_binned)], pch=21, cex=1.25,bg='grey')
  legend('bottomright',title='Quantile',legend=c(signif(c(.025,.25,.5,.75,.975),digits=3),'NA'),pch=rep(21,length(breaks)),pt.bg =c(colors,'gray'),cex=.5)
  title(TITLE)
}

pdf('stability.maps.together.annual.pdf')
par(mfrow=c(3,2),mar=rep(0,4))
map.colors.sd.1yr(x = second.deriv.keep.1yr.ed, all.sd.1yr, 'STABILITY ED')
map.colors.sd.1yr(x = second.deriv.keep.1yr.lpj.wsl, all.sd.1yr,  'STABILITY LPJ WSL')
map.colors.sd.1yr(x = second.deriv.keep.1yr.lpj.guess, all.sd.1yr,  'STABILITY LPJ GUESS')
map.colors.sd.1yr(x = second.deriv.keep.1yr.jules, all.sd.1yr, 'STABILITY JULES')
map.colors.sd.1yr(x = second.deriv.keep.1yr.linkages, all.sd.1yr, 'STABILITY LINKAGES')
dev.off()



all.sd <- c(sd.refab,sd.ed,sd.triffid,sd.lpj.guess,sd.lpj.wsl,sd.linkages)

map.colors.sd <- function(x,TITLE){
  breaks <-  c(quantile(all.sd,c(.025,.25,.5,.75,.975),na.rm = TRUE),999)
  colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
  data_binned <-  cut(x, c(breaks), include.lowest = TRUE, labels = FALSE)
  
  map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1), col='gray')
  points(paleon$lon, paleon$lat, pch=21, cex=1,bg=colors[data_binned],col='gray')
  points(paleon$lon[umw.keep][is.na(data_binned)], paleon$lat[umw.keep][is.na(data_binned)], pch=21, cex=1.25,bg='grey')
  legend('bottomright',title='Quantile',legend=c(signif(c(.025,.25,.5,.75,.975),digits=3),'NA'),pch=rep(21,length(breaks)),pt.bg =c(colors,'gray'),cex=.5)
  title(TITLE)
}


x = second.deriv.unlist
breaks <-  c(quantile(all.sd,c(.025,.25,.5,.75,.975),na.rm = TRUE),999)
data_binned <-  cut(x, c(breaks), include.lowest = TRUE, labels = FALSE)

pdf('stability.maps.together.pdf')
par(mfrow=c(3,2))
map('state', xlim=range(paleon$lon[umw.keep])+c(-2, 2), ylim=range(paleon$lat[umw.keep])+c(-1, 1), col='gray')
points(refab.lat.lon[,1], refab.lat.lon[,2], pch=21, cex=1,bg=colors[data_binned],col='gray')
legend('topright',title='Quantile',legend=signif(c(.025,.25,.5,.75,.975),digits=3),pch=rep(19,length(breaks)-1),col=colors,cex=.5)
title('ReFAB Stability')

map.colors.sd(x = sd.ed, 'STABILITY ED')
map.colors.sd(x = sd.lpj.wsl, 'STABILITY LPJ WSL')
map.colors.sd(x = sd.lpj.guess, 'STABILITY LPJ GUESS')
map.colors.sd(x = sd.triffid, 'STABILITY JULES')
map.colors.sd(x = sd.linkages, 'STABILITY LINKAGES')
dev.off()

#####
##### Principle Component ######
#####

library(vegan)

na.ed.agb <- ed.agb
na.ed.agb[is.na(na.ed.agb)]<-0
pZ <- princomp(na.ed.agb)

x <- pZ$scores[1,]

breaks <-  c(quantile(x,c(.025,.25,.5,.75,.975),na.rm = TRUE),999)
colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
data_binned <-  cut(x, c(breaks), include.lowest = TRUE, labels = FALSE)

map('state', xlim=range(paleon$lon)+c(-2, 2), ylim=range(paleon$lat)+c(-1, 1), col='gray')
points(paleon$lon, paleon$lat, pch=15, cex=1.5,col=colors[data_binned])
legend('topright',title='Quantile',legend=signif(c(.025,.25,.5,.75,.975),digits=3),pch=rep(19,length(breaks)-1),col=colors,cex=.5)




