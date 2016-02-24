# -------------------------------------------------------
# Scripts to read in, combine, & map output from the PalEON Regional Runs
# Created: Christy Rollinson, crollinson@gmail.com, Feb 2016
# -------------------------------------------------------


# -----------------------------
# libraries we're going to be using
# -----------------------------
library(ncdf4)
library(raster)
library(car)
library(abind)
library(ggplot2)
# -----------------------------

# -----------------------------
# Setting paths & directories
# -----------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP2_Region")
path.figs <- "phase2_model_output/Exploratory_Figs"
if(!dir.exists(path.figs)) dir.create(path.figs)

path.guess <- "phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1"
path.jules <- "phase2_model_output/jules"
path.wsl   <- "phase2_model_output/LPJ-WSL/LPJ-WSL.v1"

sec2yr <- 60*60*24*365

mo.yr <- vector()
for(y in 850:2010){
  mo.yr <- c(mo.yr, rep(y, 12))
}
length(mo.yr)
mo.yr[1:37]

# -----------------------------


# -----------------------------
# 1. LPJ-GUESS
# -----------------------------
# Looking at a few annual things first
guess.ann <- nc_open(file.path(path.guess, "LPJ-GUESS_annual.nc"))
summary(guess.ann$var)

guess.out <- list()
guess.out$lat <- ncvar_get(guess.ann, "lat")
guess.out$lon <- ncvar_get(guess.ann, "lon")
guess.out$Year <- ncvar_get(guess.ann, "time") + 850
guess.out$PFT <- ncvar_get(guess.ann, "PFT")
for(v in names(guess.ann$var)){
  guess.out[[v]] <- ncvar_get(guess.ann, v)  
}

dim(guess.out$Fcomp)
# -----------------------------

# -----------------------------
# 2. LPJ-WSL
# -----------------------------
# Getting a list of the files we have
files.wsl <- dir(file.path(path.wsl), ".nc")

# Creating a list to store things in
wsl.out <- list()

# Adding some of the global variables
wsl.temp <- nc_open(file.path(path.wsl, files.wsl[1]))
wsl.out$lat <- ncvar_get(wsl.temp, "lat")
wsl.out$lon <- ncvar_get(wsl.temp, "lon")

wsl.out$Month <- ncvar_get(wsl.temp, "Month")
nc_close(wsl.temp)

for(i in 1:length(files.wsl)){
  print(paste0(" ====== ", "Processing File: ", files.wsl[i], " ====== "))
  wsl.temp <- nc_open(file.path(path.wsl, files.wsl[i]))
  if(i == 1){
    wsl.out$Year <- ncvar_get(wsl.temp, "Year")
    for(v in names(wsl.temp$var)){
      wsl.out[[v]] <- ncvar_get(wsl.temp, v)  
    }
  } else {
    wsl.out$Year <- c(wsl.out$Year, ncvar_get(wsl.temp, "Year"))
    for(v in names(wsl.temp$var)){
      wsl.out[[v]] <- abind(wsl.out[[v]], ncvar_get(wsl.temp, v), along=3)
    }    
  }
  nc_close(wsl.temp)
}

summary(wsl.out)
dim(wsl.out$Fcomp)
dim(wsl.out$NPP)
# -----------------------------

# -----------------------------
# 3. Jules
# -----------------------------


# -----------------------------

# -----------------------------
# Summarizing Fcomp
# -----------------------------
# -----------
# LPJ-GUESS
# -----------
for(y in 1:dim(guess.out$Fcomp)[3]){
  print(paste0(" ---- Lat: ", y, " ---- "))
  dat.evg <- stack(data.frame(apply(guess.out$Fcomp[c(1,2,7,8,9),,y,], c(2,3), FUN=sum)))
  names(dat.evg) <- c("Fcomp", "Year")
  dat.evg$Year <- as.numeric(substr(dat.evg$Year,2,nchar(paste(dat.evg$Year)))) + 849
  dat.evg$lat  <- guess.out$lat[y]
  dat.evg$lon  <- guess.out$lon
  dat.evg$PFT  <- as.factor("Evergreen")      
  
  dat.decid <- stack(data.frame(apply(guess.out$Fcomp[c(3:6,10),,y,], c(2,3), FUN=sum)))
  names(dat.decid) <- c("Fcomp", "Year")
  dat.decid$Year <- as.numeric(substr(dat.decid$Year,2,nchar(paste(dat.decid$Year)))) + 849
  dat.decid$lat  <- guess.out$lat[y]
  dat.decid$lon  <- guess.out$lon
  dat.decid$PFT  <- as.factor("Deciduous")
  
  dat.grass <- stack(data.frame(apply(guess.out$Fcomp[c(11:12),,y,], c(2,3), FUN=sum)))
  names(dat.grass) <- c("Fcomp", "Year")
  dat.grass$Year <- as.numeric(substr(dat.grass$Year,2,nchar(paste(dat.grass$Year)))) + 849
  dat.grass$lat  <- guess.out$lat[y]
  dat.grass$lon  <- guess.out$lon
  dat.grass$PFT  <- "Grass"        
  
  if(y==1) guess.fcomp <- rbind(dat.evg, dat.decid, dat.grass) else guess.fcomp <- rbind(guess.fcomp, dat.evg, dat.decid, dat.grass)
}
summary(guess.fcomp)

# Graphing Fraction PFT
ggplot(data=guess.fcomp[guess.fcomp$Year %in% c(850, 1850, 2010),]) +
  facet_grid(Year~PFT) +
  geom_raster(aes(x=lon, y=lat, fill=Fcomp)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  ggtitle("LPJ-GUESS") +
  coord_equal(ratio=1)
# -----------

# -----------
# LPJ-WSL
# -----------
for(y in 1:dim(wsl.out$Fcomp)[2]){
  print(paste0(" ---- Lat: ", y, " ---- "))
  dat.evg <- stack(data.frame(apply(wsl.out$Fcomp[,y,,c(1,3,4,6)], c(1,2), FUN=sum)))
  names(dat.evg) <- c("Fcomp", "Year")
  dat.evg$Year <- as.numeric(substr(dat.evg$Year,2,nchar(paste(dat.evg$Year)))) + 849
  dat.evg$lat  <- wsl.out$lat[y]
  dat.evg$lon  <- wsl.out$lon
  dat.evg$PFT  <- as.factor("Evergreen")      
  
  dat.decid <- stack(data.frame(apply(wsl.out$Fcomp[,y,,c(2,5,7)], c(1,2), FUN=sum)))
  names(dat.decid) <- c("Fcomp", "Year")
  dat.decid$Year <- as.numeric(substr(dat.decid$Year,2,nchar(paste(dat.decid$Year)))) + 849
  dat.decid$lat  <- wsl.out$lat[y]
  dat.decid$lon  <- wsl.out$lon
  dat.decid$PFT  <- as.factor("Deciduous")
  
  dat.grass <- stack(data.frame(apply(wsl.out$Fcomp[,y,,c(8:9)], c(1,2), FUN=sum)))
  names(dat.grass) <- c("Fcomp", "Year")
  dat.grass$Year <- as.numeric(substr(dat.grass$Year,2,nchar(paste(dat.grass$Year)))) + 849
  dat.grass$lat  <- wsl.out$lat[y]
  dat.grass$lon  <- wsl.out$lon
  dat.grass$PFT  <- "Grass"        
  
  if(y==1) wsl.fcomp <- rbind(dat.evg, dat.decid, dat.grass) else wsl.fcomp <- rbind(wsl.fcomp, dat.evg, dat.decid, dat.grass)
}

# Graphing Fraction PFT
ggplot(data=wsl.fcomp[wsl.fcomp$Year %in% c(850, 1850, 2010),]) +
  facet_grid(Year~PFT) +
  geom_raster(aes(x=lon, y=lat, fill=Fcomp)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  ggtitle("LPJ-WSL") +
  coord_equal(ratio=1)
# -----------


# -----------
# Comparing Multiple Models
# -----------
wsl.fcomp$Model   <- as.factor("LPJ-WSL")
guess.fcomp$Model <- as.factor("LPJ-GUESS")

fcomp <- rbind(guess.fcomp, wsl.fcomp)
png(file.path(path.figs, "Fcomp_Evergreen_0850_1850_2010.png"), height=8.5, width=11, units="in", res=180)
ggplot(data=fcomp[fcomp$Year %in% c(850, 1850, 2010) & fcomp$PFT=="Evergreen",]) +
  facet_grid(Year~Model) +
  geom_raster(aes(x=lon, y=lat, fill=Fcomp)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_fill_gradientn(limits=c(0,1), colours=c("black", "green3")) +
  ggtitle("Evergreen") +
  coord_equal(ratio=1)
dev.off()

png(file.path(path.figs, "Fcomp_Deciduous_0850_1850_2010.png"), height=8.5, width=11, units="in", res=180)
ggplot(data=fcomp[fcomp$Year %in% c(850, 1850, 2010) & fcomp$PFT=="Deciduous",]) +
  facet_grid(Year~Model) +
  geom_raster(aes(x=lon, y=lat, fill=Fcomp)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_fill_gradientn(limits=c(0,1), colours=c("black", "dodgerblue2")) +
  ggtitle("Deciduous") +
  coord_equal(ratio=1)
dev.off()

png(file.path(path.figs, "Fcomp_Grass_0850_1850_2010.png"), height=8.5, width=11, units="in", res=180)
ggplot(data=fcomp[fcomp$Year %in% c(850, 1850, 2010) & fcomp$PFT=="Grass",]) +
  facet_grid(Year~Model) +
  geom_raster(aes(x=lon, y=lat, fill=Fcomp)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  scale_fill_gradientn(limits=c(0,1), colours=c("black", "goldenrod2")) +
  ggtitle("Grass") +
  coord_equal(ratio=1)
dev.off()
# -----------

# -----------------------------


# -----------------------------
# Summarizing NPP
# -----------------------------
# -----------
# LPJ-GUESS
# -----------
dim(guess.out$NPP)

for(y in 1:dim(guess.out$NPP)[3]){
  print(paste0(" ---- Lat: ", y, " ---- "))
  dat.temp <- stack(data.frame(guess.out$NPP[13,,y,]))
  names(dat.temp) <- c("NPP", "Year")
  dat.temp$Year <- as.numeric(substr(dat.temp$Year,2,nchar(paste(dat.temp$Year)))) + 849
  dat.temp$lat  <- guess.out$lat[y]
  dat.temp$lon  <- guess.out$lon

  if(y==1) guess.npp <- dat.temp else guess.npp <- rbind(guess.npp, dat.temp)
}
guess.npp$NPP.yr <- guess.npp$NPP*sec2yr
summary(guess.npp)
mean(guess.npp$NPP.yr, na.rm=T)
max(guess.npp$NPP.yr, na.rm=T)

mean(guess.npp$NPP, na.rm=T)

# Graphing
ggplot(data=guess.npp[guess.npp$Year %in% c(850, 1850, 2010),]) +
  facet_grid(Year~.) +
  geom_raster(aes(x=lon, y=lat, fill=NPP.yr)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  ggtitle("LPJ-GUESS") +
  coord_equal(ratio=1)
# -----------

# -----------
# LPJ-WSL
# -----------
dim(wsl.out$NPP)

for(y in 1:dim(wsl.out$NPP)[2]){
  print(paste0(" ---- Lat: ", y, " ---- "))
  dat.temp <- stack(data.frame(wsl.out$NPP[,y,]))
  names(dat.temp) <- c("NPP", "Year")
  dat.temp$Year <- mo.yr
  dat.temp$lat  <- wsl.out$lat[y]
  dat.temp$lon  <- wsl.out$lon
  dat.temp$NPP.yr <- dat.temp$NPP*sec2yr
  
  dat.temp <- aggregate(dat.temp[,c("NPP", "NPP.yr")], by=dat.temp[,c("lat", "lon", "Year")], FUN=mean, na.rm=T)

  if(y==1) wsl.npp <- dat.temp else wsl.npp <- rbind(wsl.npp, dat.temp)
}
wsl.npp$NPP.yr <- wsl.npp$NPP*sec2yr
summary(wsl.npp)
dim(wsl.npp)

mean(wsl.npp$NPP.yr, na.rm=T)
max(wsl.npp$NPP.yr, na.rm=T)

mean(wsl.npp$NPP, na.rm=T)

# Graphing
ggplot(data=wsl.npp[wsl.npp$Year %in% c(850, 1850, 2010),]) +
  facet_grid(Year~.) +
  geom_raster(aes(x=lon, y=lat, fill=NPP.yr)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  ggtitle("LPJ-wsl") +
  coord_equal(ratio=1)
# -----------

# -----------
# Comparing Models
# -----------
wsl.npp$Model   <- as.factor("LPJ-WSL")
guess.npp$Model <- as.factor("LPJ-GUESS")

guess.npp <- guess.npp[,c("lat", "lon", "Year", "NPP", "NPP.yr", "Model")]

npp <- rbind(guess.npp, wsl.npp)

png(file.path(path.figs, "NPP_0850_1850_2010.png"), height=8.5, width=11, units="in", res=180)
ggplot(data=npp[npp$Year %in% c(850, 1850, 2010),]) +
  facet_grid(Year~Model) +
  geom_raster(aes(x=lon, y=lat, fill=NPP.yr)) +
  scale_fill_gradientn(colours=c("black", "green3")) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  coord_equal(ratio=1)
dev.off()

# -----------

# -----------------------------


# -----------------------------
# Summarizing NEE
# -----------------------------
# -----------
# LPJ-GUESS
# -----------
dim(guess.out$NEE)

for(y in 1:dim(guess.out$NEE)[2]){
  print(paste0(" ---- Lat: ", y, " ---- "))
  dat.temp <- stack(data.frame(guess.out$NEE[,y,]))
  names(dat.temp) <- c("NEE", "Year")
  dat.temp$Year <- as.numeric(substr(dat.temp$Year,2,nchar(paste(dat.temp$Year)))) + 849
  dat.temp$lat  <- guess.out$lat[y]
  dat.temp$lon  <- guess.out$lon
  
  if(y==1) guess.NEE <- dat.temp else guess.NEE <- rbind(guess.NEE, dat.temp)
}
guess.NEE$NEE.yr <- guess.NEE$NEE*sec2yr
summary(guess.NEE)
mean(guess.NEE$NEE.yr, na.rm=T)
max(guess.NEE$NEE.yr, na.rm=T)

mean(guess.NEE$NEE, na.rm=T)

summary(guess.NEE[guess.NEE$Year==2009,])
# Graphing
ggplot(data=guess.NEE[guess.NEE$Year %in% c(850, 1850, 2000),]) +
  facet_grid(Year~.) +
  geom_raster(aes(x=lon, y=lat, fill=NEE.yr)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  ggtitle("LPJ-GUESS") +
  coord_equal(ratio=1)
# -----------

# -----------
# LPJ-WSL
# -----------
dim(wsl.out$NEE)

for(y in 1:dim(wsl.out$NEE)[2]){
  print(paste0(" ---- Lat: ", y, " ---- "))
  dat.temp <- stack(data.frame(wsl.out$NEE[,y,]))
  names(dat.temp) <- c("NEE", "Year")
  dat.temp$Year <- mo.yr
  dat.temp$lat  <- wsl.out$lat[y]
  dat.temp$lon  <- wsl.out$lon
  dat.temp$NEE.yr <- dat.temp$NEE*sec2yr
  
  dat.temp <- aggregate(dat.temp[,c("NEE", "NEE.yr")], by=dat.temp[,c("lat", "lon", "Year")], FUN=mean, na.rm=T)
  
  if(y==1) wsl.NEE <- dat.temp else wsl.NEE <- rbind(wsl.NEE, dat.temp)
}
wsl.NEE$NEE.yr <- wsl.NEE$NEE*sec2yr
summary(wsl.NEE)
dim(wsl.NEE)

mean(wsl.NEE$NEE.yr, na.rm=T)
max(wsl.NEE$NEE.yr, na.rm=T)

mean(wsl.NEE$NEE, na.rm=T)

# Doing Biome Classifications
ggplot(data=wsl.NEE[wsl.NEE$Year %in% c(850, 1850, 2000),]) +
  facet_grid(Year~.) +
  geom_raster(aes(x=lon, y=lat, fill=NEE.yr)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  ggtitle("LPJ-wsl") +
  coord_equal(ratio=1)
# -----------

# -----------
# Comparing Models
# -----------
wsl.NEE$Model   <- as.factor("LPJ-WSL")
guess.NEE$Model <- as.factor("LPJ-GUESS")

guess.NEE <- guess.NEE[,c("lat", "lon", "Year", "NEE", "NEE.yr", "Model")]

NEE <- rbind(guess.NEE, wsl.NEE)

png(file.path(path.figs, "NEE_0850_1850_2000.png"), height=8.5, width=11, units="in", res=180)
ggplot(data=NEE[NEE$Year %in% c(850, 1850, 2000),]) +
  facet_grid(Year~Model) +
  geom_raster(aes(x=lon, y=lat, fill=NEE.yr)) +
  scale_fill_gradient2(low="green3", mid="gray10", high="red3") +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  coord_equal(ratio=1)
dev.off()
# -----------
# -----------------------------


# -----------------------------
# Summarizing TotLivBiom
# -----------------------------
# -----------
# LPJ-GUESS
# -----------
dim(guess.out$TotLivBiom)

for(y in 1:dim(guess.out$TotLivBiom)[3]){
  print(paste0(" ---- Lat: ", y, " ---- "))
  dat.temp <- stack(data.frame(guess.out$TotLivBiom[13,,y,]))
  names(dat.temp) <- c("TotLivBiom", "Year")
  dat.temp$Year <- as.numeric(substr(dat.temp$Year,2,nchar(paste(dat.temp$Year)))) + 849
  dat.temp$lat  <- guess.out$lat[y]
  dat.temp$lon  <- guess.out$lon
  
  if(y==1) guess.TotLivBiom <- dat.temp else guess.TotLivBiom <- rbind(guess.TotLivBiom, dat.temp)
}
summary(guess.TotLivBiom)
mean(guess.TotLivBiom$TotLivBiom.yr, na.rm=T)
max(guess.TotLivBiom$TotLivBiom.yr, na.rm=T)

mean(guess.TotLivBiom$TotLivBiom, na.rm=T)

summary(guess.TotLivBiom[guess.TotLivBiom$Year==2009,])
# Graphing
ggplot(data=guess.TotLivBiom[guess.TotLivBiom$Year %in% c(850, 1850, 2000),]) +
  facet_grid(Year~.) +
  geom_raster(aes(x=lon, y=lat, fill=TotLivBiom)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  ggtitle("LPJ-GUESS") +
  coord_equal(ratio=1)
# -----------

# -----------
# LPJ-WSL
# -----------
dim(wsl.out$TotLivBiom)

for(y in 1:dim(wsl.out$TotLivBiom)[2]){
  print(paste0(" ---- Lat: ", y, " ---- "))
  dat.temp <- stack(data.frame(wsl.out$TotLivBiom[,y,]))
  names(dat.temp) <- c("TotLivBiom", "Year")
  dat.temp$Year <- wsl.out$Year
  dat.temp$lat  <- wsl.out$lat[y]
  dat.temp$lon  <- wsl.out$lon
  
  if(y==1) wsl.TotLivBiom <- dat.temp else wsl.TotLivBiom <- rbind(wsl.TotLivBiom, dat.temp)
}
summary(wsl.TotLivBiom)
dim(wsl.TotLivBiom)

mean(wsl.TotLivBiom$TotLivBiom.yr, na.rm=T)
max(wsl.TotLivBiom$TotLivBiom.yr, na.rm=T)

mean(wsl.TotLivBiom$TotLivBiom, na.rm=T)

# Doing Biome Classifications
ggplot(data=wsl.TotLivBiom[wsl.TotLivBiom$Year %in% c(850, 1850, 2000),]) +
  facet_grid(Year~.) +
  geom_raster(aes(x=lon, y=lat, fill=TotLivBiom)) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  ggtitle("LPJ-wsl") +
  coord_equal(ratio=1)
# -----------

# -----------
# Comparing Models
# -----------
wsl.TotLivBiom$Model   <- as.factor("LPJ-WSL")
guess.TotLivBiom$Model <- as.factor("LPJ-GUESS")

TotLivBiom <- rbind(guess.TotLivBiom, wsl.TotLivBiom)

png(file.path(path.figs, "TotLivBiom_0850_1850_2000.png"), height=8.5, width=11, units="in", res=180)
ggplot(data=TotLivBiom[TotLivBiom$Year %in% c(850, 1850, 2000),]) +
  facet_grid(Year~Model) +
  geom_raster(aes(x=lon, y=lat, fill=TotLivBiom)) +
  scale_fill_gradientn(colors=c("black", "green3")) +
  scale_y_continuous(name="Latitude", expand=c(0,0)) +
  scale_x_continuous(name="Longitude", expand=c(0,0)) +
  coord_equal(ratio=1)
dev.off()
# -----------
# -----------------------------
