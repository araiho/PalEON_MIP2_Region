# -------------------------------------------------------
# Scripts to do some quick QAQC of regional output by creating maps & time series
# Created: Christy Rollinson, crollinson@gmail.com, May 2016
# -------------------------------------------------------


# -----------------------------
# libraries we're going to be using
# -----------------------------
# Set the working directory
setwd("~/Dropbox/PalEON_CR/PalEON_MIP2_Region")

library(ncdf4)
library(raster)
library(car)
library(abind)
library(ggplot2)

# My extraction utility that can be downloaded here: 
#   https://github.com/PalEON-Project/MIP_Utils/blob/master/Phase2_region/extract_output_region.R 
source("/Users/crolli/Dropbox/PalEON_CR/MIP_Utils/Phase2_region/extract_output_region.R") # Generalized Script for extraction
source("/Users/crolli/Dropbox/PalEON_CR/MIP_Utils/Phase2_region/interpolate_model.R") # Script for interpolating ED

# Set some file paths
path.figs <- paste0("QAQC_Figs")
if(! dir.exists(path.figs)) dir.create(path.figs)

path.guess <- "phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1"
path.wsl   <- "phase2_model_output/LPJ-WSL/LPJ-WSL.v1"
path.ed    <- "phase2_model_output/ED2/ED2.v1.2016-05-25"
# -----------------------------


# -----------------------------
# Extract & graph NPP
# -----------------------------
{
npp.wsl <- list()
npp.wsl[["spinup"]]  <- extract.paleon.site(model="LPJ-WSL", 
                                            model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                            vars="NPP", yrmin= 850, yrmax= 869)[["NPP"]]
npp.wsl[["pre-ind"]] <- extract.paleon.site(model="LPJ-WSL", 
                                            model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                            vars="NPP", yrmin=1830, yrmax=1849)[["NPP"]]
npp.wsl[["modern"]]  <- extract.paleon.site(model="LPJ-WSL", 
                                            model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                            vars="NPP", yrmin=1991, yrmax=2010)[["NPP"]]

npp.guess <- list()
npp.guess[["spinup"]]  <- extract.paleon.site(model="LPJ-GUESS", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                              vars="NPP", yrmin= 850, yrmax= 869)[["NPP"]]
npp.guess[["pre-ind"]] <- extract.paleon.site(model="LPJ-GUESS", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                              vars="NPP", yrmin=1830, yrmax=1849)[["NPP"]]
npp.guess[["modern"]]  <- extract.paleon.site(model="LPJ-GUESS", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                              vars="NPP", yrmin=1991, yrmax=2010)[["NPP"]]


npp.ed <- list()
npp.ed[["spinup"]]  <- extract.paleon.site(model="ED2", 
                                           model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                           vars="NPP", yrmin= 850, yrmax= 869)[["NPP"]]
npp.ed[["pre-ind"]] <- extract.paleon.site(model="ED2", 
                                           model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                           vars="NPP", yrmin=1830, yrmax=1849)[["NPP"]]
npp.ed[["modern"]]  <- extract.paleon.site(model="ED2", 
                                           model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                           vars="NPP", yrmin=1991, yrmax=2010)[["NPP"]]

# Condensing the different layers down to the means to make graphing easier
for(i in 1:length(npp.wsl)){
  npp.wsl  [[i]] <- apply(npp.wsl  [[i]], c(1,2), mean, na.rm=T)
  npp.guess[[i]] <- apply(npp.guess[[i]], c(1,2), mean, na.rm=T)
  npp.ed   [[i]] <- apply(npp.ed   [[i]], c(1,2), mean, na.rm=T)
}

# Shaping things so that they're easy to graph in ggplot
for(i in 1:length(npp.wsl)){
  wsl.tmp <- stack(data.frame(npp.wsl[[i]]))
  names(wsl.tmp) <- c("NPP", "lon")
  wsl.tmp$lon <- as.numeric(paste0("-", substr(wsl.tmp$lon, 3,7)))
  wsl.tmp$lat <- as.numeric(dimnames(npp.wsl[[i]])[[1]])
  wsl.tmp$Model <- as.factor("LPJ-WSL")
  wsl.tmp$Time  <- as.factor(names(npp.wsl)[i])

  guess.tmp <- stack(data.frame(npp.guess[[i]]))
  names(guess.tmp) <- c("NPP", "lon")
  guess.tmp$lon <- as.numeric(paste0("-", substr(guess.tmp$lon, 3,7)))
  guess.tmp$lat <- sort(as.numeric(dimnames(npp.guess[[i]])[[1]]), decreasing=T)
  guess.tmp$Model <- as.factor("LPJ-GUESS")
  guess.tmp$Time  <- as.factor(names(npp.guess)[i])

  ed.tmp <- stack(data.frame(npp.ed[[i]]))
  names(ed.tmp) <- c("NPP", "lon")
  ed.tmp$lon <- as.numeric(paste0("-", substr(ed.tmp$lon, 3,7)))
  ed.tmp$lat <- as.numeric(dimnames(npp.ed[[i]])[[1]])
  ed.tmp$Model <- as.factor("ED2")
  ed.tmp$Time  <- as.factor(names(npp.ed)[i])
  
  if(i==1){
    npp.out <- rbind(wsl.tmp, guess.tmp, ed.tmp)
  } else {
    npp.out <- rbind(npp.out, rbind(wsl.tmp, guess.tmp, ed.tmp))
  }
}
summary(npp.out)

# interpolate the model
npp.out <- interp.mod(dat.all=npp.out, mod.missing="ED2", mod.mask="LPJ-WSL", var.interp="NPP")
npp.out$flag <- as.factor(npp.out$flag)
summary(npp.out)

npp.out$NPP.corr <- ifelse(npp.out$Model=="LPJ-GUESS", npp.out$NPP*1/(24*60*60), npp.out$NPP)

# The plot
png(file.path(path.figs, "NPP.png"), height=8.5, width=11, units="in", res=180)
ggplot(npp.out[complete.cases(npp.out),]) +
  facet_grid(Model ~ Time, drop=T) +
  geom_raster(data=npp.out[is.na(npp.out$flag), ], aes(x=lon, y=lat, fill=NPP.corr)) +
  geom_raster(data=npp.out[complete.cases(npp.out) & npp.out$flag=="interpolated", ], aes(x=lon, y=lat, fill=NPP.corr), alpha=0.8) +
  coord_equal(ratio=1) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ggtitle(paste0("Date Created: ", Sys.Date())) +
  theme(panel.background=element_blank())
dev.off()
}
# -----------------------------

# -----------------------------
# Extract & graph NEE
# -----------------------------
{
  nee.wsl <- list()
  nee.wsl[["spinup"]]  <- extract.paleon.site(model="LPJ-WSL", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                              vars="NEE", yrmin= 850, yrmax= 869)[["NEE"]]
  nee.wsl[["pre-ind"]] <- extract.paleon.site(model="LPJ-WSL", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                              vars="NEE", yrmin=1830, yrmax=1849)[["NEE"]]
  nee.wsl[["modern"]]  <- extract.paleon.site(model="LPJ-WSL", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                              vars="NEE", yrmin=1991, yrmax=2010)[["NEE"]]
  
  nee.guess <- list()
  nee.guess[["spinup"]]  <- extract.paleon.site(model="LPJ-GUESS", 
                                                model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                                vars="NEE", yrmin= 850, yrmax= 869)[["NEE"]]
  nee.guess[["pre-ind"]] <- extract.paleon.site(model="LPJ-GUESS", 
                                                model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                                vars="NEE", yrmin=1830, yrmax=1849)[["NEE"]]
  nee.guess[["modern"]]  <- extract.paleon.site(model="LPJ-GUESS", 
                                                model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                                vars="NEE", yrmin=1991, yrmax=2010)[["NEE"]]
  
  
  nee.ed <- list()
  nee.ed[["spinup"]]  <- extract.paleon.site(model="ED2", 
                                             model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                             vars="NEE", yrmin= 850, yrmax= 869)[["NEE"]]
  nee.ed[["pre-ind"]] <- extract.paleon.site(model="ED2", 
                                             model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                             vars="NEE", yrmin=1830, yrmax=1849)[["NEE"]]
  nee.ed[["modern"]]  <- extract.paleon.site(model="ED2", 
                                             model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                             vars="NEE", yrmin=1991, yrmax=2010)[["NEE"]]
  
  # Condensing the different layers down to the means to make graphing easier
  for(i in 1:length(nee.wsl)){
    nee.wsl  [[i]] <- apply(nee.wsl  [[i]], c(1,2), mean, na.rm=T)
    nee.guess[[i]] <- apply(nee.guess[[i]], c(1,2), mean, na.rm=T)
    nee.ed   [[i]] <- apply(nee.ed   [[i]], c(1,2), mean, na.rm=T)
  }
  
  # Shaping things so that they're easy to graph in ggplot
  for(i in 1:length(nee.wsl)){
    wsl.tmp <- stack(data.frame(nee.wsl[[i]]))
    names(wsl.tmp) <- c("NEE", "lon")
    wsl.tmp$lon <- as.numeric(paste0("-", substr(wsl.tmp$lon, 3,7)))
    wsl.tmp$lat <- as.numeric(dimnames(nee.wsl[[i]])[[1]])
    wsl.tmp$Model <- as.factor("LPJ-WSL")
    wsl.tmp$Time  <- as.factor(names(nee.wsl)[i])
    
    guess.tmp <- stack(data.frame(nee.guess[[i]]))
    names(guess.tmp) <- c("NEE", "lon")
    guess.tmp$lon <- as.numeric(paste0("-", substr(guess.tmp$lon, 3,7)))
    guess.tmp$lat <- sort(as.numeric(dimnames(nee.guess[[i]])[[1]]), decreasing=T)
    guess.tmp$Model <- as.factor("LPJ-GUESS")
    guess.tmp$Time  <- as.factor(names(nee.guess)[i])
    
    ed.tmp <- stack(data.frame(nee.ed[[i]]))
    names(ed.tmp) <- c("NEE", "lon")
    ed.tmp$lon <- as.numeric(paste0("-", substr(ed.tmp$lon, 3,7)))
    ed.tmp$lat <- as.numeric(dimnames(nee.ed[[i]])[[1]])
    ed.tmp$Model <- as.factor("ED2")
    ed.tmp$Time  <- as.factor(names(nee.ed)[i])
    
    if(i==1){
      nee.out <- rbind(wsl.tmp, guess.tmp, ed.tmp)
    } else {
      nee.out <- rbind(nee.out, rbind(wsl.tmp, guess.tmp, ed.tmp))
    }
  }
  summary(nee.out)
  
  # interpolate the model
  nee.out <- interp.mod(dat.all=nee.out, mod.missing="ED2", mod.mask="LPJ-WSL", var.interp="NEE")
  nee.out$flag <- as.factor(nee.out$flag)
  summary(nee.out)
  
  nee.out[,"NEE.corr"] <- ifelse(nee.out$Model=="LPJ-GUESS", nee.out$NEE*1/(24*60*60), nee.out$NEE)
  
  # The plot
  png(file.path(path.figs, "NEE.png"), height=8.5, width=11, units="in", res=180)
  ggplot(nee.out[complete.cases(nee.out),]) +
    facet_grid(Model ~ Time, drop=T) +
    geom_raster(data=nee.out[is.na(nee.out$flag), ], aes(x=lon, y=lat, fill=NEE.corr)) +
    geom_raster(data=nee.out[complete.cases(nee.out) & nee.out2$flag=="interpolated", ], aes(x=lon, y=lat, fill=NEE.corr), alpha=0.8) +
    coord_equal(ratio=1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_gradient2(low="green3", mid="gray10", high="red3") +
    ggtitle(paste0("Date Created: ", Sys.Date())) +
    theme(panel.background=element_blank())
  dev.off()
}
# -----------------------------

# -----------------------------
# Extract & graph Biomass
# -----------------------------
{
agb.wsl <- list()
agb.wsl[["spinup"]]  <- extract.paleon.site(model="LPJ-WSL", 
                                            model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                            vars="TotLivBiom", yrmin= 850, yrmax= 869)[["TotLivBiom"]]
agb.wsl[["pre-ind"]] <- extract.paleon.site(model="LPJ-WSL", 
                                            model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                            vars="TotLivBiom", yrmin=1830, yrmax=1849)[["TotLivBiom"]]
agb.wsl[["modern"]]  <- extract.paleon.site(model="LPJ-WSL", 
                                            model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                            vars="TotLivBiom", yrmin=1991, yrmax=2010)[["TotLivBiom"]]

agb.guess <- list()
agb.guess[["spinup"]]  <- extract.paleon.site(model="LPJ-GUESS", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                              vars="AGB", yrmin= 850, yrmax= 869)[["AGB"]][,,,13]
agb.guess[["pre-ind"]] <- extract.paleon.site(model="LPJ-GUESS", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                              vars="AGB", yrmin=1830, yrmax=1849)[["AGB"]][,,,13]
agb.guess[["modern"]]  <- extract.paleon.site(model="LPJ-GUESS", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                              vars="AGB", yrmin=1991, yrmax=2010)[["AGB"]][,,,13]


agb.ed <- list()
agb.ed[["spinup"]]  <- extract.paleon.site(model="ED2", 
                                           model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                           vars="AGB", yrmin= 850, yrmax= 869)[["AGB"]]
agb.ed[["pre-ind"]] <- extract.paleon.site(model="ED2", 
                                           model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                           vars="AGB", yrmin=1830, yrmax=1849)[["AGB"]]
agb.ed[["modern"]]  <- extract.paleon.site(model="ED2", 
                                           model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                           vars="AGB", yrmin=1991, yrmax=2010)[["AGB"]]

# Condensing the different layers down to the means to make graphing easier
for(i in 1:length(agb.wsl)){
  agb.wsl  [[i]] <- apply(agb.wsl  [[i]], c(1,2), mean, na.rm=T)
  agb.guess[[i]] <- apply(agb.guess[[i]], c(1,2), mean, na.rm=T)
  agb.ed   [[i]] <- apply(agb.ed   [[i]], c(1,2), mean, na.rm=T)
}

# Shaping things so that they're easy to graph in ggplot
for(i in 1:length(agb.wsl)){
  wsl.tmp <- stack(data.frame(agb.wsl[[i]]))
  names(wsl.tmp) <- c("AGB", "lon")
  wsl.tmp$lon <- as.numeric(paste0("-", substr(wsl.tmp$lon, 3,7)))
  wsl.tmp$lat <- as.numeric(dimnames(agb.wsl[[i]])[[1]])
  wsl.tmp$Model <- as.factor("LPJ-WSL")
  wsl.tmp$Time  <- as.factor(names(agb.wsl)[i])
  
  guess.tmp <- stack(data.frame(agb.guess[[i]]))
  names(guess.tmp) <- c("AGB", "lon")
  guess.tmp$lon <- as.numeric(paste0("-", substr(guess.tmp$lon, 3,7)))
  guess.tmp$lat <- as.numeric(dimnames(agb.guess[[i]])[[1]])
  guess.tmp$Model <- as.factor("LPJ-GUESS")
  guess.tmp$Time  <- as.factor(names(agb.guess)[i])
  
  ed.tmp <- stack(data.frame(agb.ed[[i]]))
  names(ed.tmp) <- c("AGB", "lon")
  ed.tmp$lon <- as.numeric(paste0("-", substr(ed.tmp$lon, 3,7)))
  ed.tmp$lat <- as.numeric(dimnames(agb.ed[[i]])[[1]])
  ed.tmp$Model <- as.factor("ED2")
  ed.tmp$Time  <- as.factor(names(agb.ed)[i])
  
  if(i==1){
    agb.out <- rbind(wsl.tmp, guess.tmp, ed.tmp)
  } else {
    agb.out <- rbind(agb.out, rbind(wsl.tmp, guess.tmp, ed.tmp))
  }
}

# interpoalte ED
agb.out <- interp.mod(dat.all=agb.out, mod.missing="ED2", mod.mask="LPJ-WSL", var.interp="AGB")
agb.out$flag <- as.factor(agb.out$flag)
summary(agb.out)

png(file.path(path.figs, "AGB.png"), height=8.5, width=11, units="in", res=180)
ggplot(agb.out[complete.cases(agb.out),]) +
  facet_grid(Model ~ Time, drop=T) +
  geom_raster(data=agb.out[is.na(agb.out$flag), ], aes(x=lon, y=lat, fill=AGB)) +
  geom_raster(data=agb.out[complete.cases(agb.out) & agb.out$flag=="interpolated", ], aes(x=lon, y=lat, fill=AGB), alpha=0.6) +
  coord_equal(ratio=1) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ggtitle(paste0("Date Created: ", Sys.Date())) +
  theme(panel.background=element_blank())
dev.off()
}
# -----------------------------

# -----------------------------
# Extract & graph Fcomp
# -----------------------------
{
  fcomp.wsl <- list()
  fcomp.wsl[["spinup"]]  <- extract.paleon.site(model="LPJ-WSL", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                              vars="Fcomp", yrmin= 850, yrmax= 869)[["Fcomp"]]
  fcomp.wsl[["pre-ind"]] <- extract.paleon.site(model="LPJ-WSL", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                              vars="Fcomp", yrmin=1830, yrmax=1849)[["Fcomp"]]
  fcomp.wsl[["modern"]]  <- extract.paleon.site(model="LPJ-WSL", 
                                              model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-WSL/LPJ-WSL.v1", 
                                              vars="Fcomp", yrmin=1991, yrmax=2010)[["Fcomp"]]
  
  fcomp.guess <- list()
  fcomp.guess[["spinup"]]  <- extract.paleon.site(model="LPJ-GUESS", 
                                                model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                                vars="Fcomp", yrmin= 850, yrmax= 869)[["Fcomp"]]
  fcomp.guess[["pre-ind"]] <- extract.paleon.site(model="LPJ-GUESS", 
                                                model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                                vars="Fcomp", yrmin=1830, yrmax=1849)[["Fcomp"]]
  fcomp.guess[["modern"]]  <- extract.paleon.site(model="LPJ-GUESS", 
                                                model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/LPJ-GUESS/LPJ-GUESS.v1.1", 
                                                vars="Fcomp", yrmin=1991, yrmax=2010)[["Fcomp"]]
  
  
  fcomp.ed <- list()
  fcomp.ed[["spinup"]]  <- extract.paleon.site(model="ED2", 
                                             model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                             vars="Fcomp", yrmin= 850, yrmax= 869)[["Fcomp"]]
  fcomp.ed[["pre-ind"]] <- extract.paleon.site(model="ED2", 
                                             model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                             vars="Fcomp", yrmin=1830, yrmax=1849)[["Fcomp"]]
  fcomp.ed[["modern"]]  <- extract.paleon.site(model="ED2", 
                                             model.dir="/Users/crolli/Dropbox/PalEON_CR/PalEON_MIP2_Region/phase2_model_output/ED2/ED2.v1.2016-05-25", 
                                             vars="Fcomp", yrmin=1991, yrmax=2010)[["Fcomp"]]
  
  # Condensing the different layers down to the means to make graphing easier
  for(i in 1:length(fcomp.wsl)){
    fcomp.wsl  [[i]] <- apply(fcomp.wsl  [[i]], c(1,2,4), mean, na.rm=T)
    fcomp.guess[[i]] <- apply(fcomp.guess[[i]], c(1,2,4), mean, na.rm=T)
    fcomp.ed   [[i]] <- apply(fcomp.ed   [[i]], c(1,2,4), mean, na.rm=T)
  }
  
  # Shaping things so that they're easy to graph in ggplot
  for(i in 1:length(fcomp.wsl)){
    wsl.tmp <- stack(data.frame(apply(fcomp.wsl[[i]][,,c(1:2,4:5,7)], c(1,2), sum, na.rm=F)))
    names(wsl.tmp) <- c("Decid", "lon")
    wsl.tmp$Evg   <- stack(data.frame(apply(fcomp.wsl[[i]][,,c(3,6)], c(1,2), sum, na.rm=F)))[,1]
    wsl.tmp$Grass <- stack(data.frame(apply(fcomp.wsl[[i]][,,c(8,9)], c(1,2), sum, na.rm=F)))[,1]
    wsl.tmp$lon <- as.numeric(paste0("-", substr(wsl.tmp$lon, 3,7)))
    wsl.tmp$lat <- as.numeric(dimnames(fcomp.wsl[[i]])[[1]])
    wsl.tmp$Model <- as.factor("LPJ-WSL")
    wsl.tmp$Time  <- as.factor(names(fcomp.wsl)[i])
    
    guess.tmp <- stack(data.frame(apply(fcomp.guess[[i]][,,c(3:6,10)], c(1,2), sum, na.rm=F)))
    names(guess.tmp) <- c("Decid", "lon")
    guess.tmp$Evg   <- stack(data.frame(apply(fcomp.guess[[i]][,,c(1:2,7:9)], c(1,2), sum, na.rm=F)))[,1]
    guess.tmp$Grass <- stack(data.frame(apply(fcomp.guess[[i]][,,c(11:12)], c(1,2), sum, na.rm=F)))[,1]
    guess.tmp$lon <- as.numeric(paste0("-", substr(guess.tmp$lon, 3,7)))
    guess.tmp$lat <- as.numeric(dimnames(fcomp.guess[[i]])[[1]])
    guess.tmp$Model <- as.factor("LPJ-GUESS")
    guess.tmp$Time  <- as.factor(names(fcomp.guess)[i])
    
    ed.tmp <- stack(data.frame(apply(fcomp.ed[[i]][,,c(9:11)], c(1,2), sum, na.rm=F)))
    names(ed.tmp) <- c("Decid", "lon")
    ed.tmp$Evg   <- stack(data.frame(apply(fcomp.ed[[i]][,,c(7:8)], c(1,2), sum, na.rm=F)))[,1]
    ed.tmp$Grass <- stack(data.frame(apply(fcomp.ed[[i]][,,c(1,5,12:16)], c(1,2), sum, na.rm=F)))[,1]
    ed.tmp$lon <- as.numeric(paste0("-", substr(ed.tmp$lon, 3,7)))
    ed.tmp$lat <- as.numeric(dimnames(fcomp.ed[[i]])[[1]])
    ed.tmp$Model <- as.factor("ED2")
    ed.tmp$Time  <- as.factor(names(fcomp.ed)[i])
    
    if(i==1){
      fcomp.out <- rbind(wsl.tmp, guess.tmp, ed.tmp)
    } else {
      fcomp.out <- rbind(fcomp.out, rbind(wsl.tmp, guess.tmp, ed.tmp))
    }
  }
  
  # interpoalte ED
  fcomp.out <- interp.mod(dat.all=fcomp.out, mod.missing="ED2", mod.mask="LPJ-WSL", var.interp="Decid")
  fcomp.out <- interp.mod(dat.all=fcomp.out, mod.missing="ED2", mod.mask="LPJ-WSL", var.interp="Evg")
  fcomp.out <- interp.mod(dat.all=fcomp.out, mod.missing="ED2", mod.mask="LPJ-WSL", var.interp="Grass")
  fcomp.out$flag <- as.factor(fcomp.out$flag)
  summary(fcomp.out)
  
  png(file.path(path.figs, "Fcomp_Decid.png"), height=8.5, width=11, units="in", res=180)
  ggplot(fcomp.out[complete.cases(fcomp.out),]) +
    facet_grid(Model ~ Time, drop=T) +
    geom_raster(data=fcomp.out[is.na(fcomp.out$flag), ], aes(x=lon, y=lat, fill=Decid)) +
    geom_raster(data=fcomp.out[complete.cases(fcomp.out) & fcomp.out$flag=="interpolated", ], aes(x=lon, y=lat, fill=Decid), alpha=0.8) +
    coord_equal(ratio=1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(paste0("Date Created: ", Sys.Date())) +
    theme(panel.background=element_blank())
  dev.off()
  
  png(file.path(path.figs, "Fcomp_Evergreen.png"), height=8.5, width=11, units="in", res=180)
  ggplot(fcomp.out[complete.cases(fcomp.out),]) +
    facet_grid(Model ~ Time, drop=T) +
    geom_raster(data=fcomp.out[is.na(fcomp.out$flag), ], aes(x=lon, y=lat, fill=Evg)) +
    geom_raster(data=fcomp.out[complete.cases(fcomp.out) & fcomp.out$flag=="interpolated", ], aes(x=lon, y=lat, fill=Evg), alpha=0.8) +
    coord_equal(ratio=1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(paste0("Date Created: ", Sys.Date())) +
    theme(panel.background=element_blank())
  dev.off()
  
  png(file.path(path.figs, "Fcomp_Grass.png"), height=8.5, width=11, units="in", res=180)
  ggplot(fcomp.out[complete.cases(fcomp.out),]) +
    facet_grid(Model ~ Time, drop=T) +
    geom_raster(data=fcomp.out[is.na(fcomp.out$flag), ], aes(x=lon, y=lat, fill=Grass)) +
    geom_raster(data=fcomp.out[complete.cases(fcomp.out) & fcomp.out$flag=="interpolated", ], aes(x=lon, y=lat, fill=Grass), alpha=0.8) +
    coord_equal(ratio=1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(paste0("Date Created: ", Sys.Date())) +
    theme(panel.background=element_blank())
  dev.off()
  
  # --------------------------
  # Individual models by PFT
  # --------------------------
  for(i in 1:length(fcomp.wsl)){
    for(j in 1:dim(fcomp.wsl[[i]])[3]){
      wsl.tmp <- stack(data.frame(fcomp.wsl[[i]][,,j]))
      names(wsl.tmp) <- c("Fcomp", "lon")
      wsl.tmp$lon <- as.numeric(paste0("-", substr(wsl.tmp$lon, 3,7)))
      wsl.tmp$lat <- as.numeric(dimnames(fcomp.wsl[[i]])[[1]])
      wsl.tmp$Model <- as.factor("LPJ-WSL")
      wsl.tmp$Time  <- as.factor(names(fcomp.wsl)[i])
      wsl.tmp$PFT   <- as.factor(dimnames(fcomp.wsl[[i]])[[3]][j])
      
      if(j==1){
        wsl.tmp2 <- wsl.tmp
      } else {
        wsl.tmp2 <- rbind(wsl.tmp2, wsl.tmp)
      }
    }
    wsl.tmp <- wsl.tmp2

    for(j in 1:dim(fcomp.guess[[i]])[3]){
      guess.tmp <- stack(data.frame(fcomp.guess[[i]][,,j]))
      names(guess.tmp) <- c("Fcomp", "lon")
      guess.tmp$lon <- as.numeric(paste0("-", substr(guess.tmp$lon, 3,7)))
      guess.tmp$lat <- as.numeric(dimnames(fcomp.guess[[i]])[[1]])
      guess.tmp$Model <- as.factor("LPJ-GUESS")
      guess.tmp$Time  <- as.factor(names(fcomp.guess)[i])
      guess.tmp$PFT   <- as.factor(dimnames(fcomp.guess[[i]])[[3]][j])
      
      if(j==1){
        guess.tmp2 <- guess.tmp
      } else {
        guess.tmp2 <- rbind(guess.tmp2, guess.tmp)
      }
    }
    guess.tmp <- guess.tmp2
    
    for(j in 1:dim(fcomp.ed[[i]])[3]){
      ed.tmp <- stack(data.frame(fcomp.ed[[i]][,,j]))
      names(ed.tmp) <- c("Fcomp", "lon")
      ed.tmp$lon <- as.numeric(paste0("-", substr(ed.tmp$lon, 3,7)))
      ed.tmp$lat <- as.numeric(dimnames(fcomp.ed[[i]])[[1]])
      ed.tmp$Model <- as.factor("ED2")
      ed.tmp$Time  <- as.factor(names(fcomp.ed)[i])
      ed.tmp$PFT   <- as.factor(dimnames(fcomp.ed[[i]])[[3]][j])
      
      if(j==1){
        ed.tmp2 <- ed.tmp
      } else {
        ed.tmp2 <- rbind(ed.tmp2, ed.tmp)
      }
    }
    ed.tmp <- ed.tmp2
    
    if(i==1){
      fcomp.out2 <- rbind(wsl.tmp, guess.tmp, ed.tmp)
    } else {
      fcomp.out2 <- rbind(fcomp.out2, rbind(wsl.tmp, guess.tmp, ed.tmp))
    }
  }

  # interpoalte ED
  fcomp.out2$flag <- NA
  for(p in unique(fcomp.out2[fcomp.out2$Model=="ED2", "PFT"])){
    if(!p %in% unique(fcomp.out2[fcomp.out2$Model=="ED2", "PFT"])[c(5:6,8:11)]) next
    print(paste0(" *** PFT: ", p, " *** "))
    fcomp.out2[(fcomp.out2$Model=="ED2" & fcomp.out2$PFT==p) | (fcomp.out2$Model=="LPJ-WSL" & fcomp.out2$PFT=="TeNE"),] <- interp.mod(dat.all=fcomp.out2[(fcomp.out2$Model=="ED2" & fcomp.out2$PFT==p) | (fcomp.out2$Model=="LPJ-WSL"  & fcomp.out2$PFT=="TeNE"),], mod.missing="ED2", mod.mask="LPJ-WSL", var.interp="Fcomp")
  }
  fcomp.out2$flag <- as.factor(fcomp.out2$flag)
  summary(fcomp.out2)
  
  guess.use <- c("BNE", "BINE", "BIBS", "TeBS", "TeIBS", "TeBE", "C3G")
  png(file.path(path.figs, "Fcomp_LPJ-GUESS.png"), height=11, width=8.5, units="in", res=180)
  ggplot(fcomp.out2[fcomp.out2$Model=="LPJ-GUESS" & fcomp.out2$PFT %in% guess.use & complete.cases(fcomp.out2),]) +
    facet_grid(PFT ~ Time, drop=T) +
    geom_raster(data=fcomp.out2[fcomp.out2$Model=="LPJ-GUESS" & fcomp.out2$PFT %in% guess.use & is.na(fcomp.out2$flag), ], aes(x=lon, y=lat, fill=Fcomp)) +
    geom_raster(data=fcomp.out2[fcomp.out2$Model=="LPJ-GUESS" & fcomp.out2$PFT %in% guess.use & complete.cases(fcomp.out2) & fcomp.out2$flag=="interpolated", ], aes(x=lon, y=lat, fill=Fcomp), alpha=0.6) +
    coord_equal(ratio=1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(paste0("Date Created: ", Sys.Date())) +
    theme(panel.background=element_blank())
  dev.off()
  
  wsl.use <- c("BNE", "BBS", "TeBS", "C3G")
  png(file.path(path.figs, "Fcomp_LPJ-WSL.png"), height=8.5, width=11, units="in", res=180)
  ggplot(fcomp.out2[fcomp.out2$Model=="LPJ-WSL" & fcomp.out2$PFT %in% wsl.use & complete.cases(fcomp.out2),]) +
    facet_grid(PFT ~ Time, drop=T) +
    geom_raster(data=fcomp.out2[fcomp.out2$Model=="LPJ-WSL" & fcomp.out2$PFT %in% wsl.use & is.na(fcomp.out2$flag), ], aes(x=lon, y=lat, fill=Fcomp)) +
    geom_raster(data=fcomp.out2[fcomp.out2$Model=="LPJ-WSL" & fcomp.out2$PFT %in% wsl.use & complete.cases(fcomp.out2) & fcomp.out2$flag=="interpolated", ], aes(x=lon, y=lat, fill=Fcomp), alpha=0.6) +
    coord_equal(ratio=1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(paste0("Date Created: ", Sys.Date())) +
    theme(panel.background=element_blank())
  dev.off()
  
  ed.use <- c("grass.c3.temp", "pine.north", "conifer.late", "temp.decid.early", "temp.decid.mid", "temp.decid.late")
  png(file.path(path.figs, "Fcomp_ED2.png"), height=11, width=8.5, units="in", res=180)
  ggplot(fcomp.out2[fcomp.out2$Model=="ED2" & fcomp.out2$PFT %in% ed.use & complete.cases(fcomp.out2),]) +
    facet_grid(PFT ~ Time, drop=T) +
    geom_raster(data=fcomp.out2[fcomp.out2$Model=="ED2" & fcomp.out2$PFT %in% ed.use & is.na(fcomp.out2$flag), ], aes(x=lon, y=lat, fill=Fcomp)) +
    geom_raster(data=fcomp.out2[fcomp.out2$Model=="ED2" & fcomp.out2$PFT %in% ed.use & complete.cases(fcomp.out2) & fcomp.out2$flag=="interpolated", ], aes(x=lon, y=lat, fill=Fcomp), alpha=0.8) +
    coord_equal(ratio=1) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(paste0("Date Created: ", Sys.Date())) +
    theme(panel.background=element_blank())
  dev.off()
  # --------------------------
  
}
# -----------------------------