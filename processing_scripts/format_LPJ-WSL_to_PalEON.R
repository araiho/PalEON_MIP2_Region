# -------------------------------------------------------
# Script to convert the Raw LPJ-WSL output into the PalEON Format
# Christy Rollinson, crollinson@gmail.com, Feb 2016
# -------------------------------------------------------


# -----------------------------
# Loading handy-dandy libraries
# -----------------------------
library(ncdf4)
# -----------------------------

# -----------------------------
# Filepaths etc
# -----------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP2_Region")

path.raw    <- "phase2_model_output/LPJ-WSL/V1_rawformat"
path.paleon <- "phase2_model_output/LPJ-WSL/LPJ-WSL.v1"

# If we don't have a folder for the output, create it now
if(!dir.exists(path.paleon)) dir.create(path.paleon)
# -----------------------------

# -----------------------------
# 1. Reading in the raw lpj-wsl format
# -----------------------------
files.raw <- dir(path.raw, ".nc")
files.raw

# Create a blank list to store everything in
wsl.raw <- list()

# storing some global factors
nc.yr <- nc_open(file.path(path.raw, "LPJ_landCoverFrac.nc"))
nc.mo <- nc_open(file.path(path.raw, "LPJ_NPP.nc"))

wsl.raw$lat   <- ncvar_get(nc.yr, "lat")
wsl.raw$lon   <- ncvar_get(nc.yr, "lon")
wsl.raw$Year  <- ncvar_get(nc.yr, "time") + 850
wsl.raw$Month <- ncvar_get(nc.mo, "time")
nc_close(nc.yr); nc_close(nc.mo)

for(i in 1:length(files.raw)){
  ncT <- nc_open(file.path(path.raw, files.raw[i]))
  wsl.raw[[names(ncT$var)]] <- ncvar_get(ncT, names(ncT$var))
  nc_close(ncT)
}
# Note: MRSO gets overwritten by msl, which is fine since we'd rather have by layer
# -----------------------------


# -----------------------------
# 2. Changing units to match PalEON specs
# -----------------------------
sec2yr <- (60*60*24*365)
sec2day <- (60*60*24)
months <- rep(1:12, length.out=length(wsl.raw$Month))

dpm    <- c(31,28,31,30,31,30,31,31,30,31,30,31) # days per month
summary(months)

flux.mo  <- c("gpp", "nbp", "npp", "ra", "rh", "tran", "mrro", "evapotrans", "tran")
flux.yr  <- c("fFire", "fGrazing", "fLuc")
# # H2O: 1 mm = 1 kg/m2 : 1 g/cm3 * 0.1 cm = 0.1 g/cm2 = 0.1 g/cm2 * 1e-3 g/kg  * 1e4 cm2/m2 = 1 kg/m2

# Converting fluxes in per month to per second
for(v in flux.mo){
  print(paste0(" ====== Unit Conversion: ", v, " ====== "))
  for(i in unique(months)){ # because the number of days per month varies, we need to go by each month
    month.now <- which(months==i)
    wsl.raw[[v]][,,month.now] <- wsl.raw[[v]][,,month.now]*1/(sec2day*dpm[i])
  }
}

# converting fluxes in per year to per second
for(v in flux.yr){
  print(paste0(" ====== Unit Conversion: ", v, " ====== "))
  wsl.raw[[v]] <- wsl.raw[[v]] * 1/sec2yr
}

# Coverting soil temp to K
wsl.raw$tsl <- wsl.raw$tsl + 273.15
# -----------------------------


# -----------------------------
# Saving as PalEON-Formatted .nc files
# -----------------------------
# -------------
# Sticking everythign into a paleon-labeled list
# -------------
# PFT labels:
PFTs <- data.frame(PFT.Num= 1:dim(wsl.raw$landCoverFrac)[4], PFT=c("TrBE", "TrBR", "TeNE", "TeBE", "TeBS", "BNE", "BBS", "C3G", "C4G"))

wsl.paleon <- list()
wsl.paleon$Fcomp       <- wsl.raw$landCoverFrac
wsl.paleon$TotLivBiom  <- wsl.raw$cVeg
wsl.paleon$TotSoilCarb <- wsl.raw$cSoil
wsl.paleon$GPP         <- wsl.raw$gpp
wsl.paleon$AutoResp    <- wsl.raw$ra
wsl.paleon$HeteroResp  <- wsl.raw$rh
wsl.paleon$NPP         <- wsl.raw$npp
wsl.paleon$NEE         <- wsl.raw$nbp*-1 
wsl.paleon$Fire        <- wsl.raw$fFire
wsl.paleon$Qs_sb       <- wsl.raw$mrro
wsl.paleon$Evap        <- wsl.raw$evapotrans - wsl.raw$tran
wsl.paleon$Transp      <- wsl.raw$tran
wsl.paleon$SnowDepth   <- wsl.raw$snowpack
wsl.paleon$SoilMoist   <- wsl.raw$mrso
wsl.paleon$SoilTemp    <- wsl.raw$tsl
wsl.paleon$FireArea    <- wsl.raw$burntArea
wsl.paleon$Grazing     <- wsl.raw$fGrazing
wsl.paleon$LULC        <- wsl.raw$fLuc
wsl.paleon$Litter      <- wsl.raw$cLitter
# wsl.paleon$PFT         <- PFTs$PFT
# -------------


# -------------
# Creating bins to group years by
# -------------
start.run = 850
end.run   = 2010
block.yr  = 100 # Number of years per bin

# Setting up the year bins
bins <- c(start.run, seq(from=round(start.run+block.yr-1,-2), to=end.run, by=block.yr), end.run) # Creating a vector with X year bins for the time period of interest

# Adding leading 0 to file name
yr.lab <- ifelse(nchar(bins)==3, paste0(0, bins), bins)

# A vector that helps us figure out which monthly files belong in each 100-year file
mo.yr <- vector()
for(y in unique(wsl.raw$Year)){
  mo.yr <- c(mo.yr, rep(y, 12))
}
# -------------

# -------------
# Actually writing the files!
# -------------
for(i in 1:(length(bins)-1)){  
  yr.min <- max(start.run, bins[i])
  if(i < (length(bins)-1)) yr.max <- min(end.run, bins[i+1]-1) else yr.max <- min(end.run, bins[i+1])

  print(paste0(" ======= ","Processing Years: ", yr.min, " - ", yr.max, " ======= "))
  
  # -------------
  # Setting up the variable dimensions
  # -------------
  # # NOTE: These need to change to fit each time slice
  dim.mo <- ncdim_def(name = "Month",
                      units = paste0("months since run start:", start.run),
                      vals = wsl.raw$Month[which(mo.yr>=yr.min & mo.yr<=yr.max)], # calculating the number of months in this run
                      calendar = "standard", unlim = TRUE)
  dim.yr <- ncdim_def(name = "Year",
                      units = "Calendar Year",
                      vals = wsl.raw$Year[wsl.raw$Year>=yr.min & wsl.raw$Year<=yr.max], # year labels
                      calendar = "standard", unlim = TRUE)
  
  dim.lat <- ncdim_def("lat", "degrees_east",
                       vals =  wsl.raw$lat,
                       longname = "station_latitude") 
  dim.lon <- ncdim_def("lon", "degrees_north",
                       vals = wsl.raw$lon,
                       longname = "station_longitude")
  dim.string <- ncdim_def("names", "", 1:24, create_dimvar=FALSE)
  dim.pft1 <- ncdim_def("pft", "",
                        1:length(PFTs$PFT),
                        longname = "Plant Functional Type", create_dimvar=FALSE)                 
  dim.pft <- ncdim_def("pft", "",
                       vals = PFTs$PFT.Num,
                       longname = "Plant Functional Type")                 
  dim.pft2 <- ncdim_def("pft.dims", "",
                        vals = 1:ncol(PFTs),
                        longname = "Plant Functional Type Description")                 
  dim.soil <- ncdim_def("SoilLayer", "meters",
                        vals = 1:dim(wsl.paleon$SoilMoist)[4],
                        longname = "Soil Layer")                 
  # -------------
  
  # -------------
  # Setting up the variable list
  # -------------
  var <- list() # Create a blank list for the variables
  var[[ 1]] <- ncvar_def("Fcomp", units="Fraction Area", dim=list(dim.lon, dim.lat, dim.yr, dim.pft), longname="Fractional Composition of PFTs by AREA")
  var[[ 2]] <- ncvar_def("TotLivBiom", units="kg m-2", dim=list(dim.lon, dim.lat, dim.yr), longname="Total Vegetation Biomass")
  var[[ 3]] <- ncvar_def("TotSoilCarb", units="kg m-2", dim=list(dim.lon, dim.lat, dim.yr), longname="Total Soil Carbon")
  var[[ 4]] <- ncvar_def("GPP", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.mo), longname="Gross Primary Productivity")
  var[[ 5]] <- ncvar_def("AutoResp", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.mo), longname="Autotrophic Respiration")
  var[[ 6]] <- ncvar_def("HeteroResp", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.mo), longname="Heterotrophic Respiration")
  var[[ 7]] <- ncvar_def("NPP", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.mo), longname="Net Primary Productivity") # NOTE: Not broken down by PFT
  var[[ 8]] <- ncvar_def("NEE", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.mo), longname="Net Ecosystem Exchange")
  var[[ 9]] <- ncvar_def("Fire", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.yr), longname="Fire Emissions")
  var[[10]] <- ncvar_def("Qs_sb", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.mo), longname="Total Runoff")
  var[[11]] <- ncvar_def("Evap", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.mo), longname="Total Evaporation")
  var[[12]] <- ncvar_def("Transp", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.mo), longname="Total Transpiration") # NOTE: not broken down by PFT
  var[[13]] <- ncvar_def("SnowDepth", units="mm (=kg/m2)", dim=list(dim.lon, dim.lat, dim.mo), longname="Total Snow/Water Depth") # NOTE: Units differ from the protocol sheet
  var[[14]] <- ncvar_def("SoilMoist", units="kg m-2", dim=list(dim.lon, dim.lat, dim.mo, dim.soil), longname="Soil Moisture") # NOTE: Units differ from the protocol sheet
  var[[15]] <- ncvar_def("SoilTemp", units="K", dim=list(dim.lon, dim.lat, dim.mo), longname="Soil Temperature")
  var[[16]] <- ncvar_def("FireArea", units="Area", dim=list(dim.lon, dim.lat, dim.yr), longname="Fire Emissions")
  var[[17]] <- ncvar_def("Grazing", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.yr), longname="Fire Emissions")
  var[[18]] <- ncvar_def("LULC", units="kg m-2 s-1", dim=list(dim.lon, dim.lat, dim.yr), longname="Fire Emissions")
  var[[19]] <- ncvar_def("Litter", units="kg m-2", dim=list(dim.lon, dim.lat, dim.yr), longname="Fire Emissions")
#   var[[20]] <- ncvar_def("PFT", units="", dim=list(dim.pft2, dim.string), longname="Plant Functional Type", prec="char")
  # -------------
  
  nc <- nc_create(file.path(path.paleon, paste0("LPJ-WSL_", yr.lab[i], ".nc")), var)
  for(v in 1:length(var)){
    if(names(wsl.paleon)[v]=="PFT"){
      ncvar_put(nc, var[[v]], wsl.paleon[[v]])
    } else { 
      # Figuring out which months/years to pull
      # I think in all cases time is the 3rd dimension
      if(dim(wsl.paleon[[v]])[3]>length(wsl.raw$Year)){ 
        time.use <- which(mo.yr>=yr.min & mo.yr<=yr.max) 
      } else { 
        time.use <- which(wsl.raw$Year>=yr.min & wsl.raw$Year<=yr.max) 
      }
      
      # Figuring out where we need to place the time index
      if(length(dim(wsl.paleon[[v]]))==3){
        ncvar_put(nc, var[[v]], wsl.paleon[[v]][,,time.use])
      } else {
        ncvar_put(nc, var[[v]], wsl.paleon[[v]][,,time.use,])
      }
    }
  }
  nc_close(nc)
}
# -------------


# -----------------------------

