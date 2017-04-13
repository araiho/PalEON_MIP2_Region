####
#### This script creates maps so select priority regional run sites
#### Author: Andria Dawson
#### Second Author: Ann Raiho
#### 

library(neotoma)
library(maptools)
library(ggplot2)

## map data
us.shp <- readShapeLines('~/Downloads/map_data/lat-long/statesp020.shp',
                         proj4string=CRS('+init=epsg:4326'))
us.shp@data$id <- rownames(us.shp@data)
us.fort <- fortify(us.shp, region='id') 

## grid completed
mod_coords = read.csv('~/PalEON_MIP2_Region/site_selection/data/Paleon_MIP_Phase2_ED_Order_Status.csv')

## kh sites
tr_coords = read.csv('~/PalEON_MIP2_Region/site_selection/data/KH_Treering_sites_PALEON_model_grid.csv')

## stepps sites
stepps_change    = read.csv('~/PalEON_MIP2_Region/site_selection/data/stepps_site_change.csv')[,-c(1,2)]
stepps_site_meta = read.csv('~/PalEON_MIP2_Region/site_selection/data/stepps_sites_meta.csv')[,-1]
stepps_dat       = merge(stepps_change, stepps_site_meta, by='id')

## refab
refab = read.csv('~/PalEON_MIP2_Region/site_selection/data/refab_second_deriv_with_meta.csv')

# get dataset ids for refab sites
pol_list = sapply(refab$site.id, get_dataset, datasettype="pollen")

dataset_id = unlist(sapply(pol_list, 
                           function(x) if(length(x)>1){ sapply(x, function(y) y$dataset.meta$dataset.id) } 
                           else { x[[1]]$dataset.meta$dataset.id }))

site_id = rep(refab$site.id, times = sapply(pol_list, length))
ids     = data.frame(site_id=site_id, dataset_id=dataset_id)

stepps_dat$site_id = ids[match(stepps_dat$id, ids$dataset_id), 'site_id']

sites_dat = data.frame(stepps_dat, refab[match(stepps_dat$site_id, refab$site.id),c('settlebiomass', 'second_deriv1000')])
sites_dat$refab = !is.na(sites_dat$settlebiomass)


##### Adding more ED coordinates per group disscussion 04/07/2017
new.ed <- c(which(mod_coords$lon==-88.25&mod_coords$lat==44.75),
            which(mod_coords$lon==-88.25&mod_coords$lat==46.75),
            which(mod_coords$lon==-88.25&mod_coords$lat==42.25),
            which(mod_coords$lon==-91.25&mod_coords$lat==43.75),
            which(mod_coords$lon==-91.25&mod_coords$lat==45.25),
            which(mod_coords$lon==-93.25&mod_coords$lat==44.75),
            which(mod_coords$lon==-95.75&mod_coords$lat==46.25),
            which(mod_coords$lon==-92.75&mod_coords$lat==38.75),
            which(mod_coords$lon==-88.75&mod_coords$lat==38.75),
            which(mod_coords$lon==-75.25&mod_coords$lat==43.75),
            which(mod_coords$lon==-73.75&mod_coords$lat==43.75),
            which(mod_coords$lon==-72.25&mod_coords$lat==43.75),
            which(mod_coords$lon==-69.25&mod_coords$lat==46.75))

new.ed.mat <- mod_coords[new.ed,]
write.csv(new.ed.mat,file='new.ed.mat.csv')

mod_coords$done = !sapply(mod_coords$runs, function(x) x=="")
mod_coords$done[new.ed] <- "DO"
mod_coords$done[grep('ERROR',mod_coords$spininital)] <- "ERROR"


sites_dat$bc_binned = cut(stepps_dat$bc,
                          breaks=3,
                          labels=FALSE) # this is sorted low medium high bre-curtis distance to show where the biggest changes occured in the STEPPS2 data product

sites_dat$d2_binned = cut(sites_dat$second_deriv1000,
                          breaks=3,
                          labels=FALSE)

### Plot Pollen Sites
p <- ggplot() +  geom_path(data=us.fort, aes(x=long, y=lat, group=group),  colour='grey55') + xlim(c(-99,-66)) + ylim(c(36, 50)) 
p <- p + geom_point(data=sites_dat, aes(x=long, y=lat, colour='black'))
p <- p + geom_point(data=tr_coords, aes(x=longitude, y=latitude, colour='brown'))
p <- p + theme_bw()
print(p)

#bre-curtis distance to measure change in STEPPS estimates
hist(stepps_dat$bc, breaks=20) 


#### Plot tree ring and stepps sites (stepps sites colored by bre-curtis, yellow little change, red high change)

p <- ggplot() +  geom_path(data=us.fort, aes(x=long, y=lat, group=group),  colour='grey55') + xlim(c(-99,-66)) + ylim(c(36, 50)) 
p <- p + geom_point(data=sites_dat, aes(x=long, y=lat, colour=factor(bc_binned)))
p <- p + scale_colour_manual(values=c('#fed976', '#fd8d3c', '#800026'))
p <- p + geom_point(data=tr_coords, aes(x=longitude, y=latitude), shape=8, colour='black')
p <- p + theme_bw()
print(p)
ggsave('figures/paleon_sites_stepps_change.pdf')+coord_equal()
dev.off()


#### Plot tree ring and refab sites (stepps sites colored by bre-curtis, yellow little change, red high change)
#### We can run stepps sites with refab if we want. This is just the subset used for whole Holocene runs

p <- ggplot() +  geom_path(data=us.fort, aes(x=long, y=lat, group=group),  colour='grey55') + xlim(c(-99,-66)) + ylim(c(36, 50)) 
p <- p + geom_point(data=subset(sites_dat, !is.na(d2_binned)), aes(x=long, y=lat, colour=factor(d2_binned)))
p <- p + scale_colour_manual(values=c('#fed976', '#fd8d3c', '#800026'))
p <- p + geom_point(data=tr_coords, aes(x=longitude, y=latitude), shape=8, colour='black')
p <- p + theme_bw()
print(p)
ggsave('figures/paleon_sites_refab_change.pdf')+coord_equal()
dev.off()

#### Plot Ed runs completed v. not completed
which((sites_dat$d2_binned > 1))

p <- ggplot() +  geom_path(data=us.fort, aes(x=long, y=lat, group=group),  colour='grey55') + xlim(c(-99,-66)) + ylim(c(36, 50)) 
p <- p + geom_raster(data=mod_coords, aes(x=lon, y=lat, fill=done), alpha=0.7)
p <- p + scale_fill_manual(values=c('magenta','#dfc27d', '#80cdc1'))
p <- p + geom_point(data=tr_coords, aes(x=longitude, y=latitude), colour='black')
p <- p + geom_point(data=sites_dat, aes(x=long, y=lat, colour=factor(bc_binned)), alpha=0.7)
p <- p + theme_bw()
print(p)
dev.off()

#### Plot another plot of ed runs completed v. not completed
q <- ggplot() +  geom_path(data=us.fort, aes(x=long, y=lat, group=group),  colour='grey55') + xlim(c(-99,-66)) + ylim(c(36, 50)) + 
  geom_tile(data=mod_coords, aes(x=lon, y=lat, fill=done), alpha=0.7) + 
  scale_fill_manual(values=c('magenta','red','#dfc27d', '#80cdc1')) + 
  geom_point(data=tr_coords, aes(x=longitude, y=latitude), colour='black', shape=8) + 
  geom_point(data=mod_coords[1:6,], aes(x=lon, y=lat), colour='blue', shape=19,cex=3) +
  geom_point(data=sites_dat, aes(x=long, y=lat, colour=factor(bc_binned)), alpha=0.7) + 
  scale_colour_manual(values=c('#fed976', '#fd8d3c', '#800026'),name='BreCurtisFactor') +
  theme_bw()
print(q)
ggsave('figures/paleon_sites_ed_grid.pdf')+coord_equal()
ggsave('figures/paleon_sites_ed_grid.png')+coord_equal()
dev.off()

#### Plot STEPPS and Tree Ring data with Ed Runs
#### Need to do q plot above to print this one
q <- ggplot() +  geom_path(data=us.fort, aes(x=long, y=lat, group=group),  colour='grey55') + 
  xlim(c(-98,-82)) +  
  ylim(c(40, 50)) + 
  coord_fixed() + 
  geom_tile(data=mod_coords, aes(x=lon, y=lat, fill=done), alpha=0.7) + 
  scale_fill_manual(values=c('magenta','red','#dfc27d', '#80cdc1')) + 
  geom_point(data=tr_coords, aes(x=longitude, y=latitude), colour='black', shape=8) + 
  geom_point(data=mod_coords[1:6,], aes(x=lon, y=lat), colour='blue', shape=19,cex=3) +
  geom_point(data=sites_dat, aes(x=long, y=lat, colour=factor(bc_binned)), alpha=0.7) +
  p <- p + scale_colour_manual(values=c('#fed976', '#fd8d3c', '#800026'),name='BreCurtisFactor')
  theme_bw()
print(q)
ggsave('paleon_sites_ed_grid_zoom.pdf')+coord_equal()

#### Save cleaned dataframe
sites_dat_cleaned = sites_dat[, c('id', 'stat_id', 'sitename', 'lat', 'long', 'x', 'y', 'mag', 'bc', 
                                  'settlebiomass', 'second_deriv1000', 'min_age', 'max_age', 'n_ages')]
colnames(sites_dat_cleaned) = c('id', 'stat_id', 'sitename', 'lat', 'long', 'x', 'y', 'stepps-abs_mag', 'stepps-BC', 'refab-settlebiomass', 'refab-second_deriv1000', 'pol-min_age', 'pol-max_age', 'pol-n_ages')

write.csv(sites_dat_cleaned, 'data/sites_meta.csv')
