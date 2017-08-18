# Alexander ZIzka & Daniele Silvestro 2017-04-03

# libraries
library(raster)
library(maptools)
library(tidyverse)
library(ggplot2)
library(rgeos)
library(geiger)
library(phytools)
library(phangorn)
library(biogeo)
library(stringr)

PlotTraitTime <- function(dat, tipdata, branchingtimes) {
  
  plo <- dat %>% 
    mutate(time_new = cut(time, breaks = seq(-47, 0, by = 1), labels = seq(-46.5, -0.5, by = 1))) %>% 
    mutate(time_new = parse_number(time_new)) %>% 
    mutate(time_new = ifelse(time_new <= -37.5, -43.5, time_new)) %>% 
    group_by(time_new) %>% 
    summarize(lwr = quantile(latitude, probs = 0.025), upr = quantile(latitude, probs = 0.975))
  
  # create format for the polygon
  plo95 <- data.frame(xvals = c(plo$time_new, rev(plo$time_new)), yval = c(plo$lwr, rev(plo$upr)))
  
  # extract tips and fossils
  pts <- data.frame(lat = tipdata, time = branchingtimes) %>% 
    mutate(ID = ifelse(time > -1e-04, "Recent", "Fossil"))
  
  # The plot
  ggplot() + 
    geom_polygon(data = plo95, aes(x = xvals, y = yval), fill = "grey50") +
    geom_point(data = pts,aes(x = time, y = lat, col = ID, shape = ID)) + 
    xlab("Time") + 
    ylab("Latitude") + 
    theme_bw() + 
    theme(legend.title = element_blank())
}

# Plotting ranges on Map Thats a plotting function for the range maps, based on a data table as output from
# the plotpointsal script and a time interval

MapRanges <- function(dat, min.time, max.time, fos.loc, lon.min = -120, lon.max = -35) {
  # The map load data
  data(wrld_simpl)
  e <- extent(-110, -30, -57, 30)
  
  map.plo <- dat %>% 
    filter(time > min.time & time < max.time) %>%
    summarize(minage = quantile(latitude, probs = 0.025), 
              maxage = quantile(latitude, probs = 0.975))
  
  # this is your min max input data
  lats.in <- map.plo %>% 
    mutate(ID = parse_character(seq(1, nrow(map.plo)))) %>% 
    mutate(colo = seq(1, nrow(map.plo)))
  
  # background map
  data(wrld_simpl)
  sa <- wrld_simpl %>% 
    gBuffer(width = 1e-05) %>% 
    aggregate() %>% 
    crop(e)
  
  # create the stupid polygons
  lats <- split(lats.in, f = lats.in$ID)
  
  pp <- lapply(lats, function(k) {
    dums <- data.frame(lon = c(lon.min, lon.min, lon.max, lon.max, lon.min), 
                       lat = unlist(c(k[1], k[2], k[2], k[1], k[1])))
    Polygons(list(Polygon(dums)), ID = k$ID)
  })
  
  pp <- SpatialPolygons(pp) %>% gIntersection(sa)
  
  #pick relevant fossils
  fos.p <- filter(fos.loc, mid.age> min.time & mid.age < max.time) %>%
    mutate(lab =  paste(species, mid.age, sep = " "))
  
  # plot
  sa <- fortify(sa)
  poly <- fortify(pp)
  
  ggplot() + 
    geom_polygon(data = sa, aes(x = long, y = lat, group = group), fill = "green") + 
    geom_polygon(data = poly, aes(x = long, y = lat, group = group), fill = "yellow", col = "pink") + 
    geom_polygon(data = sa, aes(x = long, y = lat, group = group), colour = "grey20", fill = NA) + 
    geom_point(data = fos.p, aes(x = decimallongitude, y = decimallatitude), col = "purple", size = 3) +
    # geom_label(data = fos.p, aes(x = decimallongitude, y = decimallatitude, label = lab), nudge_x = +15) +
    coord_fixed() + 
    ggtitle(sprintf("%s - %s Ma", abs(min.time), abs(max.time))) + 
    xlim(-110, -30) + 
    ylim(-57, 30) + 
    xlab("Longitude") + 
    ylab("Latitude") + 
    theme_bw()
}


# run the analyses data
dat <- read_delim("ancestral_lat_all_samples.txt", delim = "\t")

# Traits through time
load("platyrrhines0221corrected_names.trees1LAT.rda")

bt = nodeHeights(tree)
print(length(bt))
root_age <- max(bt[, 2])
names(root_age) <- ntips + 1
bt = as.numeric(bt[, 2] - max(bt[, 2]), digits = 4)
names(bt) <- tree$edge[, 2]
bt <- append(bt, -root_age)
bt = bt[order(as.numeric(names(bt)))]
bt = bt[1:ntips]
names(bt) <- tree$tip.label
bt_ord = bt[order(names(bt))]

data = data[order(names(data))]/0.05

PlotTraitTime(dat = dat, tipdata = data, branchingtimes = bt_ord)

# maps
##Load the fossil coordinates and ages
setwd("/Users/danielesilvestro/Dropbox_Personal/fossilizedBM/code/plot_empirical_results")
fos <- read_csv("Fossil platyrrhine latitude-longitude.csv")
fos.age <- read_csv("platy_fossil_ages.csv")
fos <- fos%>%
  bind_cols(data.frame(lat = str_split_fixed(fos$Latitude, " ", n = 4)))%>%
  bind_cols(data.frame(lon = str_split_fixed(fos$Longitude, " ", n = 4)))%>%
  mutate(lat.1 = parse_number(lat.1))%>%
  mutate(lat.2 = parse_number(lat.2))%>%
  mutate(lat.3 = parse_number(lat.3))%>%
  mutate(lon.1 = parse_number(lon.1))%>%
  mutate(lon.2 = parse_number(lon.2))%>%
  mutate(lon.3 = parse_number(lon.3)) %>%
  mutate(decimallongitude = dms2dd(dd = lon.1, mm = lon.2, ss = lon.3, ns = lon.4))%>%
  mutate(decimallatitude = dms2dd(dd = lat.1, mm = lat.2, ss = lat.3, ns = lat.4))%>%
  dplyr::select(species, decimallongitude, decimallatitude)%>%
  inner_join(fos.age, by = "species")%>%
  mutate(mid.age = mid.age.ma * -1)
  
write.table(fos,file="fossil_coordinates.txt",sep="\t",quote=F)


#Plot maps
MapRanges(dat, min.time = -18, max.time = -15, fos.loc = fos, lon.min = -120, lon.max = -35)
MapRanges(dat, min.time = -2, max.time = 0, fos.loc = fos, lon.min = -120, lon.max = -35)
MapRanges(dat, min.time = -5, max.time = -2, fos.loc = fos, lon.min = -120, lon.max = -35)
MapRanges(dat, min.time = -4.6, max.time = -4.5, fos.loc = fos, lon.min = -120, lon.max = -35)
