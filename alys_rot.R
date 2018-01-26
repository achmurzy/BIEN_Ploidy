library(BIEN)
library(gdata)
library(ggplot2)
library(maps)
library(plyr)
library(raster)

#flowering: get only flowering plants
#stebbins: use stebbins fraction (greater than mode of the genus) to classift ploidy
clean_raw_barker <- function(flowering=TRUE, stebbins=TRUE)
{
  print("Computing ploidies")
  ploidy_list <- read.csv("barker_list.csv")
  print("")
  if(flowering)
  {
    ploidy_list <- ploidy_list[which(ploidy_list$Plant.type == "Angiosperm"),]
  }
  
  #Weighted average across duplicate names by num. observed records ?? Or just take mode directly??
  ploidy_list <- ddply(ploidy_list, .(resolved_binomial), 
                       function(x) weighted.mean(x$inferred_n, x$num_records))
  
  colnames(ploidy_list) <- c("resolved_binomial", "ploidy")
  
  #Convert to string
  ploidy_list$resolved_binomial <- 
    as.character(ploidy_list$resolved_binomial)
  
  if(stebbins)
  {
    #split list by genus
    ploidy_list$genus <- strsplit(ploidy_list$resolved_binomial, "_")
    ploidy_list$genus <- unlist(lapply(ploidy_list$genus, '[[', 1))
    
    #compute mode on inferred_n - this function returns the first value if all are equally frequent -- avg instead?
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    ploidy_key <- aggregate(ploidy_list$inferred_n, by=list(Genus=ploidy_list$genus), FUN="Mode")
    
    #take Stebbins fraction to assign ploidy
    ploidy_list$key_ind <- match(ploidy_list$genus, cg$Genus)
    ploidy_list$ploidy <- 0
    
    pp <- which(ploidy_list$inferred_n > ploidy_key$x[ploidy_list$key_ind])
    dp <- which(ploidy_list$inferred_n <= ploidy_key$x[ploidy_list$key_ind])
    
    #find integer multiples
    ploidy_list[pp,]$ploidy <- ploidy_list[pp,]$inferred_n / cg$x[ploidy_list[pp,]$key_ind]
    ploidy_list[dp,]$ploidy <- 1
    ploidy_list$ploidy <- ploidy_list$ploidy + 1
    
    ploidy_list$polyploid <- ploidy_list$ploidy > 1
  }
  else
    ploidy_list$ploidy <- ploidy_list$inferred_n
  
  #remove tissues with no chromosomes sampled
  #valid <- which(ploidy_list$inferred_n > 0)
  #ploidy_list <- ploidy_list[valid,]
  #remove("valid")
  
  #Change underscores to spaces
  ploidy_list$resolved_binomial <- 
    gsub("_", " ", ploidy_list$resolved_binomial) 
  
  write.csv(ploidy_list, file="clean_barker.csv")
  return(ploidy_list)
}

#ploidy_list: cleaned barker data
resolve_names <- function(ploidy_list)
{
  print("Resolving names")
  #Resolved names from TNRS
  names_resolved <- read.csv("names_full.csv", sep="\t")
  names_resolved$Name_matched <- 
    as.character(names_resolved$Name_matched)
  names_resolved$Name_submitted <- 
    as.character(names_resolved$Name_submitted)
  
  #Merge ploidy data with TNRS names - should be same length as original file
  plants <- merge(ploidy_list, names_resolved, by.x="resolved_binomial",by.y="Name_submitted")
  
  print(paste("Original name count: ", nrow(ploidy_list)))
  print(paste("Resolved name count: ", nrow(plants)))
  
  res <- which(names_resolved$Name_submitted != names_resolved$Name_matched)
  print(paste("Number of names resolved: ", length(res)))
  
  #How many names have occurences in BIEN 
  bien <- data.frame(species=BIEN_list_all())
  bien_plants <- merge(bien, plants, by.x="species", by.y="Name_matched")
  return(bien_plants)
}

#plants:species list from TNRS output
#new_world: only include new world species
#natives: only include native species at localities
get_bien_occurences <- function(plants, new_world = FALSE, natives = FALSE)
{
  print("Getting BIEN occurences")
  library(data.table)
  
  #Just grab species names and ploidies
  plants <- data.table(plants[c("species", "ploidy")])
  
  #Scrub invalid lat/lon observations
  occ <- BIEN_occurrence_species(plants$species, only.new.world=new_world, native.status=natives, natives.only=natives)
  occ <- occ[which(!is.na(occ$latitude) | !is.na(occ$longitude)),]
  occ <- occ[which(occ$latitude>-90 & occ$latitude<90 & occ$longitude>-180 & occ$longitude<180),]
  occ <- data.table(occ[c("scrubbed_species_binomial", "latitude", "longitude")])

  #Retrieve native and invasive list subsets
  if(natives)
  {
    inv <- occ[which(occ$native_status == "I" | occ$native_status == "Ie"),]
    write.csv(inv, "invasives.csv")
    nat <- occ[which(occ$native_status == "N" | occ$native_status == "Ne"),]
    write.csv(nat, "natives.csv")
  }
  
  #merge with barker data and extract geo coords and ploidy
  occ<-merge(occ, plants, by.x="scrubbed_species_binomial",
            by.y="species")[c("latitude", "longitude", "ploidy")]
  print(paste("Total occurence records", nrow(occ)))
  write.csv(occ, "occurences.csv")
  return(occ)
}

#occ: bien occurrences with ploidies
#res: grid cell resolution ; 100 res = approx 600km^2 area
#FUN: function for downsampling a vector of ploidies for a locality e.g. mean
grid_downsample <- function(occ, res = 100, FUN = mean)
{
  print("Downsampling occurences")
  
  #insert two fake observations at geographic extremes to perform cut
  occ <- rbind(occ, data.frame(latitude=-90, longitude=-180, ploidy=0))
  occ <- rbind(occ, data.frame(latitude=90, longitude=180, ploidy=0))
  
  #Overlay grid of resolution 'res'
  occ$y <- cut(occ$latitude, breaks=res, labels=FALSE)
  occ$x <- cut(occ$longitude, breaks=res, labels=FALSE)
  occ$index <- ((occ$x-1)*res)+(occ$y-1)
  
  #Downsample within grid cells by index
  ploidies <- ddply(occ, .(index), FUN(x$ploidy)) #Expect errors here
  colnames(mean_ploidies) <- c("index", "ploidy")
  
  #Convert grid cells back to lat/long
  retrieve_longitude <- function(x, res) {long_offset = 180/(res*2) 
  return ((((floor((x$index)/res)+1)/res)-0.5)*360) + long_offset}
  retrieve_latitude <- function(x, res) {lat_offset = 90/(res*2) 
  return ((((((x$index) %% res)+1)/res)-0.5)*180) + lat_offset}
  ploidies$longitude <- retrieve_longitude(ploidies, res)
  ploidies$latitude <- retrieve_latitude(ploidies, res)
  
  print(paste("Downsampled records", nrow(ploidies)))
  
  return(ploidies)
}

#Assumes located in a data folder
plot_earth <- function(dat, filename = "world.tiff")
{
  world <- map_data("world")
  wm <- ggplot(world, aes(x=long, y=lat, group=group)) + 
    geom_polygon(fill="white", colour="black", size=0.05) + coord_quickmap()
  
  #Plot occurences by ploidy
  pm <- wm + geom_point(data=dat, aes(x=longitude, y=latitude, group=NULL, color=ploidy), size=0.1) +
    scale_color_gradientn(colors=c('#dcdcdc','#b4b4b4','#8d8d8d','#686868','#444444','#252525','#000000'))
  
  ggsave(file=filename, width=48, height=24, units="cm")
}

#Plotting functions take aggregated data such as 'mean_ploidies'
plot_latitude <- function(ploidies, filename = "lat_plot.png")
{
  res <- lm(ploidy ~ abs(latitude), data=mean_ploidies)
  vv <- paste("R^2: ", summary(res)$r.squared)
  vv <- paste(vv, "\nP < 0.001")
  ggplot(mean_ploidies, aes(x=abs(latitude), y=ploidy)) + geom_point(alpha=0.1) + 
    geom_smooth(method="lm") + annotate("text", x=75, y=75, label=vv)
  ggsave(file=filename)
}

plot_bioclim_mean_annual_temp <- function(ploidies, filename = "temp_plot.png")
{
  mat <- getData('worldclim', var='bio', res=10)
  spatial_points <- SpatialPoints(ploidies[,3:4])
  ploidies$temp <- extract(mat$bio1, spatial_points) * 0.1
  res <- lm(ploidy ~ temp, data=ploidies)
  vv <- paste("R^2: ", summary(res)$r.squared)
  vv <- paste(vv, "\nP < 0.001")
  gj <- ggplot(ploidies, aes(x=temp, y=ploidy)) + geom_point(alpha=0.1) +
    #scale_x_continuous(limits = c(-20, 40)) + 
    geom_smooth(method="lm") + annotate("text", x=-5, y=75, label=vv)
  ggsave(file=filename)
}

#WARNING: setting directory only works if file is sourced, then location of script is used as root directory
rotation_analysis <- function()
{
  this.dir <- getSrcDirectory(rotation_analysis)
  setwd(this.dir)
  setwd("data/rotation_test")
  ploidies <- clean_raw_barker(flowering = TRUE, stebbins=FALSE)
  ploidies <- resolve_names(ploidies)
  ploidies <- get_bien_occurences(ploidies)
  ploidies <- grid_downsample(ploidies, res = 1000)
  setwd("../../figures/working")
  plot_earth(ploidies)
  plot_latitude(ploidies)
  plot_bioclim_mean_annual_temp(ploidies)
}

stebbins_analysis <- function()
{
  this.dir <- getSrcDirectory(rotation_analysis)
  setwd(this.dir)
  setwd("data/working")
  ploidies <- clean_raw_barker(flowering = TRUE, stebbins=TRUE)
  
  ploidies <- ploidies[which(ploidies$polyploid),]
  
  ploidies <- resolve_names(ploidies)
  ploidies <- get_bien_occurences(ploidies)
  setwd("../../figures/working")
  plot_earth(ploidies)
  plot_latitude(ploidies)
  plot_bioclim_mean_annual_temp(ploidies)
}