library(BIEN)
library(gdata)
library(ggplot2)
library(maps)
library(plyr)
library(raster)

do_analysis <- function()
{
  clean_barker=FALSE
  resolve_occurrences = FALSE
  get_occurrences = FALSE
  get_invasives=FALSE
  get_flowering=TRUE
  
  if(clean_barker)
  {
    #Barker's original file
    ploidy_list <- read.csv("barker_list.csv")
    
    #remove tissues with no chromosomes sampled
    valid <- which(ploidy_list$inferred_n > 0)
    ploidy_list <- ploidy_list[valid,]
    remove("valid")
    
    if(get_flowering)
    {
      ploidy_list <- ploidy_list[which(ploidy_list$Plant.type == "Angiosperm"),]
    }
    
    #Weighted average across duplicate names by num. observed records
    ploidy_list <- ddply(ploidy_list, .(resolved_binomial), 
                         function(x) weighted.mean(x$inferred_n, x$num_records))
    
    colnames(ploidy_list) <- c("resolved_binomial", "ploidy")
    
    #Convert to string
    ploidy_list$resolved_binomial <- 
      as.character(ploidy_list$resolved_binomial)
    #Change underscores to spaces
    ploidy_list$resolved_binomial <- 
      gsub("_", " ", ploidy_list$resolved_binomial) 
    
    if(get_flowering)
    {
      write.csv(ploidy_list, file="clean_barker_flowering.csv")
    }
    else
    {
      write.csv(ploidy_list, file="clean_barker.csv")  
    }
    
  }else{
    if(get_flowering)
    {
      ploidy_list <- read.csv("clean_barker_flowering.csv")
    }
    else
    {
      ploidy_list <- read.csv("clean_barker.csv") 
    }
    
    #Resolved names from TNRS
    names_resolved <- read.csv("names_full.csv", sep="\t")
    names_resolved$Name_matched <- 
      as.character(names_resolved$Name_matched)
    names_resolved$Name_submitted <- 
      as.character(names_resolved$Name_submitted)
    
    #Merge ploidy data with TNRS names
    plants <- merge(ploidy_list, names_resolved, by.x="resolved_binomial",by.y="Name_submitted")
    
    #remove("ploidy_list")
    #remove("names_resolved")
    gc()
    #1731 / ~100,000 matched names resolved
    res <- which(names_resolved$Name_submitted != names_resolved$Name_matched)
    
    if(resolve_occurrences)
    {
      #How many names have occurences in BIEN 
      bien <- data.frame(species=BIEN_list_all())
      bien_plants <- merge(bien, plants, by.x="species", by.y="Name_matched")
      
      remove("bien")
      remove("plants")
      gc()
      #Use resolved names to retrieve all occurences (non-natives, global)
      
      #Get occurrences for matched species, remove occurences without lat/long, merge with ploidy data
      number <- seq(1, 1000)
      number <- seq(1, 5000)
      number <- seq(1, nrow(bien_plants))
      
      if(get_occurrences)
      {
        print("Getting BIEN occurences")
        #Kinda effed up - BIEN contains numerous NA's and latitudes > 90
        occ <- BIEN_occurrence_species(bien_plants$species[number], only.new.world=FALSE,
                                       native.status=TRUE, natives.only=FALSE)
        occ <- occ[which(!is.na(occ$latitude) | !is.na(occ$longitude)),]
        occ <- occ[which(occ$latitude>-90 & occ$latitude<90 & occ$longitude>-180 & occ$longitude<180),]
        write.csv(occ, "occurences.csv")
        #quit()
      }else
      {
        occ <- read.csv("occurences.csv")
      }
      
      if(get_invasives)
      {
        inv <- occ[which(occ$native_status == "I" | occ$native_status == "Ie"),]
        write.csv(inv, "invasives.csv")
        
        nat <- occ[which(occ$native_status == "N" | occ$native_status == "Ne"),]
        write.csv(nat, "natives.csv")
        quit()
      }
      if(get_flowering)
      {
        occ <- occ[which(occ$higher_plant_group == "flowering plants"),]
        mm<-merge(occ, bien_plants[number,], by.x="scrubbed_species_binomial",
                  by.y="resolved_binomial")[c("latitude", "longitude", "ploidy")]
        write.csv(mm, "resolved_occurences_angiosperms.csv")
        quit()
      }
      mm<-merge(occ, bien_plants[number,], by.x="scrubbed_species_binomial",
                by.y="resolved_binomial")[c("latitude", "longitude", "ploidy")]
      write.csv(mm, "resolved_occurences.csv")
    }else
    {
      mm<-read.csv("resolved_occurences_all_taxa.csv")[,-1]
      #mm <- read.csv("resolved_invasive_occurences.csv")[,-1]
      #mm<-read.csv("resolved_native_occurences.csv")[,-1]
      
      print("Total occurence records")
      print(nrow(mm))
      
      #Apply grid-based downsampling
      #Problem: grid cells undersampled at higher latitudes
      #1000 resolution grid makes cells approx. .1 lat x .2 long = ~2000km^2
      res = 10000
      
      #insert values to scale cut interval to lat=[-90,90], long=[-180,180]
      mm <- rbind(mm, data.frame(latitude=-90, longitude=-180, ploidy=0))
      mm <- rbind(mm, data.frame(latitude=90, longitude=180, ploidy=0))
      
      #Overlay grid of resolution 'res'
      mm$y <- cut(mm$latitude, breaks=res, labels=FALSE)
      mm$x <- cut(mm$longitude, breaks=res, labels=FALSE)
      mm$index <- ((mm$x-1)*res)+(mm$y-1)
      mm <- mm[which(mm$ploidy < 120),]
      
      occ_mean <- mean(mm$ploidy)
      
      #Downsample within grid cells by index
      mean_ploidies <- ddply(mm, .(index), function(x) mean(x$ploidy))
      max_ploidies <- ddply(mm, .(index), function(x) max(x$ploidy))
      colnames(mean_ploidies) <- c("index", "ploidy")
      
      #Convert grid cells back to lat/long
      retrieve_longitude <- function(x, res) {long_offset = 180/(res*2) 
      return ((((floor((x$index)/res)+1)/res)-0.5)*360) + long_offset}
      retrieve_latitude <- function(x, res) {lat_offset = 90/(res*2) 
      return ((((((x$index) %% res)+1)/res)-0.5)*180) + lat_offset}
      mean_ploidies$longitude <- retrieve_longitude(mean_ploidies, res)
      mean_ploidies$latitude <- retrieve_latitude(mean_ploidies, res)
      #mean_ploidies$ploidy <- mean_ploidies$ploidy
      
      
      print("Downsampled records")
      print(nrow(mean_ploidies))
      
      #up <- subset(dp, ploidy > 20)
      #down <- subset(dp, ploidy <= 20)
      
      #Get world map data and draw it
      world <- map_data("world")
      wm <- ggplot(world, aes(x=long, y=lat, group=group)) + 
        geom_polygon(fill="white", colour="black", size=0.05) + coord_quickmap()
      
      #Plot occurences by ploidy
      #pm <- wm + geom_point(data=mean_ploidies, aes(x=longitude, y=latitude, group=NULL, color=ploidy), size=0.1) +
      #  scale_color_gradientn(colors=c('#dcdcdc','#b4b4b4','#8d8d8d','#686868','#444444','#252525','#000000'))
      
      #ggsave(file="world.tiff", width=48, height=24, units="cm") 
    }
  }
}

#Plotting functions take aggregated data such as 'mean_ploidies'
plot_latitude <- function(ploidies)
{
  res <- lm(ploidy ~ abs(latitude), data=mean_ploidies)
  vv <- paste("R^2: ", summary(res)$r.squared)
  vv <- paste(vv, "\nP < 0.001")
  ggplot(mean_ploidies, aes(x=abs(latitude), y=ploidy)) + geom_point(alpha=0.1) + 
    geom_smooth(method="lm") + annotate("text", x=75, y=75, label=vv)
}

plot_bioclim_mean_annual_temp <- function(ploidies)
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
  print(gj)
}
