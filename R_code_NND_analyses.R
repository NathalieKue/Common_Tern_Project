# Packages ---------------------------------------------------------------

### Load the package you will need for this R session (do it every time you open this session)
require(dplyr)
require(lubridate)
require(geosphere)
require(rworldmap)


# here specify where your working folder is: where you save the R script, where your data file is and where the plots will be downloaded
wd <- "your_working_directory"
setwd(wd)
getwd()

# DATA FILE ---------------------------------------------------------------

d <- read.csv("positions_per_day.csv",h=T,sep=";")
d$departure_colony <- as.POSIXct(d$departure_colony, format = '%d.%m.%Y')
d$arrival_wa <- as.POSIXct(d$arrival_wa, format = '%d.%m.%Y')
d$departure_wa <- as.POSIXct(d$departure_wa, format = '%d.%m.%Y')
d$arrival_date_colony <- as.POSIXct(d$arrival_date_colony, format = '%d.%m.%Y')
d$position_date_time  <- as.POSIXct(d$position_date_time , format = '%d.%m.%Y')

################################################# Within-individual randomisation autumn ######################################################################


###Autumn migration

###calculate the NND of a track to other tracks of the same bird

d$pointID <- seq(1,nrow(d),1)

d_autumn <- d[which(d$position_date_time > d$departure_colony &   #select all point between departure colony and arrival wintering area
                      d$position_date_time < d$arrival_wa),]

d_winter <- d[which(d$position_date_time > d$arrival_wa & #select all point at wintering area
                      d$position_date_time < d$departure_wa),]


# identify minimum wintering latitude to cut tracks to same end point
autsum <- d_winter %>% group_by(track_id) %>% summarise(end_point = min(position_lat))

d_autumn <- d_autumn[which(d_autumn$position_lat > max(autsum$end_point)),]


no_points <- d_autumn %>% group_by(track_id) %>% summarise(n = n())

d_autumn$no_points <- no_points$n[match(d_autumn$track_id,no_points$track_id)]


dist_mat <- matrix(nrow = nrow(d_autumn),   #build a matrix
                   ncol = nrow(d_autumn))

colnames(dist_mat) <- d_autumn$pointID
rownames(dist_mat) <- d_autumn$pointID

#loop -> calculate the distance between all points and store in matrix

for(i in 1:nrow(dist_mat)){
  
  dists <- distHaversine(data.frame(d_autumn$position_lon[i],
                                    d_autumn$position_lat[i]),
                         data.frame(d_autumn$position_lon,
                                    d_autumn$position_lat))/1000
  dist_mat[i,] <- dists
  
  print(i)
  
}



tracks <- unique(d_autumn$track_id)  #give a unique name to each point

# for each track, calculate an NND using this loop
track_output <- data.frame(track_id = tracks,
                           bird = NA,
                           NND = NA,
                           n_tracks = NA)

d_autumn$NND <- NA
for(i in 1:length(tracks)){
  
  this_bird <- d_autumn[which(d_autumn$track_id == tracks[i]),]   #focal track of focal bird
  
  next_tracks <- unique(d_autumn$track_id[which(d_autumn$individual == this_bird$individual[1] &  # select the next tracks from the focal bird
                                  d_autumn$track_id != tracks[i])])
  
  track_output$bird[i] <- as.character(this_bird$individual[1])
  track_output$n_tracks[i] <- length(next_tracks)
  
  
  # if the number of other tracks is not 0, we can calculate a NND
  if(length(next_tracks) > 0){
    
    NND_per_track <- c()
    for(x in 1:length(next_tracks)){
      
      #select comparison track
      comparison_track <- d_autumn[which(d_autumn$track_id == next_tracks[x]),]
      
      #IMPORTANT: THE TRACK WITH THE LOWER RESOLUTION IS THEN COMPARED TO THE ONE WITH HIGHER RESOLUTION (SEE METHODS)
      if(nrow(comparison_track) < nrow(this_bird)){
        
        int <- comparison_track
        comparison_track <- this_bird
        this_bird <- int
        
        
      }
      
      # For each point, calculate a NND
      track_NND <- c()
      for(j in 1:nrow(this_bird)){    
        
        distances <- dist_mat[which(rownames(dist_mat) == this_bird$pointID[j]),]
        distances_ct <- distances[which(colnames(dist_mat) %in% comparison_track$pointID)]
        
        track_NND[j] <- min(distances_ct, na.rm = T)
        
        
      }
      
      # Take a mean per track
      NND_per_track[x] <- mean(track_NND)
      
    }
  
    # To avoid pseudoreplication, take a mean per bird too
    track_output$NND[i] <- mean(track_NND)
    
  }
  
  
  
  print(tracks[i])
  
}


NND_sum <- track_output[which(is.na(track_output$NND) == F),] %>% group_by(bird) %>% summarise(mn_NND = mean(NND))

grand_mean <- mean(NND_sum$mn_NND)
grand_mean


# calculate the NND of the focal track to a randomly selected track from another bird
# essentially, this is the exact process above, but with randomly selected birds

n <- 1000    # and repeat that 1.000 times

rand_NND <- c()
for(k in 1:n){
  
  track_output <- data.frame(track_id = tracks,
                             bird = NA,
                             NND = NA,
                             n_tracks = NA)
  
  for(i in 1:length(tracks)){
    
    this_bird <- d_autumn[which(d_autumn$track_id == tracks[i]),]   #focal track of focal bird
    
    next_tracks <- unique(d_autumn$track_id[which(d_autumn$individual == this_bird$individual[1] &  # select the next tracks from the focal bird
                                                    d_autumn$track_id != tracks[i])])
    
    next_tracks <- sample(tracks, length(next_tracks), replace = T)
    
    track_output$bird[i] <- as.character(this_bird$individual[1])
    track_output$n_tracks[i] <- length(next_tracks)
    
    
    
    if(length(next_tracks) > 0){
      
      NND_per_track <- c()
      for(x in 1:length(next_tracks)){
        
        comparison_track <- d_autumn[which(d_autumn$track_id == next_tracks[x]),]
        
        if(nrow(comparison_track) < nrow(this_bird)){
          
          int <- comparison_track
          comparison_track <- this_bird
          this_bird <- int
          
          
        }
        
        track_NND <- c()
        for(j in 1:nrow(this_bird)){    
          
          distances <- dist_mat[which(rownames(dist_mat) == this_bird$pointID[j]),]
          distances_ct <- distances[which(colnames(dist_mat) %in% comparison_track$pointID)]
          
          track_NND[j] <- min(distances_ct, na.rm = T)
          
          
        }
        
        NND_per_track[x] <- mean(track_NND)
        
      }
      
      track_output$NND[i] <- mean(track_NND)
      
    }
    
    
    
    
  }
  
  
  NND_sum <- track_output[which(is.na(track_output$NND) == F),] %>% group_by(bird) %>% summarise(mn_NND = mean(NND))
  
  rand_NND[k] <- mean(NND_sum$mn_NND)
  
  
  print(rand_NND[k])
  
  
  
}

# make histogram of randomisation output
hist(rand_NND, xlab='mean NND (km)')
# add line for real results
abline(v = grand_mean, lwd = 2, lty = 2, col = 'red')

# then calculate p value (double it to account for 2-tailed test)
length(which(rand_NND < grand_mean))/length(rand_NND) * 2

# and save output (change file as appropriate)
save(rand_NND, file ="...")

load(file ="...")



######################################################## Within-individual randomization spring #########################################################################
###Spring migration

### SEE ABOVE FOR ANNOTATION ###

###calculate the NND of a track to other tracks of the same bird

d$pointID <- seq(1,nrow(d),1)

d_spring <- d[which(d$position_date_time > d$departure_wa &   #select all point between departure colony and arrival wintering area
                      d$position_date_time < d$arrival_date_colony),]

d_winter <- d[which(d$position_date_time > d$arrival_wa &
                      d$position_date_time < d$departure_wa),]


sprisum <- d_winter %>% group_by(track_id) %>% summarise(start_point = min(position_lat))

d_spring <- d_spring[which(d_spring$position_lat > max(sprisum$start_point)),] 


no_points <- d_spring %>% group_by(track_id) %>% summarise(n = n())

d_spring$no_points <- no_points$n[match(d_spring$track_id,no_points$track_id)]

dist_mat <- matrix(nrow = nrow(d_spring),   #build a matrix
                   ncol = nrow(d_spring))

colnames(dist_mat) <- d_spring$pointID
rownames(dist_mat) <- d_spring$pointID

#loop -> calculate the distance between all points

for(i in 1:nrow(dist_mat)){
  
  dists <- distHaversine(data.frame(d_spring$position_lon[i],
                                    d_spring$position_lat[i]),
                         data.frame(d_spring$position_lon,
                                    d_spring$position_lat))/1000
  dist_mat[i,] <- dists
  
  print(i)
  
}


tracks <- unique(d_spring$track_id)  #give a unique name to each point

track_output <- data.frame(track_id = tracks,
                           bird = NA,
                           NND = NA,
                           n_tracks = NA)

d_spring$NND <- NA
for(i in 1:length(tracks)){
  
  this_bird <- d_spring[which(d_spring$track_id == tracks[i]),]   #focal track of focal bird
  
  next_tracks <- unique(d_spring$track_id[which(d_spring$individual == this_bird$individual[1] &  # select the next tracks from the focal bird
                                                  d_spring$track_id != tracks[i])])
  
  track_output$bird[i] <- as.character(this_bird$individual[1])
  track_output$n_tracks[i] <- length(next_tracks)
  
  
  
  if(length(next_tracks) > 0){
    
    NND_per_track <- c()
    for(x in 1:length(next_tracks)){
      
      comparison_track <- d_spring[which(d_spring$track_id == next_tracks[x]),]
      
      if(nrow(comparison_track) < nrow(this_bird)){
        
        int <- comparison_track
        comparison_track <- this_bird
        this_bird <- int
        
        
      }
      
      track_NND <- c()
      for(j in 1:nrow(this_bird)){    
        
        distances <- dist_mat[which(rownames(dist_mat) == this_bird$pointID[j]),]
        distances_ct <- distances[which(colnames(dist_mat) %in% comparison_track$pointID)]
        
        track_NND[j] <- min(distances_ct, na.rm = T)
        
        
      }
      
      NND_per_track[x] <- mean(track_NND)
      
    }
    
    track_output$NND[i] <- mean(track_NND)
    
  }
  
  
  
  print(tracks[i])
  
}


NND_sum <- track_output[which(is.na(track_output$NND) == F),] %>% group_by(bird) %>% summarise(mn_NND = mean(NND))

grand_mean <- mean(NND_sum$mn_NND)
grand_mean


### calculate the NND of the focal track to a randomly selected track from another bird

n <- 1000    # and repeat that 1.000 times

rand_NND <- c()
for(k in 1:n){
  
  track_output <- data.frame(track_id = tracks,
                             bird = NA,
                             NND = NA,
                             n_tracks = NA)
  
  for(i in 1:length(tracks)){
    
    this_bird <- d_spring[which(d_spring$track_id == tracks[i]),]   #focal track of focal bird
    
    next_tracks <- unique(d_spring$track_id[which(d_spring$individual == this_bird$individual[1] &  # select the next tracks from the focal bird
                                                    d_spring$track_id != tracks[i])])
    
    next_tracks <- sample(tracks, length(next_tracks), replace = T)
    
    track_output$bird[i] <- as.character(this_bird$individual[1])
    track_output$n_tracks[i] <- length(next_tracks)
    
    
    
    if(length(next_tracks) > 0){
      
      NND_per_track <- c()
      for(x in 1:length(next_tracks)){
        
        comparison_track <- d_spring[which(d_spring$track_id == next_tracks[x]),]
        
        if(nrow(comparison_track) < nrow(this_bird)){
          
          int <- comparison_track
          comparison_track <- this_bird
          this_bird <- int
          
          
        }
        
        track_NND <- c()
        for(j in 1:nrow(this_bird)){    
          
          distances <- dist_mat[which(rownames(dist_mat) == this_bird$pointID[j]),]
          distances_ct <- distances[which(colnames(dist_mat) %in% comparison_track$pointID)]
          
          track_NND[j] <- min(distances_ct, na.rm = T)
          
          
        }
        
        NND_per_track[x] <- mean(track_NND)
        
      }
      
      track_output$NND[i] <- mean(track_NND)
      
    }
    
    
    
    
  }
  
  
  NND_sum <- track_output[which(is.na(track_output$NND) == F),] %>% group_by(bird) %>% summarise(mn_NND = mean(NND))
  
  rand_NND[k] <- mean(NND_sum$mn_NND)
  
  
  print(rand_NND[k])
  
  
  
}

hist(rand_NND, xlab='mean NND (km)')
abline(v = grand_mean, lwd = 2, lty = 2, col = 'red')

length(which(rand_NND < grand_mean))/length(rand_NND) * 2 

save(rand_NND, file ="...")

load(file ="...")


################################################# Temporal randomisation autumn ######################################################################


###Autumn migration


### SEE ABOVE FOR GENERAL ANNOTATION ###
# Randomisation-specific annotation is included below.

###calculate the NND of a track to other tracks of birds leaving at the same time window

d$pointID <- seq(1,nrow(d),1)

d_autumn <- d[which(d$position_date_time > d$departure_colony &   #select all point between departure colony and arrival wintering area
                      d$position_date_time < d$arrival_wa),]

d_winter <- d[which(d$position_date_time > d$arrival_wa &
                      d$position_date_time < d$departure_wa),]


autsum <- d_winter %>% group_by(track_id) %>% summarise(end_point = min(position_lat))

d_autumn <- d_autumn[which(d_autumn$position_lat > max(autsum$end_point)),]


no_points <- d_autumn %>% group_by(track_id) %>% summarise(n = n())

d_autumn$no_points <- no_points$n[match(d_autumn$track_id,no_points$track_id)]



dist_mat <- matrix(nrow = nrow(d_autumn),   #build a matrix
                   ncol = nrow(d_autumn))

colnames(dist_mat) <- d_autumn$pointID
rownames(dist_mat) <- d_autumn$pointID

#loop -> calculate the distance between all points

for(i in 1:nrow(dist_mat)){
  
  dists <- distHaversine(data.frame(d_autumn$position_lon[i],
                                    d_autumn$position_lat[i]),
                         data.frame(d_autumn$position_lon,
                                    d_autumn$position_lat))/1000
  dist_mat[i,] <- dists
  
  print(i)
  
}


tracks <- unique(d_autumn$track_id)  #give a unique name to each point

d_autumn$NND <- NA

track_output <- data.frame(track_id = tracks,
                           bird = NA,
                           NND = NA,
                           n_tracks = NA)

for(i in 1:length(tracks)){
  
  this_bird <- d_autumn[which(d_autumn$track_id == tracks[i]),]   #focal track of focal bird
  
  # Select tracks leaving around the same time as the focal track
  # Calculate NND between focal bird and each comparison track
  next_tracks <- unique(d_autumn$track_id[which(abs(d_autumn$departure_colony-this_bird$departure_colony[1]) < 3*24*60*60 &
                                                d_autumn$track_id != tracks[i])]) # each individual is compared to an individual of the same time window (+/- 3 days) (hyp 2)
  
  track_output$bird[i] <- as.character(this_bird$individual[1])
  track_output$n_tracks[i] <- length(next_tracks)
  
  if(length(next_tracks) > 0){
    
    NND_per_track <- c()
    for(x in 1:length(next_tracks)){
      
      comparison_track <- d_autumn[which(d_autumn$track_id == next_tracks[x]),]
      
      if(nrow(comparison_track) < nrow(this_bird)){
        
        int <- comparison_track
        comparison_track <- this_bird
        this_bird <- int
        
        
      }
      
      track_NND <- c()
      for(j in 1:nrow(this_bird)){    
        
        distances <- dist_mat[which(rownames(dist_mat) == this_bird$pointID[j]),]
        distances_ct <- distances[which(colnames(dist_mat) %in% comparison_track$pointID)]
        
        track_NND[j] <- min(distances_ct, na.rm = T)
        
        
      }
      
      NND_per_track[x] <- mean(track_NND)
      
    }
    
    track_output$NND[i] <- mean(track_NND)
    
  }
  
  
  
  #print(tracks[i])
    
}


#hist(d_autumn$NND)

NND_sum <- track_output[which(is.na(track_output$NND) == F),] %>% group_by(bird) %>% summarise(mn_NND = mean(NND))

grand_mean <- mean(NND_sum$mn_NND)
grand_mean

### calculate the NND of the focal track to a randomly selected track from another bird independent of the time window

n <- 1000    # and repeat that 1.000 times

rand_NND <- c()
for(k in 1:n){
  
  d_autumn$NND <- NA
  
  track_output <- data.frame(track_id = tracks,
                             bird = NA,
                             NND = NA,
                             n_tracks = NA)
  
  for(i in 1:length(tracks)){
    
    this_bird <- d_autumn[which(d_autumn$track_id == tracks[i]),]   #focal track of focal bird
    
    track_output$bird[i] <- as.character(this_bird$individual[1])
    
    # Select tracks leaving around the same time as the focal track
    next_tracks <- unique(d_autumn$track_id[which(abs(d_autumn$departure_colony-this_bird$departure_colony[1]) < 3*24*60*60 &
                                                    d_autumn$track_id != tracks[i])]) 
    # Based on the number of real comparisons, select random comparisons. 
    next_tracks <- sample(tracks, length(next_tracks), replace = T)
    
    track_output$n_tracks[i] <- length(next_tracks)
    
    if(length(next_tracks) > 0){
      
      NND_per_track <- c()
      for(x in 1:length(next_tracks)){
        
        comparison_track <- d_autumn[which(d_autumn$track_id == next_tracks[x]),]  
        
        if(nrow(comparison_track) < nrow(this_bird)){
          
          int <- comparison_track
          comparison_track <- this_bird
          this_bird <- int
          
          
        }
        
        track_NND <- c()
        for(j in 1:nrow(this_bird)){    
          
          distances <- dist_mat[which(rownames(dist_mat) == this_bird$pointID[j]),]
          distances_ct <- distances[which(colnames(dist_mat) %in% comparison_track$pointID)]
          
          track_NND[j] <- min(distances_ct, na.rm = T)
          
          
        }
        
        NND_per_track[x] <- mean(track_NND)
        
      }
      
      track_output$NND[i] <- mean(track_NND)
      
    }
    
    
    
    #print(tracks[i])
    
  }
  
  
  NND_sum <- track_output[which(is.na(track_output$NND) == F),] %>% group_by(bird) %>% summarise(mn_NND = mean(NND))
  
  rand_NND[k] <- mean(NND_sum$mn_NND)
  
  
  print(rand_NND[k])
  
  
  
}

hist(rand_NND, ylim=c(0,250), xlab='mean NND (km)')
abline(v = grand_mean, lwd = 2, lty = 2, col = 'red')

length(which(rand_NND < grand_mean))/length(rand_NND) *2

save(rand_NND, file ="...")

load(file ="...")

################################################# Temporal randomisation spring ######################################################################


### SEE ABOVE FOR GENERAL ANNOTATION ###
# Randomisation-specific annotation is included below.

###spring migration

###calculate the NND of a track to other tracks of birds leaving at the same time window

d$pointID <- seq(1,nrow(d),1)

d_spring<- d[which(d$position_date_time > d$departure_wa &   #select all point between departure colony and arrival wintering area
                      d$position_date_time < d$arrival_date_colony),]

d_winter <- d[which(d$position_date_time > d$arrival_wa &
                      d$position_date_time < d$departure_wa),]


springsum <- d_winter %>% group_by(track_id) %>% summarise(start_point = min(position_lat))

d_spring <- d_spring[which(d_spring$position_lat > max(springsum$start_point)),]


no_points <- d_spring %>% group_by(track_id) %>% summarise(n = n())

d_spring$no_points <- no_points$n[match(d_spring$track_id,no_points$track_id)]


dist_mat <- matrix(nrow = nrow(d_spring),   #build a matrix
                   ncol = nrow(d_spring))

colnames(dist_mat) <- d_spring$pointID
rownames(dist_mat) <- d_spring$pointID

#loop -> calculate the distance between all points

for(i in 1:nrow(dist_mat)){
  
  dists <- distHaversine(data.frame(d_spring$position_lon[i],
                                    d_spring$position_lat[i]),
                         data.frame(d_spring$position_lon,
                                    d_spring$position_lat))/1000
  dist_mat[i,] <- dists
  
  print(i)
  
}


tracks <- unique(d_spring$track_id)  #give a unique name to each point

d_spring$NND <- NA

track_output <- data.frame(track_id = tracks,
                           bird = NA,
                           NND = NA,
                           n_tracks = NA)

for(i in 1:length(tracks)){
  
  this_bird <- d_spring[which(d_spring$track_id == tracks[i]),]   #focal track of focal bird
  
  

  next_tracks <- unique(d_spring$track_id[which(abs(d_spring$departure_colony-this_bird$departure_colony[1]) < 3*24*60*60 &
                                                  d_spring$track_id != tracks[i])]) 
  
  track_output$bird[i] <- as.character(this_bird$individual[1])
  track_output$n_tracks[i] <- length(next_tracks)
  
  if(length(next_tracks) > 0){
    
    NND_per_track <- c()
    for(x in 1:length(next_tracks)){
      
      comparison_track <- d_spring[which(d_spring$track_id == next_tracks[x]),]
      
      if(nrow(comparison_track) < nrow(this_bird)){
        
        int <- comparison_track
        comparison_track <- this_bird
        this_bird <- int
        
        
      }
      
      track_NND <- c()
      for(j in 1:nrow(this_bird)){    
        
        distances <- dist_mat[which(rownames(dist_mat) == this_bird$pointID[j]),]
        distances_ct <- distances[which(colnames(dist_mat) %in% comparison_track$pointID)]
        
        track_NND[j] <- min(distances_ct, na.rm = T)
        
        
      }
      
      NND_per_track[x] <- mean(track_NND)
      
    }
    
    track_output$NND[i] <- mean(track_NND)
    
  }
  
  
  
  print(tracks[i])
  
}


#hist(d_spring$NND)

NND_sum <- track_output[which(is.na(track_output$NND) == F),] %>% group_by(bird) %>% summarise(mn_NND = mean(NND))

grand_mean <- mean(NND_sum$mn_NND)
grand_mean

### calculate the NND of the focal track to a randomly selected track from another bird independent of the time window

n <- 1000    # and repeat that 1.000 times

rand_NND <- c()
for(k in 1:n){
  
  d_spring$NND <- NA
  
  track_output <- data.frame(track_id = tracks,
                             bird = NA,
                             NND = NA,
                             n_tracks = NA)
  
  for(i in 1:length(tracks)){
    
    this_bird <- d_spring[which(d_spring$track_id == tracks[i]),]   #focal track of focal bird
    
    track_output$bird[i] <- as.character(this_bird$individual[1])
    

    next_tracks <- unique(d_spring$track_id[which(abs(d_spring$departure_colony-this_bird$departure_colony[1]) < 3*24*60*60 &
                                                    d_spring$track_id != tracks[i])]) 
    
    next_tracks <- sample(tracks, length(next_tracks), replace = T)
    
    track_output$n_tracks[i] <- length(next_tracks)
    
    if(length(next_tracks) > 0){
      
      NND_per_track <- c()
      for(x in 1:length(next_tracks)){
        
        comparison_track <- d_spring[which(d_spring$track_id == next_tracks[x]),]  
        
        if(nrow(comparison_track) < nrow(this_bird)){
          
          int <- comparison_track
          comparison_track <- this_bird
          this_bird <- int
          
          
        }
        
        track_NND <- c()
        for(j in 1:nrow(this_bird)){    
          
          distances <- dist_mat[which(rownames(dist_mat) == this_bird$pointID[j]),]
          distances_ct <- distances[which(colnames(dist_mat) %in% comparison_track$pointID)]
          
          track_NND[j] <- min(distances_ct, na.rm = T)
          
          
        }
        
        NND_per_track[x] <- mean(track_NND)
        
      }
      
      track_output$NND[i] <- mean(track_NND)
      
    }
    
    
    
    #print(tracks[i])
    
  }
  
  
  NND_sum <- track_output[which(is.na(track_output$NND) == F),] %>% group_by(bird) %>% summarise(mn_NND = mean(NND))
  
  rand_NND[k] <- mean(NND_sum$mn_NND)
  
  
  print(rand_NND[k])
  
  
  
}

hist(rand_NND, xlab='mean NND (km)')
abline(v = grand_mean, lwd = 2, lty = 2, col = 'red')

length(which(rand_NND < grand_mean))/length(rand_NND) *2 

save(rand_NND, file ="...")

load(file ="...")

