require(lubridate);require(raster);require(geosphere);require(circular); require(dplyr);require(rworldxtra)

# Function for calculating difference in bearing between 2 angles

diff_bear <- function(x,y){  #calculates the angle between the bird and wind
  z <- ifelse(x > y,
              ifelse(x-y < 180,
                     ifelse(x > y,
                            -(x - y),
                            (x - y)),
                     ifelse(x > y,
                            360 - (x - y),
                            -(360 - (x - y)))),
              ifelse(y - x < 180,
                     ifelse(y > x,
                            (y - x),
                            -(y - x)),
                     ifelse(y > x,
                            -(360 - (y - x)),
                            (360 - (y - x)))))
  return(z)
}


# load position data
comics <- read.csv('...postitions_per_day.csv', sep = ';')

# load wind data as raster stacks of u and v wind
uwind <- stack('...uwind_terns.grib')
vwind <- stack('...vwind_terns.grib')

# date of each layer is implicit from stack position
# so we can label each layer as a different date like this:
layer_names <- seq(as.POSIXct('2016-01-01'), as.POSIXct('2021-12-31'), by = 'days')
layer_years <- year(layer_names)
layer_jdate <- yday(layer_names)
layer_match <- paste(layer_years,layer_jdate)

# format dates for use later
comics$departure_colony <- as.POSIXct(comics$departure_colony, format = '%d.%m.%Y', tz = 'CET')
comics$departure_wa <- as.POSIXct(comics$departure_wa, format = '%d.%m.%Y', tz = 'CET')
comics$arrival_date_colony <- as.POSIXct(comics$arrival_date_colony, format = '%d.%m.%Y', tz = 'CET')
comics$arrival_wa <- as.POSIXct(comics$arrival_wa, format = '%d.%m.%Y', tz = 'CET')
comics$position_date_time <- as.POSIXct(comics$position_date_time, format = '%d.%m.%Y %H:%M', tz = 'CET')

# we have added this toggle to make it easier to switch between seasons. 
aut <-  F  #false = spring, true = autumn


# subset for points that are above the maximum wintering latitude 
comics <- comics[which(comics$position_date_time > comics$departure_wa &
                                                  comics$position_lat > 21.15),]

if(aut == F){
  
  #subset for spring points
  comics <- comics[which(comics$position_date_time > comics$departure_wa &
                           comics$position_date_time < comics$arrival_date_colony),]
  
  
} else{
  
  # subset for autumn points
  comics <- comics[which(comics$position_date_time > comics$departure_colony &
                           comics$position_date_time < comics$arrival_wa),]
  
}

### CALCULATING THE BEARING BETWEEN COMMON TERN GLS POINTS ###
comics$bearing <- NA
for(i in 1:(nrow(comics)-1)){
  
  
  
  if(comics$track_id[i] != comics$track_id[i + 1]){
    
    print('First point, do not calculate bearing!')  
    
  } else{
    # bearing function from geosphere
    comics$bearing[i] <- (bearing(data.frame(comics$position_lon[i],
                                             comics$position_lat[i]),
                                  data.frame(comics$position_lon[i+1],
                                             comics$position_lat[i+1])))
    
    comics$dist_next[i] <- (distHaversine(data.frame(comics$position_lon[i],
                                                     comics$position_lat[i]),
                                          data.frame(comics$position_lon[i+1],
                                                     comics$position_lat[i+1])))/1000
    
  }
  
  
  
}

comics$bearing <- ifelse(comics$bearing < 0,    #we changes the scale from 0 to 360
                         360 + comics$bearing,
                         comics$bearing)

comics <- comics[which(comics$dist_next > 100),] #removes "stopovers"

comics$pointID <- seq(1,nrow(comics),1) #making a matrix -> each data point a unique ID

departure_dates <- unique(comics$departure_colony) #every potential date a bird left on

# Essentially, all we want to do is calculate the wind a bird would 
# experience at every point along its trajectory if it left on every conceivable departure date.
# The way we do this is using a matrix, where column = departure date and row = GLS fix identity:

# ---- departure dates ---->
# |
# |
# | GLS fixes 
# |
# |

dat_mat <- matrix(ncol = length(departure_dates),
                  nrow = nrow(comics))


tracks <- unique(comics$track_id)




for(i in 1:length(tracks)){
  
  # select a track from our dataset
  bird <- comics[which(comics$track_id == tracks[i]),]
  
  for(j in 1:length(departure_dates)){
    
    sim_bird <- bird
    
    # For each possible departure date, calculate the time at which the bird *would have* reached this point 
    # had it left on *that* date. 
    sim_bird$position_date_time_relative <- as.numeric(difftime(sim_bird$position_date_time, sim_bird$departure_colony, unit = 'secs'))

    # Store that date in the matrix we created above
    # later we can use this date to reference our raster stack and extract relevant wind information.
    dat_mat[sim_bird$pointID,j] <- paste(year(as.POSIXct(departure_dates[j] + sim_bird$position_date_time_relative,origin = '1970')), yday(as.POSIXct(departure_dates[j] + sim_bird$position_date_time_relative,origin = '1970')))
    
  }
  
  print(i)
  
  
}

# For the headwind component, we make a separate matrix.
# This calculates the wind the bird *would have* encountered had it reached the point in question
# having started migration on the date in question.

hw_mat <- matrix(ncol = length(departure_dates), #headwind
                 nrow = nrow(comics))

ws_mat <- matrix(ncol = length(departure_dates), #wind speed
                 nrow = nrow(comics))


date_combos <- unique(as.vector(as.matrix(dat_mat[-1]))) #every possible date we sample the wind for


for(i in 1:length(date_combos)){
  
  # select points from the day in question
  pointsz <- data.frame(which(dat_mat == date_combos[i], arr.ind=TRUE))
  
  # select the layers for the day in question
  uw <- subset(uwind,which(layer_match == date_combos[i]))
  vw <- subset(vwind,which(layer_match == date_combos[i]))
  
  # extract the uwind and vwind for the points in question
  positions <- comics[match(pointsz$row,comics$pointID),]
  positions$u_wind  <- extract(uw, SpatialPoints(data.frame(positions$position_lon,positions$position_lat)))
  positions$v_wind  <- extract(vw, SpatialPoints(data.frame(positions$position_lon,positions$position_lat)))
  
  positions$wind_speed <- sqrt(positions$u_wind^2 + positions$v_wind^2)
  positions$wind_direction <- (90 - (atan2(positions$v_wind,positions$u_wind) * 180/pi)) %% 360
  
  # use function to give a bearing between bird direction and wind direction
  positions$wind_to_course <- diff_bear(positions$wind_direction,
                                        positions$bearing)
  
  # calculate headwind component (need to use abs before the sine!)
  
  
  positions$headwind <- ifelse(abs(positions$wind_to_course) < 90,
                                    cos(rad(abs(positions$wind_to_course)))*positions$wind_speed,
                                    cos(rad((180-abs(positions$wind_to_course))))*positions$wind_speed)
  
  
  positions$headwind <- ifelse(abs(positions$wind_to_course) < 90,
                                    -positions$headwind,
                                    positions$headwind)
  
  
  # save headwind to the relevant point in the matrix
  hw_mat[which(dat_mat == date_combos[i], arr.ind=TRUE)] <- positions$headwind
  ws_mat[which(dat_mat == date_combos[i], arr.ind=TRUE)] <- positions$wind_speed
  print(i)
  
}

# exhaustive list of departure dates (with frequency distribution maintained)
dep_dates <- comics %>% group_by(track_id) %>% summarise(departure_colony = departure_colony[1]) 


tracks <- unique(comics$track_id)

n <- 10000 #number of randomisation
rand_hw <- c()
for(i in 1:n){
  
  simbird <- comics
  
  # select a random departure date for each bird
  rand_dep_date <- sample(dep_dates$departure_colony,length(tracks), replace = T)
  
  # So essentially all this does is replicate the departure date for each bird by the number of fixes that are in the track in question.
  runout <- rle(comics$track_id)
  fix_times <- rep(rand_dep_date, runout$lengths)
  
  # These are the times that the GLS fixes *would* have been taken *if* the bird had left on the day in question.
  # *NOTE* the relative difference in time between each fix and its departure date is accounted for above.
  # This means we just replicate each departure date by the number of fixes and the time difference is implicit.
  fix_times_jdate <- paste(year(as.POSIXct(fix_times)), yday(as.POSIXct(fix_times)))
  
  # Extract randomised headwind from matrix
  # Calculate mean per track
  # Calculate mean per individual
  simbird$hw_rand <- hw_mat[cbind(simbird$pointID,match(fix_times,departure_dates))]/ws_mat[cbind(simbird$pointID,match(fix_times,departure_dates))]
  op1 <- simbird %>% group_by(track_id, individual) %>% summarise(hw_rand = mean(hw_rand, na.rm = T))
  op2 <- op1 %>% group_by(individual) %>% summarise(hw_rand = mean(hw_rand, na.rm = T))
  
  # calculate a grand mean from all the individuals
  rand_hw[i] <- mean(op2$hw_rand)
  
  print(i)
  
}




fix_times <- comics$departure_colony

# calculates the real head wind proportion for the real track
fix_times_jdate <- paste(year(as.POSIXct(fix_times)), yday(as.POSIXct(fix_times)))
comics$hw <- hw_mat[cbind(comics$pointID,match(fix_times,departure_dates))]/ws_mat[cbind(comics$pointID,match(fix_times,departure_dates))]
comics$headwind <- hw_mat[cbind(comics$pointID,match(fix_times,departure_dates))]/ws_mat[cbind(comics$pointID,match(fix_times,departure_dates))]

# Then calculate a mean per track...
# ...and a mean per individual...
# ...and a grand mean!!!
op1 <- comics %>% group_by(track_id, individual) %>% summarise(hw = mean(hw, na.rm = T))
op2 <- op1 %>% group_by(individual) %>% summarise(hw = mean(hw, na.rm = T))

true_hw <- mean(op2$hw)

hist(rand_hw, 100, xlab='mean headwind proportion', ylim = c(0,355))
abline(v = true_hw, col = 'red', lwd = 2, lty = 2)


length(which(rand_hw < true_hw))/n  * 2 

# save output
save(rand_hw, file ="...")
save(true_hw, file = "...")

#load(file ="...")

