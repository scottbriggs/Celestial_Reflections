# Calculate the position of Saturn
# using the JPL DE 441 Ephemeris data

positionSaturnSSB <- function(jd)
{
  # Read data file
  df_sat_raw <- read_parquet(here("data", "processed", "saturn_441.parquet"))
  
  # Subset the data for the julian day number of interest
  df_jd <- df_sat_raw[df_sat_raw$Julian_Day_Start <= jd & 
                        df_sat_raw$Julian_Day_End > jd,]
  
  # Calculate the length of the subinterval
  length_of_subinterval <- DE441DAYSPERBLOCK / DE441NUMSUBINTSATURN
  
  subinterval <- floor(as.integer(jd - df_jd[1,1]) / length_of_subinterval)
  
  # Add 1 to get the right subinterval. The above algorithm assumes the
  # subinterval begins with 0, but the subinterval begins with 1 in the 
  # database
  subinterval <- subinterval + 1
  
  # Subset data for the interval of interest
  df_saturn <- df_jd[df_jd$INTERVAL == subinterval,]
  df_saturn <- subset(df_saturn, 
                     select = -c(Julian_Day_Start, Julian_Day_End, INTERVAL))
  
  # Normalize the Julian Day
  valid_start <- df_jd[1,1] + ((subinterval - 1) * length_of_subinterval)
  valid_end <- valid_start + length_of_subinterval
  temp <- jd - valid_start
  x <- (temp / length_of_subinterval * 2.0) - 1.0
  
  # Calculate the Chebyshev polynomials for position and velocity. The velocity
  # is the first derivative of the position polynomial
  chebyshev <- data.frame(matrix(0.0, nrow=DE441NUMCOEFFSATURN, ncol=2))
  chebyshev[1,1] <- 1.0
  chebyshev[2,1] <- x
  chebyshev[1,2] <- 0.0
  chebyshev[2,2] <- 1.0
  
  # Calculate the position coefficients
  for (i in seq(from = 3, to = DE441NUMCOEFFSATURN, by = 1)){
    chebyshev[i,1] <- (2 * x * chebyshev[i-1,1]) - chebyshev[i-2,1]
  }
  
  # Calculate the velocity coefficients
  for (i in seq(from = 3, to = DE441NUMCOEFFSATURN, by = 1)){
    chebyshev[i,2] <- (2 * x * chebyshev[i-1,2]) - chebyshev[i-2,2] + (2 * chebyshev[i-1,1])
  }
  
  # Calculate the position in kilometers and the velocity in kilometers/sec
  pos_vel <- data.frame(matrix(0.0, nrow=3, ncol=2))
  v <- 0
  for (v in seq(from = DE441NUMCOEFFSATURN, to = 1, by = -1)){
    pos_vel[1,1] <- pos_vel[1,1] + (chebyshev[v,1] * df_saturn[1,v])
    pos_vel[2,1] <- pos_vel[2,1] + (chebyshev[v,1] * df_saturn[1,v+DE441NUMCOEFFSATURN])
    pos_vel[3,1] <- pos_vel[3,1] + (chebyshev[v,1] * df_saturn[1,v+2*DE441NUMCOEFFSATURN])
    
    pos_vel[1,2] <- pos_vel[1,2] + (chebyshev[v,2] * df_saturn[1,v])
    pos_vel[2,2] <- pos_vel[2,2] + (chebyshev[v,2] * df_saturn[1,v+DE441NUMCOEFFSATURN])
    pos_vel[3,2] <- pos_vel[3,2] + (chebyshev[v,2] * df_saturn[1,v+2*DE441NUMCOEFFSATURN])
  }
  
  # Scale the velocity
  scale_value <- 2 * DE441NUMSUBINTSATURN / DE441DAYSPERBLOCK
  pos_vel[1,2] = pos_vel[1,2] * scale_value
  pos_vel[2,2] = pos_vel[2,2] * scale_value
  pos_vel[3,2] = pos_vel[3,2] * scale_value
  
  colnames(pos_vel) <- c('Position Vector', 'Velocity Vector')
  
  # Return the data
  return(pos_vel)
}