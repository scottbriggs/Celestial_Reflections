# Calculate the position of the nutation angles (Longitude, Obliquity)
# using the JPL DE 441 Ephemeris data

nutationAngles <- function(jd)
{
  # Read data file
  df_nut_raw <- read_parquet(here("data", "processed", "nutation_441.parquet"))
  
  # Subset the data for the julian day number of interest
  df_jd <- df_nut_raw[df_nut_raw$Julian_Day_Start <= jd & 
                        df_nut_raw$Julian_Day_End > jd,]
  
  # Calculate the length of the subinterval
  length_of_subinterval <- DE441DAYSPERBLOCK / DE441NUMSUBINTNUTATION
  
  subinterval <- floor(as.integer(jd - df_jd[1,1]) / length_of_subinterval)
  
  # Add 1 to get the right subinterval. The above algorithm assumes the
  # subinterval begins with 0, but the subinterval begins with 1 in the 
  # database
  subinterval <- subinterval + 1
  
  # Subset data for the interval of interest
  df_nutation <- df_jd[df_jd$INTERVAL == subinterval,]
  df_nutation <- subset(df_nutation, 
                     select = -c(Julian_Day_Start, Julian_Day_End, INTERVAL))
  
  # Normalize the Julian Day
  valid_start <- df_jd[1,1] + ((subinterval - 1) * length_of_subinterval)
  valid_end <- valid_start + length_of_subinterval
  temp <- jd - valid_start
  x <- (temp / length_of_subinterval * 2.0) - 1.0
  
  # Calculate the Chebyshev polynomials for the nutations
  chebyshev_nut <- data.frame(matrix(0.0, nrow=DE441NUMCOEFFNUTATION, ncol=2))
  chebyshev_nut[1,1] <- 1.0
  chebyshev_nut[2,1] <- x

  for (i in seq(from = 3, to = DE441NUMCOEFFNUTATION, by = 1)){
    chebyshev_nut[i,1] <- (2 * x * chebyshev_nut[i-1,1]) - chebyshev_nut[i-2,1]
  }
  
  # Calculate the nutation (longitude, obliquity) in radians
  nut_ang <- data.frame(matrix(0.0, nrow=1, ncol=2))
  v <- 0
  for (v in seq(from = DE441NUMCOEFFNUTATION, to = 1, by = -1)){
    nut_ang[1,1] <- nut_ang[1,1] + (chebyshev_nut[v,1] * df_nutation[1,v])
    nut_ang[1,2] <- nut_ang[1,2] + (chebyshev_nut[v,1] * df_nutation[1,v+DE441NUMCOEFFNUTATION])
  }
  
  colnames(nut_ang) <- c('Nutation in Longitude', 'Nutation in Obliquity')

  # Return the data
  return(nut_ang)
}