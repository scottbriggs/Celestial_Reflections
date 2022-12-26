# Calculate the position of the nutation angles (Longitude, Obliquity)
# using the JPL DE 441 Ephemeris data

# The nutation in longitude and obliquity angles are in units of radians
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
  
  colnames(nut_ang) <- c('Nut_Long', 'Nut_Oblq')

  return(nut_ang)
}

# Create the nutation matrix
nutation_matrix <- function(jd)
{
  nut_angles <- nutationAngles(jd)
  obliq <- obliquity(jd, nut_angles)
  
  # Calculate the nutation matrix elements
  cos_long <- cos(nut_angles[[1]])
  sin_long <- sin(nut_angles[[1]])
  sin_mean_obliquity <- sin(obliq[[1]])
  cos_mean_obliquity <- cos(obliq[[1]])
  sin_true_obliquity <- sin(obliq[[2]])
  cos_true_obliquity <- cos(obliq[[2]])
  
  nut_mat <- matrix(0.0, nrow=3, ncol=3)
  nut_mat[1,1] <- cos_long
  nut_mat[1,2] <- -sin_long * cos_mean_obliquity
  nut_mat[1,3] <- -sin_long * sin_mean_obliquity
  nut_mat[2,1] <- sin_long * cos_true_obliquity
  nut_mat[2,2] <- cos_long * cos_true_obliquity * cos_mean_obliquity + 
    sin_true_obliquity * sin_mean_obliquity
  nut_mat[2,3] <- cos_long * cos_true_obliquity * sin_mean_obliquity -
    sin_true_obliquity * cos_mean_obliquity
  nut_mat[3,1] <- sin_long * sin_true_obliquity
  nut_mat[3,2] <- cos_long * sin_true_obliquity * cos_mean_obliquity -
    cos_true_obliquity * sin_mean_obliquity
  nut_mat[3,3] <- cos_long * sin_true_obliquity * sin_mean_obliquity +
    cos_true_obliquity * cos_mean_obliquity
  
  return(nut_mat)
}