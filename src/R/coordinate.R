
# Functions to calculate polar coordinates from position vectors
# and perform coordinate transformations

# Convert a rectangular position vector to polar coordinates
# x - x component of vector
# y - y component of vector
# z - z component of vector
# Returns the polar distance r
# phi the latitudinal component in radians
# theta the longitudinal component in radians
rectToPolar <- function(x, y, z){
  rho_sqr <- x * x + y * y
  m_r <- sqrt(rho_sqr + z * z)
  
  # Calculate the longitudinal component of the vector
  m_theta <- 0.0
  if (x == 0 & y == 0) {
    m_theta <- 0
  }
  else {
    m_theta <- atan2(y, x)
  }
  if (m_theta < 0) {
    m_theta <- m_theta + PI2
  }
  
  # Calculate the latitudinal component of the vector
  m_phi <- 0
  rho <- sqrt(rho_sqr)
  if (z == 0 & rho == 0) {
    m_phi = 0
  }
  else {
    m_phi <- atan2(z, rho)
  }
  
  z <- c(m_r, m_phi, m_theta)
  
  return (z)
}

# Returns a rotation matrix based on the axis (x, y, or z) and the angle phi
# The angle phi is expected to be in radians
rotationMatrix <- function(axis, phi)
{
  mat <- matrix(0.0, nrow = 3, ncol = 3)
  cosphi <- cos(phi)
  sinphi <- sin(phi)
  
  if (axis == 1) {
    mat[1,1] <- 1.0
    mat[2,2] <- cosphi
    mat[3,3] <- cosphi
    mat[2,3] <- sinphi
    mat[3,2] <- -sinphi
  } else if (axis == 2) {
    mat[1,1] <- cosphi
    mat[1,3] <- -sinphi
    mat[2,2] <- 1.0
    mat[3,1] <- sinphi
    mat[3,3] <- cosphi
  } else if (axis == 3){
    mat[1,1] <- cosphi
    mat[2,2] <- cosphi
    mat[3,3] <- 1.0
    mat[2,1] <- -sinphi
    mat[1,2] <- sinphi
  } else {print("Axis is wrong value")}
  
  return (mat)
}

# Convert equatorial coordinates - right ascension and declination to ecliptic
# latitude and longitude
# ra - right ascension in decimal hours
# dec - declination in decimal degrees
# obliquity in decimal degrees
# eclipitic latitude in decimal degrees
# ecliptic longitude in decimal degrees
equatorialToEcliptic <- function(ra, dec, obliquity)
{
  sinRA <- sin(ra*15*DEG2RAD)
  cosRA <- cos(ra*15*DEG2RAD)
  sinObliq <- sin(obliquity*DEG2RAD)
  cosObliq <- cos(obliquity*DEG2RAD)
  sinDec <- sin(dec*DEG2RAD)
  cosDec <- cos(dec*DEG2RAD)
  tanDec <- tan(dec*DEG2RAD)
  
  sinBeta <- sinDec * cosObliq - cosDec * sinObliq * sinRA
  beta <- asin(sinBeta) * RAD2DEG
  
  numerator <- sinRA * cosObliq + tanDec * sinObliq
  lambda <- atan(numerator/cosRA)*RAD2DEG
  
  if (numerator > 0 & cosRA < 0) {
    lambda <- lambda + 180
  } else if (numerator < 0 & cosRA > 0) {
    lambda <- lambda + 360
  } else if (numerator < 0 & cosRA < 0) {
    lambda <- lambda + 180
  }
  
  z <- c(beta, lambda)
  names(z) <- c("Ecliptic Latitude", "Ecliptic Longitude")
  
  return (z)
}

# Convert ecliptical coordinates to equatorial coordinates
# ecliptic latitude - beta in decimal degrees, [-90 - + 90] degrees
# ecliptic longitude - lambda in decimal degrees, [0 360] degrees
# obliquity in decimal degrees
eclipticToEquatorial <- function(beta, lambda, obliquity)
{
  sinBeta <- sin(beta*DEG2RAD)
  cosBeta <- cos(beta*DEG2RAD)
  tanBeta <- tan(beta*DEG2RAD)
  sinObliq <- sin(obliquity*DEG2RAD)
  cosObilq <- cos(obliquity*DEG2RAD)
  sinLambda <- sin(lambda*DEG2RAD)
  cosLambda <- cos(lambda*DEG2RAD)
  
  sinDec <- sinBeta * cosObilq + cosBeta * sinObliq * sinLambda
  dec <- asin(sinDec) * RAD2DEG
  
  numerator <- sinLambda * cosObilq - tanBeta * sinObliq
  ra <- atan(numerator/cosLambda) * RAD2DEG
  
  if (numerator > 0 & cosLambda < 0) {
    ra <- ra + 180
  } else if (numerator < 0 & cosLambda > 0) {
    ra <- ra + 360
  } else if (numerator < 0 & cosLambda < 0) {
    ra <- ra + 180
  }
  
  z <- c(ra/15, dec)
  names(z) <- c("RA", "DEC")
  
  return (z)
}

# Convert equatorial coordinates to horizon coordinates
# azimuth in decimal degrees, [0 - 360] degrees
# altitude in decimal degrees, [-90 - +90] degrees
# hour angle in decimal hours, [0 - 24] hours
# declination in decimal degrees, [-90 - +90] degrees
# observer latitude in decimal degrees, [-90, +90] degrees
equatorialToHorizon <- function(hourAngle, dec, obs_lat)
{
  cosH <- cos(hourAngle*15*DEG2RAD)
  cosLat <- cos(obs_lat*DEG2RAD)
  sinLat <- sin(obs_lat*DEG2RAD)
  cosDec <- cos(dec*DEG2RAD)
  sinDec <- sin(dec*DEG2RAD)

  sinAlt <- sinDec * sinLat + cosDec * cosLat * cosH
  h <- asin(sinAlt) * RAD2DEG
  
  cosAz <- (sinDec - sinLat * sinAlt) / (cosLat * cos(h*DEG2RAD))
  az <- acos(cosAz) * RAD2DEG
  
  sinH <- sin(hourAngle*15*DEG2RAD)
  
  if (sinH > 0) {
    az <- 360 - az
  }
  
  z <- c(az, h)
  names(z) <- c("Azimuth", "Altitude")
  
  return (z)
}

# Convert horizon coordinates to equatorial coordinates
# azimuth in decimal degrees, [0 - 360] degrees
# altitude in decimal degrees, [-90 - +90] degrees
# hour angle in decimal hours, [0 - 24] hours
# declination in decimal degrees, [-90 - +90] degrees
horizonToEquatorial <- function(azimuth, altitude, obs_lat)
{
  sinAlt <- sin(altitude*DEG2RAD)
  cosAlt <- cos(altitude*DEG2RAD)
  cosAz <- cos(azimuth*DEG2RAD)
  sinLat <- sin(obs_lat*DEG2RAD)
  cosLat <- cos(obs_lat*DEG2RAD)
  
 sinDec <- sinAlt * sinLat + cosAlt * cosLat * cosAz
 dec <- asin(sinDec) * RAD2DEG
 
 numerator <- sinAlt - sinLat * sinDec
 cosH <- numerator / (cosLat * cos(dec*DEG2RAD))
 H <- acos(cosH) * RAD2DEG
 
 # Convert the hour angle from a range of [0 - 180] degrees to a range of
 # [0 - 360] degrees or [0 - 24] hours
 sinAz <- sin(azimuth*DEG2RAD)
 if (sinAz > 0) {
   H <- 360 - H
 }
 
 # Convert hour angle in degrees to hours
 H <- H /15
 
 z <- c(H, dec)
 names(z) <- c("Hour Angle", "Declination")
  
  return(z)
}

findHorizonCoordSun <- function(year, month, day, hour, tz, dst, obsLat, obsLong)
{
  # Convert local civil time to universal time
  hr <- hour
  dy <- day
  
  # Adjust for daylight savings time
  if (dst == TRUE) {
    hr <- hr - 1
  }
   # Adjust for the time zone
  hr <- hr - tz
  
  # Check boundary conditions
  if (hr > 24) {
    hr <- hr - 24
    # Time is for the next day
    dy <- dy + 1
  } else if (hr < 0) {
    hr <- hr + 24
    # The time is for the previous day
    dy <- dy - 1
  }
  
  # Convert universal time to Greenwich sidereal time
  # Calculate the julian day number for universal time
  jd_ut <- julianDayNumber(year, month, dy + hr/24)
  # Calculate delta t in seconds
  delta_t <- deltaT(year, month)
  # Calculate the julian day number for dynamical time
  jd_td <- jd_ut + delta_t / SEC2DAY
  # Calculate Greenwich sidereal time
  gSidTime <- gst(jd_ut, delta_t)
  
  # Convert Greenwich sidereal time to local sideral time
  lst <- gstToLst(gSidTime[[1]][1], obsLong)
  
  # Calculate the apparent place of the Sun in RA and Dec
  pos_sun <- apparentPlaceSun(jd_td)
  polar <- rectToPolar(pos_sun[[1]][[1]], pos_sun[[1]][[2]], pos_sun[[1]][[3]])
  RA <- polar[3]
  Dec <- polar[2]
  
  # Convert the right ascension to the hour angle
  ha <- lst - RA
  
  # Convert equatorial coordinates to horizon coordinates
  hor_coord <- equatorialToHorizon(ha*RAD2HR, Dec*RAD2DEG, obsLat)
  
  res <- list(year, month, dy, hr, jd_ut, hr, hor_coord)
  names(res) <- c("Year", "Month", "Day", "Hour", "JD_UT", "Horizon Coord")
  
  return(res)
}

# Calculate the altitude and azimuth of the Sun over a single day
findAltAzmSunDay <- function(year, month, day, tz, dst, obsLat, obsLong)
{
  # Create data frame to store data
  df <- data.frame(matrix(0.0, nrow=24, ncol=6))
  
  for (i in seq(from = 1, to = 24, by = 1))
  {
    df[i,1] <- year
    df[i,2] <- month
    df[i,3] <- day
    hr <- i-1
    df[i,4] <- hr
    hor_coord <- findHorizonCoordSun(year, month, day, 
                                     hr, tz, dst, obsLat, obsLong)
    df[i,5] <- hor_coord[[7]][1]
    df[i,6] <- hor_coord[[7]][2]
  }
  
  return(df)
}
