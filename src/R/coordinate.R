
# Functions to calculate right ascension (RA) and declination (Dec) from position vectors
# and perform coordinate transformations

# Convert a rectangular position vector to polar coordinates
# pos is a list with the position vector being a member of the list
# Returns the polar distance r
# phi the longitudinal component in radians
# theta the latitudinal component in radians
rectToPolar <- function(pos){
  rho_sqr <- pos[["Position Vector"]][1] * pos[["Position Vector"]][1] + 
    pos[["Position Vector"]][2] * pos[["Position Vector"]][2]
  m_r <- sqrt(rho_sqr + pos[["Position Vector"]][3] * pos[["Position Vector"]][3])
  m_phi <- 0.0
  
  if (pos[["Position Vector"]][1] == 0.0 & pos[["Position Vector"]][2] == 0.0) {
    m_phi <- 0.0
  } else {
    m_phi <- atan2(pos[["Position Vector"]][2], pos[["Position Vector"]][1])
  }
  
  if (m_phi < 0.0) {m_phi <- m_phi + PI2}
  
  rho <- sqrt(rho_sqr)
  m_theta <- 0.0
  
  if (pos[["Position Vector"]][3] == 0.0 & rho == 0.0) {
    m_theta <- 0.0
  } else {
    m_theta <- atan2(pos[["Position Vector"]][3], rho)
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
equatorialToEcliptical <- function(ra, dec, obliquity)
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
  names(z) <- c("Eclipt_Long", "Eclipt_Lat")
  
  return (z)
}

# Convert ecliptical coordinates to equatorial coordinates
# ecliptic latitude - beta in decimal degrees, [-90 - + 90] degrees
# ecliptic longitude - lambda in decimal degrees, [0 360] degrees
# obliquity in decimal degrees
eclipticalToEquatorial <- function(beta, lambda, obliquity)
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
