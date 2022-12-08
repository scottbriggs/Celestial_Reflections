
# Functions to calculate right ascension (RA) and declination (Dec) from position vectors
# and perform coordinate transformations

# Convert a rectangular position vector to polar coordinates
# pos is a position vector (x, y, z)
# Returns the polar distance r
# phi the longitudinal component
# theta the latitudinal component
rectToPolar <- function(pos){
  rho_sqr <- pos[1] * pos[1] + pos[2] * pos[2]
  m_r <- sqrt(rho_sqr + pos[3] * pos[3])
  m_phi <- 0.0
  
  if (pos[1] == 0.0 & pos[2] == 0.0) {
    m_phi <- 0.0
  } else {
    m_phi <- atan2(pos[2], pos[1])
  }
  
  if (m_phi < 0.0) {m_phi <- m_phi + PI2}
  
  rho <- sqrt(rho_sqr)
  m_theta <- 0.0
  
  if (pos[3] == 0.0 & rho == 0.0) {
    m_theta <- 0.0
  } else {
    m_theta <- atan2(pos[3], rho)
  }
  
  z <- c(m_r, m_phi, m_theta)
  names(z) <- c("Norm", "Polar_Long", "Polar_Lat")
  
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

# Convert an equatorial position vector to an ecliptical position vector
# Obliquity is in radians
equatorialToEcliptical <- function(ra, dec, obliquity)
{
  sinra <- sin(ra)
  cosra <- cos(ra)
  sine <- sin(obliquity)
  cose <- cos(obliquity)
  sindec <- sin(dec)
  cosdec <- cos(dec)
  
  numerator <- sinra * cose + tan(dec) * sine
  eclipt_long <- atan2(numerator, cosra)
  eclipt_lat <- asin(sindec * cose - cosdec * sine * sinra)
  
  z <- c(eclipt_long, eclipt_lat)
  names(z) <- c("Eclipt_Long", "Eclipt_Lat")
  
  return (z)
}

# Convert an ecliptical position vector to an equatorial position vector
# Obliquity is in radians
eclipticalToEquatorial <- function(elong, elat, obliquity)
{
  sinelat <- sin(elat)
  coselat <- cos(elat)
  sine <- sin(obliquity)
  cose <- cos(obliquity)
  sinelong <- sin(elong)
  coselong <- cos(elong)
  
  numerator <- sinelong * cose - tan(elat) * sine
  ra <- atan2(numerator, coselong)
  dec <- asin(sinelat * cose + coselat * sine * sinelong)
  
  z <- c(ra, dec)
  names(z) <- c("RA", "DEC")
  
  return (z)
}

# Convert equatorial coordinates (hour angle, declination) in radians to
# horizon coordinates (azimuth, altitude) in radians.
equatorialToHorizon <- function(hourAngle, dec, obs_lat)
{
  sinH <- sin(hourAngle)
  cosH <- cos(hourAngle)
  cosLat <- cos(obs_lat)
  sinLat <- sin(obs_lat)
  cosDec <- cos(dec)
  sinDec <- sin(dec)
  tanDec <- tan(dec)
  
  azimuth <- atan2(-sinH, -sinLat * cosH + cosLat * tanDec)
  altitude <- arcsin(cosLat * cosDec * cosH + sinLat * sinDec)
  
  z <- c(azimuth, altitude)
  names(z) <- c("Azimuth", "Altitude")
  
  return (z)
}

# Convert horizon coordinates (azimuth, altitude) in radians to
# equatorial coordinates (hour angle, declination) in radians.
horizonToEquatorial <- function(azimuth, altitude, obs_lat)
{
  sinAlt <- sin(altitude)
  cosAlt <- cos(altitude)
  tanAlt <- tan(altitude)
  sinAz <- sin(azimuth)
  cosAz <- cos(azimuth)
  sinLat <- sin(obs_lat)
  cosLat <- cos(obs_lat)
  
  hourAngle <- atan2(-sinAz, -cosAz * sinLat + tanAlt * sinLat)
  dec <- arcsin(cosAlt * cosAz * cosLat + sinAlt * sinLat)
  
  z <- c(hourAngle, dec)
  names(z) <- c("Hour Angle", "DEC")
  
  return(z)
}
