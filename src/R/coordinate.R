
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
  
  return (c(m_r, m_phi, m_theta))
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
equatorialToEcliptical <- function(pos, obliquity)
{
  pos1 <- rotationMatrix(1, obliquity) %*% pos
  
  return (pos1)
}

# Convert an ecliptical position vector to an equatorial position vector
# Obliquity is in radians
eclipticalToEquatorial <- function(pos, obliquity)
{
  pos1 <- rotationMatrix(1, -obliquity) %*% pos
  
  return(pos1)
}

# Convert equatorial coordinates (hour angle, declination) in radians to
# horizon coordinates (azimuth, altitude) in radians.
equatorialToHorizon <- function(hour_angle, dec, obs_lat)
{
  vec1 <- c(hour_angle, dec, 1)
  vec2 <- rotationMatrix(2, (PI/2 - obs_lat)) %*% vec1
  
  azimuth <- vec2[1]
  altitude <- vec2[2]
  
  return(c(azimuth, altitude))
}

# Convert horizon coordinates (azimuth, altitude) in radians to
# equatorial coordinates (hour angle, declination) in radians.
horizonToEquatorial <- function(azimuth, altitude, obs_lat)
{
  vec1 <- c(azimuth, altitude, 1)
  vec2 <- rotationMatrix(2, -(PI/2 - obs_lat)) %*% vec1
  
  hour_angle <- vec2[1]
  dec <- vec2[2]
  
  return(c(hour_angle, dec))
}