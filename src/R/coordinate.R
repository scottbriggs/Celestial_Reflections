
# Functions to calculate right ascension (RA) and declination (Dec) from position vectors
# and perform coordinate transformations

# Convert a position vector to RA and Dec
# pos is a position vector (x, y, z)
# Returns R, phi, and theta
pos2RaDec <- function(pos){
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