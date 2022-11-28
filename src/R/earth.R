# Calculate the position of the Earth
# using the JPL DE 441 Ephemeris data

positionEarthSSB <- function(jd)
{
  # Get the position and velocity of the earth-moon barycenter and moon
  emb_pos_vel <- positionEMBSSB(jd)
  moon_pos_vel <- positionMoonGEO(jd)
  
  # Calculate the position and velocity of the earth
  earth_pos_vel <- data.frame(matrix(0.0, nrow=2, ncol=3))
  earth_pos_vel <- emb_pos_vel - (moon_pos_vel / (1 + EMRAT))
  
  colnames(earth_pos_vel) <- c('Position Vector', 'Velocity Vector')
  
  return(earth_pos_vel)
}