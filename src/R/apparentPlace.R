
# Functions to calculate the apparent place of a solar system body

# Epoch of observation
# jd = julian day number
epochOfObs <- function(jd)
{
  t_prime <- (jd - EPOCHJ2000) / DAYSJULCENT
  
  mean_anomaly <- (357.528 + 35999.050 * t_prime) * PI2 / 360
  
  s <- 0.001658 * sin(mean_anomaly + 0.01671 * sin(mean_anomaly))
  
  t <- jd + (s / SEC2DAY)
  
  return(t)
}

# Geometric distance of a solar system body from the Earth
# Body is a position vector of length 3
# Earth is a position vector of length 3
geometricDistance <- function(body, earth)
{
  return(vecNorm(body - earth))
}

# body is a position vector of length 3
# earth is a velocity vector of length 3
aberrationOfLight <- function(body, earth)
{
  V <- earth / CAUD
  
  return(vecNorm(body) * V)
}

