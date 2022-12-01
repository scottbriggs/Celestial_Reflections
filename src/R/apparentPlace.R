
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

# Body is a position vector of length 3
# Earth is a position vector of length 3
# Sun is a position vector of length 3
relativisticDeflectionOfLight <- function(body, earth, sun)
{
  body_geo <- body - earth
  body_helio <- body - sun
  earth_helio <- earth - sun
  
  u <- unitVector(body_geo)
  q <- unitVector(body_helio)
  e <- uitVector(earth_helio)
  
  g1 <- MUC / vecNorm(earth_helio)
  g2 <- 1 + dotProduct(q,e)
  
  tmp1 <- dotProduct(u,q) * e
  tmp2 <- dotProduct(e,u) * q
  
  u1 <- vecNorm(u) * (u + g1/g2 * (tmp1 - tmp2))
  
  return(u1)
}

# body is a position vector of length 3
# earth is a velocity vector of length 3
aberrationOfLight <- function(body, earth)
{
  V <- earth / CAUD
  
  return(vecNorm(body) * V)
}

