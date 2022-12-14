
# Functions to calculate the apparent place solar system bodies

# Calculate the apparent place of the Sun
# jd = the julian day number of interest
apparentPlaceSun <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of the Sun and convert from
  # KM and KM/sec to AU and AU/day
  sun_ssb <- positionSunSSB(t)
  sun_ssb <- sun_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the sun and the Earth
  u <- sun_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the sun and Earth and update the
  # sun's position and velocity
  tau <- geom_dist / CAUD
  sun_ssb <- positionSunSSB(t - tau)
  sun_ssb <- sun_ssb / KM2AU
  u1 <- sun_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Calculate the apparent place of the Moon
# jd = the julian day number of interest
apparentPlaceMoon <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/Day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU

  # Extract the geocentric position and velocity of the moon and convert from
  # KM and KM/sec to AU and AU/Day
  moon_geo <- positionMoonGEO(t)
  moon_geo <- moon_geo / KM2AU

  # Calculate the geometric distance of the moon in AU
  geom_dist <- vecNorm(moon_geo[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  moon_geo <- positionMoonGEO(t - tau)
  moon_geo <- moon_geo / KM2AU

  # Adjust for the abberation of light
  abberation <- aberrationOfLight(moon_geo[,1], earth_ssb[,2])
  u1 <- moon_geo[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u2 <- prec %*% u1
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u3 <- nut_matrix %*% u2
  
  sinpi = EARTHRADKM / (geom_dist * KM2AU)
  hor_parallax <- asin(sinpi)
  semi_diameter <- asin(0.272493 * sinpi)
  
  # Calculate the horizontal parallax in radians
  hor_parallax <- asin(EARTHRADKM/(geom_dist*KM2AU))
  
  # Return the moon's position vector, geocentric distance in earth radii,
  # horizontal parallax in radians, and the semi-diameter in radians
  z <- list(c(u3[1,1], u3[2,1], u3[3,1]), geom_dist*KM2AU/EARTHRADKM,
            hor_parallax, semi_diameter)
  names(z) <- c("Position Vector", "Geometric Distance",
                "Horizontal Parallax", "Semi-Diameter")
  
  return (z)
}

# Calculate the apparent place of Mercury
# jd = the julian day number of interest
apparentPlaceMercury <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of Mercury and convert from
  # KM and KM/sec to AU and AU/day
  planet_ssb <- positionMercurySSB(t)
  planet_ssb <- planet_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the planet and the Earth
  u <- planet_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  planet_ssb <- positionMercurySSB(t - tau)
  planet_ssb <- planet_ssb / KM2AU
  u1 <- planet_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Calculate the apparent place of Venus
# jd = the julian day number of interest
apparentPlaceVenus <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of Venus and convert from
  # KM and KM/sec to AU and AU/day
  planet_ssb <- positionVenusSSB(t)
  planet_ssb <- planet_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the planet and the Earth
  u <- planet_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  planet_ssb <- positionVenusSSB(t - tau)
  planet_ssb <- planet_ssb / KM2AU
  u1 <- planet_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Calculate the apparent place of Mars
# jd = the julian day number of interest
apparentPlaceMars <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of Mars and convert from
  # KM and KM/sec to AU and AU/day
  planet_ssb <- positionMarsSSB(t)
  planet_ssb <- planet_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the planet and the Earth
  u <- planet_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  planet_ssb <- positionMarsSSB(t - tau)
  planet_ssb <- planet_ssb / KM2AU
  u1 <- planet_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Calculate the apparent place of Jupiter
# jd = the julian day number of interest
apparentPlaceJupiter <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of Jupiter and convert from
  # KM and KM/sec to AU and AU/day
  planet_ssb <- positionJupiterSSB(t)
  planet_ssb <- planet_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the planet and the Earth
  u <- planet_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  planet_ssb <- positionJupiterSSB(t - tau)
  planet_ssb <- planet_ssb / KM2AU
  u1 <- planet_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Calculate the apparent place of Saturn
# jd = the julian day number of interest
apparentPlaceSaturn <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of Saturn and convert from
  # KM and KM/sec to AU and AU/day
  planet_ssb <- positionSaturnSSB(t)
  planet_ssb <- planet_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the planet and the Earth
  u <- planet_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  planet_ssb <- positionSaturnSSB(t - tau)
  planet_ssb <- planet_ssb / KM2AU
  u1 <- planet_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Calculate the apparent place of Uranus
# jd = the julian day number of interest
apparentPlaceUranus <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of Uranus and convert from
  # KM and KM/sec to AU and AU/day
  planet_ssb <- positionUranusSSB(t)
  planet_ssb <- planet_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the planet and the Earth
  u <- planet_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  planet_ssb <- positionUranusSSB(t - tau)
  planet_ssb <- planet_ssb / KM2AU
  u1 <- planet_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Calculate the apparent place of Neptune
# jd = the julian day number of interest
apparentPlaceNeptune <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of Neptune and convert from
  # KM and KM/sec to AU and AU/day
  planet_ssb <- positionNeptuneSSB(t)
  planet_ssb <- planet_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the planet and the Earth
  u <- planet_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  planet_ssb <- positionNeptuneSSB(t - tau)
  planet_ssb <- planet_ssb / KM2AU
  u1 <- planet_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Calculate the apparent place of Pluto
# jd = the julian day number of interest
apparentPlacePluto <- function(jd)
{
  # Calculate the epoch of observation
  t <- epochOfObs(jd)
  
  # Extract the barycentric position and velocity of the Earth and convert from
  # KM and KM/sec to AU and AU/day
  earth_ssb <- positionEarthSSB(t)
  earth_ssb <- earth_ssb / KM2AU
  
  # Extract the barycentric position and velocity of Pluto and convert from
  # KM and KM/sec to AU and AU/day
  planet_ssb <- positionPlutoSSB(t)
  planet_ssb <- planet_ssb / KM2AU
  
  # Calculate the geometric distance between the positions of the center of mass
  # of the planet and the Earth
  u <- planet_ssb - earth_ssb
  geom_dist <- vecNorm(u[,1])
  
  # Calculate the light travel time between the planet and Earth and update the
  # planet position and velocity
  tau <- geom_dist / CAUD
  planet_ssb <- positionPlutoSSB(t - tau)
  planet_ssb <- planet_ssb / KM2AU
  u1 <- planet_ssb - earth_ssb
  
  # Adjust for the abberation of light
  abberation <- aberrationOfLight(u1[,1], earth_ssb[,2])
  u2 <- u1[,1] + abberation
  
  # Adjust for Precession
  prec <- precessionMatrix(jd)
  u3 <- prec %*% u2
  
  # Adjust for Nutation
  nut_matrix <- nutation_matrix(jd)
  u4 <- nut_matrix %*% u3
  
  z <- list(c(u4[1,1], u4[2,1], u4[3,1]), geom_dist)
  names(z) <- c("Position Vector", "Geometric Distance")
  
  return (z)
}

# Epoch of observation
# jd = julian day number
epochOfObs <- function(jd)
{
  t_prime <- (jd - EPOCHJ2000) / DAYSJULCENT
  
  mean_anomaly <- (357.528 + 35999.050 * t_prime) * PI2 / 360
  
  s <- 0.001658 * sin(mean_anomaly + 0.01671 * sin(mean_anomaly))
  
  t <- jd + (s / SEC2DAY)
  
  return (t)
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
  
  return (u1)
}

# body is a position vector of length 3
# earth is a velocity vector of length 3
aberrationOfLight <- function(body, earth)
{
  V <- earth / CAUD
  
  return (vecNorm(body) * V)
}
