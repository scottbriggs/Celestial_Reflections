
# Functions to convert between degrees, arcminutes, and arcseconds and
# decimal degrees, and hours, minutes, and seconds to decimal hours

# Convert degrees, arcminutes, and arcseconds to decimal degrees
# dms is a vector consisting of degrees, arcminutes, and arcseconds
dmsToDeg <- function(dms)
{
  sgn <- 1
  deg <- dms[1]
  amin <- dms[2]
  asec <- dms[3]
  
  if (deg < 0 | amin < 0 | asec < 0){
    sgn <- -1
  }
  
  decDeg <- sgn * (abs(deg) + abs(amin)/60 + abs(asec)/3600)
  
  return (decDeg)
}

# Convert decimal degrees to degrees, arcminutes, and arcseconds
# decDeg represents decimal degrees
degToDMS <- function(decDeg)
{
  sgn <- "+"
  if (decDeg <0){
    sgn <- "-"
  }
  
  x <- abs(decDeg)
  deg <- as.integer(x)
  x <- (x - deg) * 60
  amin <- as.integer(x)
  asec <- (x - amin) * 60
  
  dms <- c(deg, amin, asec)
  
  return (list(sgn, dms))
}

# Convert hours, minutes, and seconds to decimal hours
# hms is a vector consisting of hours, minutes, and seconds
hmsToHour <- function(hms)
{
  hr <- hms[1]
  min <- hms[2]
  sec <- hms[3]
  
  decHr <- hr + min/60 + sec / 3600
  
  return (decHr)
}

# Convert decimal hours to hours, minutes, and seconds
# decHr represents decimal hours
hourToHMS <- function(decHr)
{
  x <- decHr
  hr <- as.integer(x)
  x <- (x - hr) * 60
  min <- as.integer(x)
  sec <- (x - min) * 60
  
  hms <- c(hr, min, sec)
  
  return (hms)
}

# Convert decimal hours to hours and minutes
# hm is a vector consisting of hours and minutes
hourToHM <- function(decHr)
{
  x <- decHr
  hr <- as.integer(x)
  x <- (x - hr) * 60
  min <- as.integer(x)
  x <- (x - min) * 60
  if (x >= 30) {
    min <- min + 1
    if (min == 60) {
      hr <- hr + 1
      min = 0
    }
  }
  
  hm <- c(hr, min)
  
  return (hm)
}

# Return string for degrees, arcminutes, and arcseconds
dmsString <- function(dms)
{
  str <- sprintf("%s %d\u00B0 %d' %.2f",
                 dms[[1]], dms[[2]][1], dms[[2]][2], dms[[2]][3])
  
  return (str)
}

# Return string for hours, minutes, and seconds
hmsString <- function(hms)
{
  str <- sprintf("%dh %dm %.2fs", hms[1], hms[2], hms[3])
  
  return(str)
}

# Return string for hours and minutes
hmString <- function(hm)
{
  str <- sprintf("%dh %dm ", hm[1], hm[2])
  
  return(str)
}

# Calculate delta-t in seconds to get an estimate of the difference between
# universal and dynamical time for any year in the past or future
deltaT <- function(year, month)
{
  delta_t <- 0.0
  
  y <- year + (month - 0.5) / 12
  
  if (year < -500) {
    u <- (year - 1820) / 100
    delta_t <- -20 + (32 * u * u)
  } else if (year >= -500 & year < 500) {
    u <- y/100
    delta_t <- 10583.6 + u*(-1014.41 + u*(33.78311 + u*(-5.952053 + 
                                                          u*(-0.1798452 + u*(0.022174192 + u*(0.0090316521))))))
  } else if (year >= 500 & year < 1600) {
    u = (y - 1000)/100
    delta_t <- 1574.2 + u*(-556.01 + u*(71.23472 + u*(0.319781 + 
                                                        u*(-0.8503463 + u*(-0.005050998 + u*(0.0083572073))))))
  } else if (year >= 1600 & year < 1700) {
    t <- y - 1600
    delta_t <- 120 + t*(-0.9808 + t*(-0.01532 + t*(1/7129)))
  } else if (year >= 1700 & year < 1800) {
    t <- y - 1700
    delta_t <- 8.83 + t*(0.1603 + t*(-0.0059285 + t*(0.00013336 +
                                                       t*(-1/1174000))))
  } else if (year >= 1800 & year < 1860) {
    t <- y - 1800
    delta_t <- 13.72 + t*(-0.332447 + t*(0.0068612 + t*(0.00411116 +
                                                          t*(-0.00037436 + t*(0.0000121272 + t*(-0.0000001699 +
                                                                                                  t*(0.000000000875)))))))
  } else if (year >= 1860 & year < 1900) {
    t <- y - 1860
    delta_t <- 7.62 + t*(0.5737 + t*(-0.251754 + t*(0.01680668 +
                                                      t*(-0.0004473624 + t*(1/233174)))))
  } else if (year >= 1900 & year < 1920) {
    t <- y - 1900
    delta_t <- -2.79 + t*(1.494119 + t*(-0.0598939 + t*(0.0061966 +
                                                          t*(-0.000197))))
  } else if (year >= 1920 & year < 1941) {
    t <- y - 1920
    delta_t <- 21.20 + t*(0.84493 + t*(-0.076100 + t*(0.0020936)))
  } else if (year >= 1941 & year < 1961) {
    t <- y - 1950
    delta_t <- 29.07 + t*(0.407 + t*(-2/233)) + t*(1/2547)
  } else if (year >= 1961 & year < 1986) {
    t <- y - 1975
    delta_t <- 45.45 + t*(1.067 + t*(-1/260)) + t*(-1/718)
  } else if (year >= 1986 & year < 2005) {
    t <- y - 2000
    delta_t <- 63.86 + t*(0.3345 + t*(-0.060374 + t*(0.0017275 +
                                                       t*(0.000651814 + t*(0.00002373599)))))
  } else if (year >= 2005 & year < 2050) {
    t <- y - 2000
    delta_t <- 62.92 + t*(0.32217 + t*(0.005589))
  } else if (year >= 2050 & year < 2150) {
    t <- y - 2000
    delta_t <- -20 + 32 * ((y - 1829)/100) * ((y - 1829)/100) - 
      0.5628 * (2150 - y)
  } else {
    u <- (year - 1820) / 100
    delta_t <- -20 + (32 * u * u)
  }
  
  return (delta_t)
}

# Calculate the mean and apparent sidereal time in radians
# jd is the julian day number in universal time at the time of observation
siderealTime <- function(jd_ut, deltaT)
{
  T <- (jd_ut - EPOCHJ2000) / DAYSJULCENT
  
  # Calculate mean sidereal time in degrees
  gmst <- 280.46061837 + 360.98564736629 * (jd_ut - EPOCHJ2000) +
    (0.000387933 - (1/38710000) * T) * T * T
  
  # Convert result to the range 0 - 360 degrees
  gmst <- amodulo(gmst, 360)
  
  # Convert degrees to radians
  gmst <- gmst * DEG2RAD
  
  # Calculate dynamical time
  jd_td <- jd_ut + deltaT/24
  
  # Get the nutation angles in longitude and obliquity
  nut_angles <- nutationAngles(jd_td)
  
  # Get the mean and true obliquity of the ecliptic
  ob <- obliquity(jd_td, nut_angles)
  
  # Calculate apparent sidereal time
  gast <- gmst + nut_angles["Nut_Long"] * cos(ob["True_Obliquity"])
  
  z <- c(gmst, gast)
  names(z) <- c("GMST", "GAST")
  
  return (z)
}

# Julian day number corrected for local civil time including
# daylight savings time and time zone
# year, month, and day of interest
# dst - daylight savings time, 1 for dst, otherwise 0
# tz - time zone correction, negative if time zone is west of Greenwich,
# positive if east of Greenwich
timeZoneCorrection <- function(year, month, day, dst, tz)
{
  
}
