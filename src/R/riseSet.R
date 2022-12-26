
# Functions to calculate rise, transit, and set times of the planets, Sun,
# and Moon

# jd_ut is the universal time for the event of interest
# deltaT is the different in seconds between universal time and dynamical time
# obslat is the observer latitude in degrees, +N, -S
# obslong is the observer longitude in degrees, +E, - W
# event is the event of interest - sun, moon, or planet/star
# func1 is a function passed in for planets or stars
riseSetEvent <- function(jd_ut, deltaT, obsLat, obsLong, event, func1)
{
  # Initialization
  std_alt <- 0
  riseFlag <- FALSE
  setFlag <- FALSE
  culminationFlag <- FALSE
  riseString <- ""
  setString <- ""

  latRad <- obsLat * DEG2RAD
  sinLat <- sin(latRad)
  cosLat <- cos(latRad)
  
  # RA and Dec for the event of interest at day-1, day, and day+1
  ra1 <- 0
  ra2 <- 0
  ra3 <- 0
  dec1 <- 0
  dec2 <- 0
  dec3 <- 0
  
  # Calculate dynamical time
  jd_td <- jd_ut + deltaT/86400
  
  # Calculate the Greenwich Mean and Apparent Sidereal Time
  sid_time <- siderealTime(jd_ut, deltaT)
  
  if (event == "Sun") {
    
    # standard altitude for the sun
    std_alt <- -50/60*DEG2RAD
    
    # apparent place of the sun on day-1
    pos <- apparentPlaceSun(jd_ut -1)
    polar <- rectToPolar(pos)
    ra1 <- polar[2]
    dec1 <- polar[3]
    
    # apparent place of the sun on day
    pos <- apparentPlaceSun(jd_ut)
    polar <- rectToPolar(pos)
    ra2 <- polar[2]
    dec2 <- polar[3]
    
    # apparent place of the sun on day+1
    pos <- apparentPlaceSun(jd_ut +1)
    polar <- rectToPolar(pos)
    ra3 <- polar[2]
    dec3 <- polar[3]
    
  } else if (event == "Moon") {
    
    # apparent place of the moon on day-1
    pos <- apparentPlaceMoon(jd_ut -1)
    polar <- rectToPolar(pos)
    ra1 <- polar[2]
    dec1 <- polar[3]
    
    # apparent place of the moon on day
    pos <- apparentPlaceMoon(jd_ut)
    polar <- rectToPolar(pos)
    ra2 <- polar[2]
    dec2 <- polar[3]
    
    # standard altitude for the moon
    std_alt = -34/60*DEG2RAD - pos[["Semi-Diameter"]][1] + 
      pos[["Horizontal Parallax"]][1]
    
    # apparent place of the Moon on day+1
    pos <- apparentPlaceMoon(jd_ut +1)
    polar <- rectToPolar(pos)
    ra3 <- polar[2]
    dec3 <- polar[3]
    
  } else {
    
    # Planets and stars
    # Standard Altitude
    std_alt <- -34/60*DEG2RAD
    
    # apparent place of the planet on day-1
    pos <- func1(jd_ut - 1)
    polar <- rectToPolar(pos)
    ra1 <- polar[2]
    dec1 <- polar[3]
    
    # apparent place of the planet on day
    pos <- func1(jd_ut)
    polar <- rectToPolar(pos)
    ra2 <- polar[2]
    dec2 <- polar[3]
    
    # apparent place of the planet on day+1
    pos <- func1(jd_ut + 1)
    polar <- rectToPolar(pos)
    ra3 <- polar[2]
    dec3 <- polar[3]
  }
  
  # Check RA values on 24 hour boundary
  if ((ra2 < ra1) & (ra3 > ra2)){
    ra2 <- ra2 + PI2
    ra3 <- ra3 + PI2
  }else if ((ra2 > ra1) & (ra3 < ra2)){
    ra3 <- ra3 + PI2
  }
  
  # Calculate cosH0
  cosH0 <- (sin(std_alt) - sinLat * sin(dec2)) / cosLat * cos(dec2)
  
  if (cosH0 > 1) {
    riseString <- "Never Rises"
    setString <- "Never Rises"
    riseFlag <- FALSE
    setFlag <- FALSE
    culminationFlag <- FALSE
  } else if (cosH0 < -1) {
    riseString <- "Circumpolar"
    setString <- "Circumpolar"
    riseFlag <- FALSE
    setFlag <- FALSE
    culminationFlag <- TRUE
  } else {
    riseString <- "Rise"
    setString <- "Set"
    riseFlag <- TRUE
    setFlag <- TRUE
    culminationFlag <- TRUE
  }
  
  H0 <- acos(cosH0) * RAD2DEG
  H0 <- amodulo(H0, 180)
  
  m0 <- (ra2 * RAD2DEG + obsLong - sid_time[2] * RAD2DEG) / 360
  m1 <- m0 - (H0 / 360)
  m2 <- m0 + (H0 / 360)
  
  if (m0 < 0) {
    m0 <- m0 + 1
  } else if (m0 > 1) {
    m0 <- m0 - 1
  }
  
  if (m1 < 0) {
    m1 <- m1 + 1
  } else if (m1 > 1) {
    m1 <- m1 - 1
  }
  
  if (m2 < 0) {
    m2 <- m2 + 1
  } else if (m2 > 1) {
    m2 <- m2 - 1
  }
  
  thetaM0 <- sid_time[2] * RAD2DEG + 360.985647 * m0
  thetaM1 <- sid_time[2] * RAD2DEG + 360.985647 * m1
  thetaM2 <- sid_time[2] * RAD2DEG + 360.985647 * m2
  
  nM0 <- m0 + deltaT/86400
  nM1 <- m1 + deltaT/86400
  nM2 <- m2 + deltaT/86400
  
  raM0 <- interpolate(c(ra1, ra2, ra3), nM0)
  raM1 <- interpolate(c(ra1, ra2, ra3), nM1)
  raM2 <- interpolate(c(ra1, ra2, ra3), nM2)
  decM1 <- interpolate(c(dec1, dec2, dec3), nM1)
  decM2 <- interpolate(c(dec1, dec2, dec3), nM2)
  
  
  
}

# Function to calculate roots of a quadratic function based on three
# equidistant values of the function
# yMinus - value of the function at x = -1
# y - value of the function at x = 0
# yPlus - value of the function at x = 1
# xe - abscisa of the extreme value
# ye - ordinate of the extreme value
# root1 - first root of the function
# roo2t - second root of the function
# nRoots - number of roots found in the interval [-1, 1]
quadraticInterpolation <- function(yMinus, y, yPlus)
{
  dx <- 0
  root1 <- 0
  root2 <- 0
  nRoot <- 0
  
  # Coefficients of the interpolating parabola, y = ax2 + bx + c
  a <- 0.5 * (yPlus + yMinus) - y
  b <- 0.5 * (yPlus - yMinus)
  c <- y
  
  # Find extreme value
  xe <- -b / (2.0 * a)
  ye <- (a * xe  + b) * xe +c
  
  # Discriminant
  dis <- b * b - 4.0 * a * c
  
  # If the discriminant is >= 0, the function has roots
  if (dis >= 0){
    dx <- 0.5 * sqrt(dis) / abs(a)
    root1 <- xe - dx
    root2 <- xe + dx
    
    if (abs(root1) <= 1.0){
      nRoot <- nRoot + 1
    }
    
    if (abs(root2) <= 1) {
      nRoot <- nRoot + 1
    }
    
    if (root1 < -1) {
      root1 <- root2
    }
    
    z <- c(xe, ye, root1, root2, nRoot)
    names(z) <- c("Abscissa", "Ordinate", "Root 1", "Root2", "Num Roots")
    
    return (z)
  }
}

# Rise - Set for the Sun using the method outlined in Astronomy on the 
# Personal Computer
# Calculates the rise and/or set of the Sun for a given day
# jd_ut is the julian day number for universal time
# delta T is the correction in seconds for dynamical time
# ObsLat is the observer latitude in degrees, + North, - South
# ObsLong is the observer longitude in degrees, + East, - West
riseSetEventSun <- function(jd_ut, deltaT, obsLat, obsLong)
{
  # Take the integer part of the julian day number so that the rise set calculations
  # start at midnight
  jd_ut_start <- as.integer(jd_ut)
  
  # Calculate the julian day number in dynamical time
  jd_td <- jd_ut_start + deltaT/SEC2DAY
  
  # Calculate the standard altitude for rise and set of the sun
  std_alt <- -50/60 * DEG2RAD
  
  # Convert latitude and longitude to radians
  latRad <- obsLat * DEG2RAD
  cosLat <- cos(latRad)
  sinLat <- sin(latRad)
  longRad <- obsLong * DEG2RAD
  
  # Search for rise and set times for the day of interest
  ap_sun <- apparentPlaceSun(jd_td)
  pos <- c(ap_sun[["Position Vector"]][1], ap_sun[["Position Vector"]][2], 
           ap_sun[["Position Vector"]][3])
  polar <- rectToPolar(pos)
  ra <- polar[2]
  dec <- polar[3]
  hourAngle <- meanSiderealTime(jd_ut_start) * HR2RAD + longRad - ra
  if (hourAngle < 0.0) {
    hourAngle <- hourAngle + PI2
  }
  sinAlt <- (sinLat * sin(dec) + cosLat * cos(dec) * cos(hourAngle)) - sin(std_alt)
  yMinus <- sinAlt
  aboveFlag <- FALSE
  aboveFlag <- (yMinus > 0)
  riseFlag <- FALSE
  setFlag <- FALSE
  localTimeRise <- 0
  localTimeSet <- 0
  hour <- 1
  
  while ( !((hour == 25) | (riseFlag & setFlag)) ) {
    ap_sun <- apparentPlaceSun(jd_td + hour/24)
    pos <- c(ap_sun[["Position Vector"]][1], ap_sun[["Position Vector"]][2], 
             ap_sun[["Position Vector"]][3])
    polar <- rectToPolar(pos)
    ra <- polar[2]
    dec <- polar[3]
    hourAngle <- meanSiderealTime(jd_ut_start + hour/24) * HR2RAD + longRad - ra
    if (hourAngle < 0.0) {
      hourAngle <- hourAngle + PI2
    }
    sinAlt <- (sinLat * sin(dec) + cosLat * cos(dec) * cos(hourAngle)) - sin(std_alt)
    y0 <- sinAlt
    
    ap_sun <- apparentPlaceSun(jd_td + (hour + 1)/24)
    pos <- c(ap_sun[["Position Vector"]][1], ap_sun[["Position Vector"]][2], 
             ap_sun[["Position Vector"]][3])
    polar <- rectToPolar(pos)
    ra <- polar[2]
    dec <- polar[3]
    hourAngle <- meanSiderealTime(jd_ut_start + (hour+1)/24) * HR2RAD + longRad - ra
    if (hourAngle < 0.0) {
      hourAngle <- hourAngle + PI2
    }
    sinAlt <- (sinLat * sin(dec) + cosLat * cos(dec) * cos(hourAngle)) - sin(std_alt)
    yPlus <- sinAlt
    
    vect <- quadraticInterpolation(yMinus, y0, yPlus)
    
    if (vect["Num Roots"] == 1) {
      if (yMinus < 0.0) {
        localTimeRise <- hour + vect["Root 1"]
        riseFlag <- TRUE
      } else {
        localTimeSet <- hour + vect["Root 1"]
        setFlag <- TRUE
      }
    }
    
    if (vect["Num Roots"] == 2) {
      if (vect["Ordinate"] < 0.0) {
        localTimeRise <- hour + vect["Root 2"]
        localTimeSet <- hour + vect["Root 1"]
      } else {
        localTimeRise <- hour + vect["Root 1"]
        localTimeSet <- hour + vect["Root 2"]
      }
      riseFlag <- TRUE
      setFlag <- TRUE
    }
    yMinus <- yPlus
    hour <- hour + 2
  }
  
  z <- list(riseFlag, localTimeRise, setFlag, localTimeSet, aboveFlag)
  names(z) <- c("Rise Flag", "Local Time Rise", "Set Flag", 
                "Local Time Set", "Above Flag")
  
  return(z)
}

