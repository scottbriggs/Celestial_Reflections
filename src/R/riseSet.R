
# Functions to calculate rise, transit, and set times of the planets, Sun,
# and Moon

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

# Rise - Set Events for the Sun
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
  std_alt <- -0.50/60 * DEG2RAD
  
  # Convert latitude and longitude to radians
  latRad <- obsLat * DEG2RAD
  cosLat <- cos(latRad)
  sinLat <- sin(latRad)
  longRad <- obsLong * DEG2RAD
  
  # Search for rise and set times for the day of interest
  ap_sun <- apparentPlaceSun(jd_td)
  pos <- c(ap[["Position Vector"]][1], ap[["Position Vector"]][2], 
           ap[["Position Vector"]][3])
  polar <- rectToPolar(pos)
  ra <- polar[2]
  dec <- polar[3]
  hourAngle <- meanSiderealTime(jd_ut_start) * HR2RAD + longRad - ra
  sinAlt <- (sinLat * sin(dec) + cosLat * cos(dec) * cos(hourAngle)) - sin(std_alt)
  yMinus <- sinAlt
  aboveFlag <- FALSE
  aboveFlag <- (y_minus > 0)
  riseFlag <- FALSE
  setFlag <- FALSE
  localTimeRise <- 0
  localTimeSet <- 0
  hour <- 1
  
  while ( !((hour == 25) | (rises & sets)) ) {
    ap_sun <- apparentPlaceSun(jd_td + hour/24)
    pos <- c(ap[["Position Vector"]][1], ap[["Position Vector"]][2], 
             ap[["Position Vector"]][3])
    polar <- rectToPolar(pos)
    ra <- polar[2]
    dec <- polar[3]
    hourAngle <- meanSiderealTime(jd_ut_start + hour/24) * HR2RAD + longRad - ra
    sinAlt <- (sinLat * sin(dec) + cosLat * cos(dec) * cos(hourAngle)) - sin(std_alt)
    y0 <- sinAlt
    
    ap_sun <- apparentPlaceSun(jd_td + (hour + 1)/24)
    pos <- c(ap[["Position Vector"]][1], ap[["Position Vector"]][2], 
             ap[["Position Vector"]][3])
    polar <- rectToPolar(pos)
    ra <- polar[2]
    dec <- polar[3]
    hourAngle <- meanSiderealTime(jd_ut_start + (hour+1)/24) * HR2RAD + longRad - ra
    sinAlt <- (sinLat * sin(dec) + cosLat * cos(dec) * cos(hourAngle)) - sin(std_alt)
    yPlus <- sinAlt
    
    vect <- quadraticInterpolation(yMinus, y0, yPlus)
    
    if (vect["Num Roots"] == 1) {
      if (yMinus < 0) {
        localTimeRise <- hour + vect["Root 1"]
        riseFlag <- TRUE
      } else {
        localTimeSet <- hour + vect["Root 1"]
        setFlag <- TRUE
      }
    }
    
    if (vect["Num Roots"] == 2) {
      if (vect["Ordinate"] < 0) {
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
  
  z <- list(riseFlag, localRiseTime, setFlag, localSetTime, aboveFlag)
  names(z) <- c("Rise Flag", "Local Rise Time", "Set Flag", 
                "local Set Time", "Above Flag")
  
  return(z)
}

