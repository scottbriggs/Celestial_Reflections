
# Functions to convert between degrees, arcminutes, and arcseconds and
# decimal degrees, and hours, minutes, and seconds to decimal hours

# dms is a vector consisting of degrees, arcminutes, and arcseconds
dmsToDeg <- function(dms)
{
  sgn <- 1
  deg <- dms[1]
  amin <- dms[2]
  asec <- dms[3]
  
  if (deg < 0){
    sgn <- -1
  }
  
  decDeg <- sgn * (abs(deg) + amin/60 + asec/3600)
  
  return (decDeg)
}

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

# hms is a vector consisting of hours, minutes, and seconds
hmsToHour <- function(hms)
{
  hr <- hms[1]
  min <- hms[2]
  sec <- hms[3]
  
  decHr <- hr + min/60 + sec / 3600
  
  return (decHr)
}

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

dmsString <- function(dms)
{
  str <- sprintf("%s %d\u00B0 %d' %.2f",
                 dms[[1]], dms[[2]][1], dms[[2]][2], dms[[2]][3])
  
  return (str)
}

hmsString <- function(hms)
{
  str <- sprintf("%dh %dm %.2fs", hms[1], hms[2], hms[3])
  
  return(str)
}
