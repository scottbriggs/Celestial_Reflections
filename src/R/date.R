
# Calculate the julian day number allowing fractional parts of days. In other,
# words, the day is not limited to being an integer.

julianDayNumber <- function(year, month, day)
{
  int_day <- trunc(day)
  frac_day <- day - int_day
  
  jd <- julianDayInt(year, month, int_day)
  
  if (frac_day > 0.5){
    frac_day <- frac_day - 0.5
    jd <- jd + frac_day}
  else if (frac_day < 0.5){
    frac_day <- 0.5 - frac_day
    jd <- jd - frac_day
    }
    
    return(jd)
}

# Calculates the Julian Day given the month, day, and year. The algorithm works
# for any date in the common era (CE) or before the common era (BCE).
# The Julian Day Number is calculated for a calendar date at 12 noon.
# The year, month, and day are integers
julianDayInt <- function(year, month, day)
{
  # Gregorian calendar adopted October 15, 1582
  IGREG <- (15 + 31 * (10 + 12 * 1582))
  jm <- 0
  
  jy <- year
  
  # This code causes the julian day to be in error for years
  # less than zero
  #if (jy < 0) { jy <- jy + 1}
  
  if (month > 2)
  { 
    jm <- month + 1
  } else {
    jy <- jy - 1
    jm <- month + 13
  }
  
  julday <- floor(365.25 * jy) + floor(30.6001 * jm) + day + 1720995
  
  if (day + 31 * (month + 12 * year) >= IGREG)
  {
    ja <- trunc(0.01 * jy)
    julday <- julday + 2 - ja + trunc(0.25 * ja)
  }
  
  return (julday)
}