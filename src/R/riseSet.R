
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

# Rise - Set Events
riseSetEvent <- function()
{
  
}

