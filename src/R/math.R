# Math routines that are useful for computations

# x is a vector of length 3
# y is a vector of length 3
dotProduct <- function(x, y)
{
  return(x%*%y)
}

# x is a vector of length 3
# y is a vector of length 3
crossProduct <- function(x,y)
{
  result <- c(0.0, 0.0, 0.0)
  for (i in 1:length(x)){
    result[1] <- x[2] * y[3] - x[3] * y[2]
    result[2] <- x[3] * y[1] - x[1] * y[3]
    result[3] <- x[1] * y[2] - x[2] * y[1]
  }
  
  return(result)
}

# Calculate the magnitude of a vector x of length 3
vecNorm <- function(x)
{
  return(sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]))
}

# x is a vector of length 3
unitVector <- function(x)
{
  mag <- vecNorm(x)
  result <- c(0.0, 0.0, 0.0)
  result[1] <- x[1] / mag
  result[2] <- x[2] / mag
  result[3] <- x[3] / mag
  
  return (result)
}

# Reduce a to the range 0 <= a < b
amodulo <- function(a, b)
{
  return(x <- a - b * floor(a/b))
}

# Returns a rotation matrix based on the axis (x, y, or z) and the angle phi
# The angle phi is in radians
rotationMatrix <- function(axis, phi)
{
  mat <- matrix(0.0, nrow = 3, ncol = 3)
  cosphi <- cos(phi)
  sinphi <- sin(phi)
  
  if (axis == 1) {
    mat[1,1] <- 1.0
    mat[2,2] <- cosphi
    mat[3,3] <- cosphi
    mat[2,3] <- sinphi
    mat[3,2] <- -sinphi
  } else if (axis == 2) {
    mat[1,1] <- cosphi
    mat[1,3] <- -sinphi
    mat[2,2] <- 1.0
    mat[3,1] <- sinphi
    mat[3,3] <- cosphi
  } else if (axis == 3){
    mat[1,1] <- cosphi
    mat[2,2] <- cosphi
    mat[3,3] <- 1.0
    mat[2,1] <- -sinphi
    mat[1,2] <- sinphi
  } else {print("Axis is wrong value")}
  
  return (mat)
}

# Interpolate three values using the interpolating value n
interpolate <- function(vec, n)
{
  # Take differences
  a <- vec[2] - vec[1]
  b <- vec[3] - vec[2]
  c <- b - a
  
  y <- vec[2] + n/2 * (a + b + n * c)
  
  return(y)
}

