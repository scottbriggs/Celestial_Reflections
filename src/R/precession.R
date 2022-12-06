
# Calculate the precession matrix
precessionMatrix <- function(jd)
{
  # Periodic terms in PA and QA - Table 1
  PA_C <- c(-5486.751211, -17.127623, -617.517403, 413.442940, 78.614193,
            -180.732815, -87.676083, 46.140315)
  
  PA_S <- c(667.666730, -2354.886252, -428.152441, 376.202861, 184.778874,
            335.321713, -185.138669, -120.972830)
  
  QA_C <- c(-684.661560, 2446.283880, 399.671049, -356.652376, -186.387003,
            -316.800070, 198.296071, 101.135679)
  
  QA_S <- c(-5523.863691, -549.747450, -310.998056, 421.535876, -36.776172,
            -145.278396, -34.744450, 22.885731)
  
  P1 <- c(708.15, 2309.0, 1620.0, 492.20, 1183.0, 622.0, 882.0, 547.0)
  
  # Periodic terms in XA and YA - Table 2
  XA_C <- c(-819.940624, -8444.676815, 2600.009459, 2755.175630, -167.659835,
            871.855056, 44.769698, -512.313065, -819.415595, -538.071099,
            -189.793622, -402.922932, 179.516345, -9.814756)
  
  XA_S <- c(81491.287984, 787.163481, 1251.296102, -1257.950837, -2966.799730,
            639.744522, 131.600209, -445.040117, 584.522874, -89.756563,
            524.429630, -13.549067, -210.157124, -44.919798)
  
  YA_C <- c(75004.344875, 624.033993, 1251.136893, -1102.212834, -2660.664980,
            699.291817, 153.167220, -950.865637, 499.754645, -145.188210,
            558.116553, -23.923029, -165.405086, 9.344131)
  
  YA_S <- c(1558.515853, 7774.939698, -2219.534038, -2523.969396, 247.850422,
            -846.485643, -1393.124055, 368.526116, 749.045012, 444.704518,
            235.934465, 374.049623, -171.330180, -22.899655)
  
  P2 <- c(256.75, 708.15, 274.20, 241.45, 2309.0, 492.20, 396.10, 288.90,
          231.10, 1610.0, 620.0, 157.87, 220.30, 1200.0)
  
  # Obliquity at J2000.0
  eps0 <- 84381.406 * ARCSEC2RAD
  
  # Julian centuries since J2000.0
  T <- (jd - 2451545.0) / 36525
  T2 <- T * T
  T3 <- T2 * T
  
  # Sum the periodic terms for PA and QA
  sum_pa <- 0.0
  sum_qa <- 0.0
  for(i in 1:8) {
    sum_pa <- sum_pa + PA_C[i]*cos(PI2*T/P1[i]) + PA_S[i]*sin(PI2*T/P1[i])
    sum_qa <- sum_qa + QA_C[i]*cos(PI2*T/P1[i]) + QA_S[i]*sin(PI2*T/P1[i])
  }
  
  # Calculate PA and QA
  PA <- 5851.607687 - 0.1189000 * T - 0.00028913 * T2 + 101E-9 * T3 + sum_pa
  QA <- -1600.886300 + 1.1689818 * T - 0.00000020 * T2 - 437E-9 * T3 + sum_qa
  
  # Convert to radians
  P <- PA * ARCSEC2RAD
  Q <- QA * ARCSEC2RAD
  
  # Calculate the ecliptic pole vector
  temp <- c(sqrt(1 - P*P - Q*Q), 0.0)
  Z <- max(temp)
  S <- sin(eps0)
  C <- cos(eps0) 
  
  ecliptic_pole <- c(0.0, 0.0, 0.0)
  ecliptic_pole[1] <- P
  ecliptic_pole[2] <- -Q*C - Z*S
  ecliptic_pole[3] <- -Q*S + Z*C
  
  # Sum the periodic terms for XA and YA
  sum_xa <- 0.0
  sum_ya <- 0.0
  for(i in 1:14) {
    sum_xa <- sum_xa + XA_C[i]*cos(PI2*T/P2[i]) + XA_S[i]*sin(PI2*T/P2[i])
    sum_ya <- sum_ya + YA_C[i]*cos(PI2*T/P2[i]) + YA_S[i]*sin(PI2*T/P2[i])
  }
  
  # Calculate XA and YA
  XA <- 5453.282155 + 0.4252841 * T - 0.00037173 * T2 - 152E-9 * T3 + sum_xa
  YA <- -73750.930350 - 0.7675452 * T - 0.00018725 * T2 + 231E-9 * T3 + sum_ya
  
  # Convert to radians
  X <- XA * ARCSEC2RAD
  Y <- YA * ARCSEC2RAD
  
  # Calculate the equator pole vector
  equator_pole <- c(0.0, 0.0, 0.0)
  equator_pole[1] <- X
  equator_pole[2] <- Y
  W <- X*X + Y*Y
  
  if (W < 1.0){
    equator_pole[3] <- sqrt(1.0 - W)
  } else {
    equator_pole[3] <- 0.0
  }
  
  # Calculate the precession matrix
  prec_mat <- matrix(0.0, nrow=3, ncol=3)
  
  v <- crossProduct(equator_pole, ecliptic_pole)
  normal_vec <- unitVector(v)
  a <- crossProduct(equator_pole, normal_vec)
  
  prec_mat[1,1] <- normal_vec[1]
  prec_mat[1,2] <- normal_vec[2]
  prec_mat[1,3] <- normal_vec[3]
  prec_mat[2,1] <- a[1]
  prec_mat[2,2] <- a[2]
  prec_mat[2,3] <- a[3]
  prec_mat[3,1] <- equator_pole[1]
  prec_mat[3,2] <- equator_pole[2]
  prec_mat[3,3] <- equator_pole[3]
  
  return(prec_mat)
}