# Automatic Bandwidth selection for Bartlett and Parzen Kernel
# Describe in A Two-Stage Plug-In Bandwidth Selection And Its Implementation For Covariance Estimation

.f.hiruk.bandwidth.solve <- function(y, kernel, prewhite){

  bandwidth.solve = function(y, alpha ,q, kq, intk2, intx2k2,st0, kernel) {
    n   = length(y)
    bt  = ((alpha^2)*intk2*st0^(2*q+1)/((2*q+1)*intx2k2))^(1/(4*q+1)) # bandwidth
    rn  = kernel(y,bt,q) # estimate of spectral density at origin
    rd  = kernel(y,bt,0) # estimate of q'th generalized derivative of spectral density at origin
    r2  = (rn/rd)^2 # squared normalized curvature
    st1 = (q*(kq^2)*r2*n/intk2)^(1/(2*q+1)) # calculate new bandwidth

    #if st1 = st0 optimal bandwidth is found
    out = abs(st0 - st1)
    return(out)
  }

  n = length(y)

  if (isTRUE(prewhite)) {
  # prewithening of time serie
    y = ar(y, aic = FALSE, order.max = 1)$resid[2: n]
  }

  #ar coefficient to cacluclate Alpha(q)
  arCoefficient = ar(y, aic = FALSE, order.max = 1)$ar[1]

  if (kernel == "Bartlett") {
    q       = 1 #characteristic exponent of Bartlett kernel
    k1      = 1 #generalized derivative of Bartlett kernel at origin
    intk2   = 2/3 #squared integral of Bartlett kernel
    intx2k2 = 1/15 #2th-order moment" of Bartlett kernel
    kernel  = f.bartlett
  } else if (kernel == "Parzen") {
    q       = 2 #characteristic exponent of Parzen kernel
    k1      = 6 #generalized derivative of Parzen kernel at origin
    intk2   = 151/280 #squared integral of Parzen kernel
    intx2k2 = 929/295680 #4th-order moment" of Bartlett kernel
    kernel  = f.parzen

  } else {
    stop("only Bartlett, Parzen are supported")
  }

  #proxy for Alpha(q)
  if (q == 1) {
    alpha = (arCoefficient ^ 2  + 1) / (arCoefficient ^ 2 - 1)

  } else if (q == 2) {
    alpha = -(arCoefficient ^ 2 + 8 * arCoefficient + 1) / (arCoefficient - 1) ^ 2
  }

  #Function to optimize
  FUN = function(x) bandwidth.solve(y, alpha ,q, k1, intk2, intx2k2, x, kernel)
  #optimization
  optimalBandwidth = optimize(f = FUN, lower = 4, upper =  n)
  out = optimalBandwidth$minimum
  return(out)
}
f.hiruk.bandwidth.solve = compiler::cmpfun(.f.hiruk.bandwidth.solve)
