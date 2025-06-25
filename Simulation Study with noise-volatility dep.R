##########################################################################################
##########################################################################################
##########################################################################################
#################Replication Code for the article: Part with leverage##############ä##########################
###"Dynamically Consistent Analysis of Realized Covariations in Term Structure Models"####
############by Dennis Schroers############################################################
##########################################################################################
##########################################################################################




###required packages
library(MASS)
library(matrixcalc)
library(pracma)
library(npreg)
library(splines)


########################################################################
########################################################################
########################functions#######################################
########################################################################
########################################################################

{
  #The next functions are needed to replicate the truncated variation
  #function from the FDATSM package (https://github.com/dschroers/FDATSM).
  #Instead of bond prices, the data input here are log prices, as in the
  #simulation study in the original article. Moreover, no artificial date
  #column is attached to the simulated data and the normalization of the
  # estimator is by (1/Delta)^2=100^2 instead of (1/Delta)^2=365^2 like
  # in the FDATSM package.
  #
  #For further descriptions of these functions, consult the documentation
  #of the FDATSM package.

  tc<-function (x, tq = 0.75, l = 3, sparse = TRUE)
  {
    n <- nrow(x)
    m <- ncol(x)

    adj_inc <- matrix(0, n - 1, m - 2)
    for (i in 1:(n - 1)) {
      adj_inc[i, ] <- diff(x[i + 1, 1:(m - 1)]) - diff(x[i, 2:m])
    }
    if (sparse) {
      bspline_result <- bspline_smooth(adj_inc)
      adj_data <- bspline_result$coeffs
      B <- bspline_result$B_grid
      x_grid = bspline_result$x_grid
    }
    else {
      adj_data <- adj_inc
      B <- NULL
    }

    if(l<1000){
      C.Prel <- preliminary_covariation_estimator(data = adj_data,tq)
      ft <- functional_data_truncation(d = d_star(C = C.Prel, tq),C = C.Prel, data = adj_data, Delta = 1/n, sd = l)
      locs <- ft$locations
      if (length(locs) == 0) {
        C.final <- variation(adj_data)
      }
      else {
        C.final <- variation(adj_data[-locs, ])
      }
      if (!is.null(B)) {
        Truncated.variation <- B %*% C.final %*% t(B)
      }
      else {
        Truncated.variation <- C.final
      }}else{
        C.final <- variation(adj_data)
        if (!is.null(B)) {
          Truncated.variation <- B %*% C.final %*% t(B)
        }
        else {
          Truncated.variation <- C.final
        }
      }#these cases are not necessary, but they reduce the simulation time
    Truncated.variation <- Truncated.variation * (100^2)
    eig <- eigen(Truncated.variation, only.values = TRUE)
    expl.var <- cumsum(eig$values)/sum(eig$values)
    return(list(IV = Truncated.variation, expl.var = expl.var))
  }

  d_star<-function(C,rho =.9){
    eigvals <- eigen(C, only.values = TRUE)$values
    scores <- cumsum(eigvals) / sum(eigvals)
    return(min(which(scores >= rho)))
  }

  functional_data_truncation<-function(d, ### number of inverted eigenvalues
                                       C, ### stat. covariance matrix
                                       data, ### data matrix (rows = time)
                                       Delta, ###increment size
                                       sd ### number of 'standard deviations' truncated
  ){
    E<-eigen(C)

    n<-nrow(data)

    #Calculate the values of the truncation function
    gn<-numeric(n)
    for (i in 1:n) {
      x <- data[i, ]

      # First part of g_n (sum of first 'd' components)
      g1 <- sum((x %*% E$vectors[, 1:d])^2 / E$values[1:d])

      # Second part of g_n (remaining components, if any)
      idx <- seq(d + 1, length(E$values))
      g2_num <- sum((x %*% E$vectors[, idx])^2)
      g2_den <- sum(E$values[idx])
      g2 <- g2_num / g2_den

      #take the squareroot of their sum
      gn[i] <- sqrt(g1 + g2)
    }

    #Identify truncation locations
    truncation.locations<-which(gn > (sd*sqrt(d+1)*Delta^(0.49)))
    return(list("gn" = gn , "locations" = truncation.locations))
  }

  preliminary_covariation_estimator <- function(data, tq) {
    n= nrow(data)
    rough.locs <- quantile_truncation(data, tq)
    C.Prel <- variation(data[-rough.locs, ,drop = FALSE])
    EG <- eigen(C.Prel)

    values <- sapply(1:(n - 1), function(i) t(EG$vectors[, 1]) %*% data[i, ])
    q.75 <- quantile(values, 0.75)
    q.25 <- quantile(values, 0.25)

    rho.star <- ((q.75 - q.25)^2) / ((4 * (qnorm(0.75)^2) * EG$values[1]) / n)
    C.Prel <- rho.star * C.Prel

    return(C.Prel = C.Prel)
  }

  quantile_truncation<-function(data,q){

    # Compute Euclidean norm for each row
    norms.data <- apply(data, 1, function(row) sqrt(sum(row^2)))

    # Determine quantile threshold
    quant<- quantile(norms.data, probs = q)

    # Identify indices of points exceeding the threshold
    truncation.locations<-which(norms.data>=quant)
    return(truncation.locations)
  }

  variation <- function(Incr){
    return(t(Incr) %*% Incr)
  }

  bspline_smooth <- function(adj_increments) {
    m=ncol(adj_increments)+2
    x_grid <- seq(1/365, m / 365, length.out = m - 1)
    x_centers <- (x_grid[-1] + x_grid[-length(x_grid)]) / 2

    n_knots <- min(50, m - 2)

    # Define boundaries explicitly to avoid extrapolation issues
    boundary_knots <- range(x_grid)

    # Construct spline basis at x_centers
    B_centers <- bs(x_centers, df = n_knots, degree = 3, intercept = TRUE,
                    Boundary.knots = boundary_knots)
    B_centers_mat <- as.matrix(B_centers)

    # Fit coefficients
    coeffs <- t(apply(adj_increments, 1, function(row) lm.fit(B_centers_mat, row)$coefficients))

    # Evaluate the same basis at x_grid for reconstruction
    B_grid <- as.matrix(predict(B_centers, x_grid))

    return(list(coeffs = coeffs, B_centers = B_centers_mat, B_grid = B_grid, x_grid = x_grid, x_centers = x_centers))
  }

  L2_HS_norm <- function(M, from = 0, to = 1) {
    delta <- (to - from) / nrow(M)
    hilbert.schmidt.norm(M)* delta
  }

  L2_norm <- function(v, from = 0, to = 1) {
    n <- length(v)
    delta <- (to - from) / n
    sqrt(sum(v^2) * delta)
  }


}#from the FDATSM package
{
  CIR.SIM.EULER <- function(kappa = 1.5, mu = 0.58, b = 0.5, x0 = 1,START = 0, END = 1, DELTA = 0.001) {

    time <- seq(START, END, by = DELTA)
    N <- length(time)
    X <- numeric(N)
    X[1] <- x0

    dW <- rnorm(N - 1, mean = 0, sd = sqrt(DELTA))  # Precompute Brownian increments

    for (i in 1:(N - 1)) {
      X[i + 1] <- X[i] + kappa * (mu - X[i]) * DELTA + b * sqrt(pmax(X[i], 0)) * dW[i]
      X[i + 1] <- pmax(X[i + 1], 0)  # Enforce non-negativity
    }

    return(list("Y" = X, "T" = time, "dW" = dW))
  }
  #simulates from a CIR process by an Euler-Maruyama scheme from START to END with resolution DELTA


  ###Simulator--main simulation function
  Simulator_with_leverage<-function(kappa = 1.5, #Reversion rate for scalar volatility process
                      mu= 0.058, #Level for scalar volatility process
                      b= 0.05, #Volatility of scalar volatility process
                      x0= 0.058, #Initial volatility for scalar volatility process
                      lambda.1=1, #Jump intensity, first jumps process
                      lambda.2=4, #Jump intensity, second jumps process
                      Q=q.100, #Noise covariance of the continuous driver
                      rho.1=0.02, #Jump hight, first jump process
                      rho.2 =0.01, #Jump hight, second jumps process
                      f0=numeric(1000), #Initial forward curve
                      corr = -.5 #leverage
  ){
    {
      Vol.driver.integrals<-numeric(n)
      DW<-numeric(n)
      CIR <- CIR.SIM.EULER(kappa = kappa, mu= mu, b= b, x0= x0, START = 0, END = 1, DELTA=1/(n*100))




      for (i in 2:n) {
        Vol.driver.integrals[i] <- sum(CIR$Y[((i-1)*100):(i*100)]^2)/(n*100)
        DW[i]<- sum(CIR$dW[((i-1)*100):(i*100)])
      }#calculate integrated volatilities

      kernel.samples <- mvrnorm(n = n, numeric(M+n), Q, tol = 1e-3)# simulate discretized increments of a Q-Wiener process
      #observe that we have to sample 1100 instead of 1000 noise innovations,
      #to avoid boundary problems at maximal maturity.
      for (i in 1:n) {
        kernel.samples[i,]<-(corr*kernel.samples[i,])+(sqrt(1-corr^2)*DW[i]*10*(numeric(M+n)+1/10)) #/10 since we want to normalize the new covariance of the Wiener process
      }




      samples <- matrix(0,n,(M+n))
      samples[1,] <- f0
      for (i in 2:n) {
        squared.volatility <- Vol.driver.integrals[i]*Q
        samples[i,1:(M+n-i)] <- samples[i-1,2:(M+1+n-i)]
        samples[i,] <- samples[i,]+(sqrt(Vol.driver.integrals[i]))*kernel.samples[i,]
      }
    }#simulate continuous part
    {
      Int.Arr.times.1 <- rexp(n=100,rate=lambda.1)
      Arr.times.1 <- numeric(100)
      for (i in 1:100) {
        Arr.times.1[i] <- sum(Int.Arr.times.1[1:i])
      }
      Arr.times.1 <- Arr.times.1[which(Arr.times.1<1)]


      Int.Arr.times.2 <- rexp(n=100,rate=lambda.2)
      Arr.times.2 <- numeric(100)
      for (i in 1:100) {
        Arr.times.2[i] <- sum(Int.Arr.times.2[1:i])
      }
      Arr.times.2 <- Arr.times.2[which(Arr.times.2<1)]


      jump.locations.1 <- trunc(Arr.times.1*n)
      jump.locations.2 <- trunc(Arr.times.2*n)

      Jump.number.1 <- length(Arr.times.1)
      Jump.number.2 <- length(Arr.times.2)

      ##simulate stationary CPP
      CPP.1 <- matrix(0,n,M+n)
      if(length(Arr.times.1) > 0){
        Chi.1<-matrix(0,length(Arr.times.1),(M+n))
        for (i in 1:length(Arr.times.1)) {
          Chi.1[i,]<-as.numeric(mvrnorm(n = 1, numeric(M+n), rho.1*q.0.01, tol = 1e-3))
        }



        if(jump.locations.1[1]==1){
          CPP.1[1,]<-Chi.1[1,]
        }
        for (i in 2:n) {
          CPP.1[i,]<-CPP.1[i-1,]
          if(is.element(i,jump.locations.1)){
            for (j in which(i==jump.locations.1)) {
              CPP.1[i,]<- CPP.1[i,]+Chi.1[j,]
            }
          }
        }
      }

      ###simulate nonstationary CPP
      CPP.2 <- matrix(0,n,M+n)
      if(length(Arr.times.2) > 0){
        Chi.2<-matrix(0,length(Arr.times.2),(M+n))
        for (i in 1:length(Arr.times.2)) {
          Chi.2[i,]<-as.numeric(mvrnorm(n = 1, numeric(M+n), rho.2*k, tol = 1e-3))
        }


        if(jump.locations.2[1]==1){
          CPP.2[1,]<-Chi.2[1,]
        }
        for (i in 2:n) {
          CPP.2[i,]<-CPP.2[i-1,]
          if(is.element(i,jump.locations.2)){
            for (j in which(i==jump.locations.2)) {
              CPP.2[i,]<- CPP.2[i,]+exp(-10*(((jump.locations.2[j]+1)/n)-Arr.times.2[j]))*Chi.2[j,]
            }
          }
        }
      }


      jump.samples <- CPP.1+CPP.2

    }##simulate the jump part

    cont.data <- samples[,1:(M+1)] #difference return data without jumps
    data <- samples[,1:(M+1)]+jump.samples[,1:(M+1)] #difference return data with jumps


    #transform from difference return data to log bond price data
    Cont.Price.Data<-matrix(0,n,(M+1))
    Price.Data<-matrix(0,n,(M+1))
    for (i in 1:n) {
      Cont.Price.Data[i,]<- -cumsum(cont.data[i,])[1:(M+1)]
      Price.Data[i,]<- -cumsum(data[i,])[1:(M+1)]
    }

    #Calculate also the simulated quadratic variation/integrated volatility
    v<- (numeric(999)+1/10)
    target_cov <- (corr^2) * Q[1:999, 1:999] + (1 - corr^2) * (v%*%t(v))


    IV <- sum(Vol.driver.integrals)* target_cov


        return(list("Prices" = Price.Data, "Prices.cont" = Cont.Price.Data, "IV"= IV ))
  }

}#for the simulation


##Metaparameters & initializations
n<-100 #Number of simulated days
M<- 1000 #Number of simulated maturity grid points
noise.level<-.01 #Variance of the noise
sparseness<-100 #Number of observed supsamples in the sparse sampling case
Initial.fw<-numeric(M+n) #Create a discretized initial forward curve
probs<-(numeric(M)+1)/M #Select probabilities for the irregular sampling points
sparsing = FALSE #make TRUE to reduce simulation time, by making enforcing ex-post smoothing of difference return curves

K=500 #Number of Monte-Carlo runs


#Calculate the covariance kernels for volatility and jumps sizes
{
  q.50<-matrix(0, M+n,M+n)
  q.0.01<-matrix(0, M+n,M+n)
  k<-matrix(0,M+n,M+n)


  integral.Gaussian.kernel<- function(i,j,alpha = 0.1,Delta = 1/n){# Calculates < Q_2 1_{[(j-1)Delta,jDelta],1_{[(j-1)Delta,jDelta]>}
    fun <- function(x,y){
      exp(-alpha*(x-y)^{2})
    }
    integral2(fun, (i-1)*Delta, i*Delta,  (j-1)*Delta, j*Delta, reltol = 1e-10)$Q
  }
  #calculates the integral of f(x,y)=exp(-alpha*(x-y)^{2}) over the interval
  #[(i-1)*Delta,i*Delta]x[(j-1)*Delta,j*Delta]


  for (i in 1:(M+n)) {
    for (j in 1:i) {
      q.50[i,j]<-integral.Gaussian.kernel(i,j, alpha = 50, Delta= 1/M)
      q.50[j,i]<- q.50[i,j]
    }
    print(i/(3*(M+n)))
  }#Calculate the local averaged Gaussian kernel with parameter alpha=50
  q.50<-q.50/ L2_HS_norm(q.50, from = 0, to = 11) #normalize the local averaged kernel

  for (i in 1:(M+n)) {
    for (j in 1:i) {
      q.0.01[i,j]<-integral.Gaussian.kernel(i,j, alpha = .01, Delta= 1/M)
      q.0.01[j,i]<- q.0.01[i,j]
    }
    print(i/(3*(M+n))+(1/3))
  }#Calculate the local averaged Gaussian kernel with parameter alpha=0.01
  q.0.01<-q.0.01/ L2_HS_norm(q.0.01, from = 0, to = 11) #normalize the local averaged kernel


  for (i in 1:(M+n)) {
    for (j in 1:i) {
      k[i,j]<-((1-exp(-1/n))^2/100)*exp(-(j+i-2)/n)
      k[j,i]<- k[i,j]
    }
    print(i/(3*(M+n))+2/3)
  }##Calculate the local average-discretized kernel k(x,y)=exp(y)exp(x)
  k<-k/ L2_HS_norm(k, from = 0, to = 11) #normalize the local averaged kernel
}



########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################



{
  #norms of the true quadratic variations
  normIV<-numeric(K)

  #S1
  {



    M1.L00.Err<-numeric(K)
    M2.L00.Err<-numeric(K)
    M3.L3.Err<-numeric(K)
    M3.L4.Err<-numeric(K)
    M3.L5.Err<-numeric(K)
    M4.L3.Err<-numeric(K)
    M4.L4.Err<-numeric(K)
    M4.L5.Err<-numeric(K)

    M1.L00.dimensionality.85<-numeric(K)
    M2.L00.dimensionality.85<-numeric(K)
    M3.L3.dimensionality.85<-numeric(K)
    M3.L4.dimensionality.85<-numeric(K)
    M3.L5.dimensionality.85<-numeric(K)
    M4.L3.dimensionality.85<-numeric(K)
    M4.L4.dimensionality.85<-numeric(K)
    M4.L5.dimensionality.85<-numeric(K)

    M1.L00.dimensionality.90<-numeric(K)
    M2.L00.dimensionality.90<-numeric(K)
    M3.L3.dimensionality.90<-numeric(K)
    M3.L4.dimensionality.90<-numeric(K)
    M3.L5.dimensionality.90<-numeric(K)
    M4.L3.dimensionality.90<-numeric(K)
    M4.L4.dimensionality.90<-numeric(K)
    M4.L5.dimensionality.90<-numeric(K)

    M1.L00.dimensionality.95<-numeric(K)
    M2.L00.dimensionality.95<-numeric(K)
    M3.L3.dimensionality.95<-numeric(K)
    M3.L4.dimensionality.95<-numeric(K)
    M3.L5.dimensionality.95<-numeric(K)
    M4.L3.dimensionality.95<-numeric(K)
    M4.L4.dimensionality.95<-numeric(K)
    M4.L5.dimensionality.95<-numeric(K)

    M1.L00.dimensionality.99<-numeric(K)
    M2.L00.dimensionality.99<-numeric(K)
    M3.L3.dimensionality.99<-numeric(K)
    M3.L4.dimensionality.99<-numeric(K)
    M3.L5.dimensionality.99<-numeric(K)
    M4.L3.dimensionality.99<-numeric(K)
    M4.L4.dimensionality.99<-numeric(K)
    M4.L5.dimensionality.99<-numeric(K)

  }

  #S2
  {
    Proj.M1.L00.Err<-numeric(K)
    Proj.M2.L00.Err<-numeric(K)
    Proj.M3.L3.Err<-numeric(K)
    Proj.M3.L4.Err<-numeric(K)
    Proj.M3.L5.Err<-numeric(K)
    Proj.M4.L3.Err<-numeric(K)
    Proj.M4.L4.Err<-numeric(K)
    Proj.M4.L5.Err<-numeric(K)


    Proj.M1.L00.dimensionality.85<-numeric(K)
    Proj.M2.L00.dimensionality.85<-numeric(K)
    Proj.M3.L3.dimensionality.85<-numeric(K)
    Proj.M3.L4.dimensionality.85<-numeric(K)
    Proj.M3.L5.dimensionality.85<-numeric(K)
    Proj.M4.L3.dimensionality.85<-numeric(K)
    Proj.M4.L4.dimensionality.85<-numeric(K)
    Proj.M4.L5.dimensionality.85<-numeric(K)

    Proj.M1.L00.dimensionality.90<-numeric(K)
    Proj.M2.L00.dimensionality.90<-numeric(K)
    Proj.M3.L3.dimensionality.90<-numeric(K)
    Proj.M3.L4.dimensionality.90<-numeric(K)
    Proj.M3.L5.dimensionality.90<-numeric(K)
    Proj.M4.L3.dimensionality.90<-numeric(K)
    Proj.M4.L4.dimensionality.90<-numeric(K)
    Proj.M4.L5.dimensionality.90<-numeric(K)

    Proj.M1.L00.dimensionality.95<-numeric(K)
    Proj.M2.L00.dimensionality.95<-numeric(K)
    Proj.M3.L3.dimensionality.95<-numeric(K)
    Proj.M3.L4.dimensionality.95<-numeric(K)
    Proj.M3.L5.dimensionality.95<-numeric(K)
    Proj.M4.L3.dimensionality.95<-numeric(K)
    Proj.M4.L4.dimensionality.95<-numeric(K)
    Proj.M4.L5.dimensionality.95<-numeric(K)

    Proj.M1.L00.dimensionality.99<-numeric(K)
    Proj.M2.L00.dimensionality.99<-numeric(K)
    Proj.M3.L3.dimensionality.99<-numeric(K)
    Proj.M3.L4.dimensionality.99<-numeric(K)
    Proj.M3.L5.dimensionality.99<-numeric(K)
    Proj.M4.L3.dimensionality.99<-numeric(K)
    Proj.M4.L4.dimensionality.99<-numeric(K)
    Proj.M4.L5.dimensionality.99<-numeric(K)

    logdiff.explavar1 <- matrix(0,K,10)
    logdiff.explavar2 <- matrix(0,K,10)
    logdiff.explavar3 <- matrix(0,K,10)
    logdiff.explavar4 <- matrix(0,K,10)

  }

  #S3
  {
    Cov.Proj.M1.L00.Err<-numeric(K)
    Cov.Proj.M2.L00.Err<-numeric(K)
    Cov.Proj.M3.L3.Err<-numeric(K)
    Cov.Proj.M3.L4.Err<-numeric(K)
    Cov.Proj.M3.L5.Err<-numeric(K)
    Cov.Proj.M4.L3.Err<-numeric(K)
    Cov.Proj.M4.L4.Err<-numeric(K)
    Cov.Proj.M4.L5.Err<-numeric(K)


    Cov.Proj.M1.L00.dimensionality.85<-numeric(K)
    Cov.Proj.M2.L00.dimensionality.85<-numeric(K)
    Cov.Proj.M3.L3.dimensionality.85<-numeric(K)
    Cov.Proj.M3.L4.dimensionality.85<-numeric(K)
    Cov.Proj.M3.L5.dimensionality.85<-numeric(K)
    Cov.Proj.M4.L3.dimensionality.85<-numeric(K)
    Cov.Proj.M4.L4.dimensionality.85<-numeric(K)
    Cov.Proj.M4.L5.dimensionality.85<-numeric(K)

    Cov.Proj.M1.L00.dimensionality.90<-numeric(K)
    Cov.Proj.M2.L00.dimensionality.90<-numeric(K)
    Cov.Proj.M3.L3.dimensionality.90<-numeric(K)
    Cov.Proj.M3.L4.dimensionality.90<-numeric(K)
    Cov.Proj.M3.L5.dimensionality.90<-numeric(K)
    Cov.Proj.M4.L3.dimensionality.90<-numeric(K)
    Cov.Proj.M4.L4.dimensionality.90<-numeric(K)
    Cov.Proj.M4.L5.dimensionality.90<-numeric(K)

    Cov.Proj.M1.L00.dimensionality.95<-numeric(K)
    Cov.Proj.M2.L00.dimensionality.95<-numeric(K)
    Cov.Proj.M3.L3.dimensionality.95<-numeric(K)
    Cov.Proj.M3.L4.dimensionality.95<-numeric(K)
    Cov.Proj.M3.L5.dimensionality.95<-numeric(K)
    Cov.Proj.M4.L3.dimensionality.95<-numeric(K)
    Cov.Proj.M4.L4.dimensionality.95<-numeric(K)
    Cov.Proj.M4.L5.dimensionality.95<-numeric(K)

    Cov.Proj.M1.L00.dimensionality.99<-numeric(K)
    Cov.Proj.M2.L00.dimensionality.99<-numeric(K)
    Cov.Proj.M3.L3.dimensionality.99<-numeric(K)
    Cov.Proj.M3.L4.dimensionality.99<-numeric(K)
    Cov.Proj.M3.L5.dimensionality.99<-numeric(K)
    Cov.Proj.M4.L3.dimensionality.99<-numeric(K)
    Cov.Proj.M4.L4.dimensionality.99<-numeric(K)
    Cov.Proj.M4.L5.dimensionality.99<-numeric(K)

    Cov.explavar1 <- matrix(0,K,10)
    Cov.explavar2 <- matrix(0,K,10)
    Cov.explavar3 <- matrix(0,K,10)
    Cov.explavar4 <- matrix(0,K,10)

  }

  #S4
  {
    excess.Proj.M1.L00.Err<-numeric(K)
    excess.Proj.M2.L00.Err<-numeric(K)
    excess.Proj.M3.L3.Err<-numeric(K)
    excess.Proj.M3.L4.Err<-numeric(K)
    excess.Proj.M3.L5.Err<-numeric(K)
    excess.Proj.M4.L3.Err<-numeric(K)
    excess.Proj.M4.L4.Err<-numeric(K)
    excess.Proj.M4.L5.Err<-numeric(K)


    excess.Proj.M1.L00.dimensionality.85<-numeric(K)
    excess.Proj.M2.L00.dimensionality.85<-numeric(K)
    excess.Proj.M3.L3.dimensionality.85<-numeric(K)
    excess.Proj.M3.L4.dimensionality.85<-numeric(K)
    excess.Proj.M3.L5.dimensionality.85<-numeric(K)
    excess.Proj.M4.L3.dimensionality.85<-numeric(K)
    excess.Proj.M4.L4.dimensionality.85<-numeric(K)
    excess.Proj.M4.L5.dimensionality.85<-numeric(K)

    excess.Proj.M1.L00.dimensionality.90<-numeric(K)
    excess.Proj.M2.L00.dimensionality.90<-numeric(K)
    excess.Proj.M3.L3.dimensionality.90<-numeric(K)
    excess.Proj.M3.L4.dimensionality.90<-numeric(K)
    excess.Proj.M3.L5.dimensionality.90<-numeric(K)
    excess.Proj.M4.L3.dimensionality.90<-numeric(K)
    excess.Proj.M4.L4.dimensionality.90<-numeric(K)
    excess.Proj.M4.L5.dimensionality.90<-numeric(K)

    excess.Proj.M1.L00.dimensionality.95<-numeric(K)
    excess.Proj.M2.L00.dimensionality.95<-numeric(K)
    excess.Proj.M3.L3.dimensionality.95<-numeric(K)
    excess.Proj.M3.L4.dimensionality.95<-numeric(K)
    excess.Proj.M3.L5.dimensionality.95<-numeric(K)
    excess.Proj.M4.L3.dimensionality.95<-numeric(K)
    excess.Proj.M4.L4.dimensionality.95<-numeric(K)
    excess.Proj.M4.L5.dimensionality.95<-numeric(K)

    excess.Proj.M1.L00.dimensionality.99<-numeric(K)
    excess.Proj.M2.L00.dimensionality.99<-numeric(K)
    excess.Proj.M3.L3.dimensionality.99<-numeric(K)
    excess.Proj.M3.L4.dimensionality.99<-numeric(K)
    excess.Proj.M3.L5.dimensionality.99<-numeric(K)
    excess.Proj.M4.L3.dimensionality.99<-numeric(K)
    excess.Proj.M4.L4.dimensionality.99<-numeric(K)
    excess.Proj.M4.L5.dimensionality.99<-numeric(K)


    excess.explavar1 <- matrix(0, K , 10)
    excess.explavar2 <- matrix(0, K , 10)
    excess.explavar3 <- matrix(0, K , 10)
    excess.explavar4 <- matrix(0, K , 10)
  }
}#define the Monte-Carlo sample vectors

set.seed(123)#to reproduce the results
ptm <-proc.time()
for (step in 1:K) {
  DATA <- Simulator_with_leverage(kappa = 1.5, mu= .058, b= .05, x0= 0.058,Q=q.50, lambda.1=1, lambda.2=4, rho.1=0.0116, rho.2 =0.0029, f0=Initial.fw)



  M1.price.data <- DATA$Prices.cont
  M2.price.data <- matrix(0,n,M+1)
  M3.price.data <- DATA$Prices
  M4.price.data <- matrix(0,n,M+1)


  {
    sparseness = 100

    X.Sparse <- matrix(0,n,sparseness)
    M1.price.data.Sparse <- matrix(0,n,sparseness)
    M3.price.data.Sparse <- matrix(0,n,sparseness)
    for (z in 1:n) {
      Errors <- rnorm(sparseness, 0, noise.level)
      X.Sparse[z,] <- sort(sample(1:M, size = sparseness, prob = probs))
      M1.price.data.Sparse[z,] <- M1.price.data[z,X.Sparse[z,]]+Errors
      M3.price.data.Sparse[z,] <- M3.price.data[z,X.Sparse[z,]]+Errors
    }###Derive sparse and noisy samples (100 out of 1000 points)

  }#Calculate the sparse & noisy price data
  for (i in 1:n) {
    SS.1 <- ss(x = X.Sparse[i,], y = M1.price.data.Sparse[i,], method = "BIC", m=3)
    M2.price.data[i,] <- predict(SS.1, 1:(M+1))$y
    SS.2 <- ss(x = X.Sparse[i,], y = M3.price.data.Sparse[i,], method = "BIC", m=3)
    M4.price.data[i,] <- predict(SS.2, 1:(M+1))$y
  }#Calculate the smoothed log price data

  Integrated.Volatility <- DATA$IV*(100^2)#Calculate the true quadratic covariation
  normIV[step] <- L2_HS_norm(Integrated.Volatility, from = 0, to = 10)#Calculate the norm of  true quadratic covariation
  # For Scenario S1: (unprojected data)
  {
    {
      M1.L00 <- tc(x= M1.price.data, sparse = sparsing, tq = 0.75,l = 10000)

      M2.L00 <- tc(x= M2.price.data, sparse = sparsing, tq = 0.75,l = 10000)

      M3.L3 <- tc(x= M3.price.data, sparse = sparsing, tq = 0.75,l = 3)
      M3.L4 <- tc(x= M3.price.data, sparse = sparsing, tq = 0.75,l = 4)
      M3.L5 <- tc(x= M3.price.data, sparse = sparsing, tq = 0.75,l = 5)

      M4.L3 <- tc(x= M4.price.data, sparse = sparsing, tq = 0.75,l = 3)
      M4.L4 <- tc(x= M4.price.data, sparse = sparsing, tq = 0.75,l = 4)
      M4.L5 <- tc(x= M4.price.data, sparse = sparsing, tq = 0.75,l = 5)
    }#calculate the estimators


    {
      {
        M1.L00.loadings <- M1.L00$expl.var[1:20]
        M2.L00.loadings <- M2.L00$expl.var[1:20]
        M3.L3.loadings <- M3.L3$expl.var[1:20]
        M3.L4.loadings <- M3.L4$expl.var[1:20]
        M3.L5.loadings <- M3.L5$expl.var[1:20]
        M4.L3.loadings <- M4.L3$expl.var[1:20]
        M4.L4.loadings <- M4.L4$expl.var[1:20]
        M4.L5.loadings <- M4.L5$expl.var[1:20]

        M1.L00.dimensionality.85[step] <- min(which(M1.L00.loadings > .85))
        M2.L00.dimensionality.85[step] <- min(which(M2.L00.loadings > .85))
        M3.L3.dimensionality.85[step] <- min(which(M3.L3.loadings > .85))
        M3.L4.dimensionality.85[step] <- min(which(M3.L4.loadings > .85))
        M3.L5.dimensionality.85[step] <- min(which(M3.L5.loadings > .85))
        M4.L3.dimensionality.85[step] <- min(which(M4.L3.loadings > .85))
        M4.L4.dimensionality.85[step] <- min(which(M4.L4.loadings > .85))
        M4.L5.dimensionality.85[step] <- min(which(M4.L5.loadings > .85))

        M1.L00.dimensionality.90[step] <- min(which(M1.L00.loadings > .9))
        M2.L00.dimensionality.90[step] <- min(which(M2.L00.loadings > .9))
        M3.L3.dimensionality.90[step] <- min(which(M3.L3.loadings > .9))
        M3.L4.dimensionality.90[step] <- min(which(M3.L4.loadings > .9))
        M3.L5.dimensionality.90[step] <- min(which(M3.L5.loadings > .9))
        M4.L3.dimensionality.90[step] <- min(which(M4.L3.loadings > .9))
        M4.L4.dimensionality.90[step] <- min(which(M4.L4.loadings > .9))
        M4.L5.dimensionality.90[step] <- min(which(M4.L5.loadings > .9))

        M1.L00.dimensionality.95[step] <- min(which(M1.L00.loadings > .95))
        M2.L00.dimensionality.95[step] <- min(which(M2.L00.loadings > .95))
        M3.L3.dimensionality.95[step] <- min(which(M3.L3.loadings > .95))
        M3.L4.dimensionality.95[step] <- min(which(M3.L4.loadings > .95))
        M3.L5.dimensionality.95[step] <- min(which(M3.L5.loadings > .95))
        M4.L3.dimensionality.95[step] <- min(which(M4.L3.loadings > .95))
        M4.L4.dimensionality.95[step] <- min(which(M4.L4.loadings > .95))
        M4.L5.dimensionality.95[step] <- min(which(M4.L5.loadings > .95))

        M1.L00.dimensionality.99[step] <- min(which(M1.L00.loadings > .99))
        M2.L00.dimensionality.99[step] <- min(which(M2.L00.loadings > .99))
        M3.L3.dimensionality.99[step] <- min(which(M3.L3.loadings > .99))
        M3.L4.dimensionality.99[step] <- min(which(M3.L4.loadings > .99))
        M3.L5.dimensionality.99[step] <- min(which(M3.L5.loadings > .99))
        M4.L3.dimensionality.99[step] <- min(which(M4.L3.loadings > .99))
        M4.L4.dimensionality.99[step] <- min(which(M4.L4.loadings > .99))
        M4.L5.dimensionality.99[step] <- min(which(M4.L5.loadings > .99))
      }#dimensionalities



      {
        u<-min(ncol(M1.L00$IV),ncol(Integrated.Volatility))

        M1.L00.Err[step]<-L2_HS_norm(M1.L00$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        M2.L00.Err[step]<-L2_HS_norm(M2.L00$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        M3.L3.Err[step]<-L2_HS_norm(M3.L3$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        M3.L4.Err[step]<-L2_HS_norm(M3.L4$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        M3.L5.Err[step]<-L2_HS_norm(M3.L5$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        M4.L3.Err[step]<-L2_HS_norm(M4.L3$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        M4.L4.Err[step]<-L2_HS_norm(M4.L4$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        M4.L5.Err[step]<-L2_HS_norm(M4.L5$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
      }#errors

    }#Calculate the errors and dimensions
  }
  # For Scenario S2: (data projected onto time-differenced log-price PCAs)
  {
    #make projections
    {
      E1 <- eigen(t(diff(M1.price.data))%*%diff(M1.price.data))
      E2 <- eigen(t(diff(M2.price.data))%*%diff(M2.price.data))
      E3 <- eigen(t(diff(M3.price.data))%*%diff(M3.price.data))
      E4 <- eigen(t(diff(M4.price.data))%*%diff(M4.price.data))

      logdiff.explavar1[step,] <- (cumsum(E1$values)/sum(E1$values))[1:10]
      logdiff.explavar2[step,] <- (cumsum(E2$values)/sum(E2$values))[1:10]
      logdiff.explavar3[step,] <- (cumsum(E3$values)/sum(E3$values))[1:10]
      logdiff.explavar4[step,] <- (cumsum(E4$values)/sum(E4$values))[1:10]

      d1 <- which.max(logdiff.explavar1[step,]>.99)
      d2 <- which.max(logdiff.explavar2[step,]>.99)
      d3 <- which.max(logdiff.explavar3[step,]>.99)
      d4 <- which.max(logdiff.explavar4[step,]>.99)

      Proj.M1.price.data <- M1.price.data %*% E1$vectors[,1:d1] %*% t(E1$vectors[,1:d1])
      Proj.M2.price.data <- M2.price.data %*% E2$vectors[,1:d2] %*% t(E2$vectors[,1:d2])
      Proj.M3.price.data <- M3.price.data %*% E3$vectors[,1:d3] %*% t(E3$vectors[,1:d3])
      Proj.M4.price.data <- M4.price.data %*% E4$vectors[,1:d4] %*% t(E4$vectors[,1:d4])
    }


    {
      Proj.M1.L00 <- tc(x= Proj.M1.price.data, sparse = sparsing, tq = 0.75,l = 10000)

      Proj.M2.L00 <- tc(x= Proj.M2.price.data, sparse = sparsing, tq = 0.75,l = 10000)

      Proj.M3.L3 <- tc(x= Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 3)
      Proj.M3.L4 <- tc(x= Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 4)
      Proj.M3.L5 <- tc(x= Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 5)

      Proj.M4.L3 <- tc(x= Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 3)
      Proj.M4.L4 <- tc(x= Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 4)
      Proj.M4.L5 <- tc(x= Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 5)

    }#calculate the estimators


    {
      {
        Proj.M1.L00.loadings <- Proj.M1.L00$expl.var[1:20]
        Proj.M2.L00.loadings <- Proj.M2.L00$expl.var[1:20]
        Proj.M3.L3.loadings <- Proj.M3.L3$expl.var[1:20]
        Proj.M3.L4.loadings <- Proj.M3.L4$expl.var[1:20]
        Proj.M3.L5.loadings <- Proj.M3.L5$expl.var[1:20]
        Proj.M4.L3.loadings <- Proj.M4.L3$expl.var[1:20]
        Proj.M4.L4.loadings <- Proj.M4.L4$expl.var[1:20]
        Proj.M4.L5.loadings <- Proj.M4.L5$expl.var[1:20]

        Proj.M1.L00.dimensionality.85[step] <- min(which(Proj.M1.L00.loadings > .85))
        Proj.M2.L00.dimensionality.85[step] <- min(which(Proj.M2.L00.loadings > .85))
        Proj.M3.L3.dimensionality.85[step] <- min(which(Proj.M3.L3.loadings > .85))
        Proj.M3.L4.dimensionality.85[step] <- min(which(Proj.M3.L4.loadings > .85))
        Proj.M3.L5.dimensionality.85[step] <- min(which(Proj.M3.L5.loadings > .85))
        Proj.M4.L3.dimensionality.85[step] <- min(which(Proj.M4.L3.loadings > .85))
        Proj.M4.L4.dimensionality.85[step] <- min(which(Proj.M4.L4.loadings > .85))
        Proj.M4.L5.dimensionality.85[step] <- min(which(Proj.M4.L5.loadings > .85))

        Proj.M1.L00.dimensionality.90[step] <- min(which(Proj.M1.L00.loadings > .9))
        Proj.M2.L00.dimensionality.90[step] <- min(which(Proj.M2.L00.loadings > .9))
        Proj.M3.L3.dimensionality.90[step] <- min(which(Proj.M3.L3.loadings > .9))
        Proj.M3.L4.dimensionality.90[step] <- min(which(Proj.M3.L4.loadings > .9))
        Proj.M3.L5.dimensionality.90[step] <- min(which(Proj.M3.L5.loadings > .9))
        Proj.M4.L3.dimensionality.90[step] <- min(which(Proj.M4.L3.loadings > .9))
        Proj.M4.L4.dimensionality.90[step] <- min(which(Proj.M4.L4.loadings > .9))
        Proj.M4.L5.dimensionality.90[step] <- min(which(Proj.M4.L5.loadings > .9))

        Proj.M1.L00.dimensionality.95[step] <- min(which(Proj.M1.L00.loadings > .95))
        Proj.M2.L00.dimensionality.95[step] <- min(which(Proj.M2.L00.loadings > .95))
        Proj.M3.L3.dimensionality.95[step] <- min(which(Proj.M3.L3.loadings > .95))
        Proj.M3.L4.dimensionality.95[step] <- min(which(Proj.M3.L4.loadings > .95))
        Proj.M3.L5.dimensionality.95[step] <- min(which(Proj.M3.L5.loadings > .95))
        Proj.M4.L3.dimensionality.95[step] <- min(which(Proj.M4.L3.loadings > .95))
        Proj.M4.L4.dimensionality.95[step] <- min(which(Proj.M4.L4.loadings > .95))
        Proj.M4.L5.dimensionality.95[step] <- min(which(Proj.M4.L5.loadings > .95))

        Proj.M1.L00.dimensionality.99[step] <- min(which(Proj.M1.L00.loadings > .99))
        Proj.M2.L00.dimensionality.99[step] <- min(which(Proj.M2.L00.loadings > .99))
        Proj.M3.L3.dimensionality.99[step] <- min(which(Proj.M3.L3.loadings > .99))
        Proj.M3.L4.dimensionality.99[step] <- min(which(Proj.M3.L4.loadings > .99))
        Proj.M3.L5.dimensionality.99[step] <- min(which(Proj.M3.L5.loadings > .99))
        Proj.M4.L3.dimensionality.99[step] <- min(which(Proj.M4.L3.loadings > .99))
        Proj.M4.L4.dimensionality.99[step] <- min(which(Proj.M4.L4.loadings > .99))
        Proj.M4.L5.dimensionality.99[step] <- min(which(Proj.M4.L5.loadings > .99))
      }#dimensionalities


      {

        Proj.M1.L00.Err[step] <- L2_HS_norm(Proj.M1.L00$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        Proj.M2.L00.Err[step] <- L2_HS_norm(Proj.M2.L00$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        Proj.M3.L3.Err[step] <- L2_HS_norm(Proj.M3.L3$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        Proj.M3.L4.Err[step] <- L2_HS_norm(Proj.M3.L4$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        Proj.M3.L5.Err[step] <- L2_HS_norm(Proj.M3.L5$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        Proj.M4.L3.Err[step] <- L2_HS_norm(Proj.M4.L3$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        Proj.M4.L4.Err[step] <- L2_HS_norm(Proj.M4.L4$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        Proj.M4.L5.Err[step] <- L2_HS_norm(Proj.M4.L5$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
      }#errors

    }#Calculate the errors and dimensions
  }
  # For Scenario S3: (data projected onto log-price PCAs)
  {
    #make projections
    {
      E1 <- eigen(t(M1.price.data)%*%M1.price.data)
      E2 <- eigen(t(M2.price.data)%*%M2.price.data)
      E3 <- eigen(t(M3.price.data)%*%M3.price.data)
      E4 <- eigen(t(M4.price.data)%*%M4.price.data)

      Cov.explavar1[step,] <- (cumsum(E1$values)/sum(E1$values))[1:10]
      Cov.explavar2[step,] <- (cumsum(E2$values)/sum(E2$values))[1:10]
      Cov.explavar3[step,] <- (cumsum(E3$values)/sum(E3$values))[1:10]
      Cov.explavar4[step,] <- (cumsum(E4$values)/sum(E4$values))[1:10]

      d1 <- which.max(Cov.explavar1[step,]>.99)
      d2 <- which.max(Cov.explavar2[step,]>.99)
      d3 <- which.max(Cov.explavar3[step,]>.99)
      d4 <- which.max(Cov.explavar4[step,]>.99)

      Cov.Proj.M1.price.data <- M1.price.data %*% E1$vectors[,1:d1] %*% t(E1$vectors[,1:d1])
      Cov.Proj.M2.price.data <- M2.price.data %*% E2$vectors[,1:d2] %*% t(E2$vectors[,1:d2])
      Cov.Proj.M3.price.data <- M3.price.data %*% E3$vectors[,1:d3] %*% t(E3$vectors[,1:d3])
      Cov.Proj.M4.price.data <- M4.price.data %*% E4$vectors[,1:d4] %*% t(E4$vectors[,1:d4])
    }

    {
      Cov.Proj.M1.L00 <- tc(x= Cov.Proj.M1.price.data, sparse = sparsing, tq = 0.75,l = 10000)

      Cov.Proj.M2.L00 <- tc(x= Cov.Proj.M2.price.data, sparse = sparsing, tq = 0.75,l = 10000)

      Cov.Proj.M3.L3 <- tc(x= Cov.Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 3)
      Cov.Proj.M3.L4 <- tc(x= Cov.Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 4)
      Cov.Proj.M3.L5 <- tc(x= Cov.Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 5)

      Cov.Proj.M4.L3 <- tc(x= Cov.Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 3)
      Cov.Proj.M4.L4 <- tc(x= Cov.Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 4)
      Cov.Proj.M4.L5 <- tc(x= Cov.Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 5)

    }#calculate the estimators


    {
      {
        Cov.Proj.M1.L00.loadings <- Cov.Proj.M1.L00$expl.var[1:20]
        Cov.Proj.M2.L00.loadings <- Cov.Proj.M2.L00$expl.var[1:20]
        Cov.Proj.M3.L3.loadings <- Cov.Proj.M3.L3$expl.var[1:20]
        Cov.Proj.M3.L4.loadings <- Cov.Proj.M3.L4$expl.var[1:20]
        Cov.Proj.M3.L5.loadings <- Cov.Proj.M3.L5$expl.var[1:20]
        Cov.Proj.M4.L3.loadings <- Cov.Proj.M4.L3$expl.var[1:20]
        Cov.Proj.M4.L4.loadings <- Cov.Proj.M4.L4$expl.var[1:20]
        Cov.Proj.M4.L5.loadings <- Cov.Proj.M4.L5$expl.var[1:20]

        Cov.Proj.M1.L00.dimensionality.85[step] <- min(which(Cov.Proj.M1.L00.loadings > .85))
        Cov.Proj.M2.L00.dimensionality.85[step] <- min(which(Cov.Proj.M2.L00.loadings > .85))
        Cov.Proj.M3.L3.dimensionality.85[step] <- min(which(Cov.Proj.M3.L3.loadings > .85))
        Cov.Proj.M3.L4.dimensionality.85[step] <- min(which(Cov.Proj.M3.L4.loadings > .85))
        Cov.Proj.M3.L5.dimensionality.85[step] <- min(which(Cov.Proj.M3.L5.loadings > .85))
        Cov.Proj.M4.L3.dimensionality.85[step] <- min(which(Cov.Proj.M4.L3.loadings > .85))
        Cov.Proj.M4.L4.dimensionality.85[step] <- min(which(Cov.Proj.M4.L4.loadings > .85))
        Cov.Proj.M4.L5.dimensionality.85[step] <- min(which(Cov.Proj.M4.L5.loadings > .85))

        Cov.Proj.M1.L00.dimensionality.90[step] <- min(which(Cov.Proj.M1.L00.loadings > .9))
        Cov.Proj.M2.L00.dimensionality.90[step] <- min(which(Cov.Proj.M2.L00.loadings > .9))
        Cov.Proj.M3.L3.dimensionality.90[step] <- min(which(Cov.Proj.M3.L3.loadings > .9))
        Cov.Proj.M3.L4.dimensionality.90[step] <- min(which(Cov.Proj.M3.L4.loadings > .9))
        Cov.Proj.M3.L5.dimensionality.90[step] <- min(which(Cov.Proj.M3.L5.loadings > .9))
        Cov.Proj.M4.L3.dimensionality.90[step] <- min(which(Cov.Proj.M4.L3.loadings > .9))
        Cov.Proj.M4.L4.dimensionality.90[step] <- min(which(Cov.Proj.M4.L4.loadings > .9))
        Cov.Proj.M4.L5.dimensionality.90[step] <- min(which(Cov.Proj.M4.L5.loadings > .9))

        Cov.Proj.M1.L00.dimensionality.95[step] <- min(which(Cov.Proj.M1.L00.loadings > .95))
        Cov.Proj.M2.L00.dimensionality.95[step] <- min(which(Cov.Proj.M2.L00.loadings > .95))
        Cov.Proj.M3.L3.dimensionality.95[step] <- min(which(Cov.Proj.M3.L3.loadings > .95))
        Cov.Proj.M3.L4.dimensionality.95[step] <- min(which(Cov.Proj.M3.L4.loadings > .95))
        Cov.Proj.M3.L5.dimensionality.95[step] <- min(which(Cov.Proj.M3.L5.loadings > .95))
        Cov.Proj.M4.L3.dimensionality.95[step] <- min(which(Cov.Proj.M4.L3.loadings > .95))
        Cov.Proj.M4.L4.dimensionality.95[step] <- min(which(Cov.Proj.M4.L4.loadings > .95))
        Cov.Proj.M4.L5.dimensionality.95[step] <- min(which(Cov.Proj.M4.L5.loadings > .95))

        Cov.Proj.M1.L00.dimensionality.99[step] <- min(which(Cov.Proj.M1.L00.loadings > .99))
        Cov.Proj.M2.L00.dimensionality.99[step] <- min(which(Cov.Proj.M2.L00.loadings > .99))
        Cov.Proj.M3.L3.dimensionality.99[step] <- min(which(Cov.Proj.M3.L3.loadings > .99))
        Cov.Proj.M3.L4.dimensionality.99[step] <- min(which(Cov.Proj.M3.L4.loadings > .99))
        Cov.Proj.M3.L5.dimensionality.99[step] <- min(which(Cov.Proj.M3.L5.loadings > .99))
        Cov.Proj.M4.L3.dimensionality.99[step] <- min(which(Cov.Proj.M4.L3.loadings > .99))
        Cov.Proj.M4.L4.dimensionality.99[step] <- min(which(Cov.Proj.M4.L4.loadings > .99))
        Cov.Proj.M4.L5.dimensionality.99[step] <- min(which(Cov.Proj.M4.L5.loadings > .99))
      }#dimensionalities


      {

        Cov.Proj.M1.L00.Err[step] <- L2_HS_norm(Cov.Proj.M1.L00$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        Cov.Proj.M2.L00.Err[step] <- L2_HS_norm(Cov.Proj.M2.L00$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        Cov.Proj.M3.L3.Err[step] <- L2_HS_norm(Cov.Proj.M3.L3$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        Cov.Proj.M3.L4.Err[step] <- L2_HS_norm(Cov.Proj.M3.L4$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        Cov.Proj.M3.L5.Err[step] <- L2_HS_norm(Cov.Proj.M3.L5$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        Cov.Proj.M4.L3.Err[step] <- L2_HS_norm(Cov.Proj.M4.L3$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        Cov.Proj.M4.L4.Err[step] <- L2_HS_norm(Cov.Proj.M4.L4$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        Cov.Proj.M4.L5.Err[step] <- L2_HS_norm(Cov.Proj.M4.L5$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
      }#errors

    }#Calculate the errors and dimensions

  }

  # For Scenario S4: (data projected onto excess return PCAs)
  {
    #make projections
    {
      excess_ret1<-matrix(0,n-1,M)
      excess_ret2<-matrix(0,n-1,M)
      excess_ret3<-matrix(0,n-1,M)
      excess_ret4<-matrix(0,n-1,M)

      for (i in 1:(n-1)) {
        excess_ret1[i,]<-(M1.price.data[i+1,1:M]-M1.price.data[i,2:(M+1)])-M1.price.data[i,1]
        excess_ret2[i,]<-(M2.price.data[i+1,1:M]-M2.price.data[i,2:(M+1)])-M2.price.data[i,1]
        excess_ret3[i,]<-(M3.price.data[i+1,1:M]-M3.price.data[i,2:(M+1)])-M3.price.data[i,1]
        excess_ret4[i,]<-(M4.price.data[i+1,1:M]-M4.price.data[i,2:(M+1)])-M4.price.data[i,1]
      }

      E1 <- eigen(t(excess_ret1)%*%excess_ret1)
      E2 <- eigen(t(excess_ret2)%*%excess_ret2)
      E3 <- eigen(t(excess_ret3)%*%excess_ret3)
      E4 <- eigen(t(excess_ret4)%*%excess_ret4)

      excess.explavar1[step,] <- (cumsum(E1$values)/sum(E1$values))[1:10]
      excess.explavar2[step,] <- (cumsum(E2$values)/sum(E2$values))[1:10]
      excess.explavar3[step,] <- (cumsum(E3$values)/sum(E3$values))[1:10]
      excess.explavar4[step,] <- (cumsum(E4$values)/sum(E4$values))[1:10]

      d1 <- which.max(excess.explavar1[step,]>.99)
      d2 <- which.max(excess.explavar2[step,]>.99)
      d3 <- which.max(excess.explavar3[step,]>.99)
      d4 <- which.max(excess.explavar4[step,]>.99)

      excess.projections1 <- excess_ret1 %*% E1$vectors[,1:d1] %*% t(E1$vectors[,1:d1])
      excess.projections2 <- excess_ret2 %*% E2$vectors[,1:d2] %*% t(E2$vectors[,1:d2])
      excess.projections3 <- excess_ret3 %*% E3$vectors[,1:d3] %*% t(E3$vectors[,1:d3])
      excess.projections4 <- excess_ret4 %*% E4$vectors[,1:d4] %*% t(E4$vectors[,1:d4])


      excess.Proj.M1.price.data<-matrix(0,n,M)
      excess.Proj.M2.price.data<-matrix(0,n,M)
      excess.Proj.M3.price.data<-matrix(0,n,M)
      excess.Proj.M4.price.data<-matrix(0,n,M)


      excess.Proj.M1.price.data[1,]<-M1.price.data[1,1:M]
      excess.Proj.M2.price.data[1,]<-M2.price.data[1,1:M]
      excess.Proj.M3.price.data[1,]<-M3.price.data[1,1:M]
      excess.Proj.M4.price.data[1,]<-M4.price.data[1,1:M]

      # Recursively reconstruct
      for (i in 1:(n - 1)) {
        for (j in 1:(M-1)) {
          excess.Proj.M1.price.data[i+1,j]<-excess.Proj.M1.price.data[i,j+1]+excess.projections1[i,j] +excess.projections1[i,1]
          excess.Proj.M2.price.data[i+1,j]<-excess.Proj.M2.price.data[i,j+1]+excess.projections2[i,j] +excess.projections2[i,1]
          excess.Proj.M3.price.data[i+1,j]<-excess.Proj.M3.price.data[i,j+1]+excess.projections3[i,j] +excess.projections3[i,1]
          excess.Proj.M4.price.data[i+1,j]<-excess.Proj.M4.price.data[i,j+1]+excess.projections4[i,j] +excess.projections4[i,1]
        }
      }

    }

    {
      excess.Proj.M1.L00 <- tc(x= excess.Proj.M1.price.data, sparse = sparsing, tq = 0.75,l = 10000)

      excess.Proj.M2.L00 <- tc(x= excess.Proj.M2.price.data, sparse = sparsing, tq = 0.75,l = 10000)

      excess.Proj.M3.L3 <- tc(x= excess.Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 3)
      excess.Proj.M3.L4 <- tc(x= excess.Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 4)
      excess.Proj.M3.L5 <- tc(x= excess.Proj.M3.price.data, sparse = sparsing, tq = 0.75,l = 5)

      excess.Proj.M4.L3 <- tc(x= excess.Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 3)
      excess.Proj.M4.L4 <- tc(x= excess.Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 4)
      excess.Proj.M4.L5 <- tc(x= excess.Proj.M4.price.data, sparse = sparsing, tq = 0.75,l = 5)

    }#calculate the estimators


    {
      {
        excess.Proj.M1.L00.loadings <- excess.Proj.M1.L00$expl.var[1:20]
        excess.Proj.M2.L00.loadings <- excess.Proj.M2.L00$expl.var[1:20]
        excess.Proj.M3.L3.loadings <- excess.Proj.M3.L3$expl.var[1:20]
        excess.Proj.M3.L4.loadings <- excess.Proj.M3.L4$expl.var[1:20]
        excess.Proj.M3.L5.loadings <- excess.Proj.M3.L5$expl.var[1:20]
        excess.Proj.M4.L3.loadings <- excess.Proj.M4.L3$expl.var[1:20]
        excess.Proj.M4.L4.loadings <- excess.Proj.M4.L4$expl.var[1:20]
        excess.Proj.M4.L5.loadings <- excess.Proj.M4.L5$expl.var[1:20]

        excess.Proj.M1.L00.dimensionality.85[step] <- min(which(excess.Proj.M1.L00.loadings > .85))
        excess.Proj.M2.L00.dimensionality.85[step] <- min(which(excess.Proj.M2.L00.loadings > .85))
        excess.Proj.M3.L3.dimensionality.85[step] <- min(which(excess.Proj.M3.L3.loadings > .85))
        excess.Proj.M3.L4.dimensionality.85[step] <- min(which(excess.Proj.M3.L4.loadings > .85))
        excess.Proj.M3.L5.dimensionality.85[step] <- min(which(excess.Proj.M3.L5.loadings > .85))
        excess.Proj.M4.L3.dimensionality.85[step] <- min(which(excess.Proj.M4.L3.loadings > .85))
        excess.Proj.M4.L4.dimensionality.85[step] <- min(which(excess.Proj.M4.L4.loadings > .85))
        excess.Proj.M4.L5.dimensionality.85[step] <- min(which(excess.Proj.M4.L5.loadings > .85))

        excess.Proj.M1.L00.dimensionality.90[step] <- min(which(excess.Proj.M1.L00.loadings > .9))
        excess.Proj.M2.L00.dimensionality.90[step] <- min(which(excess.Proj.M2.L00.loadings > .9))
        excess.Proj.M3.L3.dimensionality.90[step] <- min(which(excess.Proj.M3.L3.loadings > .9))
        excess.Proj.M3.L4.dimensionality.90[step] <- min(which(excess.Proj.M3.L4.loadings > .9))
        excess.Proj.M3.L5.dimensionality.90[step] <- min(which(excess.Proj.M3.L5.loadings > .9))
        excess.Proj.M4.L3.dimensionality.90[step] <- min(which(excess.Proj.M4.L3.loadings > .9))
        excess.Proj.M4.L4.dimensionality.90[step] <- min(which(excess.Proj.M4.L4.loadings > .9))
        excess.Proj.M4.L5.dimensionality.90[step] <- min(which(excess.Proj.M4.L5.loadings > .9))

        excess.Proj.M1.L00.dimensionality.95[step] <- min(which(excess.Proj.M1.L00.loadings > .95))
        excess.Proj.M2.L00.dimensionality.95[step] <- min(which(excess.Proj.M2.L00.loadings > .95))
        excess.Proj.M3.L3.dimensionality.95[step] <- min(which(excess.Proj.M3.L3.loadings > .95))
        excess.Proj.M3.L4.dimensionality.95[step] <- min(which(excess.Proj.M3.L4.loadings > .95))
        excess.Proj.M3.L5.dimensionality.95[step] <- min(which(excess.Proj.M3.L5.loadings > .95))
        excess.Proj.M4.L3.dimensionality.95[step] <- min(which(excess.Proj.M4.L3.loadings > .95))
        excess.Proj.M4.L4.dimensionality.95[step] <- min(which(excess.Proj.M4.L4.loadings > .95))
        excess.Proj.M4.L5.dimensionality.95[step] <- min(which(excess.Proj.M4.L5.loadings > .95))

        excess.Proj.M1.L00.dimensionality.99[step] <- min(which(excess.Proj.M1.L00.loadings > .99))
        excess.Proj.M2.L00.dimensionality.99[step] <- min(which(excess.Proj.M2.L00.loadings > .99))
        excess.Proj.M3.L3.dimensionality.99[step] <- min(which(excess.Proj.M3.L3.loadings > .99))
        excess.Proj.M3.L4.dimensionality.99[step] <- min(which(excess.Proj.M3.L4.loadings > .99))
        excess.Proj.M3.L5.dimensionality.99[step] <- min(which(excess.Proj.M3.L5.loadings > .99))
        excess.Proj.M4.L3.dimensionality.99[step] <- min(which(excess.Proj.M4.L3.loadings > .99))
        excess.Proj.M4.L4.dimensionality.99[step] <- min(which(excess.Proj.M4.L4.loadings > .99))
        excess.Proj.M4.L5.dimensionality.99[step] <- min(which(excess.Proj.M4.L5.loadings > .99))
      }#dimensionalities


      {u <- min(nrow(excess.Proj.M1.L00$IV),nrow(excess.Proj.M1.L00$IV))

        excess.Proj.M1.L00.Err[step] <- L2_HS_norm(excess.Proj.M1.L00$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        excess.Proj.M2.L00.Err[step] <- L2_HS_norm(excess.Proj.M2.L00$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        excess.Proj.M3.L3.Err[step] <- L2_HS_norm(excess.Proj.M3.L3$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        excess.Proj.M3.L4.Err[step] <- L2_HS_norm(excess.Proj.M3.L4$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        excess.Proj.M3.L5.Err[step] <- L2_HS_norm(excess.Proj.M3.L5$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)

        excess.Proj.M4.L3.Err[step] <- L2_HS_norm(excess.Proj.M4.L3$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        excess.Proj.M4.L4.Err[step] <- L2_HS_norm(excess.Proj.M4.L4$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
        excess.Proj.M4.L5.Err[step] <- L2_HS_norm(excess.Proj.M4.L5$IV[1:u,1:u]-Integrated.Volatility[1:u,1:u], from = 0, to = 10)
      }#errors

    }#Calculate the errors and dimensions

  }
  print(step/K)
}#Run the Monte-Carlo-Simulation
proc.time()-ptm




# Create a dataframe containing the results
{

  {
    # Define LaTeX-style row labels as a column
    latex_labels <- c(
      "S.1-rel. err.",
      "S.1-Dimension(0.85)",
      "S.1-Dimension(0.90)",
      "S.1-Dimension(0.95)",
      "S.1-Dimension(0.99)",
      "S.2-rel. err.",
      "S.2-Dimension(0.85)",
      "S.2-Dimension(0.90)",
      "S.2-Dimension(0.95)",
      "S.2-Dimension(0.99)",
      "S.3-rel. err.",
      "S.3-Dimension(0.85)",
      "S.3-Dimension(0.90)",
      "S.3-Dimension(0.95)",
      "S.3-Dimension(0.99)",
      "S.4-rel. err.",
      "S.4-Dimension(0.85)",
      "S.4-Dimension(0.90)",
      "S.4-Dimension(0.95)",
      "S.4-Dimension(0.99)"
    )
    col_names <- c( "Model 1, l = ∞", "Model 2, l = ∞",
                    "Model 3, l = 3", "Model 3, l = 4", "Model 3, l = 5",
                    "Model 4, l = 3", "Model 4, l = 4", "Model 4, l = 5")

    # Fill in with NA values
    results <- data.frame(matrix(NA, nrow = length(latex_labels), ncol = length(col_names)))
    colnames(results) <- col_names

    # Add LaTeX row labels as a column
    results <- cbind("Model & truncation level (columns) \\errors and meas.dims.(rows)" = latex_labels, results)
  }#create a Dataframe containing the Monte-Carlo results

  #S1
  {
    ###Calculate relative errors

    med.M1.L00<-median(M1.L00.Err/normIV)
    quarts.M1.L00<-quantile(M1.L00.Err/normIV, c(.25,.75))

    results$`Model 1, l = ∞`[1]<- sprintf("%.3f (%.3f, %.3f)", med.M1.L00,quarts.M1.L00[1],quarts.M1.L00[2])


    med.M2.L00<-median(M2.L00.Err/normIV)
    quarts.M2.L00<-quantile(M2.L00.Err/normIV, c(.25,.75))

    results$`Model 2, l = ∞`[1]<- sprintf("%.3f (%.3f, %.3f)", med.M2.L00,quarts.M2.L00[1],quarts.M2.L00[2])


    med.M3.L3<-median(M3.L3.Err/normIV)
    quarts.M3.L3<-quantile(M3.L3.Err/normIV, c(.25,.75))

    results$`Model 3, l = 3`[1]<- sprintf("%.3f (%.3f, %.3f)", med.M3.L3,quarts.M3.L3[1],quarts.M3.L3[2])


    med.M3.L4<-median(M3.L4.Err/normIV)
    quarts.M3.L4<-quantile(M3.L4.Err/normIV, c(.25,.75))

    results$`Model 3, l = 4`[1]<- sprintf("%.3f (%.3f, %.3f)", med.M3.L4,quarts.M3.L4[1],quarts.M3.L4[2])


    med.M3.L5<-median(M3.L5.Err/normIV)
    quarts.M3.L5<-quantile(M3.L5.Err/normIV, c(.25,.75))

    results$`Model 3, l = 5`[1]<- sprintf("%.3f (%.3f, %.3f)", med.M3.L5,quarts.M3.L5[1],quarts.M3.L5[2])


    med.M4.L3<-median(M4.L3.Err/normIV)
    quarts.M4.L3<-quantile(M4.L3.Err/normIV, c(.25,.75))

    results$`Model 4, l = 3`[1]<- sprintf("%.3f (%.3f, %.3f)", med.M4.L3,quarts.M4.L3[1],quarts.M4.L3[2])


    med.M4.L4<-median(M4.L4.Err/normIV)
    quarts.M4.L4<-quantile(M4.L4.Err/normIV, c(.25,.75))

    results$`Model 4, l = 4`[1]<- sprintf("%.3f (%.3f, %.3f)", med.M4.L4,quarts.M4.L4[1],quarts.M4.L4[2])



    med.M4.L5<-median(M4.L5.Err/normIV)
    quarts.M4.L5<-quantile(M4.L5.Err/normIV, c(.25,.75))

    results$`Model 4, l = 5`[1]<- sprintf("%.3f (%.3f, %.3f)", med.M4.L5,quarts.M4.L5[1],quarts.M4.L5[2])


    #####Dimensions for 85%

    med.dim.M1.L00.85 <- median(M1.L00.dimensionality.85)
    quarts.dim.M1.L00.85<-quantile(M1.L00.dimensionality.85, c(.25,.75))

    results$`Model 1, l = ∞`[2]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M1.L00.85,quarts.dim.M1.L00.85[1],quarts.dim.M1.L00.85[2])


    med.dim.M2.L00.85 <- median(M2.L00.dimensionality.85)
    quarts.dim.M2.L00.85 <-quantile(M2.L00.dimensionality.85, c(.25,.75))

    results$`Model 2, l = ∞`[2]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M2.L00.85,quarts.dim.M2.L00.85[1],quarts.dim.M2.L00.85[2])


    med.dim.M3.L3.85 <- median(M3.L3.dimensionality.85)
    quarts.dim.M3.L3.85 <-quantile(M3.L3.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 3`[2]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L3.85,quarts.dim.M3.L3.85[1],quarts.dim.M3.L3.85[2])


    med.dim.M3.L4.85 <- median(M3.L4.dimensionality.85)
    quarts.dim.M3.L4.85<-quantile(M3.L4.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 4`[2]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L4.85,quarts.dim.M3.L4.85[1],quarts.dim.M3.L4.85[2])


    med.dim.M3.L5.85 <- median(M3.L5.dimensionality.85)
    quarts.dim.M3.L5.85 <-quantile(M3.L5.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 5`[2]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L5.85,quarts.dim.M3.L5.85[1],quarts.dim.M3.L5.85[2])


    med.dim.M4.L3.85 <- median(M4.L3.dimensionality.85)
    quarts.dim.M4.L3.85<-quantile(M4.L3.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 3`[2]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L3.85,quarts.dim.M4.L3.85[1],quarts.dim.M4.L3.85[2])


    med.dim.M4.L4.85 <- median(M4.L4.dimensionality.85)
    quarts.dim.M4.L4.85<-quantile(M4.L4.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 4`[2]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L4.85,quarts.dim.M4.L4.85[1],quarts.dim.M4.L4.85[2])


    med.dim.M4.L5.85 <- median(M4.L5.dimensionality.85)
    quarts.dim.M4.L5.85<-quantile(M4.L5.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 5`[2]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L5.85,quarts.dim.M4.L5.85[1],quarts.dim.M4.L5.85[2])



    #####Dimensions for 90%

    med.dim.M1.L00.90 <- median(M1.L00.dimensionality.90)
    quarts.dim.M1.L00.90<-quantile(M1.L00.dimensionality.90, c(.25,.75))

    results$`Model 1, l = ∞`[3]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M1.L00.90,quarts.dim.M1.L00.90[1],quarts.dim.M1.L00.90[2])


    med.dim.M2.L00.90 <- median(M2.L00.dimensionality.90)
    quarts.dim.M2.L00.90 <-quantile(M2.L00.dimensionality.90, c(.25,.75))

    results$`Model 2, l = ∞`[3]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M2.L00.90,quarts.dim.M2.L00.90[1],quarts.dim.M2.L00.90[2])


    med.dim.M3.L3.90 <- median(M3.L3.dimensionality.90)
    quarts.dim.M3.L3.90 <-quantile(M3.L3.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 3`[3]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L3.90,quarts.dim.M3.L3.90[1],quarts.dim.M3.L3.90[2])


    med.dim.M3.L4.90 <- median(M3.L4.dimensionality.90)
    quarts.dim.M3.L4.90<-quantile(M3.L4.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 4`[3]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L4.90,quarts.dim.M3.L4.90[1],quarts.dim.M3.L4.90[2])


    med.dim.M3.L5.90 <- median(M3.L5.dimensionality.90)
    quarts.dim.M3.L5.90 <-quantile(M3.L5.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 5`[3]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L5.90,quarts.dim.M3.L5.90[1],quarts.dim.M3.L5.90[2])


    med.dim.M4.L3.90 <- median(M4.L3.dimensionality.90)
    quarts.dim.M4.L3.90<-quantile(M4.L3.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 3`[3]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L3.90,quarts.dim.M4.L3.90[1],quarts.dim.M4.L3.90[2])


    med.dim.M4.L4.90 <- median(M4.L4.dimensionality.90)
    quarts.dim.M4.L4.90<-quantile(M4.L4.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 4`[3]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L4.90,quarts.dim.M4.L4.90[1],quarts.dim.M4.L4.90[2])


    med.dim.M4.L5.90 <- median(M4.L5.dimensionality.90)
    quarts.dim.M4.L5.90<-quantile(M4.L5.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 5`[3]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L5.90,quarts.dim.M4.L5.90[1],quarts.dim.M4.L5.90[2])



    #####Dimensions for 95%

    med.dim.M1.L00.95 <- median(M1.L00.dimensionality.95)
    quarts.dim.M1.L00.95<-quantile(M1.L00.dimensionality.95, c(.25,.75))

    results$`Model 1, l = ∞`[4]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M1.L00.95,quarts.dim.M1.L00.95[1],quarts.dim.M1.L00.95[2])


    med.dim.M2.L00.95 <- median(M2.L00.dimensionality.95)
    quarts.dim.M2.L00.95 <-quantile(M2.L00.dimensionality.95, c(.25,.75))

    results$`Model 2, l = ∞`[4]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M2.L00.95,quarts.dim.M2.L00.95[1],quarts.dim.M2.L00.95[2])


    med.dim.M3.L3.95 <- median(M3.L3.dimensionality.95)
    quarts.dim.M3.L3.95 <-quantile(M3.L3.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 3`[4]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L3.95,quarts.dim.M3.L3.95[1],quarts.dim.M3.L3.95[2])


    med.dim.M3.L4.95 <- median(M3.L4.dimensionality.95)
    quarts.dim.M3.L4.95<-quantile(M3.L4.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 4`[4]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L4.95,quarts.dim.M3.L4.95[1],quarts.dim.M3.L4.95[2])


    med.dim.M3.L5.95 <- median(M3.L5.dimensionality.95)
    quarts.dim.M3.L5.95 <-quantile(M3.L5.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 5`[4]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L5.95,quarts.dim.M3.L5.95[1],quarts.dim.M3.L5.95[2])


    med.dim.M4.L3.95 <- median(M4.L3.dimensionality.95)
    quarts.dim.M4.L3.95<-quantile(M4.L3.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 3`[4]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L3.95,quarts.dim.M4.L3.95[1],quarts.dim.M4.L3.95[2])


    med.dim.M4.L4.95 <- median(M4.L4.dimensionality.95)
    quarts.dim.M4.L4.95<-quantile(M4.L4.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 4`[4]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L4.95,quarts.dim.M4.L4.95[1],quarts.dim.M4.L4.95[2])


    med.dim.M4.L5.95 <- median(M4.L5.dimensionality.95)
    quarts.dim.M4.L5.95<-quantile(M4.L5.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 5`[4]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L5.95,quarts.dim.M4.L5.95[1],quarts.dim.M4.L5.95[2])


    #####Dimensions for 99%

    med.dim.M1.L00.99 <- median(M1.L00.dimensionality.99)
    quarts.dim.M1.L00.99<-quantile(M1.L00.dimensionality.99, c(.25,.75))

    results$`Model 1, l = ∞`[5]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M1.L00.99,quarts.dim.M1.L00.99[1],quarts.dim.M1.L00.99[2])


    med.dim.M2.L00.99 <- median(M2.L00.dimensionality.99)
    quarts.dim.M2.L00.99 <-quantile(M2.L00.dimensionality.99, c(.25,.75))

    results$`Model 2, l = ∞`[5]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M2.L00.99,quarts.dim.M2.L00.99[1],quarts.dim.M2.L00.99[2])


    med.dim.M3.L3.99 <- median(M3.L3.dimensionality.99)
    quarts.dim.M3.L3.99 <-quantile(M3.L3.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 3`[5]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L3.99,quarts.dim.M3.L3.99[1],quarts.dim.M3.L3.99[2])


    med.dim.M3.L4.99 <- median(M3.L4.dimensionality.99)
    quarts.dim.M3.L4.99<-quantile(M3.L4.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 4`[5]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L4.99,quarts.dim.M3.L4.99[1],quarts.dim.M3.L4.99[2])


    med.dim.M3.L5.99 <- median(M3.L5.dimensionality.99)
    quarts.dim.M3.L5.99 <-quantile(M3.L5.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 5`[5]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M3.L5.99,quarts.dim.M3.L5.99[1],quarts.dim.M3.L5.99[2])


    med.dim.M4.L3.99 <- median(M4.L3.dimensionality.99)
    quarts.dim.M4.L3.99<-quantile(M4.L3.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 3`[5]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L3.99,quarts.dim.M4.L3.99[1],quarts.dim.M4.L3.99[2])


    med.dim.M4.L4.99 <- median(M4.L4.dimensionality.99)
    quarts.dim.M4.L4.99<-quantile(M4.L4.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 4`[5]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L4.99,quarts.dim.M4.L4.99[1],quarts.dim.M4.L4.99[2])


    med.dim.M4.L5.99 <- median(M4.L5.dimensionality.99)
    quarts.dim.M4.L5.99<-quantile(M4.L5.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 5`[5]<- sprintf("%.3f (%.3f, %.3f)", med.dim.M4.L5.99,quarts.dim.M4.L5.99[1],quarts.dim.M4.L5.99[2])
  }

  #S2
  {
    ###Calculate relative errors

    Proj.med.M1.L00<-median(Proj.M1.L00.Err/normIV)
    Proj.quarts.M1.L00<-quantile(Proj.M1.L00.Err/normIV, c(.25,.75))

    results$`Model 1, l = ∞`[6]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.M1.L00,Proj.quarts.M1.L00[1],Proj.quarts.M1.L00[2])


    Proj.med.M2.L00<-median(Proj.M2.L00.Err/normIV)
    Proj.quarts.M2.L00<-quantile(Proj.M2.L00.Err/normIV, c(.25,.75))

    results$`Model 2, l = ∞`[6]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.M2.L00,Proj.quarts.M2.L00[1],Proj.quarts.M2.L00[2])


    Proj.med.M3.L3<-median(Proj.M3.L3.Err/normIV)
    Proj.quarts.M3.L3<-quantile(Proj.M3.L3.Err/normIV, c(.25,.75))

    results$`Model 3, l = 3`[6]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.M3.L3,Proj.quarts.M3.L3[1],Proj.quarts.M3.L3[2])


    Proj.med.M3.L4<-median(Proj.M3.L4.Err/normIV)
    Proj.quarts.M3.L4<-quantile(Proj.M3.L4.Err/normIV, c(.25,.75))

    results$`Model 3, l = 4`[6]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.M3.L4,Proj.quarts.M3.L4[1],Proj.quarts.M3.L4[2])


    Proj.med.M3.L5<-median(Proj.M3.L5.Err/normIV)
    Proj.quarts.M3.L5<-quantile(Proj.M3.L5.Err/normIV, c(.25,.75))

    results$`Model 3, l = 5`[6]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.M3.L5,Proj.quarts.M3.L5[1],Proj.quarts.M3.L5[2])


    Proj.med.M4.L3<-median(Proj.M4.L3.Err/normIV)
    Proj.quarts.M4.L3<-quantile(Proj.M4.L3.Err/normIV, c(.25,.75))

    results$`Model 4, l = 3`[6]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.M4.L3,Proj.quarts.M4.L3[1],Proj.quarts.M4.L3[2])


    Proj.med.M4.L4<-median(Proj.M4.L4.Err/normIV)
    Proj.quarts.M4.L4<-quantile(Proj.M4.L4.Err/normIV, c(.25,.75))

    results$`Model 4, l = 4`[6]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.M4.L4,Proj.quarts.M4.L4[1],Proj.quarts.M4.L4[2])



    Proj.med.M4.L5<-median(Proj.M4.L5.Err/normIV)
    Proj.quarts.M4.L5<-quantile(Proj.M4.L5.Err/normIV, c(.25,.75))

    results$`Model 4, l = 5`[6]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.M4.L5,Proj.quarts.M4.L5[1],Proj.quarts.M4.L5[2])


    #####Dimensions for 85%

    Proj.med.dim.M1.L00.85 <- median(Proj.M1.L00.dimensionality.85)
    Proj.quarts.dim.M1.L00.85<-quantile(Proj.M1.L00.dimensionality.85, c(.25,.75))

    results$`Model 1, l = ∞`[7]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M1.L00.85,Proj.quarts.dim.M1.L00.85[1],Proj.quarts.dim.M1.L00.85[2])


    Proj.med.dim.M2.L00.85 <- median(Proj.M2.L00.dimensionality.85)
    Proj.quarts.dim.M2.L00.85 <-quantile(Proj.M2.L00.dimensionality.85, c(.25,.75))

    results$`Model 2, l = ∞`[7]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M2.L00.85,Proj.quarts.dim.M2.L00.85[1],Proj.quarts.dim.M2.L00.85[2])


    Proj.med.dim.M3.L3.85 <- median(Proj.M3.L3.dimensionality.85)
    Proj.quarts.dim.M3.L3.85 <-quantile(Proj.M3.L3.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 3`[7]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L3.85,Proj.quarts.dim.M3.L3.85[1],Proj.quarts.dim.M3.L3.85[2])


    Proj.med.dim.M3.L4.85 <- median(Proj.M3.L4.dimensionality.85)
    Proj.quarts.dim.M3.L4.85<-quantile(Proj.M3.L4.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 4`[7]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L4.85,Proj.quarts.dim.M3.L4.85[1],Proj.quarts.dim.M3.L4.85[2])


    Proj.med.dim.M3.L5.85 <- median(Proj.M3.L5.dimensionality.85)
    Proj.quarts.dim.M3.L5.85 <-quantile(Proj.M3.L5.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 5`[7]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L5.85,Proj.quarts.dim.M3.L5.85[1],Proj.quarts.dim.M3.L5.85[2])


    Proj.med.dim.M4.L3.85 <- median(Proj.M4.L3.dimensionality.85)
    Proj.quarts.dim.M4.L3.85<-quantile(Proj.M4.L3.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 3`[7]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L3.85,Proj.quarts.dim.M4.L3.85[1],Proj.quarts.dim.M4.L3.85[2])


    Proj.med.dim.M4.L4.85 <- median(Proj.M4.L4.dimensionality.85)
    Proj.quarts.dim.M4.L4.85<-quantile(Proj.M4.L4.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 4`[7]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L4.85,Proj.quarts.dim.M4.L4.85[1],Proj.quarts.dim.M4.L4.85[2])


    Proj.med.dim.M4.L5.85 <- median(Proj.M4.L5.dimensionality.85)
    Proj.quarts.dim.M4.L5.85<-quantile(Proj.M4.L5.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 5`[7]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L5.85,Proj.quarts.dim.M4.L5.85[1],Proj.quarts.dim.M4.L5.85[2])



    #####Dimensions for 90%

    Proj.med.dim.M1.L00.90 <- median(Proj.M1.L00.dimensionality.90)
    Proj.quarts.dim.M1.L00.90<-quantile(Proj.M1.L00.dimensionality.90, c(.25,.75))

    results$`Model 1, l = ∞`[8]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M1.L00.90,Proj.quarts.dim.M1.L00.90[1],Proj.quarts.dim.M1.L00.90[2])


    Proj.med.dim.M2.L00.90 <- median(Proj.M2.L00.dimensionality.90)
    Proj.quarts.dim.M2.L00.90 <-quantile(Proj.M2.L00.dimensionality.90, c(.25,.75))

    results$`Model 2, l = ∞`[8]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M2.L00.90,Proj.quarts.dim.M2.L00.90[1],Proj.quarts.dim.M2.L00.90[2])


    Proj.med.dim.M3.L3.90 <- median(Proj.M3.L3.dimensionality.90)
    Proj.quarts.dim.M3.L3.90 <-quantile(Proj.M3.L3.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 3`[8]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L3.90,Proj.quarts.dim.M3.L3.90[1],Proj.quarts.dim.M3.L3.90[2])


    Proj.med.dim.M3.L4.90 <- median(Proj.M3.L4.dimensionality.90)
    Proj.quarts.dim.M3.L4.90<-quantile(Proj.M3.L4.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 4`[8]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L4.90,Proj.quarts.dim.M3.L4.90[1],Proj.quarts.dim.M3.L4.90[2])


    Proj.med.dim.M3.L5.90 <- median(Proj.M3.L5.dimensionality.90)
    Proj.quarts.dim.M3.L5.90 <-quantile(Proj.M3.L5.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 5`[8]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L5.90,Proj.quarts.dim.M3.L5.90[1],Proj.quarts.dim.M3.L5.90[2])


    Proj.med.dim.M4.L3.90 <- median(Proj.M4.L3.dimensionality.90)
    Proj.quarts.dim.M4.L3.90<-quantile(Proj.M4.L3.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 3`[8]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L3.90,Proj.quarts.dim.M4.L3.90[1],Proj.quarts.dim.M4.L3.90[2])


    Proj.med.dim.M4.L4.90 <- median(Proj.M4.L4.dimensionality.90)
    Proj.quarts.dim.M4.L4.90<-quantile(Proj.M4.L4.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 4`[8]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L4.90,Proj.quarts.dim.M4.L4.90[1],Proj.quarts.dim.M4.L4.90[2])


    Proj.med.dim.M4.L5.90 <- median(Proj.M4.L5.dimensionality.90)
    Proj.quarts.dim.M4.L5.90<-quantile(Proj.M4.L5.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 5`[8]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L5.90,Proj.quarts.dim.M4.L5.90[1],Proj.quarts.dim.M4.L5.90[2])



    #####Dimensions for 95%

    Proj.med.dim.M1.L00.95 <- median(Proj.M1.L00.dimensionality.95)
    Proj.quarts.dim.M1.L00.95<-quantile(Proj.M1.L00.dimensionality.95, c(.25,.75))

    results$`Model 1, l = ∞`[9]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M1.L00.95,Proj.quarts.dim.M1.L00.95[1],Proj.quarts.dim.M1.L00.95[2])


    Proj.med.dim.M2.L00.95 <- median(Proj.M2.L00.dimensionality.95)
    Proj.quarts.dim.M2.L00.95 <-quantile(Proj.M2.L00.dimensionality.95, c(.25,.75))

    results$`Model 2, l = ∞`[9]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M2.L00.95,Proj.quarts.dim.M2.L00.95[1],Proj.quarts.dim.M2.L00.95[2])


    Proj.med.dim.M3.L3.95 <- median(Proj.M3.L3.dimensionality.95)
    Proj.quarts.dim.M3.L3.95 <-quantile(Proj.M3.L3.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 3`[9]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L3.95,Proj.quarts.dim.M3.L3.95[1],Proj.quarts.dim.M3.L3.95[2])


    Proj.med.dim.M3.L4.95 <- median(Proj.M3.L4.dimensionality.95)
    Proj.quarts.dim.M3.L4.95<-quantile(Proj.M3.L4.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 4`[9]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L4.95,Proj.quarts.dim.M3.L4.95[1],Proj.quarts.dim.M3.L4.95[2])


    Proj.med.dim.M3.L5.95 <- median(Proj.M3.L5.dimensionality.95)
    Proj.quarts.dim.M3.L5.95 <-quantile(Proj.M3.L5.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 5`[9]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L5.95, Proj.quarts.dim.M3.L5.95[1], Proj.quarts.dim.M3.L5.95[2])


    Proj.med.dim.M4.L3.95 <- median(Proj.M4.L3.dimensionality.95)
    Proj.quarts.dim.M4.L3.95<-quantile(Proj.M4.L3.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 3`[9]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L3.95,Proj.quarts.dim.M4.L3.95[1],Proj.quarts.dim.M4.L3.95[2])


    Proj.med.dim.M4.L4.95 <- median(Proj.M4.L4.dimensionality.95)
    Proj.quarts.dim.M4.L4.95<-quantile(Proj.M4.L4.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 4`[9]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L4.95,Proj.quarts.dim.M4.L4.95[1],Proj.quarts.dim.M4.L4.95[2])


    Proj.med.dim.M4.L5.95 <- median(Proj.M4.L5.dimensionality.95)
    Proj.quarts.dim.M4.L5.95<-quantile(Proj.M4.L5.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 5`[9]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L5.95,Proj.quarts.dim.M4.L5.95[1],Proj.quarts.dim.M4.L5.95[2])


    #####Dimensions for 99%

    Proj.med.dim.M1.L00.99 <- median(Proj.M1.L00.dimensionality.99)
    Proj.quarts.dim.M1.L00.99<-quantile(Proj.M1.L00.dimensionality.99, c(.25,.75))

    results$`Model 1, l = ∞`[10]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M1.L00.99,Proj.quarts.dim.M1.L00.99[1],Proj.quarts.dim.M1.L00.99[2])


    Proj.med.dim.M2.L00.99 <- median(Proj.M2.L00.dimensionality.99)
    Proj.quarts.dim.M2.L00.99 <-quantile(Proj.M2.L00.dimensionality.99, c(.25,.75))

    results$`Model 2, l = ∞`[10]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M2.L00.99, Proj.quarts.dim.M2.L00.99[1], Proj.quarts.dim.M2.L00.99[2])


    Proj.med.dim.M3.L3.99 <- median(Proj.M3.L3.dimensionality.99)
    Proj.quarts.dim.M3.L3.99 <-quantile(Proj.M3.L3.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 3`[10]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L3.99, Proj.quarts.dim.M3.L3.99[1], Proj.quarts.dim.M3.L3.99[2])


    Proj.med.dim.M3.L4.99 <- median(Proj.M3.L4.dimensionality.99)
    Proj.quarts.dim.M3.L4.99<-quantile(Proj.M3.L4.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 4`[10]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L4.99, Proj.quarts.dim.M3.L4.99[1], Proj.quarts.dim.M3.L4.99[2])


    Proj.med.dim.M3.L5.99 <- median(Proj.M3.L5.dimensionality.99)
    Proj.quarts.dim.M3.L5.99 <-quantile(Proj.M3.L5.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 5`[10]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M3.L5.99, Proj.quarts.dim.M3.L5.99[1], Proj.quarts.dim.M3.L5.99[2])


    Proj.med.dim.M4.L3.99 <- median(Proj.M4.L3.dimensionality.99)
    Proj.quarts.dim.M4.L3.99<-quantile(Proj.M4.L3.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 3`[10]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L3.99, Proj.quarts.dim.M4.L3.99[1], Proj.quarts.dim.M4.L3.99[2])


    Proj.med.dim.M4.L4.99 <- median(Proj.M4.L4.dimensionality.99)
    Proj.quarts.dim.M4.L4.99<-quantile(Proj.M4.L4.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 4`[10]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L4.99, Proj.quarts.dim.M4.L4.99[1], Proj.quarts.dim.M4.L4.99[2])


    Proj.med.dim.M4.L5.99 <- median(Proj.M4.L5.dimensionality.99)
    Proj.quarts.dim.M4.L5.99<-quantile(Proj.M4.L5.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 5`[10]<- sprintf("%.3f (%.3f, %.3f)", Proj.med.dim.M4.L5.99, Proj.quarts.dim.M4.L5.99[1], Proj.quarts.dim.M4.L5.99[2])
  }

  #S3
  {
    ###Calculate relative errors

    Cov.Proj.med.M1.L00<-median(Cov.Proj.M1.L00.Err/normIV)
    Cov.Proj.quarts.M1.L00<-quantile(Cov.Proj.M1.L00.Err/normIV, c(.25,.75))

    results$`Model 1, l = ∞`[11]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.M1.L00,Cov.Proj.quarts.M1.L00[1],Cov.Proj.quarts.M1.L00[2])


    Cov.Proj.med.M2.L00<-median(Cov.Proj.M2.L00.Err/normIV)
    Cov.Proj.quarts.M2.L00<-quantile(Cov.Proj.M2.L00.Err/normIV, c(.25,.75))

    results$`Model 2, l = ∞`[11]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.M2.L00,Cov.Proj.quarts.M2.L00[1],Cov.Proj.quarts.M2.L00[2])


    Cov.Proj.med.M3.L3<-median(Cov.Proj.M3.L3.Err/normIV)
    Cov.Proj.quarts.M3.L3<-quantile(Cov.Proj.M3.L3.Err/normIV, c(.25,.75))

    results$`Model 3, l = 3`[11]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.M3.L3,Cov.Proj.quarts.M3.L3[1],Cov.Proj.quarts.M3.L3[2])


    Cov.Proj.med.M3.L4<-median(Cov.Proj.M3.L4.Err/normIV)
    Cov.Proj.quarts.M3.L4<-quantile(Cov.Proj.M3.L4.Err/normIV, c(.25,.75))

    results$`Model 3, l = 4`[11]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.M3.L4,Cov.Proj.quarts.M3.L4[1],Cov.Proj.quarts.M3.L4[2])


    Cov.Proj.med.M3.L5<-median(Cov.Proj.M3.L5.Err/normIV)
    Cov.Proj.quarts.M3.L5<-quantile(Cov.Proj.M3.L5.Err/normIV, c(.25,.75))

    results$`Model 3, l = 5`[11]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.M3.L5,Cov.Proj.quarts.M3.L5[1],Cov.Proj.quarts.M3.L5[2])


    Cov.Proj.med.M4.L3<-median(Cov.Proj.M4.L3.Err/normIV)
    Cov.Proj.quarts.M4.L3<-quantile(Cov.Proj.M4.L3.Err/normIV, c(.25,.75))

    results$`Model 4, l = 3`[11]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.M4.L3,Cov.Proj.quarts.M4.L3[1],Cov.Proj.quarts.M4.L3[2])


    Cov.Proj.med.M4.L4<-median(Cov.Proj.M4.L4.Err/normIV)
    Cov.Proj.quarts.M4.L4<-quantile(Cov.Proj.M4.L4.Err/normIV, c(.25,.75))

    results$`Model 4, l = 4`[11]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.M4.L4,Cov.Proj.quarts.M4.L4[1],Cov.Proj.quarts.M4.L4[2])



    Cov.Proj.med.M4.L5<-median(Cov.Proj.M4.L5.Err/normIV)
    Cov.Proj.quarts.M4.L5<-quantile(Cov.Proj.M4.L5.Err/normIV, c(.25,.75))

    results$`Model 4, l = 5`[11]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.M4.L5,Cov.Proj.quarts.M4.L5[1],Cov.Proj.quarts.M4.L5[2])


    #####Dimensions for 85%

    Cov.Proj.med.dim.M1.L00.85 <- median(Cov.Proj.M1.L00.dimensionality.85)
    Cov.Proj.quarts.dim.M1.L00.85<-quantile(Cov.Proj.M1.L00.dimensionality.85, c(.25,.75))

    results$`Model 1, l = ∞`[12]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M1.L00.85,Cov.Proj.quarts.dim.M1.L00.85[1],Cov.Proj.quarts.dim.M1.L00.85[2])


    Cov.Proj.med.dim.M2.L00.85 <- median(Cov.Proj.M2.L00.dimensionality.85)
    Cov.Proj.quarts.dim.M2.L00.85 <-quantile(Cov.Proj.M2.L00.dimensionality.85, c(.25,.75))

    results$`Model 2, l = ∞`[12]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M2.L00.85,Cov.Proj.quarts.dim.M2.L00.85[1],Cov.Proj.quarts.dim.M2.L00.85[2])


    Cov.Proj.med.dim.M3.L3.85 <- median(Cov.Proj.M3.L3.dimensionality.85)
    Cov.Proj.quarts.dim.M3.L3.85 <-quantile(Cov.Proj.M3.L3.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 3`[12]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L3.85,Cov.Proj.quarts.dim.M3.L3.85[1],Cov.Proj.quarts.dim.M3.L3.85[2])


    Cov.Proj.med.dim.M3.L4.85 <- median(Cov.Proj.M3.L4.dimensionality.85)
    Cov.Proj.quarts.dim.M3.L4.85<-quantile(Cov.Proj.M3.L4.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 4`[12]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L4.85,Cov.Proj.quarts.dim.M3.L4.85[1],Cov.Proj.quarts.dim.M3.L4.85[2])


    Cov.Proj.med.dim.M3.L5.85 <- median(Cov.Proj.M3.L5.dimensionality.85)
    Cov.Proj.quarts.dim.M3.L5.85 <-quantile(Cov.Proj.M3.L5.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 5`[12]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L5.85,Cov.Proj.quarts.dim.M3.L5.85[1],Cov.Proj.quarts.dim.M3.L5.85[2])


    Cov.Proj.med.dim.M4.L3.85 <- median(Cov.Proj.M4.L3.dimensionality.85)
    Cov.Proj.quarts.dim.M4.L3.85<-quantile(Cov.Proj.M4.L3.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 3`[12]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L3.85,Cov.Proj.quarts.dim.M4.L3.85[1],Cov.Proj.quarts.dim.M4.L3.85[2])


    Cov.Proj.med.dim.M4.L4.85 <- median(Cov.Proj.M4.L4.dimensionality.85)
    Cov.Proj.quarts.dim.M4.L4.85<-quantile(Cov.Proj.M4.L4.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 4`[12]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L4.85,Cov.Proj.quarts.dim.M4.L4.85[1],Cov.Proj.quarts.dim.M4.L4.85[2])


    Cov.Proj.med.dim.M4.L5.85 <- median(Cov.Proj.M4.L5.dimensionality.85)
    Cov.Proj.quarts.dim.M4.L5.85<-quantile(Cov.Proj.M4.L5.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 5`[12]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L5.85,Cov.Proj.quarts.dim.M4.L5.85[1],Cov.Proj.quarts.dim.M4.L5.85[2])



    #####Dimensions for 90%

    Cov.Proj.med.dim.M1.L00.90 <- median(Cov.Proj.M1.L00.dimensionality.90)
    Cov.Proj.quarts.dim.M1.L00.90<-quantile(Cov.Proj.M1.L00.dimensionality.90, c(.25,.75))

    results$`Model 1, l = ∞`[13]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M1.L00.90,Cov.Proj.quarts.dim.M1.L00.90[1],Cov.Proj.quarts.dim.M1.L00.90[2])


    Cov.Proj.med.dim.M2.L00.90 <- median(Cov.Proj.M2.L00.dimensionality.90)
    Cov.Proj.quarts.dim.M2.L00.90 <-quantile(Cov.Proj.M2.L00.dimensionality.90, c(.25,.75))

    results$`Model 2, l = ∞`[13]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M2.L00.90,Cov.Proj.quarts.dim.M2.L00.90[1],Cov.Proj.quarts.dim.M2.L00.90[2])


    Cov.Proj.med.dim.M3.L3.90 <- median(Cov.Proj.M3.L3.dimensionality.90)
    Cov.Proj.quarts.dim.M3.L3.90 <-quantile(Cov.Proj.M3.L3.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 3`[13]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L3.90,Cov.Proj.quarts.dim.M3.L3.90[1],Cov.Proj.quarts.dim.M3.L3.90[2])


    Cov.Proj.med.dim.M3.L4.90 <- median(Cov.Proj.M3.L4.dimensionality.90)
    Cov.Proj.quarts.dim.M3.L4.90<-quantile(Cov.Proj.M3.L4.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 4`[13]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L4.90,Cov.Proj.quarts.dim.M3.L4.90[1],Cov.Proj.quarts.dim.M3.L4.90[2])


    Cov.Proj.med.dim.M3.L5.90 <- median(Cov.Proj.M3.L5.dimensionality.90)
    Cov.Proj.quarts.dim.M3.L5.90 <-quantile(Cov.Proj.M3.L5.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 5`[13]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L5.90,Cov.Proj.quarts.dim.M3.L5.90[1],Cov.Proj.quarts.dim.M3.L5.90[2])


    Cov.Proj.med.dim.M4.L3.90 <- median(Cov.Proj.M4.L3.dimensionality.90)
    Cov.Proj.quarts.dim.M4.L3.90<-quantile(Cov.Proj.M4.L3.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 3`[13]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L3.90,Cov.Proj.quarts.dim.M4.L3.90[1],Cov.Proj.quarts.dim.M4.L3.90[2])


    Cov.Proj.med.dim.M4.L4.90 <- median(Cov.Proj.M4.L4.dimensionality.90)
    Cov.Proj.quarts.dim.M4.L4.90<-quantile(Cov.Proj.M4.L4.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 4`[13]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L4.90,Cov.Proj.quarts.dim.M4.L4.90[1],Cov.Proj.quarts.dim.M4.L4.90[2])


    Cov.Proj.med.dim.M4.L5.90 <- median(Cov.Proj.M4.L5.dimensionality.90)
    Cov.Proj.quarts.dim.M4.L5.90<-quantile(Cov.Proj.M4.L5.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 5`[13]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L5.90,Cov.Proj.quarts.dim.M4.L5.90[1],Cov.Proj.quarts.dim.M4.L5.90[2])



    #####Dimensions for 95%

    Cov.Proj.med.dim.M1.L00.95 <- median(Cov.Proj.M1.L00.dimensionality.95)
    Cov.Proj.quarts.dim.M1.L00.95<-quantile(Cov.Proj.M1.L00.dimensionality.95, c(.25,.75))

    results$`Model 1, l = ∞`[14]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M1.L00.95,Cov.Proj.quarts.dim.M1.L00.95[1],Cov.Proj.quarts.dim.M1.L00.95[2])


    Cov.Proj.med.dim.M2.L00.95 <- median(Cov.Proj.M2.L00.dimensionality.95)
    Cov.Proj.quarts.dim.M2.L00.95 <-quantile(Cov.Proj.M2.L00.dimensionality.95, c(.25,.75))

    results$`Model 2, l = ∞`[14]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M2.L00.95,Cov.Proj.quarts.dim.M2.L00.95[1],Cov.Proj.quarts.dim.M2.L00.95[2])


    Cov.Proj.med.dim.M3.L3.95 <- median(Cov.Proj.M3.L3.dimensionality.95)
    Cov.Proj.quarts.dim.M3.L3.95 <-quantile(Cov.Proj.M3.L3.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 3`[14]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L3.95,Cov.Proj.quarts.dim.M3.L3.95[1],Cov.Proj.quarts.dim.M3.L3.95[2])


    Cov.Proj.med.dim.M3.L4.95 <- median(Cov.Proj.M3.L4.dimensionality.95)
    Cov.Proj.quarts.dim.M3.L4.95<-quantile(Cov.Proj.M3.L4.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 4`[14]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L4.95,Cov.Proj.quarts.dim.M3.L4.95[1],Cov.Proj.quarts.dim.M3.L4.95[2])


    Cov.Proj.med.dim.M3.L5.95 <- median(Cov.Proj.M3.L5.dimensionality.95)
    Cov.Proj.quarts.dim.M3.L5.95 <-quantile(Cov.Proj.M3.L5.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 5`[14]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L5.95, Cov.Proj.quarts.dim.M3.L5.95[1], Cov.Proj.quarts.dim.M3.L5.95[2])


    Cov.Proj.med.dim.M4.L3.95 <- median(Cov.Proj.M4.L3.dimensionality.95)
    Cov.Proj.quarts.dim.M4.L3.95<-quantile(Cov.Proj.M4.L3.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 3`[14]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L3.95,Cov.Proj.quarts.dim.M4.L3.95[1],Cov.Proj.quarts.dim.M4.L3.95[2])


    Cov.Proj.med.dim.M4.L4.95 <- median(Cov.Proj.M4.L4.dimensionality.95)
    Cov.Proj.quarts.dim.M4.L4.95<-quantile(Cov.Proj.M4.L4.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 4`[14]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L4.95,Cov.Proj.quarts.dim.M4.L4.95[1],Cov.Proj.quarts.dim.M4.L4.95[2])


    Cov.Proj.med.dim.M4.L5.95 <- median(Cov.Proj.M4.L5.dimensionality.95)
    Cov.Proj.quarts.dim.M4.L5.95<-quantile(Cov.Proj.M4.L5.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 5`[14]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L5.95,Cov.Proj.quarts.dim.M4.L5.95[1],Cov.Proj.quarts.dim.M4.L5.95[2])


    #####Dimensions for 99%

    Cov.Proj.med.dim.M1.L00.99 <- median(Cov.Proj.M1.L00.dimensionality.99)
    Cov.Proj.quarts.dim.M1.L00.99<-quantile(Cov.Proj.M1.L00.dimensionality.99, c(.25,.75))

    results$`Model 1, l = ∞`[15]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M1.L00.99,Cov.Proj.quarts.dim.M1.L00.99[1],Cov.Proj.quarts.dim.M1.L00.99[2])


    Cov.Proj.med.dim.M2.L00.99 <- median(Cov.Proj.M2.L00.dimensionality.99)
    Cov.Proj.quarts.dim.M2.L00.99 <-quantile(Cov.Proj.M2.L00.dimensionality.99, c(.25,.75))

    results$`Model 2, l = ∞`[15]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M2.L00.99, Cov.Proj.quarts.dim.M2.L00.99[1], Cov.Proj.quarts.dim.M2.L00.99[2])


    Cov.Proj.med.dim.M3.L3.99 <- median(Cov.Proj.M3.L3.dimensionality.99)
    Cov.Proj.quarts.dim.M3.L3.99 <-quantile(Cov.Proj.M3.L3.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 3`[15]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L3.99, Cov.Proj.quarts.dim.M3.L3.99[1], Cov.Proj.quarts.dim.M3.L3.99[2])


    Cov.Proj.med.dim.M3.L4.99 <- median(Cov.Proj.M3.L4.dimensionality.99)
    Cov.Proj.quarts.dim.M3.L4.99<-quantile(Cov.Proj.M3.L4.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 4`[15]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L4.99, Cov.Proj.quarts.dim.M3.L4.99[1], Cov.Proj.quarts.dim.M3.L4.99[2])


    Cov.Proj.med.dim.M3.L5.99 <- median(Cov.Proj.M3.L5.dimensionality.99)
    Cov.Proj.quarts.dim.M3.L5.99 <-quantile(Cov.Proj.M3.L5.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 5`[15]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M3.L5.99, Cov.Proj.quarts.dim.M3.L5.99[1], Cov.Proj.quarts.dim.M3.L5.99[2])


    Cov.Proj.med.dim.M4.L3.99 <- median(Cov.Proj.M4.L3.dimensionality.99)
    Cov.Proj.quarts.dim.M4.L3.99<-quantile(Cov.Proj.M4.L3.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 3`[15]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L3.99, Cov.Proj.quarts.dim.M4.L3.99[1], Cov.Proj.quarts.dim.M4.L3.99[2])


    Cov.Proj.med.dim.M4.L4.99 <- median(Cov.Proj.M4.L4.dimensionality.99)
    Cov.Proj.quarts.dim.M4.L4.99<-quantile(Cov.Proj.M4.L4.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 4`[15]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L4.99, Cov.Proj.quarts.dim.M4.L4.99[1], Cov.Proj.quarts.dim.M4.L4.99[2])


    Cov.Proj.med.dim.M4.L5.99 <- median(Cov.Proj.M4.L5.dimensionality.99)
    Cov.Proj.quarts.dim.M4.L5.99<-quantile(Cov.Proj.M4.L5.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 5`[15]<- sprintf("%.3f (%.3f, %.3f)", Cov.Proj.med.dim.M4.L5.99, Cov.Proj.quarts.dim.M4.L5.99[1], Cov.Proj.quarts.dim.M4.L5.99[2])
  }

  #S4
  {
    ###Calculate relative errors

    excess.Proj.med.M1.L00<-median(excess.Proj.M1.L00.Err/normIV)
    excess.Proj.quarts.M1.L00<-quantile(excess.Proj.M1.L00.Err/normIV, c(.25,.75))

    results$`Model 1, l = ∞`[16]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.M1.L00,excess.Proj.quarts.M1.L00[1],excess.Proj.quarts.M1.L00[2])


    excess.Proj.med.M2.L00<-median(excess.Proj.M2.L00.Err/normIV)
    excess.Proj.quarts.M2.L00<-quantile(excess.Proj.M2.L00.Err/normIV, c(.25,.75))

    results$`Model 2, l = ∞`[16]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.M2.L00,excess.Proj.quarts.M2.L00[1],excess.Proj.quarts.M2.L00[2])


    excess.Proj.med.M3.L3<-median(excess.Proj.M3.L3.Err/normIV)
    excess.Proj.quarts.M3.L3<-quantile(excess.Proj.M3.L3.Err/normIV, c(.25,.75))

    results$`Model 3, l = 3`[16]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.M3.L3,excess.Proj.quarts.M3.L3[1],excess.Proj.quarts.M3.L3[2])


    excess.Proj.med.M3.L4<-median(excess.Proj.M3.L4.Err/normIV)
    excess.Proj.quarts.M3.L4<-quantile(excess.Proj.M3.L4.Err/normIV, c(.25,.75))

    results$`Model 3, l = 4`[16]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.M3.L4,excess.Proj.quarts.M3.L4[1],excess.Proj.quarts.M3.L4[2])


    excess.Proj.med.M3.L5<-median(excess.Proj.M3.L5.Err/normIV)
    excess.Proj.quarts.M3.L5<-quantile(excess.Proj.M3.L5.Err/normIV, c(.25,.75))

    results$`Model 3, l = 5`[16]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.M3.L5,excess.Proj.quarts.M3.L5[1],excess.Proj.quarts.M3.L5[2])


    excess.Proj.med.M4.L3<-median(excess.Proj.M4.L3.Err/normIV)
    excess.Proj.quarts.M4.L3<-quantile(excess.Proj.M4.L3.Err/normIV, c(.25,.75))

    results$`Model 4, l = 3`[16]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.M4.L3,excess.Proj.quarts.M4.L3[1],excess.Proj.quarts.M4.L3[2])


    excess.Proj.med.M4.L4<-median(excess.Proj.M4.L4.Err/normIV)
    excess.Proj.quarts.M4.L4<-quantile(excess.Proj.M4.L4.Err/normIV, c(.25,.75))

    results$`Model 4, l = 4`[16]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.M4.L4,excess.Proj.quarts.M4.L4[1],excess.Proj.quarts.M4.L4[2])



    excess.Proj.med.M4.L5<-median(excess.Proj.M4.L5.Err/normIV)
    excess.Proj.quarts.M4.L5<-quantile(excess.Proj.M4.L5.Err/normIV, c(.25,.75))

    results$`Model 4, l = 5`[16]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.M4.L5,excess.Proj.quarts.M4.L5[1],excess.Proj.quarts.M4.L5[2])


    #####Dimensions for 85%

    excess.Proj.med.dim.M1.L00.85 <- median(excess.Proj.M1.L00.dimensionality.85)
    excess.Proj.quarts.dim.M1.L00.85<-quantile(excess.Proj.M1.L00.dimensionality.85, c(.25,.75))

    results$`Model 1, l = ∞`[17]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M1.L00.85,excess.Proj.quarts.dim.M1.L00.85[1],excess.Proj.quarts.dim.M1.L00.85[2])


    excess.Proj.med.dim.M2.L00.85 <- median(excess.Proj.M2.L00.dimensionality.85)
    excess.Proj.quarts.dim.M2.L00.85 <-quantile(excess.Proj.M2.L00.dimensionality.85, c(.25,.75))

    results$`Model 2, l = ∞`[17]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M2.L00.85,excess.Proj.quarts.dim.M2.L00.85[1],excess.Proj.quarts.dim.M2.L00.85[2])


    excess.Proj.med.dim.M3.L3.85 <- median(excess.Proj.M3.L3.dimensionality.85)
    excess.Proj.quarts.dim.M3.L3.85 <-quantile(excess.Proj.M3.L3.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 3`[17]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L3.85,excess.Proj.quarts.dim.M3.L3.85[1],excess.Proj.quarts.dim.M3.L3.85[2])


    excess.Proj.med.dim.M3.L4.85 <- median(excess.Proj.M3.L4.dimensionality.85)
    excess.Proj.quarts.dim.M3.L4.85<-quantile(excess.Proj.M3.L4.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 4`[17]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L4.85,excess.Proj.quarts.dim.M3.L4.85[1],excess.Proj.quarts.dim.M3.L4.85[2])


    excess.Proj.med.dim.M3.L5.85 <- median(excess.Proj.M3.L5.dimensionality.85)
    excess.Proj.quarts.dim.M3.L5.85 <-quantile(excess.Proj.M3.L5.dimensionality.85, c(.25,.75))

    results$`Model 3, l = 5`[17]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L5.85,excess.Proj.quarts.dim.M3.L5.85[1],excess.Proj.quarts.dim.M3.L5.85[2])


    excess.Proj.med.dim.M4.L3.85 <- median(excess.Proj.M4.L3.dimensionality.85)
    excess.Proj.quarts.dim.M4.L3.85<-quantile(excess.Proj.M4.L3.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 3`[17]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L3.85,excess.Proj.quarts.dim.M4.L3.85[1],excess.Proj.quarts.dim.M4.L3.85[2])


    excess.Proj.med.dim.M4.L4.85 <- median(excess.Proj.M4.L4.dimensionality.85)
    excess.Proj.quarts.dim.M4.L4.85<-quantile(excess.Proj.M4.L4.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 4`[17]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L4.85,excess.Proj.quarts.dim.M4.L4.85[1],excess.Proj.quarts.dim.M4.L4.85[2])


    excess.Proj.med.dim.M4.L5.85 <- median(excess.Proj.M4.L5.dimensionality.85)
    excess.Proj.quarts.dim.M4.L5.85<-quantile(excess.Proj.M4.L5.dimensionality.85, c(.25,.75))

    results$`Model 4, l = 5`[17]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L5.85,excess.Proj.quarts.dim.M4.L5.85[1],excess.Proj.quarts.dim.M4.L5.85[2])



    #####Dimensions for 90%

    excess.Proj.med.dim.M1.L00.90 <- median(excess.Proj.M1.L00.dimensionality.90)
    excess.Proj.quarts.dim.M1.L00.90<-quantile(excess.Proj.M1.L00.dimensionality.90, c(.25,.75))

    results$`Model 1, l = ∞`[18]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M1.L00.90,excess.Proj.quarts.dim.M1.L00.90[1],excess.Proj.quarts.dim.M1.L00.90[2])


    excess.Proj.med.dim.M2.L00.90 <- median(excess.Proj.M2.L00.dimensionality.90)
    excess.Proj.quarts.dim.M2.L00.90 <-quantile(excess.Proj.M2.L00.dimensionality.90, c(.25,.75))

    results$`Model 2, l = ∞`[18]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M2.L00.90,excess.Proj.quarts.dim.M2.L00.90[1],excess.Proj.quarts.dim.M2.L00.90[2])


    excess.Proj.med.dim.M3.L3.90 <- median(excess.Proj.M3.L3.dimensionality.90)
    excess.Proj.quarts.dim.M3.L3.90 <-quantile(excess.Proj.M3.L3.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 3`[18]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L3.90,excess.Proj.quarts.dim.M3.L3.90[1],excess.Proj.quarts.dim.M3.L3.90[2])


    excess.Proj.med.dim.M3.L4.90 <- median(excess.Proj.M3.L4.dimensionality.90)
    excess.Proj.quarts.dim.M3.L4.90<-quantile(excess.Proj.M3.L4.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 4`[18]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L4.90,excess.Proj.quarts.dim.M3.L4.90[1],excess.Proj.quarts.dim.M3.L4.90[2])


    excess.Proj.med.dim.M3.L5.90 <- median(excess.Proj.M3.L5.dimensionality.90)
    excess.Proj.quarts.dim.M3.L5.90 <-quantile(excess.Proj.M3.L5.dimensionality.90, c(.25,.75))

    results$`Model 3, l = 5`[18]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L5.90,excess.Proj.quarts.dim.M3.L5.90[1],excess.Proj.quarts.dim.M3.L5.90[2])


    excess.Proj.med.dim.M4.L3.90 <- median(excess.Proj.M4.L3.dimensionality.90)
    excess.Proj.quarts.dim.M4.L3.90<-quantile(excess.Proj.M4.L3.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 3`[18]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L3.90,excess.Proj.quarts.dim.M4.L3.90[1],excess.Proj.quarts.dim.M4.L3.90[2])


    excess.Proj.med.dim.M4.L4.90 <- median(excess.Proj.M4.L4.dimensionality.90)
    excess.Proj.quarts.dim.M4.L4.90<-quantile(excess.Proj.M4.L4.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 4`[18]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L4.90,excess.Proj.quarts.dim.M4.L4.90[1],excess.Proj.quarts.dim.M4.L4.90[2])


    excess.Proj.med.dim.M4.L5.90 <- median(excess.Proj.M4.L5.dimensionality.90)
    excess.Proj.quarts.dim.M4.L5.90<-quantile(excess.Proj.M4.L5.dimensionality.90, c(.25,.75))

    results$`Model 4, l = 5`[18]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L5.90,excess.Proj.quarts.dim.M4.L5.90[1],excess.Proj.quarts.dim.M4.L5.90[2])



    #####Dimensions for 95%

    excess.Proj.med.dim.M1.L00.95 <- median(excess.Proj.M1.L00.dimensionality.95)
    excess.Proj.quarts.dim.M1.L00.95<-quantile(excess.Proj.M1.L00.dimensionality.95, c(.25,.75))

    results$`Model 1, l = ∞`[19]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M1.L00.95,excess.Proj.quarts.dim.M1.L00.95[1],excess.Proj.quarts.dim.M1.L00.95[2])


    excess.Proj.med.dim.M2.L00.95 <- median(excess.Proj.M2.L00.dimensionality.95)
    excess.Proj.quarts.dim.M2.L00.95 <-quantile(excess.Proj.M2.L00.dimensionality.95, c(.25,.75))

    results$`Model 2, l = ∞`[19]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M2.L00.95,excess.Proj.quarts.dim.M2.L00.95[1],excess.Proj.quarts.dim.M2.L00.95[2])


    excess.Proj.med.dim.M3.L3.95 <- median(excess.Proj.M3.L3.dimensionality.95)
    excess.Proj.quarts.dim.M3.L3.95 <-quantile(excess.Proj.M3.L3.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 3`[19]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L3.95,excess.Proj.quarts.dim.M3.L3.95[1],excess.Proj.quarts.dim.M3.L3.95[2])


    excess.Proj.med.dim.M3.L4.95 <- median(excess.Proj.M3.L4.dimensionality.95)
    excess.Proj.quarts.dim.M3.L4.95<-quantile(excess.Proj.M3.L4.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 4`[19]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L4.95,excess.Proj.quarts.dim.M3.L4.95[1],excess.Proj.quarts.dim.M3.L4.95[2])


    excess.Proj.med.dim.M3.L5.95 <- median(excess.Proj.M3.L5.dimensionality.95)
    excess.Proj.quarts.dim.M3.L5.95 <-quantile(excess.Proj.M3.L5.dimensionality.95, c(.25,.75))

    results$`Model 3, l = 5`[19]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L5.95, excess.Proj.quarts.dim.M3.L5.95[1], excess.Proj.quarts.dim.M3.L5.95[2])


    excess.Proj.med.dim.M4.L3.95 <- median(excess.Proj.M4.L3.dimensionality.95)
    excess.Proj.quarts.dim.M4.L3.95<-quantile(excess.Proj.M4.L3.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 3`[19]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L3.95,excess.Proj.quarts.dim.M4.L3.95[1],excess.Proj.quarts.dim.M4.L3.95[2])


    excess.Proj.med.dim.M4.L4.95 <- median(excess.Proj.M4.L4.dimensionality.95)
    excess.Proj.quarts.dim.M4.L4.95<-quantile(excess.Proj.M4.L4.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 4`[19]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L4.95,excess.Proj.quarts.dim.M4.L4.95[1],excess.Proj.quarts.dim.M4.L4.95[2])


    excess.Proj.med.dim.M4.L5.95 <- median(excess.Proj.M4.L5.dimensionality.95)
    excess.Proj.quarts.dim.M4.L5.95<-quantile(excess.Proj.M4.L5.dimensionality.95, c(.25,.75))

    results$`Model 4, l = 5`[19]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L5.95,excess.Proj.quarts.dim.M4.L5.95[1],excess.Proj.quarts.dim.M4.L5.95[2])


    #####Dimensions for 99%

    excess.Proj.med.dim.M1.L00.99 <- median(excess.Proj.M1.L00.dimensionality.99)
    excess.Proj.quarts.dim.M1.L00.99<-quantile(excess.Proj.M1.L00.dimensionality.99, c(.25,.75))

    results$`Model 1, l = ∞`[20]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M1.L00.99,excess.Proj.quarts.dim.M1.L00.99[1],excess.Proj.quarts.dim.M1.L00.99[2])


    excess.Proj.med.dim.M2.L00.99 <- median(excess.Proj.M2.L00.dimensionality.99)
    excess.Proj.quarts.dim.M2.L00.99 <-quantile(excess.Proj.M2.L00.dimensionality.99, c(.25,.75))

    results$`Model 2, l = ∞`[20]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M2.L00.99, excess.Proj.quarts.dim.M2.L00.99[1], excess.Proj.quarts.dim.M2.L00.99[2])


    excess.Proj.med.dim.M3.L3.99 <- median(excess.Proj.M3.L3.dimensionality.99)
    excess.Proj.quarts.dim.M3.L3.99 <-quantile(excess.Proj.M3.L3.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 3`[20]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L3.99, excess.Proj.quarts.dim.M3.L3.99[1], excess.Proj.quarts.dim.M3.L3.99[2])


    excess.Proj.med.dim.M3.L4.99 <- median(excess.Proj.M3.L4.dimensionality.99)
    excess.Proj.quarts.dim.M3.L4.99<-quantile(excess.Proj.M3.L4.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 4`[20]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L4.99, excess.Proj.quarts.dim.M3.L4.99[1], excess.Proj.quarts.dim.M3.L4.99[2])


    excess.Proj.med.dim.M3.L5.99 <- median(excess.Proj.M3.L5.dimensionality.99)
    excess.Proj.quarts.dim.M3.L5.99 <-quantile(excess.Proj.M3.L5.dimensionality.99, c(.25,.75))

    results$`Model 3, l = 5`[20]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M3.L5.99, excess.Proj.quarts.dim.M3.L5.99[1], excess.Proj.quarts.dim.M3.L5.99[2])


    excess.Proj.med.dim.M4.L3.99 <- median(excess.Proj.M4.L3.dimensionality.99)
    excess.Proj.quarts.dim.M4.L3.99<-quantile(excess.Proj.M4.L3.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 3`[20]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L3.99, excess.Proj.quarts.dim.M4.L3.99[1], excess.Proj.quarts.dim.M4.L3.99[2])


    excess.Proj.med.dim.M4.L4.99 <- median(excess.Proj.M4.L4.dimensionality.99)
    excess.Proj.quarts.dim.M4.L4.99<-quantile(excess.Proj.M4.L4.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 4`[20]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L4.99, excess.Proj.quarts.dim.M4.L4.99[1], excess.Proj.quarts.dim.M4.L4.99[2])


    excess.Proj.med.dim.M4.L5.99 <- median(excess.Proj.M4.L5.dimensionality.99)
    excess.Proj.quarts.dim.M4.L5.99<-quantile(excess.Proj.M4.L5.dimensionality.99, c(.25,.75))

    results$`Model 4, l = 5`[20]<- sprintf("%.3f (%.3f, %.3f)", excess.Proj.med.dim.M4.L5.99, excess.Proj.quarts.dim.M4.L5.99[1], excess.Proj.quarts.dim.M4.L5.99[2])
  }

  #for the screeplots
  {
    median.logdiff.explavar1 <- numeric(10)
    median.logdiff.explavar2 <- numeric(10)
    median.logdiff.explavar3 <- numeric(10)
    median.logdiff.explavar4 <- numeric(10)

    for (d in 1:10) {
      median.logdiff.explavar1[d]<-median(logdiff.explavar1[,d])
      median.logdiff.explavar2[d]<-median(logdiff.explavar2[,d])
      median.logdiff.explavar3[d]<-median(logdiff.explavar3[,d])
      median.logdiff.explavar4[d]<-median(logdiff.explavar4[,d])
    }

    median.cov.explavar1 <- numeric(10)
    median.cov.explavar2 <- numeric(10)
    median.cov.explavar3 <- numeric(10)
    median.cov.explavar4 <- numeric(10)

    for (d in 1:10) {
      median.cov.explavar1[d]<-median(Cov.explavar1[,d])
      median.cov.explavar2[d]<-median(Cov.explavar2[,d])
      median.cov.explavar3[d]<-median(Cov.explavar3[,d])
      median.cov.explavar4[d]<-median(Cov.explavar4[,d])
    }

    median.excess.explavar1 <- numeric(10)
    median.excess.explavar2 <- numeric(10)
    median.excess.explavar3 <- numeric(10)
    median.excess.explavar4 <- numeric(10)

    for (d in 1:10) {
      median.excess.explavar1[d]<-median(excess.explavar1[,d])
      median.excess.explavar2[d]<-median(excess.explavar2[,d])
      median.excess.explavar3[d]<-median(excess.explavar3[,d])
      median.excess.explavar4[d]<-median(excess.explavar4[,d])
    }
  }
}



View(results)

{
  plot(x=1:10, y=median.logdiff.explavar1, type = "b" ,main = "expl.var. of log-bond diff. pcs" , ylab = "explained var.")
  plot(x=1:10, y=median.cov.explavar1, type = "b" ,main = "expl.var. of log-bond pcs" , ylab = "explained var.")
  plot(x=1:10, y=median.excess.explavar1, type = "b" ,main = "expl.var. of excess ret. pcs" , ylab = "explained var.")

  E.correct<-eigen(q.50, only.values = TRUE)
  factorexpl<-(cumsum(E.correct$values)/sum(E.correct$values))[1:10]
  plot(x=1:10, y=factorexpl, type = "b" ,main = "expl.var.: Q vs log-price-diff. ret. pcs" , ylab = "explained var.")
  points(x=1:10, y=median.logdiff.explavar1, type = "b" , pch = 24)

  plot(x=1:10, y=factorexpl, type = "b" ,main = "expl.var.: Q vs log-price pcs" , ylab = "explained var.")
  points(x=1:10, y=median.cov.explavar1, type = "b" , pch = 24)

  plot(x=1:10, y=factorexpl, type = "b" ,main = "expl.var.: Q vs excess ret. pcs" , ylab = "explained var.")
  points(x=1:10, y=median.excess.explavar1, type = "b" , pch = 24)
}#Create Screeplots




write.csv(results, "results/simulation_results_leverage.csv", row.names = FALSE)





###For comparing the eigenvalue decays of q.50 and the r^2 q.50 + (1-r^2) v\%*\%t(v)
EXPL.VAR.TRUE<- cumsum(eigen(DATA$IV, only.values = TRUE)$values)/sum(eigen(DATA$IV, only.values = TRUE)$values)
EXPL.VAR.Q<- cumsum(eigen(q.50, only.values = TRUE)$values)/sum(eigen(q.50, only.values = TRUE)$values)

plot(EXPL.VAR.Q[1:10],
     ylab = "Explained Variation",
     xlab = "Factors",
     type = "b",
     pch = 20,
     bg = "black",
     col = "black",
     lwd = 1.5,
     cex = 1.1,
     cex.lab = 0.9,
     cex.axis = 0.8)
lines(EXPL.VAR.TRUE[1:10],
      type = "b",
      pch = 18,
      lwd = 1.5,
      cex = 1.1)

# Legend
legend("bottomright",
       legend = c(expression(r == 0), expression(r == -0.5)),
       pch = c(20, 18),
       pt.bg = c("black", "blue"),
       col = c("black", "black"),
       lty = 1,
       lwd = 1.5,
       pt.cex = 1.2,
       cex = 0.8,        # text size in legend
       bty = "n")        # remove legend box
#For checking the true threshholds of explained variation for r^2 q.50 + (1-r^2) v\%*\%t(v)

which.max(EXPL.VAR.TRUE>.85)
which.max(EXPL.VAR.TRUE>.90)
which.max(EXPL.VAR.TRUE>.95)
which.max(EXPL.VAR.TRUE>.99)
