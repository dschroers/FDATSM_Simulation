################################################################################################
################################Replication Code for the Empirical Study########################
################################################ in ############################################
#######Dynamically Consistent Analysis of Realized Covariations in Term Structure Models########
#######################################by Dennis Schroers#######################################
################################################################################################
################################################################################################
################################################################################################


#if not already installed:
# install.packages("devtools")
#devtools::install_github("dschroers/FDATSM")
#install.packages("lubridate")
#install.packages("readr")

library(FDATSM)
library(lubridate)
library(readr)


################################################################################################
################################################################################################
#################################Functions######################################################
################################################################################################
################################################################################################
{
  d_star<-function(C,rho =.9){
    eigvals <- eigen(C, only.values = TRUE)$values
    scores <- cumsum(eigvals) / sum(eigvals)
    return(min(which(scores >= rho)))
  }
  variation <- function(Incr){
    return(t(Incr) %*% Incr)
  }

  quantile_truncation<-function(data,q){
    # Input checks
    if (!is.data.frame(data) && !is.matrix(data)) {
      stop("`data` must be a matrix or data frame.")
    }

    if (!all(sapply(data, is.numeric))) {
      stop("All columns in `data` must be numeric.")
    }

    if (!is.numeric(q) || length(q) != 1 || q <= 0 || q >= 1) {
      stop("`q` must be a single numeric value strictly between 0 and 1.")
    }

    if (anyNA(data)) {
      warning("Missing values detected; rows with NA values will be removed.")
      data <- data[complete.cases(data), , drop = FALSE]
    }


    # Compute Euclidean norm for each row
    norms.data <- apply(data, 1, function(row) sqrt(sum(row^2)))

    # Determine quantile threshold
    quant<- quantile(norms.data, probs = q)

    # Identify indices of points exceeding the threshold
    truncation.locations<-which(norms.data>=quant)
    return(truncation.locations)
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

  functional_data_truncation<-function(d, ### number of inverted eigenvalues
                                       C, ### stat. covariance matrix
                                       data, ### data matrix (rows = time)
                                       Delta, ###increment size
                                       sd ### number of 'standard deviations' truncated
  ){
    E<-eigen(C)

    # Input validation
    if ((d + 1) > length(E$values)) {
      stop("Error: The number of eigenvalues 'd' should be at most min(100, ncol(C)).")
    }

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
  truncated_covariation_diffretbased<-function (adj_data, tq = 0.75, l = 3, sparse = TRUE, sumplot = FALSE)
  {
    C.Prel <- preliminary_covariation_estimator(data = adj_data,
                                                tq)
    ft <- functional_data_truncation(d = d_star(C = C.Prel, tq),
                                     C = C.Prel, data = adj_data, Delta = 1/n, sd = l)
    locs <- ft$locations
    if (length(locs) == 0) {
      Truncated.variation <- variation(adj_data)
    }
    else {
      Truncated.variation <- variation(adj_data[-locs, ])
    }
    Truncated.variation <- Truncated.variation * (365^2)

    return(Truncated.variation)
  }
}#FDATSM equivalents


{

  FACTOR.APPROX <- function(FACTORS,  data) {
    #if(dim(FACTORS)[2]>1){
    scores <- data %*% FACTORS                    # Project onto factors (matrix multiplication)
    approximations <- scores %*% t(FACTORS)
    return(approximations)
  }




  random_sampler_of_dates <- function(DATES, size_per_year) {
    # Ensure DATES are of Date type
    DATES <- as.Date(DATES)

    # Extract year from each date
    years <- format(DATES, "%Y")

    # Split dates by year
    dates_by_year <- split(DATES, years)

    # Sample from each year
    sampled_dates <- lapply(dates_by_year, function(dates_in_year) {
      if (length(dates_in_year) < size_per_year) {
        stop(paste("Not enough dates in year", format(dates_in_year[1], "%Y")))
      }
      sample(dates_in_year, size_per_year)
    })

    # Combine sampled dates into a single vector
    return(do.call(c, sampled_dates))
  }

  MRISE<-function(dimen = 6, lag = 180, r.dates, DIFF.MEAN, LOG.MEAN){





    APPROX_DIFF_RET<-FACTOR.APPROX(FACTORS=E.DIFF, d= dimen, data =DIFF.RET, mvmean = DIFF.MEAN)
    APPROX.LOG.RET<-FACTOR.APPROX(FACTORS=E.LOG, d= dimen, data =LOG.RET, mvmean = LOG.MEAN)




    APPROX.LAGDIFF.LOG.RET<-matrix(0,nrow(LAGDIFF.RET),grid.length)
    APPROX.LAGDIFF.DIFF.RET<-matrix(0,nrow(LAGDIFF.RET),grid.length)
    for (l in 1:nrow(LAGDIFF.RET)) {
      for (p in 1:grid.length) {
        APPROX.LAGDIFF.DIFF.RET[l,p]<-(-sum(APPROX.DIFF.RET[l,grid[p]:(grid[p]+lag)]))
        APPROX.LAGDIFF.LOG.RET[l,p]<-APPROX.LOG.RET[l,grid[p]]-APPROX.LOG.RET[l,(grid[p]+(lag+1))]
      }
    }


    MRISE.LOG.ERRS<-numeric(nrow(LAGDIFF.RET))
    MRISE.DIFF.ERRS<-numeric(nrow(LAGDIFF.RET))

    for(l in 1:nrow(LAGDIFF.RET)){
      norm.lag<-L2.norm(LAGDIFF.RET[l,], from = 0, to = 10)
        MRISE.LOG.ERRS[l] <- L2.norm(APPROX.LAGDIFF.LOG.RET[l,]-LAGDIFF.RET[l,], from = 0, to = 10)/norm.lag
        MRISE.DIFF.ERRS[l] <- L2.norm(APPROX.LAGDIFF.DIFF.RET[l,]-LAGDIFF.RET[l,], from = 0, to = 10)/norm.lag
    }


    MeanRISE.LOG <-mean(MRISE.LOG.ERRS)
    MeanRISE.DIFF <-mean(MRISE.DIFF.ERRS)


    return(list("mean.rel.mrise.log" = MeanRISE.LOG, "mean.rel.mrise.diff" = MeanRISE.DIFF))
  }##function that calculates the integrated squared error

  }# aux functions





################################################################################################
################################################################################################
#################################Data Preprocessing#############################################
################################################################################################
################################################################################################

#Download the daily resolution yield Data from https://www.discount-bond-data.org and save as YIELDDATA.csv

df.fpy <- read.csv("YIELDDATA.csv", TRUE, ",")


#delete the metadata columns and convert yields to prices
price.df.fpy <- df.fpy[,-2]
for (i in 1:(ncol(price.df.fpy)-1)) {
  price.df.fpy[,1+i] <- exp(-(i/365)*price.df.fpy[,1+i])
}


#Reduce to the dates from "1990-01-01" and "2022-12-31" and maturities of maximally 3650 days (approx 10 years)
price.df.fpy<-price.df.fpy[as.Date(price.df.fpy[,1]) >= as.Date("1990-01-01") &
                             as.Date(price.df.fpy[,1]) <= as.Date("2022-12-31"),]
price.df.fpy<-price.df.fpy[,1:(3651+182)]

# calculate log prices
Log_Prices <- as.matrix(log(price.df.fpy[,-1]))

m=ncol(Log_Prices)
n=nrow(Log_Prices)

# calculate log price differences and difference returns

LOG_RET <- Log_Prices[2:n, 1:(m - 2)] - Log_Prices[1:(n - 1), 2:(m-1)]

DIFF_RET <- matrix(0, n - 1, m - 2)
for (i in 1:(n - 1)) {
  DIFF_RET[i, ] <- diff(Log_Prices[i + 1, 1:(m - 1)]) - diff(Log_Prices[i, 2:m])
}

################################################################################################
################################################################################################
###Study on "Importance of higher-order factors for short term trading strategies (Sec. 5.3)####
################################################################################################
################################################################################################


#Sample random dates
set.seed(123)
r.sample<-random_sampler_of_dates(DATES = as.Date(price.df.fpy[,1]), size_per_year = 25)

#define samples with deleted dates
LOG_RET_del <- LOG_RET[!(price.df.fpy[,1]%in%r.sample)[-8250],]
DIFF_RET_del <- DIFF_RET[!(price.df.fpy[,1]%in%r.sample)[-8250],]

#define data corresponding to out-of-sample dates
LOG_RET_oos <- LOG_RET[(as.Date(price.df.fpy[,1])%in%r.sample)[-8250],]
DIFF_RET_oos <- DIFF_RET[(as.Date(price.df.fpy[,1])%in%r.sample)[-8250],]






#Calculate a new version of the long term covariation leaving out the random dates
for(i  in 1: 33){
  start_date <- ymd("1989-01-01") + years(i)
  end_date <- ymd("1989-12-31") + years(i)

  temp <- DIFF_RET_del[(as.Date(price.df.fpy[,1]) >= as.Date(start_date) &
                         as.Date(price.df.fpy[,1]) <= as.Date(end_date))[-8250],]

C <-   truncated_covariation_diffretbased(adj_data = temp, sparse = FALSE)

  if(i==1){
    Clong_del<- C
  }else{
  Clong_del<-Clong_del+C
  }

  print(i/33)
}

E_DIFF<-eigen(Clong_del)$vectors

#Calculate the covariation of log price differences
COV<-t(LOG_RET_del)%*%LOG_RET_del

E_LOG<-eigen(COV)$vectors






#Calculate the higher order difference return at the randomly sampled dates
LAG_7_DIFF_RET_oos <- LOG_RET_oos[, (1 + 6 + 1):(3650 + 6 + 1)] - LOG_RET_oos[, 1:3650]
LAG_30_DIFF_RET_oos <- LOG_RET_oos[, (1 + 29 + 1):(3650 + 29 + 1)] - LOG_RET_oos[, 1:3650]
LAG_90_DIFF_RET_oos <- LOG_RET_oos[, (1 + 89 + 1):(3650 + 89 + 1)] - LOG_RET_oos[, 1:3650]
LAG_180_DIFF_RET_oos <- LOG_RET_oos[, (1 + 179 + 1):(3650 + 179 + 1)] - LOG_RET_oos[, 1:3650]


{
  # Define LaTeX-style row labels as a column
  latex_labels <- c(
    "Lag = 7, S1","Lag = 7, S2", "Lag = 30, S1", "Lag = 30, S2", "Lag = 90, S1","Lag = 90, S2", "Lag = 180, S1", "Lag = 180, S2"
  )

  col_names <- c( "d=1", "d=2","d=3","d=4","d=5","d=6","d=7","d=8","d=9","d=10","d=11","d=12","d=13","d=14","d=15","d=16")

  # Fill in with NA values
  results.53 <- data.frame(matrix(NA, nrow = length(latex_labels), ncol = length(col_names)))
  colnames(results.53) <- col_names

  # Add LaTeX row labels as a column
  results.53 <- cbind("Year" = latex_labels, results.53)
}#create a Dataframe containing the  results


for(d in 1:16){
  APPROX_DIFF_RET <- FACTOR.APPROX(FACTORS=E_DIFF[,1:d, drop = FALSE], data = DIFF_RET_oos)
APPROX_LOG_RET <- FACTOR.APPROX(FACTORS=E_LOG[,1:d, drop = FALSE], data = LOG_RET_oos)









#proxies for lag 7
APPROX_LAG_7_DIFF_LOG_BASIS<-matrix(0,825,3650)
APPROX_LAG_7_DIFF_DIFF_BASIS<-matrix(0,825,3650)

for(x in 1:825){
  for(l in 1:8){
    APPROX_LAG_7_DIFF_DIFF_BASIS[x,]<- APPROX_LAG_7_DIFF_DIFF_BASIS[x,]+APPROX_DIFF_RET[x,l:(3650+l-1)]
  }
}
APPROX_LAG_7_DIFF_LOG_BASIS<-APPROX_LOG_RET[,7:(3650+7-1)]-APPROX_LOG_RET[,1:3650]



#proxies for lag 30
APPROX_LAG_30_DIFF_LOG_BASIS<-matrix(0,825,3650)
APPROX_LAG_30_DIFF_DIFF_BASIS<-matrix(0,825,3650)

for(x in 1:825){
  for(l in 1:31){
    APPROX_LAG_30_DIFF_DIFF_BASIS[x,]<- APPROX_LAG_30_DIFF_DIFF_BASIS[x,]+APPROX_DIFF_RET[x,l:(3650+l-1)]
  }
}
APPROX_LAG_30_DIFF_LOG_BASIS<-APPROX_LOG_RET[,30:(3650+30-1)]-APPROX_LOG_RET[,1:3650]



#proxies for lag 90
APPROX_LAG_90_DIFF_LOG_BASIS<-matrix(0,825,3650)
APPROX_LAG_90_DIFF_DIFF_BASIS<-matrix(0,825,3650)

for(x in 1:825){
  for(l in 1:91){
    APPROX_LAG_90_DIFF_DIFF_BASIS[x,]<- APPROX_LAG_90_DIFF_DIFF_BASIS[x,]+APPROX_DIFF_RET[x,l:(3650+l-1)]
  }
}
APPROX_LAG_90_DIFF_LOG_BASIS<-APPROX_LOG_RET[,90:(3650+90-1)]-APPROX_LOG_RET[,1:3650]






#proxies for lag 180
APPROX_LAG_180_DIFF_LOG_BASIS<-matrix(0,825,3650)
APPROX_LAG_180_DIFF_DIFF_BASIS<-matrix(0,825,3650)

for(x in 1:825){
  for(l in 1:181){
    APPROX_LAG_180_DIFF_DIFF_BASIS[x,]<- APPROX_LAG_180_DIFF_DIFF_BASIS[x,]+APPROX_DIFF_RET[x,l:(3650+l-1)]
  }
}
APPROX_LAG_180_DIFF_LOG_BASIS<-APPROX_LOG_RET[,180:(3650+180-1)]-APPROX_LOG_RET[,1:3650]







# Calculate the RMAE


norms_7<-numeric(825)
log_err_7<-numeric(825)
diff_err_7<-numeric(825)

norms_30<-numeric(825)
log_err_30<-numeric(825)
diff_err_30<-numeric(825)

norms_90<-numeric(825)
log_err_90<-numeric(825)
diff_err_90<-numeric(825)

norms_180<-numeric(825)
log_err_180<-numeric(825)
diff_err_180<-numeric(825)

for (x in 1:825) {
  diff_err_7[x] <- L2_norm(APPROX_LAG_7_DIFF_DIFF_BASIS[x,]-LAG_7_DIFF_RET_oos[x,])
  log_err_7[x] <- L2_norm(APPROX_LAG_7_DIFF_LOG_BASIS[x,]-LAG_7_DIFF_RET_oos[x,])
  norms_7[x] <- L2_norm(LAG_7_DIFF_RET_oos[x,])

  diff_err_30[x] <- L2_norm(APPROX_LAG_30_DIFF_DIFF_BASIS[x,]-LAG_30_DIFF_RET_oos[x,])
  log_err_30[x] <- L2_norm(APPROX_LAG_30_DIFF_LOG_BASIS[x,]-LAG_30_DIFF_RET_oos[x,])
  norms_30[x] <- L2_norm(LAG_30_DIFF_RET_oos[x,])

  diff_err_90[x] <- L2_norm(APPROX_LAG_90_DIFF_DIFF_BASIS[x,]-LAG_90_DIFF_RET_oos[x,])
  log_err_90[x] <- L2_norm(APPROX_LAG_90_DIFF_LOG_BASIS[x,]-LAG_90_DIFF_RET_oos[x,])
  norms_90[x] <- L2_norm(LAG_90_DIFF_RET_oos[x,])

  diff_err_180[x] <- L2_norm(APPROX_LAG_180_DIFF_DIFF_BASIS[x,]-LAG_180_DIFF_RET_oos[x,])
  log_err_180[x] <- L2_norm(APPROX_LAG_180_DIFF_LOG_BASIS[x,]-LAG_180_DIFF_RET_oos[x,])
  norms_180[x] <- L2_norm(LAG_180_DIFF_RET_oos[x,])

}


rel_log_err_7 <-log_err_7/norms_7
rel_diff_err_7 <- diff_err_7/norms_7


rel_log_err_30 <-log_err_30/norms_30
rel_diff_err_30 <- diff_err_30/norms_30


rel_log_err_90 <-log_err_90/norms_90
rel_diff_err_90 <- diff_err_90/norms_90


rel_log_err_180 <-log_err_180/norms_180
rel_diff_err_180 <- diff_err_180/norms_180




RMAE_LAG_7_LOG <- mean(rel_log_err_7)
RMAE_LAG_7_DIFF <- mean(rel_diff_err_7)


RMAE_LAG_30_LOG <- mean(rel_log_err_30)
RMAE_LAG_30_DIFF <- mean(rel_diff_err_30)


RMAE_LAG_90_LOG <- mean(rel_log_err_90)
RMAE_LAG_90_DIFF <- mean(rel_diff_err_90)


RMAE_LAG_180_LOG <- mean(rel_log_err_180)
RMAE_LAG_180_DIFF <- mean(rel_diff_err_180)


results.53[1,d+1]<-RMAE_LAG_7_LOG
results.53[2,d+1]<-RMAE_LAG_7_DIFF

results.53[3,d+1]<-RMAE_LAG_30_LOG
results.53[4,d+1]<-RMAE_LAG_30_DIFF

results.53[5,d+1]<-RMAE_LAG_90_LOG
results.53[6,d+1]<-RMAE_LAG_90_DIFF

results.53[7,d+1]<-RMAE_LAG_180_LOG
results.53[8,d+1]<-RMAE_LAG_180_DIFF


print(d)
}

#Calculate how many factors are respectively needed to explain 99% of the variation in the data
EV_long<-eigen(Clong_del, only.values = TRUE)
EV_stat<-eigen(COV, only.values = TRUE)

crit_long <- which.max((cumsum(EV_long$values[1:16])/sum(EV_long$values))>.99)
crit_stat <- which.max((cumsum(EV_stat$values[1:16])/sum(EV_stat$values))>.99)
names(results.53)[crit_long+1] <- paste0("d=16**", names(df)[crit_long]) # mark the coloumn with the dimension chosen by long term volatility estimate
names(results.53)[crit_stat+1] <- paste0("d=2*", names(df)[crit_stat])# mark the coloumn with the dimension chosen by stationary covariance estimate



View(results.53)


write.csv(results.53, "results/empirical_results_outofsample.csv", row.names = FALSE)

