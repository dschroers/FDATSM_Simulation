################################################################################################
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
price.df.fpy<-price.df.fpy[,1:3651]



################################################################################################
################################################################################################
###########Calculate the yearly Realized Covariations and the long-time Covariation#############
################################################################################################
################################################################################################









################################################################################################
################################################################################################
#####################Study on "Impact of jumps" (Sec. 5.1)#######################################
################################################################################################
################################################################################################

{
  # Define LaTeX-style row labels as a column
  latex_labels <- c(
    "1990", "1991", "1992", "1993", "1994", "1995", "1996", "1997", "1998", "1999",
    "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009",
    "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019",
    "2020", "2021", "2022"
  )

  col_names <- c( "#j,||q-||/||q-|| (l=3)", "#j,||q-||/||q-|| (l=4)",
                  "#j,||q-||/||q-|| (l=5)", "#||q-||")

  # Fill in with NA values
  results.51 <- data.frame(matrix(NA, nrow = length(latex_labels), ncol = length(col_names)))
  colnames(results.51) <- col_names

  # Add LaTeX row labels as a column
  results.51 <- cbind("Year" = latex_labels, results.51)
}#create a Dataframe containing the  results

#The yearly truncated realized covariations for l=3, needed for calculation of the long term covariation
Estimates.XC.l3<-list()
for(i in 1:33){
  start_date <- ymd("1989-01-01") + years(i)
  end_date <- ymd("1989-12-31") + years(i)

  temp<-price.df.fpy[as.Date(price.df.fpy[,1]) >= as.Date(start_date) &
                      as.Date(price.df.fpy[,1]) <= as.Date(end_date),]


  Estimates.XC.l3[[as.character(i)]] <- truncated_covariation(x=temp)

  Estimates.XC.l4.temp<-truncated_covariation(x=temp, l=4)
  Estimates.XC.l5.temp<-truncated_covariation(x=temp, l=5)
  Estimates.X.temp <- truncated_covariation(x=temp, l=1000)


  {
    j.l3<-length(Estimates.XC.l3[[as.character(i)]]$locs)
    j.l4<-length(Estimates.XC.l4.temp$locs)
    j.l5<-length(Estimates.XC.l5.temp$locs)

    mag.X<-L2_HS_norm(Estimates.X.temp$IV, from = 0, to = 10)

    rel.magn.l3 <- L2_HS_norm(Estimates.XC.l3[[as.character(i)]]$IV, from = 0, to = 10)/mag.X
    rel.magn.l4 <- L2_HS_norm(Estimates.XC.l4.temp$IV, from = 0, to = 10)/mag.X
    rel.magn.l5 <- L2_HS_norm(Estimates.XC.l5.temp$IV, from = 0, to = 10)/mag.X

    results.51$`#j,||q-||/||q-|| (l=3)`[i]<- sprintf("%.3f, %.3f", j.l3,rel.magn.l3)
    results.51$`#j,||q-||/||q-|| (l=4)`[i]<- sprintf("%.3f, %.3f", j.l4,rel.magn.l4)
    results.51$`#j,||q-||/||q-|| (l=5)`[i]<- sprintf("%.3f, %.3f", j.l5,rel.magn.l5)
    results.51$`#||q-||`[i]<- mag.X
  }#Results for Section 5.1: Analysis of Jumps


  print(i/33)
}


View(results.51)

write.csv(results.51, "results/empirical_results_jumps.csv", row.names = FALSE)



################################################################################################
################################################################################################
#####################Study on "Dimensionality" (Sec. 5.2)#######################################
################################################################################################
################################################################################################

#Estimate the long-time Covariation
Clong<-Estimates.XC.l3[[as.character(1)]]$IV

for(i in 2:33){
  Clong<-Clong+Estimates.XC.l3[[as.character(i)]]$IV
}
Clong<-Clong/33

E <- eigen(Clong)



{
  # Define LaTeX-style row labels as a column
  latex_labels <- c(
    "1990", "1991", "1992", "1993", "1994", "1995", "1996", "1997", "1998", "1999",
    "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009",
    "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019",
    "2020", "2021", "2022"
  )

  col_names <- c( "D^ei(0.85)", "D^ei(0.90)",
                  "D^ei(0.95)", "D^ei(0.99)","D^elong(0.85)", "D^elong(0.90)",
                  "D^elong(0.95)", "D^elong(0.99)")

  # Fill in with NA values
  results.52 <- data.frame(matrix(NA, nrow = length(latex_labels), ncol = length(col_names)))
  colnames(results.52) <- col_names

  # Add LaTeX row labels as a column
  results.52 <- cbind("Year" = latex_labels, results.52)
}#create a Dataframe containing the  results

{
  for(i in 1:33){


    De.85<-which.max(Estimates.XC.l3[[as.character(i)]]$expl.var>=.85)
    De.90<-which.max(Estimates.XC.l3[[as.character(i)]]$expl.var>=.9)
    De.95<-which.max(Estimates.XC.l3[[as.character(i)]]$expl.var>=.95)
    De.99<-which.max(Estimates.XC.l3[[as.character(i)]]$expl.var>=.99)

    XC.XC<-Estimates.XC.l3[[as.character(i)]]
    projections <- diag(t(E$vectors[,1:20]) %*% XC.XC$IV %*% E$vectors[,1:20])  # gives <Q_i e^long_i, e^long_i> for each i

    cum_sums <- cumsum(projections)/sum(eigen(XC.XC$IV, only.values = TRUE)$values)

    Delong.85<-which.max(cum_sums > .85)
    Delong.90<-which.max(cum_sums > .9)
    Delong.95<-which.max(cum_sums > .95)
    Delong.99<-which.max(cum_sums > .99)

    results.52$`D^ei(0.85)`[i]<-De.85
    results.52$`D^ei(0.90)`[i]<-De.90
    results.52$`D^ei(0.95)`[i]<-De.95
    results.52$`D^ei(0.99)`[i]<-De.99

    results.52$`D^elong(0.85)`[i]<-Delong.85
    results.52$`D^elong(0.90)`[i]<-Delong.90
    results.52$`D^elong(0.95)`[i]<-Delong.95
    results.52$`D^elong(0.99)`[i]<-Delong.99

    print(i/33)
  }
}#fill in the results


View(results.52)

write.csv(results.52, "results/empirical_results_dimensions.csv", row.names = FALSE)
