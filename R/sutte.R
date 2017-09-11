#################################################################################
##
##   R package sutteForecastR by Ansari Saleh Ahmar Copyright (C) 2015-2035.
##   This file is part of the R package sutteForecastR
##
##   The R package sutteForecastR is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package sutteForecastR is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

alpha.sutte <- function(x){

  if (!requireNamespace("forecast", quietly = TRUE)) {
    stop("forecast needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("fracdiff", quietly = TRUE)) {
    stop("fracdiff needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("robets", quietly = TRUE)) {
    stop("robets needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("forecastHybrid", quietly = TRUE)) {
    stop("forecastHybrid needed for this function to work. Please install it.",
         call. = FALSE)
  }

  al<-as.numeric(x)
  sutte<-c()
  for(i in 2:length(al)){
    sutte[i] <- al[i]-al[i-1]
  }
  sutte2<-c()
  for(i in 2:length(al)){
    sutte2[i] <- (al[i]+al[i-1])/2
  }
  sutte3<-c()
  for(i in 2:length(al)){
    sutte3[i] <- (sutte[i]/sutte2[i])
  }
  sutte4<-c()
  for(i in 2:length(al)){
    sutte4[i] <- (sutte3[i]*al[i])
  }
  sutte5<-c()
  for(i in 2:length(al)){
    sutte5[i+2] <- (sutte4[i]+sutte4[i+1]+sutte4[i+2])/3
  }
  sutteoke<-c()
  for(i in 2:(length(al)-3)){
    sutteoke[i+3-4] <- sutte5[i+2]+al[i+2]
  }

  panjang <- length(al)
  sutteokes<-c()
  for(i in (panjang-9):(panjang-3)){
    sutteokes[i+10-panjang] <- sutte5[i+2]+al[i+2]
  }

  als<-c()
  for(i in 1:(length(al))){
    als[i] <- al[i]
  }

  al_mi_10<-c()
  for(i in 1:(length(al)-7)){
    al_mi_10[i] <- al[i]
  }

  abserror<-c()
  for(i in 3:(length(al)-2)){
    abserror[i-2] <- abs(als[i+2]-sutteoke[i-2])
  }
  mse<-c()
  for(i in 3:(length(al)-2)){
    mse[i-2] <- (abs(als[i+2]-sutteoke[i-2]))^2
  }

  alass<-c()
  for(i in 2:(length(al)-3)){
    alass[i-1] <- al[i+3]
  }

  panjang <- length(al)
  alassss<-c()
  for(i in (panjang-6):(panjang)){
    alassss[i+7-panjang] <- al[i]
  }

  graphics.off()
  #par(mfrow=c(2,2))
  par(mfrow=c(2,3))

  # Plot the bar chart.
  cie <- plot(sutteoke,type = "o",pch = ".", col = "red", xlab = "Data", ylab = "", main = "Forecasts from Alpha-Sutte Indicator")
  cie2 <- lines(alass, type = "o", pch = ".",col = "blue")

  #AutoARIMA
  fit <- auto.arima(al_mi_10)
  foreauto <- forecast(fit,h=7)
  plotfore <- plot(foreauto)

  #HoltWinter
  fit2 <- HoltWinters(al_mi_10,gamma=FALSE)
  foreHW <- forecast(fit2,h=7)
  plotforeHW <- plot(foreHW)

  #Neural Network Time Series Forecasts
  fit3 <- nnetar(al_mi_10)
  foreNNETAR <- forecast(fit3,h=7)
  plotforeNNETAR <- plot(foreNNETAR)

  #Robust exponential smoothing model Forecasts
  fit4 <- robets(al_mi_10)
  forerobets <- forecast(fit4,h=7)
  plotforerobets <- plot(forerobets)

  #The theta model Forecasts
  fit5 <- thetam(al_mi_10)
  forethetam <- forecast(fit5,h=7)
  plotforethetam <- plot(forethetam)

  mae <- mean(abserror)
  mse <- mean(mse)
  rmse <- sqrt(mse)
  hasil=list(Tes_Data=alassss, Forecast_AlphaSutte=sutteokes, Forecast_AutoARIMA=foreauto,Forecast_HoltWinters=foreHW,Forecast_NNETAR=foreNNETAR, Forecast_Robust_exponential_smoothing=forerobets, Forecast_Theta=forethetam, AutoARIMA=fit, HoltWinters=fit2, NNETAR=fit3, Robust_exponential_smoothing=fit4, Theta_Model=fit5)
  return(hasil)
}
