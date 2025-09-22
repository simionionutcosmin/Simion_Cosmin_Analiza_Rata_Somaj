install.packages("ggplot2")
install.packages("forecast")
install.packages("tseries")
install.packages("urca")
install.packages("lmtest") 
install.packages("TSA")
install.packages("FinTS")

library(vars)
library(tseries)
library(stats)
library(dplyr)
#library(uroot)
library(TSA)
library(FinTS)
library(urca)              
library(ggplot2)
library(forecast)
library(tidyverse)
library(tseries)
library(readr)
library(forecast)
library(lmtest)

setwd("D:/Cosmin facultate/Serii de timp/proiect_seriitimp")

###############################################################################################################

## Prelucrare + analiza initiala fisier input

# Citire fișier
somaj <- read.csv("unemployment_romania.csv")

summary(somaj$Unemployment_Romania)

# Creează seria de timp DOAR cu coloana Unemployment_Romania
ts_somaj <- ts(somaj$Unemployment_Romania, start = c(2005, 5), frequency = 12)

# Window pentru întregul interval
ts_cv <- window(ts_somaj, start = c(2005, 5), end = c(2025, 3))

# Training până la sfârșit de 2023
training_cv <- window(ts_cv, start = c(2005, 5), end = c(2023, 12))

# Test de la 2024 încolo (ultimele 15 luni)
test_cv <- window(ts_cv, start = c(2024, 1))

# Plot
autoplot(training_cv) + labs(title = "Șomaj în România (2005-2025)", y = "Șomaj", x = "An")

# Plot sezonalitate
ggseasonplot(training_cv) +
  labs(title = "Sezonalitatea șomajului în România", y = "Șomaj", x = "Lună")

# Grafic ACF training + test
autoplot(training_cv) + labs(title = "Șomaj în România (2005-2025)", y = "Șomaj", x = "An")
ggseasonplot(training_cv) +
  labs(title = "Sezonalitatea șomajului în România", y = "Șomaj", x = "Lună")
ggAcf(training_cv)


# Descompunerea seriei in componente
descompunere_stl <- stl(ts_cv, s.window = "periodic")

autoplot(ts_cv, series = "Data") +
  autolayer(trendcycle(descompunere_stl), series = "Trend") +
  autolayer(seasadj(descompunere_stl), series = "Seasonally Adjusted") +
  xlab("An") + ylab("%") +
  ggtitle("Șomaj") +
  scale_colour_manual(values = c("green", "blue", "red"),
                      breaks = c("Data", "Seasonally Adjusted", "Trend"))


descompunere_stl %>% seasonal() %>% ggsubseriesplot() + ylab("Sezonier")
seasadj(descompunere_stl) %>% ggsubseriesplot() + ylab("Sezonier")

# Ar avea sezonalitate

ggseasonplot(ts_cv,  continuous=TRUE)

stl_decomp <- stl(training_cv, s.window="periodic")
plot(stl_decomp)


#Corelograma
ggtsdisplay(seasadj(descompunere_stl))


###############################################################################################################

## Stationaritatea seriei 

autoplot(training_cv) + 
  ggtitle('Evoluția șomajului') +
  theme_bw() 

#Test ADF
rcv_none <- ur.df(training_cv, type='none', selectlags = c("AIC"))
summary(rcv_none)

rcv_drift <- ur.df(training_cv, type='drift', selectlags = c("AIC"))
summary(rcv_drift) 

rcv_trend <- ur.df(training_cv, type='trend', selectlags = c("AIC"))
summary(rcv_trend) 

diff_training <- diff(training_cv)
ur.df(diff_training, type = "drift", selectlags = "AIC")


# => Seria este NEstationara in oricare din cele 3 variante => trebuie stationarizata

#Verificare daca nestationara DS in medie
mean_diff <- mean(diff(training_cv))
print(mean_diff) 

#=> DS în medie după o singură diferențiere (val apropiata de 0)

#Test Phillips-Perron si KPSS
PP.test(training_cv)
ts_cv %>% ur.kpss() %>% summary()


#ADF (toate 3: none, drift, trend): statistica testului > valorile critice ⇒ NU respinge H₀ ⇒ serie nestationară
#Phillips-Perron: p-value = 0.3562 ⇒ NU respinge H₀ ⇒ serie nestationară
#KPSS: test-statistic = 4.3541 > orice valoare critică ⇒ respinge H₀ ⇒ serie nestationară


ndiffs(training_cv) # nr de diferente de care ar fi nevoie

ggAcf(training_cv)
ggAcf(diff(training_cv)) # = > ar fi stationara dupa prima diferenta


#Prima diff
rcv_none2 <- ur.df(diff(training_cv), type='none', selectlags = c("AIC"))
summary(rcv_none2)

rcv_drift2 <- ur.df(diff(training_cv), type='drift', selectlags = c("AIC"))
summary(rcv_drift2)

rcv_trend2 <- ur.df(diff(training_cv), type='trend', selectlags = c("AIC"))
summary(rcv_trend2)

PP.test(diff(training_cv))
diff(training_cv) %>% ur.kpss() %>% summary()

#=> dupa prima diferenta seria devine stationara

###############################################################################################################

### APLICARE MODELE 

###############################################################################################################

## I. Tehnici de netezire - Holt Winters : sezonalitate

#Holt-Winters aditiv/multiplicativ = il alegem pe cel mai bun
fit_hw_ad <- hw(training_cv, seasonal = "additive")
fit_hw_mult <- hw(training_cv, seasonal = "multiplicative")

round(accuracy(fit_hw_ad), 2) 
round(accuracy(fit_hw_mult), 2)

#Valori un pic mai mici la HW aditiv (erori mai mici) - explicat si de sezonalitatea din grafic care pare stabila (nu creste odata cu trendul)

#Prognoza
forecast::accuracy(fit_hw_ad,h = length(test_cv))
autoplot(fit_hw_ad)


# Analiza reziduuri pt HW:

reziduri_hw_ad <- residuals(fit_hw_ad) 
autoplot(reziduri_hw_ad) + xlab("Luna") + ylab("") +
  ggtitle("Reziduri pentru HW aditiv")

# Testul J-B = Normalitate reziduuri

jarque.bera.test(reziduri_hw_ad) #=> reziduurile nu sunt normal distribuite               

ggAcf(reziduri_hw_ad) + ggtitle("ACF of residuals")

# Testul Box-Pierce = Autocorelatie reziduuri

Box.test(reziduri_hw_ad, lag=1)
Box.test(reziduri_hw_ad, lag=2)
Box.test(reziduri_hw_ad, lag=3)
Box.test(reziduri_hw_ad, lag=4)
Box.test(reziduri_hw_ad, lag=5)  
Box.test(reziduri_hw_ad, lag=10)  

# Testul Ljung-Box

Box.test(reziduri_hw_ad, lag=1,type="Lj")
Box.test(reziduri_hw_ad, lag=2, type="Lj")
Box.test(reziduri_hw_ad, lag=3, type="Lj")
Box.test(reziduri_hw_ad, lag=4, type="Lj")
Box.test(reziduri_hw_ad, lag=5, type="Lj")  
Box.test(reziduri_hw_ad, lag=10, type="Lj") 

checkresiduals(fit_hw_ad)

###############################################################################################################

## II. Modele ARMA-ARIMA-SARIMA

best_model<-auto.arima(training_cv)
summary(auto.arima(training_cv)) #=> ARIMA(1, 1, 0)
coeftest(auto.arima(training_cv))

# pe baza ACF si PACF => #=> ARIMA(1, 1, 1)
ggAcf(diff(training_cv)) # MA(1)
ggPacf(diff(training_cv)) #AR(1)


# Lista de modele (p, d=1, q) pe care le testam

ARIMA_1_1_0 = Arima(training_cv, order = c(1,1,0))
coeftest(ARIMA_1_1_0)
summary(ARIMA_1_1_0)

ARIMA_0_1_1 = Arima(training_cv, order = c(0,1,1))
coeftest(ARIMA_0_1_1)
summary(ARIMA_0_1_1)

ARIMA_1_1_1 = Arima(training_cv, order = c(1,1,1))
coeftest(ARIMA_1_1_1)
summary(ARIMA_1_1_1)

ARIMA_2_1_0 = Arima(training_cv, order = c(2,1,0))
coeftest(ARIMA_2_1_0)
summary(ARIMA_2_1_0)

ARIMA_2_1_1 = Arima(training_cv, order = c(2,1,1))
coeftest(ARIMA_2_1_1)
summary(ARIMA_2_1_1)

ARIMA_0_1_2 = Arima(training_cv, order = c(0,1,2))
coeftest(ARIMA_0_1_2)
summary(ARIMA_0_1_2)

ARIMA_3_1_0 = Arima(training_cv, order = c(3,1,0))
coeftest(ARIMA_3_1_0)
summary(ARIMA_3_1_0)

ARIMA_3_1_1 = Arima(training_cv, order = c(3,1,1))
coeftest(ARIMA_3_1_1)
summary(ARIMA_3_1_1)


#=> modelele cu cele mai bine valori: ARIMA(1,1,0) si ARIMA(3,1,1) 
# dar daca ne uitam la coeficienti, din cele doua doar ARIMA(1,1,0) are toti coef semnificativi


residuals_ARIMA_1_1_0 <- residuals(ARIMA_1_1_0)
ggtsdisplay(residuals_ARIMA_1_1_0) # nu avem autocorelatie in reziduuri conform ACF 
Box.test(residuals_ARIMA_1_1_0, lag=1) # nu avem autocorelare in reziduuri deoarece p > 0.1
Box.test(residuals_ARIMA_1_1_0, lag=8,type="Lj")
Box.test(residuals_ARIMA_1_1_0, lag=22,type="Lj")

ArchTest(residuals_ARIMA_1_1_0, lags = 1) # avem efecte ARCH
ArchTest(residuals_ARIMA_1_1_0, lags=12) 
ArchTest(residuals_ARIMA_1_1_0, lags=24) # nu avem efecte Garch
ArchTest(residuals_ARIMA_1_1_0, lags=36) # nu avem efecte Garch
ArchTest(residuals_ARIMA_1_1_0, lags=48) # nu avem efecte Garch

#Testul ARCH indică un efect marginal semnificativ doar la lag 1 (p = 0.015), 
#însă absența acestui efect pe intervale mai mari (lags = 12–48) sugerează lipsa persistenței varianței 
#condiționate, nejustificând utilizarea unui model GARCH.

jarque.bera.test(residuals_ARIMA_1_1_0) # distributie nenormala
checkresiduals(ARIMA_1_1_0) # din histograma se observa ca distributie este leptocurtică, asimetrie stângă



forecast_result <- forecast(ARIMA_1_1_0, h = length(test_cv))

autoplot(forecast_result) +
  autolayer(test_cv, series = "Test", color = "red") +
  ggtitle("Prognoză ARIMA (1,1,0)") +
  ylab("Valoare") + xlab("Timp")

accuracy(forecast_result, test_cv)



#Testare model ETS

ets_model <- ets(training_cv)
forecast_ets <- forecast(ets_model, h = length(test_cv))

autoplot(forecast_ets) + autolayer(test_cv, series="Test")
accuracy(forecast_ets, test_cv)



#Testare SARIMA
# Test Zivot-Andrews pe seria originală
# Test Zivot-Andrews pe seria cursului
za_test <- ur.za(training_cv, model = "both", lag = 4)  # "both" = intercept și trend

summary(za_test)

auto.arima(training_cv)

fit <- Arima(training_cv, order = c(1, 1, 1), seasonal = c(1, 0, 0))
summary(fit)
checkresiduals(fit)

forecast_auto <- forecast(auto.arima(training_cv),  h = length(test_cv))
autoplot(forecast_auto)

autoplot(forecast_auto) +
  autolayer(test_cv, series = "Test", color = "red") +
  labs(title = "Forecast vs Test Data",
       x = "Timp", y = "Valoare",
       subtitle = "Linie albastră = forecast, roșu = valori reale") +
  theme_minimal()

accuracy(forecast_auto, test_cv)


######################################################################################################
#Dieblo Mariano

round(forecast::accuracy(fit_hw_ad),2)
round(forecast::accuracy(ARIMA_1_1_0),2) 
round(forecast::accuracy(ets_model),2)

# Compararea modeluelor cu Diebold Mariano
dm.test(residuals(fit_hw_ad),residuals(ARIMA_1_1_0))
dm.test(residuals(fit_hw_ad),residuals(ets_model))
dm.test(residuals(ARIMA_1_1_0),residuals(ets_model)) 


# Forecasts for the three models
forecast_arima <- forecast(ARIMA_1_1_0, h = length(test_cv))
forecast_ets <- forecast(ets_model, h = length(test_cv))
#forecast_ad <- forecast::accuracy(fit_hw_ad,test_cv)

autoplot(ts_cv) + 
  # Add ARIMA forecast (blue)
  autolayer(forecast_arima$mean, series = "ARIMA(1,1,0)", color = "blue") +
  autolayer(forecast_ets$mean, series = "ETS", color = "darkgreen") +
  autolayer(test_cv, series = "Real", color = "red") +
  labs(title = "Comparative Forecast: ARIMA vs ETS vs HW Aditiv",
       x = "Time", y = "Value",
       color = "Series") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "darkgreen", "orange", "red"))

###############################################################################################################

  autoplot(ts_cv, series = "Curs real") +  # Adaugă nume și la seria de bază
  autolayer(forecast_arima$mean, series = "ARIMA(1,1,0)") +
  autolayer(forecast_ets$mean, series = "ETS") +
  autolayer(fit_hw_ad$mean, series = "HW Aditiv") +
  autolayer(test_cv, series = "Valori reale") +
  labs(
    title = "Comparative Forecast: ARIMA vs ETS vs HW Aditiv",
    x = "Time", y = "Value",
    color = "Serie"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "Curs real" = "black",
      "ARIMA(1,1,0)" = "blue",
      "ETS" = "darkgreen",
      "HW Aditiv" = "purple",
      "Valori reale" = "red"
    )
  )

  


###############################################################################################################
install.packages("ggfortify")
install.packages("tsDyn")




# --- Librării ---
library(vars)
library(tseries)
library(stats)
library(dplyr)
library(urca)              
library(ggplot2)
library(forecast)
library(tidyverse)
library(readr)
library(lmtest)
library(FinTS)
library(tsDyn)
library(ggfortify)  
library(tsDyn)

# --- Setare director lucru ---
setwd("D:/Cosmin facultate/Serii de timp/proiect_seriitimp")
  
# --- Citire fișiere ---
somaj <- read.csv("unemployment_romania.csv")


hicp <- read.csv("hicp.csv", stringsAsFactors = FALSE)

# Conversie la numeric și eliminare NA
hicp$HICP <- as.numeric(trimws(hicp$HICP))
hicp <- hicp[!is.na(hicp$HICP), ]
hicp <- data.frame(HICP = hicp)

# --- Serii de timp ---
ts_hicp <- ts(hicp$HICP, start = c(2005, 1), frequency = 12)
ts_somaj <- ts(somaj$Unemployment_Romania, start = c(2005, 1), frequency = 12)

# --- Matrice multivariată ---
df <- cbind(ts_hicp, ts_somaj)
colnames(df) <- c("HICP", "Somaj")

# --- Ploturi ---
autoplot(df) + 
  ylab('') + 
  ggtitle('Graficul seriei HICP și Șomaj') + 
  theme_bw()

# --- Test ADF pentru staționaritate ---
adf_hicp <- ur.df(ts_hicp, type = "trend", selectlags = "AIC")
summary(adf_hicp)

adf_somaj <- ur.df(ts_somaj, type = "trend", selectlags = "AIC")
summary(adf_somaj)

# --- Selectare lag optim ---
lagselect <- VARselect(df, lag.max = 8, type = 'const')
print(lagselect$selection)
lag_johansen <- max(2, lagselect$selection["AIC(n)"] - 1)


# --- Testul Johansen ---
johansen_trace <- ca.jo(df, type = 'trace', ecdet = 'const', K = lag_johansen)
summary(johansen_trace)

johansen_eigen <- ca.jo(df, type = 'eigen', ecdet = 'const', K = lag_johansen)
summary(johansen_eigen)

# --- Dacă există cointegrare: VECM ---
vecm_model <- VECM(df, lag = lag_johansen, r = 1, estim = "2OLS", LRinclude = "const")
summary(vecm_model)

vecm_var <- vec2var(johansen_trace, r = 1)

# --- Diagnostic VECM ---
resid_vecm <- residuals(vecm_var)
apply(resid_vecm, 2, function(x) Box.test(x, lag = 10, type = "Ljung-Box"))
normality.test(vecm_var, multivariate.only = TRUE)

# --- IRF ---
irf1 <- irf(vecm_var, impulse = "HICP", response = "Somaj", n.ahead = 20, boot = TRUE)
plot(irf1, ylab = "Șomaj", main = "Șoc HICP asupra șomajului")

irf2 <- irf(vecm_var, impulse = "Somaj", response = "HICP", n.ahead = 20, boot = TRUE)
plot(irf2, ylab = "HICP", main = "Șoc șomaj asupra HICP")

# --- Descompunere varianță ---
desc_var <- fevd(vecm_var, n.ahead = 10)
plot(desc_var)

# --- Prognoză ---
forecast <- predict(vecm_var, n.ahead = 4, ci = 0.90)
plot(forecast, name = "HICP")
plot(forecast, name = "Somaj")
fanchart(forecast, name = "HICP")
fanchart(forecast, name = "Somaj")



# --- Dacă nu există cointegrare: VAR pe diferențiate ---
diff_hicp <- diff(ts_hicp)
diff_somaj <- diff(ts_somaj)

# --- Test ADF pe diferențiate ---
adf_diff_hicp <- ur.df(diff_hicp, type = "drift", selectlags = "AIC")
summary(adf_diff_hicp)

adf_diff_somaj <- ur.df(diff_somaj, type = "drift", selectlags = "AIC")
summary(adf_diff_somaj)

# --- Matrice VAR ---
df_diff <- cbind(diff_hicp, diff_somaj)
colnames(df_diff) <- c("dHICP", "dSomaj")

# --- Selectare lag optim VAR ---
lagselect_diff <- VARselect(df_diff, lag.max = 8, type = 'const')
print(lagselect_diff$selection)
lag_opt <- lagselect_diff$selection["AIC(n)"]

# --- Estimare model VAR ---
model_VAR <- VAR(df_diff, p = lag_opt, type = "const")
summary(model_VAR)

# --- Diagnostic VAR ---
serial.test(model_VAR, lags.pt = 8, type = "PT.asymptotic")
arch.test(model_VAR, lags.multi = 12, multivariate.only = TRUE)
normality.test(model_VAR, multivariate.only = TRUE)
plot(stability(model_VAR, type = "OLS-CUSUM"))

# --- Cauzalitate Granger ---
Granger_HICP <- causality(model_VAR, cause = 'dHICP')
Granger_Somaj <- causality(model_VAR, cause = 'dSomaj')
Granger_HICP
Granger_Somaj

# --- IRF pe VAR ---
irf_var1 <- irf(model_VAR, impulse = "dHICP", response = "dSomaj", n.ahead = 20, boot = TRUE)
plot(irf_var1)

irf_var2 <- irf(model_VAR, impulse = "dSomaj", response = "dHICP", n.ahead = 20, boot = TRUE)
plot(irf_var2)

# --- Descompunere varianță pe VAR ---
plot(fevd(model_VAR, n.ahead = 10))







