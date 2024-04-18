######################################
# Projet Séries temporelles linéaires 
# Julia LU
######################################


############################
# Chargement des librairies
############################

library(zoo)
library(tseries)
library(forecast)
library(fUnitRoots)
library(stargazer)
library(xtable)


##########################
# Partie 1 : Les données
##########################

### Q1 - Présentation de la série

# Importation des données et renommage des colonnes
# Lien : https://www.insee.fr/fr/statistiques/serie/010767804
data <- read.csv("valeurs_mensuelles_pesticides.csv", sep=";", col.names = c('Dates', 'Indice', 'Codes'))

# On enlève les 3 premières lignes qui ne sont pas des données et on enlève la troisième colonne qui n'est pas utile
data <- data[-(1:3), 1:2]

# On réinitialise l'index du DataFrame
rownames(data) <- NULL

# On convertit les valeurs de la colonne 'Indice' en données numériques
data$Indice <- as.numeric(data$Indice)

# On définit les dates de la série
data$Dates[1] #première date : janvier 1990
tail(data$Dates, 1) #dernière date : février 2024
dates <- as.yearmon(seq(from=1990+0/12, to=2024+1/12, by=1/12))

# On créé la série temporelle associée aux valeurs prises par l'indice de production
indice.source <- zoo(data$Indice, order.by=dates)

# On supprime les 2 dernières valeurs pour la prédiction
indice <- indice.source[1:(length(indice.source)-2)]
dates2 <- dates[1:(length(dates)-2)]

# On trace la série 
png('./images/serie_initiale.png', width=600, height=450)
plot(indice, xlab = "Dates", ylab = "Indice de production industrielle", main="Fabrication de pesticides et d'autres produits agrochimiques" , col="blue")
dev.off()


### Q2 - Transformation de la série

# On régresse les valeurs de la série sur les dates et une constante pour vérifier si la série présente une tendance 
reg <- lm(indice ~ dates2)
summary(reg)

# La régression linéaire met en évidence une tendance croissante pour la série (coefficient significativement positif) et une constante significativement non nulle

# On effectue le test de racine unitaire ADF dans le cas avec constante et tendance pour déterminer si la série est stationnaire ou non
# L'hypothèse nulle est la non stationnarité (présence de racine unitaire)
adf <- adfTest(indice, lag=0, type="ct")
adf

# Cependant, on ne peut pas interpréter ce test car on ne sait pas s'il est valide
# Pour que le test soit valide, il faut que les résidus de la régression soient décorrélés

# On teste donc l’autocorrélation des résidus dans la régression .

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients)) 

# Toutes les p-valeurs sont en dessous du seuil de 5% donc l’hypothèse nulle d'absence d’autocorrélation des résidus est rejetée.
# Le test ADF avec aucun retard n’est donc pas valide. 

# On ajoute alors des lags jusqu’à obtenir des résidus décorrélés avec la fonction ci-dessous (cf TD5)

adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}
adf <- adfTest_valid(indice, 24, "ct")
adf

# On doit considérer 18 lags pour que les résidus ne soient plus autocorrélés 
# Le test ADF avec 18 lags est donc valide
# pval=0.6971 > 0.05 donc on ne rejette pas l'hypothèse nulle de non stationnarité 
# ie. la série initiale n'est pas stationnaire

# On va donc considérer la série différenciée à l'ordre 1
diff_indice <- diff(indice, 1)

# De même, on régresse la série différenciée sur les dates et une constante 
reg_diff <- lm(diff_indice ~ dates2[-1])
summary(reg_diff)

# Il n'y a pas de tendance ni de constante significatives (p-valeur > 0.05)

# On effectue le test ADF dans le cas "nc" (sans canstante et sans tendance) en contrôlant pour l'absence d'autocorrélation entre les résidus
adf_diff <- adfTest_valid(diff_indice, 24, "nc")
adf_diff

# Le test ADF avec 7 lags est valide
# pval=0.01 < 0.05 donc on rejette l'hypothèse de non stationnarité au seuil de 5%
# Conclusion : la série différenciée est stationnaire


### Q3 - Représentation graphique avant et après transformation

png('./images/serie_avant_apres.png', width=600, height=450)
plot(cbind(indice, diff_indice), xlab = "Dates", main="Série avant et après différenciation", col="blue")
dev.off()


##########################
# Partie 2 : Modèles ARMA
##########################

### Q4 - Choix du modèle ARMA(p,q) pour la série différenciée

# On détermine d'abord qmax et pmax grâce à l'acf et la pacf

# On trace l'autocorrélogramme de la série différenciée
png('./images/acf.png', width=600, height=450)
acf(diff_indice, main="Autocorrélogramme")
dev.off()

# L'ACF ne présente plus de pics significativement différents de zéro au delà de 2 retards (en comptant à partir du lag 0)
# On fait le choix d'ignorer les éventuels pics significatifs pour des retards supérieurs à 10
# On a donc qmax=2

# On trace l'autocorrélogramme partiel de la série différenciée
png('./images/pacf.png', width=600, height=450)
pacf(diff_indice, main="Autocorrélogramme partiel")
dev.off()

# La PACF ne présente plus de pics significatifs au delà de 7 retards (en comptant à partir du lag 1)
# On fait le choix d'ignorer les éventuels pics significatifs pour des retards supérieurs à 10
# On a donc pmax=7

# qmax and pmax
qmax <- 2
pmax <- 7

# On calcule le AIC/BIC pour chaque combinaison possible de p<=pmax et q<=qmax
pqs <- expand.grid(0:pmax,0:qmax) 
mat <- matrix(NA, nrow=pmax+1, ncol=qmax+1)
rownames(mat) <- paste0("p=",0:pmax) 
colnames(mat) <- paste0("q=",0:qmax) 
AICs <- mat
BICs <- mat 
for (row in 1:dim(pqs)[1]){
  p <- pqs[row,1]
  q <- pqs[row,2]
  estim <- try(arima(diff_indice,c(p,0,q), include.mean=F))
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim)
}

AICs
AICs==min(AICs)
xtable(AICs) #pour obtenir la table Latex

# Le modèle ARMA(0,2) minimise le critère AIC donc on garde ce modèle
arima012 <- arima(indice,c(0,1,2),include.mean=F)

BICs
BICs==min(BICs)
xtable(BICs) #pour obtenir la table Latex

# Le modèle ARMA(0,2) minimise également le critère BIC

# A présent, on va donc étudier l'ajustement et la validité du modèle ARMA(0,2) 

# Fonction de test des significations individuelles des coefficients (cf. TD4)
signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

signif(arima012)
# Les deux coefficients du modèle sont significatifs (pval<0.05) donc le modèle ARMA(0,2) est bien ajusté

# Tests d'absence d'autocorrélation des résidus 
qtest_arima012 <- Qtests(arima012$residuals,24,length(arima012$coef)-1)
qtest_arima012
xtable(qtest_arima012)
# Toutes les p-valeurs sont au dessus du seuil de 5% donc le modèle ARMA(0,2) est valide (résidus non autocorrélés)

# Conclusion : le modèle ARMA(0,2) est bien ajusté et valide

# On étudie également la normalité des résidus de notre modèle 
png('./images/residus_arima012.png', width=800, height=600)
checkresiduals(arima012)
dev.off()

# Diagramme Quantile-Quantile
png('./images/qqplot_residus_arima012.png', width=600, height=450)
qqnorm(arima012$residuals)
qqline(arima012$residuals)
dev.off()

# Test de normalité de Shapiro Wilk (H0 : les résidus sont gaussiens)
shapiro.test(arima012$residuals)
# pval=0.0037 < 0.01 donc on rejette H0 

# D'apprès le diagramme Q-Q et le test de Shapiro Wilk, l'hypothèse de normalité des résidus n'est pas vérifiée

# Représentation graphique de la série initiale et sa modélisation
png('./images/serie_initiale_et_modele012.png', width=600, height=450)
plot(indice, xlab="Dates", col="grey")
lines(fitted(arima012), col="red")
legend("topleft", legend = c(expression(X[t]), "ARIMA(0,1,2)"), col = c("grey", "red"), lty = 1)
dev.off()

########################
# Partie 3 : Prévision
########################

# On suppose tout de même que nos résidus sont gaussiens

# On estime l'écart-type de nos résidus
sigma <- sd(arima012$residuals)
sigma # sigma=8.6

# On prédit les valeurs de la série aux horizons T+1 et T+2 (janvier et février 2024) et les régions de confiance à 95% associées
forecast <- forecast(arima012, h=2, level=0.95)

# On représente les prédictions avec leur intervalle de confiance à 95% et la série initiale
# On limite la série à partir de l'année 2020 pour plus de lisibilité
png('./images/predictions_arima012.png', width=600, height=450)
plot(forecast, fcol="red", xlim=c(2020,2024+1/12), main="Prédictions et intervalles de confiance à 95% pour T+1 et T+2")
lines(x=c(2024, 2024+1/12), y=forecast$mean, col="red", lwd=2)
par(new=T)
plot(indice.source, xlim=c(2020,2024+1/12), xlab="Dates", ylab=expression(X[t]))
dev.off()

# Les points rouges correspondent prédictions et les zones en gris clair représentent les intervalles de confiance à 95%

