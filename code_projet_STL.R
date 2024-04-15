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

##########################
# Partie 1 : Les données
##########################

# Importation des données et renommage des colonnes
# Lien : https://www.insee.fr/fr/statistiques/serie/010767633
data <- read.csv("valeurs_mensuelles_fromage.csv", sep=";", col.names = c('Dates', 'Indice', 'Codes'))

# On enlève les 3 premières lignes et la colonne 'Codes' qui ne sont pas pertinentes
data <- data[-(1:3), 1:2]

# On réinitialise l'index du DataFrame
rownames(data) <- NULL

# On convertit les valeurs de la colonne 'Indice' en données numériques
data$Indice <- as.numeric(data$Indice)

# On définit les dates de la série
data$Dates[1] #première date : janvier 1990
tail(data$Dates, 1) #dernière date : décembre 2019
dates <- as.yearmon(seq(from=1990+0/12, to=2019+11/12, by=1/12))

# On créé la série temporelle associée aux valeurs prises par l'indice de production
prod.source <- zoo(data$Indice, order.by=dates)

# On supprime les 4 dernières valeurs pour la prédiction
prod <- prod.source[1:(length(prod.source)-4)]
plot(prod, xlab = "Année", ylab = "Indice de production industrielle")

# On créé la série différenciée
dprod <- diff(prod, 1)
plot(cbind(prod, dprod))
