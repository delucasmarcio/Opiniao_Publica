library(foreign)
library(lm.beta)
library(lavaan)
library(mediation)
library(cem)
library(mgcv)
library(randomForest)

am <- read.spss("863896541Brazil LAPOP AmericasBarometer 2014 Espanol v3.0_W.sav",
                reencode = TRUE, to.data.frame = TRUE)

am$sexo <- as.numeric(am$q1 == "Mujer")

am$idade <- 2014 - as.integer(as.character(am$q2y))

am$cor <- as.numeric(am$etid == "Blanco")

am$escola <- as.character(am$ed)
am$escola[am$escola == "Ninguno"] <- "0"
am$escola <- as.numeric(am$escola)

am$total_income <- as.numeric(am$q10new) - 1

am$renda <- am$total_income / as.integer(as.character(am$q12c))

am$conf <- 4 - as.numeric(am$it1)

am$vitima <- as.numeric(am$vic1ext == "Sí")

am$segur <- 4 - as.numeric(am$aoj11)

am$judiciario <- 4 - as.numeric(am$aoj12)

am$policia <- 4 - as.numeric(am$pole2n)

am$punitiv <- am$judiciario + am$policia

am$satisfdem <- 4 - as.numeric(am$pn4)

am <- am[is.na(am$conf) == FALSE &
         is.na(am$vitima) == FALSE &
         is.na(am$segur) == FALSE &
         is.na(am$punitiv) == FALSE &
         is.na(am$satisfdem) == FALSE, ]

am <- am[, c("sexo", "idade", "cor", "escola", "renda", "conf", "vitima",
             "segur", "punitiv", "satisfdem")]

### Completando valores omissos
am$cor[is.na(am$cor)] <- sample(am$cor[!is.na(am$cor)], sum(is.na(am$cor)), replace = TRUE)
am$escola[is.na(am$escola)] <- sample(am$escola[!is.na(am$escola)], sum(is.na(am$escola)),
                                      replace = TRUE)
am$renda[is.na(am$renda)] <- sample(am$renda[!is.na(am$renda)], sum(is.na(am$renda)),
                                    replace = TRUE)

### Modelos de Regressão

ols1 <- lm(satisfdem ~ conf + punitiv + segur + vitima + sexo +
           idade + cor + renda + escola, data = am)
summary(lm.beta(ols1))

ols2 <- lm(conf ~ punitiv + segur + vitima + sexo + idade +
           cor + renda + escola, data = am)
summary(lm.beta(ols2))

### Modelo de Equações Simultâneas

model <- "satisfdem ~ b1*conf + b2*punitiv + b3*segur + b4*vitima + sexo + idade + cor + renda + escola

          conf ~ c1*punitiv + c2*segur + c3*vitima + sexo + idade + cor + renda + escola

          ef_direto := b2 + b3 + b4

          ef_indireto := b1*(c1 + c2 + c3)

          ef_total := ef_direto + ef_indireto
         "

sem_model <- sem(model, data = am, estimator = "MLR",
                 missing = "ml", mimic = "Mplus")

summary(sem_model, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

### Modelo Quase completamente não paramétrico de Mediação Causal (Wood)
non_pa_dem  <- gam(satisfdem ~ s(conf, k = 4) + s(segur, k = 4) + s(punitiv, k = 7) + vitima,
                   fit = TRUE, data = am, method = "REML")

summary(non_pa_dem)

non_pa_tru  <- gam(conf ~ s(segur, k = 4) + s(punitiv, k = 7) + vitima,
                   fit = TRUE, data = am, method = "REML")

summary(non_pa_tru)

## Calculando F da Vitimização
f_tes <- gam(satisfdem ~ s(conf, k = 4) + s(segur, k = 4) + s(punitiv, k = 7),
             fit = TRUE, data = am, method = "REML")

vit_F <-  sum(f_tes$residuals ^ 2) - sum(non_pa_dem$residuals ^ 2)
vit_F <- vit_F * nrow(am) - 3 - 1
vit_F <- vit_F / sum(non_pa_dem$residuals ^ 2)
# F Vitimização: Modelo 1
print(vit_F)
#
f_tes <-  gam(conf ~ s(segur, k = 4) + s(punitiv, k = 7),
              fit = TRUE, data = am, method = "REML")

vit_F <-  sum(f_tes$residuals ^ 2) - sum(non_pa_tru$residuals ^ 2)
vit_F <- vit_F * nrow(am) - 3 - 1
vit_F <- vit_F / sum(non_pa_tru$residuals ^ 2)
# F Vitimização: Modelo 2
print(vit_F)

### Modelo Mediação Causal (Imai)
mediator <- non_pa_tru

general <- non_pa_dem

med_vit <- mediate(mediator, general, treat = "vitima", mediator = "conf",
                   robustSE = TRUE, sims = 100, boot = TRUE)

med_seg <- mediate(mediator, general, treat = "segur", mediator = "conf",
                   robustSE = TRUE, sims = 100, boot = TRUE)

med_pun <- mediate(mediator, general, treat = "punitiv", mediator = "conf",
                   robustSE = TRUE, sims = 100, boot = TRUE)

summary(med_vit)
summary(med_seg)
summary(med_pun)

### Modelo Coarsened Exact Matching (King)
b_conf <- c(-1, 0, 1, 2, 3)
b_punitiv <- c(-1, 0, 1, 2, 3, 4, 5, 6)
b_segur <- c(-1, 0, 1, 2, 3)

mat <- cem(treatment = "vitima", data = am, drop = "satisfdem",
           cutpoints = list(b_conf, b_punitiv, b_segur), k2k = TRUE)

est <- att(mat, satisfdem ~ conf + vitima + punitiv + segur,
           data = am, model = "rf")
est

### Técnicas de Learning
atr <- as.matrix(am[c(colnames(am)[1:9])])
atr <- (atr - colMeans(atr)) / apply(atr, 2, sd)

tar <- as.matrix(scale(am$satisfdem))
tar <- (tar - min(tar)) / (max(tar) - min(tar))
colnames(tar) <- "y"

## Random Forrest
simulations <- matrix(ncol = 9, nrow = 10)
for( i in 1:10){

  forest <- randomForest(y = tar, x = atr, importance = TRUE)
  var_imp <- importance(forest)[, 1]

  simulations[i, 1] <- var_imp[1]
  simulations[i, 2] <- var_imp[2]
  simulations[i, 3] <- var_imp[3]
  simulations[i, 4] <- var_imp[4]
  simulations[i, 5] <- var_imp[5]
  simulations[i, 6] <- var_imp[6]
  simulations[i, 7] <- var_imp[7]
  simulations[i, 8] <- var_imp[8]
  simulations[i, 9] <- var_imp[9]
}

dat_imp <- data.frame(sexo = mean(simulations[, 1]),
                      idade = mean(simulations[, 2]),
                      cor = mean(simulations[, 3]),
                      escola = mean(simulations[, 4]),
                      renda = mean(simulations[, 5]),
                      conf = mean(simulations[, 6]),
                      vitima = mean(simulations[, 7]),
                      segur = mean(simulations[, 8]),
                      punitiv = mean(simulations[, 9]))

names(dat_imp) <- sub("^sexo", "Sexo", names(dat_imp))
names(dat_imp) <- sub("^idade", "Idade", names(dat_imp))
names(dat_imp) <- sub("^cor", "Cor", names(dat_imp))
names(dat_imp) <- sub("^escola", "Escolaridade", names(dat_imp))
names(dat_imp) <- sub("^renda", "Renda", names(dat_imp))
names(dat_imp) <- sub("^conf", "Confiança", names(dat_imp))
names(dat_imp) <- sub("^vitima", "Vitimização", names(dat_imp))
names(dat_imp) <- sub("^segur", "Percep. segurança", names(dat_imp))
names(dat_imp) <- sub("^punitiv", "Conf. Punição", names(dat_imp))

dat_imp <- as.matrix(dat_imp)

dat_imp <- dat_imp[, order(dat_imp)]

par(mar = c(11, 6, 4, 2) + 0.1)

barplot(scale(dat_imp),
        names.arg = names(dat_imp),
        beside = TRUE,
        las = 2,
        ylab = "% Garson Relevância Padronizada",
        cex.names = 1.4,
        cex.axis = 1.4,
        cex.lab = 1.6)
