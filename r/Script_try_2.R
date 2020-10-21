# Title     : Script 2 con análisis comentados en la call
# Objective : Analizar el problema eliminando los valores del -1
# Created by: borja
# Created on: 17/10/2020

### Cambio de metodología - hay que buscar otra manera de hacer esto
### Quitar los -1 y hacer estudio independiente (variable unidimensional)
### Ajuste para el resto de valores según se crea conveniente.
    ## Regresión quantil
    ## Ajustes no paramétricos y paramétricos (incluir ICs)
    ## Ver resultados
### Usar estimaciones de la desviación típica condicional.
require(ggplot2)
require(pscl)
require(boot)
require(VGAM)
require(gamlss)
require(gamlss.dist)
require(locfit)
require(KernSmooth)
require(MASS)
###
data <- na.omit(read.csv('output/lsu_df.csv', sep = ',', dec = '.', header = T))
dim(data)
head(data)
data$CVTOT <- as.numeric(as.character(data$CVTOT))
data$GB.Seqs <- na.omit(as.numeric(as.character(data$GB.Seqs)))
data <- na.omit(data)
data_neg <- data[data$CVTOT == -1,]
data_filtered <- data[data$CVTOT != -1,]
# Visualización de los hist de data_filtered
sum(data_filtered$GB.Seqs == 0)
sum(data_filtered$GB.Seqs != 0)
hist(data_filtered$GB.Seqs, xlab = 'Apariciones', breaks = 'fd', freq = F)
hist(data_filtered$CVTOT[data_filtered$GB.Seqs != 0],
     breaks = 'fd', freq = F, xlab = 'CVindex Apariciones > 0')
hist(data_filtered$CVTOT[data_filtered$GB.Seqs == 0],
     breaks = 'fd', freq = F, xlab = 'CVindex Apariciones == 0')
summary(data_neg)
summary(data_filtered)
plot(data_filtered$CVTOT,data_filtered$GB.Seqs,pch=1,ylim = c(0,30),
     xlab="CVTOT",ylab="Appearances")
x <- unique(data_filtered$CVTOT)
y <- data_filtered$GB.Seqs
## Regresión quantil

## Conditional variance - modo univariante pensar en el modo multi
con_var <- function(df, residuos, col_y, col_x, new_data, FUN = glm, ...)
{
  df$con_var_up <- df[,col_y]+residuos^2
  df$con_var_down <- df[,col_y] - residuos^2
  print(length(df[,col_x]))
  print(length(df$con_var_up))
  model_adjust_up <- FUN(df[,'con_var_up'] ~ df[,col_x], ...)
  model_adjust_down <- FUN(df[,'con_var_down'] ~ df[,col_x], ...)
  return(c(sqrt(predict(model_adjust_up, newdata = new_data)),
           sqrt(predict(model_adjust_down, newdata= new_data))))
}

## Estudios paramétricos:
### GLM
## Poisson
## Geometrica
## BN
# Poisson + ZerosInflated
## Caso 1: Datos completos
## Caso 2: Datos filtrados al 1% superior
## Caso 3: Datos filtrados al 5$ superior
l_a <- c(0, 0.01, 0.05)
for (a in l_a)
{
  # Eliminamos datos fuera de lugar
  lower_bound <- quantile(y, 0.00);lower_bound
  upper_bound <- quantile(y, 1 - a);upper_bound

  outlier_ind <- which(y < lower_bound | y > upper_bound)
  if (!is.na(outlier_ind[1]))
    data_tmp <- data_filtered[-outlier_ind,]
  else
    data_tmp <- data_filtered

  # Comprobamos datos
  summary(data_tmp)
  head(data_tmp)
  dim(data_tmp)
  # Analisis visual
  windows()
  hist(data_tmp$GB.Seqs, freq = F, breaks = 'fd')

  windows()
  plot(data_tmp[,'CVTOT'],data_tmp[,'GB.Seqs'],
       pch=16,col="blue",ylim=c(-2,100),xlab="CvIndex",ylab="Appearances",main=paste("Filtro: ",a))
  ## Regresión paramétrica - poisson (sistema de recuento, con valores > 0, makes sense) - sin tocar 0s
  RegZero_s_inf <- glm(GB.Seqs ~ CVTOT , data = data_tmp,family = poisson)
  RegZero_s_inf_theta <- glm(GB.Seqs ~ CVTOT,  data = data_tmp, family = quasipoisson(link = "log"))
  summary(RegZero_s_inf_theta)
  summary(RegZero_s_inf)
  x<-seq(min(data_tmp[,'CVTOT']),max(data_tmp[,'CVTOT']),by=0.01)
  Pr_s_inf<-predict(RegZero_s_inf, newdata=data.frame("CVTOT"=x),se.fit=TRUE)
  lines(x,Pr_s_inf$fit,col="green",type="l",lwd=2)
  # El ajuste va malament
  f1 <- formula(GB.Seqs ~ CVTOT | 1 ) # Los 0's son estructurales
  f2 <- formula(GB.Seqs ~ CVTOT | CVTOT) # Los 0s dependen de CVTOT
  f3 <- formula(GB.Seqs ~ CVTOT)
  RegZero_inf_sP_c1 <- zeroinfl(f1 ,data= data_tmp, dist="poisson", link = "logit")
  RegZero_inf_sP_c2 <- zeroinfl(f2, data = data_tmp, dist="poisson", link = "logit")
  RegZero_inf_sP_c3 <- zeroinfl(f3, data = data_tmp, dist="poisson", link = "logit")
  RegZero_inf_sP_def <- zeroinfl(f3, data = data_tmp)
  RegZero_inf_sP_geo1 <- zeroinfl(f1, data = data_tmp, dist = "geometric")
  RegZero_inf_sP_geo2 <- zeroinfl(f2, data = data_tmp, dist = "geometric")
  RegZero_inf_sP_ngbin_1 <- zeroinfl(f1, data = data_tmp, dist = "negbin")
  RegZero_inf_sP_ngbin_2 <- zeroinfl(f2, data = data_tmp, dist = "negbin")
  summary(RegZero_inf_sP_c1)
  summary(RegZero_inf_sP_c2)
  summary(RegZero_inf_sP_c3)
  summary(RegZero_inf_sP_def)
  summary(RegZero_inf_sP_geo1)
  summary(RegZero_inf_sP_geo2)
  summary(RegZero_inf_sP_ngbin_1)
  summary(RegZero_inf_sP_ngbin_2)
  Pr_inf_sP_c1<-predict(RegZero_inf_sP_c1,
                        newdata=data.frame("CVTOT"=x),se.fit=TRUE, type = "response")
  Pr_inf_sP_c2<-predict(RegZero_inf_sP_c2,
                        newdata=data.frame("CVTOT"=x),se.fit=TRUE,type = "response")
  Pr_inf_sP_c3<-predict(RegZero_inf_sP_c3,
                        newdata=data.frame("CVTOT"=x),se.fit=TRUE,type = "response")
  Pr_inf_sP_def<-predict(RegZero_inf_sP_def,
                         newdata=data.frame("CVTOT"=x),se.fit=TRUE,type = "response")
  Pr_inf_sP_geo1<-predict(RegZero_inf_sP_geo1,
                         newdata=data.frame("CVTOT"=x),se.fit=TRUE,type = "response")
  Pr_inf_sP_geo2<-predict(RegZero_inf_sP_geo2,
                         newdata=data.frame("CVTOT"=x),se.fit=TRUE,type = "response")
  Pr_inf_sP_ngbin_1<-predict(RegZero_inf_sP_ngbin_1,
                         newdata=data.frame("CVTOT"=x),se.fit=TRUE,type = "response")
  Pr_inf_sP_ngbin_2<-predict(RegZero_inf_sP_ngbin_2,
                         newdata=data.frame("CVTOT"=x),se.fit=TRUE,type = "response")
  # Poisson - blue
  # Geometric - red
  # Default - pink
  # NegBin - orange
  lines(x,Pr_inf_sP_c1,col="blue",type="l",lwd=2,lty=1)
  lines(x,Pr_inf_sP_c2,col="blue",type="l",lwd=2,lty=2)
  lines(x,Pr_inf_sP_c3,col="blue",type="l",lwd=2,lty=3)
  lines(x,Pr_inf_sP_def,col="pink",type="l",lwd=2,lty=1)
  # Curiosamente estas dos no van tan mal
  lines(x,Pr_inf_sP_geo1,col="red",type="l",lwd=2,lty=1)
  lines(x,Pr_inf_sP_geo2,col="red",type="l",lwd=2,lty=2)
  # BinomNeg
  lines(x,Pr_inf_sP_ngbin_1,col="orange",type="l",lwd=2,lty=1)
  lines(x,Pr_inf_sP_ngbin_2,col="orange",type="l",lwd=2,lty=2)
  legend("topleft", c("Poisson 1","Poisson 2","Poisson 3","Default"
    ,"Geometric 1","Geometric 2", "Binom 1", "Binom 2"),
         lty=c(1,2,3,1,1,2,1,2),col=c(rep("blue",3),"pink",rep("red",2),rep("orange",2)),
         ncol=1, lwd=rep(2,7), inset = .05)
  # Comparación de los modelos
  AIC(RegZero_inf_sP_c1,RegZero_inf_sP_c2,RegZero_inf_sP_c3,
      RegZero_inf_sP_def,RegZero_inf_sP_geo1,RegZero_inf_sP_geo2,RegZero_inf_sP_ngbin_1,
      RegZero_inf_sP_ngbin_2)
  # Una de alternativa suele ser aplicar el logaritmo de los datos
  data_tmp$GB.Seqs.log <- as.numeric(log(data_tmp$GB.Seqs + 1))
  summary(data_tmp$GB.Seqs.log)
  windows()
  plot(data_tmp$CVTOT, data_tmp$GB.Seqs.log, col = "blue",
       xlab = "CvTot",ylab = "log(Appearances) + 1", main = paste("Log appearances, ",a))
  lines(x,log(Pr_inf_sP_c1+1),col="blue",type="l",lty =1, lwd = 2)
  lines(x,log(Pr_inf_sP_c2+1),col="blue",type="l",lty = 2, lwd = 2)
  lines(x,log(Pr_inf_sP_c3+1),col="blue",type="l",lty = 3, lwd = 2)
  lines(x, log(Pr_inf_sP_def+1), col="pink", type="l", lty = 1,lwd = 2)
  lines(x,log(Pr_inf_sP_geo1+1),col="red",type="l",lty = 1,lwd = 2)
  lines(x,log(Pr_inf_sP_geo1+1),col="red",type="l",lty = 2,lwd = 2)
  lines(x,log(Pr_inf_sP_ngbin_1+1),col="orange",type="l",lty = 1, lwd = 2)
  lines(x,log(Pr_inf_sP_ngbin_2+1),col="orange",type="l",lty = 2, lwd = 2)
  legend("topleft", c("Poisson 1","Poisson 2","Poisson 3","Default"
    ,"Geometric 1","Geometric 2", "Binom 1", "Binom 2"),
         lty=c(1,2,3,1,1,2,1,2),col=c(rep("blue",3),"pink",rep("red",2),rep("orange",2)),
         ncol=1, lwd=rep(2,7), inset = .05)
}

## Modelos paramétricos via GAMLSS (Generalize Aditive models for location scale and shape)
l_a <- c(0.0, 0.01)
for (a in l_a){
  # Eliminamos datos fuera de lugar
  lower_bound <- quantile(y, 0.00);lower_bound
  upper_bound <- quantile(y, 1 - a);upper_bound

  outlier_ind <- which(y < lower_bound | y > upper_bound)
  if (!is.na(outlier_ind[1])){
    data_tmp <- data_filtered[-outlier_ind,]
  }else{
    data_tmp <- data_filtered}
  ## Estudio de las distribuciones
  windows()
  op <- par(mfrow = c(2, 3))
  dPO <- histDist(data_tmp$GB.Seqs, "PO",
                  main = "PO", trace = FALSE, xlim=c(0,20))     # Poisson
  dNBI <- histDist(data_tmp$GB.Seqs, "NBI",
                   main = "NBI", trace = FALSE,xlim=c(0,20))   # Binomial negativa tipo 1
  dPIG <- histDist(data_tmp$GB.Seqs, "PIG",
                   main = "PIG", trace = FALSE,xlim=c(0,20))   # Poisson inverse gaussian
  dSI <- histDist(data_tmp$GB.Seqs, "SICHEL",
                  main = "SICHEL", trace = FALSE, xlim=c(0,20)) # Sichel? Inverse Gaussian distribution
  dZIP <- histDist(data_tmp$GB.Seqs, "ZIP",
                   main = "ZIP", trace = FALSE, xlim=c(0,20))  # Zero inflated poisson
  dZIP2 <- histDist(data_tmp$GB.Seqs, "ZIP2",
                    main = "ZIP2", trace = FALSE, xlim=c(0,20)) # Zero inf. pois. 2
  par(op)
  AIC(dPO, dNBI, dPIG, dSI, dZIP, dZIP2)

  ## En base al AIC nos quedaríamos con dPIG que es un modelo de poisson pero descartariamos claramente los de ceros
  ## inflados.

  ## A través del ajuste de gamlss ajustamos distintos modelos con diferentes familias de funciones
  ## Eliminapos PO, ZIP y ZIP2 dado que sus AIC están fuera de lugar en comparación a los otros 3
  ?gamlss
  fam.gamlss <- c("NBI", "NBII", "PIG", "DEL", "SICHEL")
  m.l <- m.q <- m.s <- m.lq <-list()
  ## Ajsute lineal con enlaces
  for (i in 1:5) {
    m.l[[fam.gamlss[i]]] <- gamlss(GB.Seqs~ CVTOT, data=na.omit(data_tmp),
                                   family = fam.gamlss[i], n.cyc = 60, trace = FALSE)$aic
  }
  ## Regresión polinómico de grado 2 con los enlaces
  for (i in 1:5) {
    m.q[[fam.gamlss[i]]] <-GAIC(  gamlss(GB.Seqs ~ poly(CVTOT,2),
                                         data=data_tmp, family = fam.gamlss[i],
                                         n.cyc = 60, trace = FALSE))
  }
  ## Smooth cubic splines (suavización via splines cúbicos)
  for (i in 1:5) {
    m.s[[fam.gamlss[i]]] <-GAIC(  gamlss(GB.Seqs ~ cs(CVTOT),
                                         data=na.omit(data_tmp), family = fam.gamlss[i],
                                         n.cyc = 60, trace = FALSE))
  }
  unlist(m.l)
  unlist(m.q)
  unlist(m.s)

  ## Igual que se visualizaba anteriormente, los mejores resultados vienen de SICHEL y de PIG (Igual
  # a como ocurre en la versión de Manuel)

  # sigma.fo permite otorgar una formula para ajustar un modelo para el parámetro sigma (desviación)
  m.ql <- list()
  for (i in 1:5) {
    m.ql[[fam.gamlss[i]]] <- GAIC(gamlss(GB.Seqs ~ poly(CVTOT, 2),
                                         data=na.omit(data_tmp), sigma.fo = ~CVTOT,
                                         family = fam.gamlss[i], n.cyc = 60, trace = FALSE))
  }
  unlist(m.ql)

  # Nuevamente apenas hay mejoría.
  ## En cuanto al "mejor" modelo:
  ##  PIG - Poisson Inverse Gaussian o SICHEL (no existen diferencias significativas en el AIC)
  ##  El ajuste vía smoothing spline - aunque el ajuste por regresión polinómico local (grado 2)
  # tampoco queda descartado. - Nuevamente lo mismo que para Manuel.

  ## En base a esto vamos a evaluar los resultados via plot
  # Lineal + ajuste de sd
  mL <- gamlss(GB.Seqs~CVTOT , data=na.omit(data_tmp), family =PIG,
               sigma.fo=~CVTOT, n.cyc = 60, trace = FALSE)
  # Cuadrático + ajuste de sd
  mC <- gamlss(GB.Seqs~CVTOT+I(CVTOT^2) , data=na.omit(data_tmp), family =PIG,
               sigma.fo=~CVTOT, n.cyc = 60, trace = FALSE)
  # Vía smoothing spline + ajuste de sd
  mS <- gamlss(GB.Seqs~cs(CVTOT) , data=na.omit(data_tmp), family =PIG,
               sigma.fo=~CVTOT, n.cyc = 60, trace = FALSE)

  # Obtenemos predicciones
  mL.p <- predict(mL, what = "mu", newdata=data.frame("CVTOT"=x),
                  type = "response", terms = NULL, se.fit = FALSE, data = na.omit(data_tmp))
  mC.p <- predict(mC, what = "mu", newdata=data.frame("CVTOT"=x),
                  type = "response", terms = NULL, se.fit = FALSE, data = na.omit(data_tmp))
  mS.p <- predict(mS, what = "mu", newdata=data.frame("CVTOT"=x),
                  type = "response", terms = NULL, se.fit = FALSE, data = na.omit(data_tmp))
  # Ploteamos
  windows()
  op <- par(mfrow = c(1,1))
  plot(data_tmp$CVTOT,data_tmp$GB.Seqs,pch=16,ylim=c(0,200), xlab="CVTOT",
      ylab="Appearances",main=paste("Filtro: ",a))
  lines(x,mL.p,col="red",type="l",lwd=2)
  lines(x,mC.p,col="blue",type="l",lwd=2)
  lines(x,mS.p,col="orange",type="l",lwd=2)
  legend("topright", c("PIG.Lineal","PIG.Cuadrático","PIG.Splines"),
         lty=c(2,2,2),col=c("red","blue","orange"), ncol=1, lwd=rep(2,3), inset = .05, cex=0.75)
}

## Modelos no paramétricos:
#### Splines.
####    * Smoothing spline
####    * B-Splines
####    * P-Splines
dfs <- seq(2, 10, by = 2)
for (df in dfs)
{
  spl <- smooth.spline(x,y, df = 6)
  summary(spl)
}
#### Modelos aditivos.
#### Modelos semiparamétricos