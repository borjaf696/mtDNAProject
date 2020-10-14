# Title     : TODO
# Objective : TODO
# Created by: borja
# Created on: 08/10/2020
### Paquetes #######
require(ggplot2)
require(pscl)
require(boot)
require(VGAM)
require(gamlss)
require(gamlss.dist)
####################
data <- read.csv('output/lsu_df.csv', sep = ',', dec = '.', header = T)
head(data)
data$CVTOT <- as.numeric(as.character(data$CVTOT))
data$GB.Seqs <- na.omit(as.numeric(as.character(data$GB.Seqs)))
plot(data[,'CVTOT'],data[,'GB.Seqs'],pch=16,col="blue",ylim=c(0,38),xlab="CvIndex",ylab="Appearances",main="")
data <- data[data$CVTOT != -1,]
data$GB.Seqs[data$CVTOT == -1] <- 0
x_ind <- is.na(data$CVTOT)
data <- data[!x_ind,]
x <- data$CVTOT
y <- data$GB.Seqs
sum(is.na(x))
sum(is.na(y))
# Comprobamos datos
summary(data)
head(data)
dim(data)
hist(data$GB.Seqs, breaks = 'fd', freq = F)

# Visualización
par(mfrow = c(1,1))
plot(data[,'CVTOT'],data[,'GB.Seqs'],pch=16,col="blue",ylim=c(-2,38),xlab="CvIndex",ylab="Appearances",main="")
sum(y == 0)
# Efectivamente hay muchos 0 lo que da cierto motivo a las técnicas de zeros inflados
sum(y != 0)
# Poisson + ZerosInflated
## Caso 1: Datos completos
## Caso 2: Datos filtrados al 1% superior
## Caso 3: Datos filtrados al 5$ superior
l_a <- c(0,0.01,0.05)
for (a in l_a)
{
  # Eliminamos datos fuera de lugar
  lower_bound <- quantile(y, 0.00);lower_bound
  upper_bound <- quantile(y, 1 - a);upper_bound

  outlier_ind <- which(y < lower_bound | y > upper_bound)
  if (!is.na(outlier_ind[1]))
    data_tmp <- data[-outlier_ind,]
  else
    data_tmp <- data

  # Comprobamos data_outlier

  # Comprobamos datos
  summary(data_tmp)
  head(data_tmp)
  dim(data_tmp)
  # Analisis visual
  windows()
  hist(data_tmp$GB.Seqs, freq = F, breaks = 'fd')

  windows()
  plot(data_tmp[,'CVTOT'],data_tmp[,'GB.Seqs'],pch=16,col="blue",ylim=c(-2,38),xlab="CvIndex",ylab="Appearances",main="")
  ## Regresión paramétrica - poisson (sistema de recuento, con valores > 0, makes sense) - sin tocar 0s
  RegZero_s_inf <- glm(GB.Seqs ~ CVTOT , data = data_tmp,family = poisson)
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
  summary(RegZero_inf_sP_c1)
  summary(RegZero_inf_sP_c2)
  summary(RegZero_inf_sP_c3)
  summary(RegZero_inf_sP_def)
  summary(RegZero_inf_sP_geo1)
  summary(RegZero_inf_sP_geo2)
  Pr_inf_sP_c1<-predict(RegZero_inf_sP_c1,
                        newdata=data.frame("CVTOT"=x),se.fit=TRUE)
  Pr_inf_sP_c2<-predict(RegZero_inf_sP_c2,
                        newdata=data.frame("CVTOT"=x),se.fit=TRUE)
  Pr_inf_sP_c3<-predict(RegZero_inf_sP_c3,
                        newdata=data.frame("CVTOT"=x),se.fit=TRUE)
  Pr_inf_sP_def<-predict(RegZero_inf_sP_def,
                         newdata=data.frame("CVTOT"=x),se.fit=TRUE)
  Pr_inf_sP_geo1<-predict(RegZero_inf_sP_geo1,
                         newdata=data.frame("CVTOT"=x),se.fit=TRUE)
  Pr_inf_sP_geo2<-predict(RegZero_inf_sP_geo2,
                         newdata=data.frame("CVTOT"=x),se.fit=TRUE)

  lines(x,Pr_inf_sP_c1,col="blue",type="l",lwd=2,lty=2)
  lines(x,Pr_inf_sP_c2,col="blue",type="l",lwd=2,lty=2)
  lines(x,Pr_inf_sP_c3,col="red",type="l",lwd=2,lty=2)
  lines(x,Pr_inf_sP_def,col="pink",type="l",lwd=2,lty=2)
  # Curiosamente estas dos no van tan mal
  lines(x,Pr_inf_sP_geo1,col="black",type="l",lwd=2,lty=2)
  lines(x,Pr_inf_sP_geo2,col="orange",type="l",lwd=2,lty=2)
  # Comparación de los modelos
  pchisq(2 * (logLik(RegZero_inf_sP_c1) - logLik(RegZero_inf_sP_c2)), df = 1, lower.tail = FALSE)
  pchisq(2 * (logLik(RegZero_inf_sP_c1) - logLik(RegZero_inf_sP_c3)), df = 1, lower.tail = FALSE)
  pchisq(2 * (logLik(RegZero_inf_sP_c2) - logLik(RegZero_inf_sP_c3)), df = 1, lower.tail = FALSE)

  # Una de alternativa suele ser aplicar el logaritmo de los datos
  data_tmp$GB.Seqs.log <- as.numeric(log(data_tmp$GB.Seqs + 1))
  summary(data_tmp$GB.Seqs.log)
  windows()
  plot(data_tmp$CVTOT, data_tmp$GB.Seqs.log, col = "blue", xlab = "CvTot",ylab = "log(Appearances) + 1", main = "Log appearances")
  lines(x,log(Pr_inf_sP_c1+1),col="red",type="l",lwd=2)
  lines(x,log(Pr_inf_sP_c2+1),col="orange",type="l",lwd=2)
  lines(x,log(Pr_inf_sP_c3+1),col="black",type="l",lwd=2)
  lines(x,log(Pr_inf_sP_geo1+1),col="yellow",type="l",lwd=2)
  legend("topleft", c("ZIsRPo","ZIcRPo","ZIsRGeom","ZIcRGeom"), lty=c(1,2,1,2),col=c("blue","blue","orange","orange"), ncol=1, lwd=rep(2,4), inset = .05)

}
# No se muy bien si se captura así el comportamiento de manera apropiada, la tendencia queda más o menos clara
# a mayor CvIndex menos ocurrencias en GenBank. Por otro lado, la presencia de valores tan elevados de Apariciones
# provoca que juegen demasiado en las estimaciones, sin embargo al eliminarlas tenemos el problema de excesiva
# carga a la izquierda. El único comportamiento reseñable es que con la regresión de Poisson básica parece como
# que se invierte la tendencia de decrecimiento a mayor CV.

## Estudio de las distribuciones
op <- par(mfrow = c(2, 3))
dPO <- histDist(data$GB.Seqs, "PO", main = "PO", trace = FALSE, xlim=c(0,20))     # Poisson
dNBI <- histDist(data$GB.Seqs, "NBI", main = "NBI", trace = FALSE,xlim=c(0,20))   # Binomial negativa tipo 1
dPIG <- histDist(data$GB.Seqs, "PIG", main = "PIG", trace = FALSE,xlim=c(0,20))   # Poisson inverse gaussian
dSI <- histDist(data$GB.Seqs, "SICHEL", main = "SICHEL", trace = FALSE, xlim=c(0,20)) # Sichel?
dZIP <- histDist(data$GB.Seqs, "ZIP", main = "ZIP", trace = FALSE, xlim=c(0,20))  # Zero inflated poisson
dZIP2 <- histDist(data$GB.Seqs, "ZIP2", main = "ZIP2", trace = FALSE, xlim=c(0,20)) # Zero inf. pois. 2
par(op)
AIC(dPO, dNBI, dPIG, dSI, dZIP, dZIP2)

## En base al AIC nos quedaríamos con dPIG que es un modelo de poisson pero descartariamos claramente los de ceros
## inflados.

## A través del ajuste de gamlss ajustamos distintos modelos con diferentes familias de funciones
fam.gamlss <- c("PO", "NBI", "NBII", "PIG", "DEL", "SICHEL","ZIP","ZIP2")
m.l <- m.q <- m.s <-list()
## Ajsute lineal con enlaces
for (i in 1:8) {
  m.l[[fam.gamlss[i]]] <- gamlss(GB.Seqs~ CVTOT , data=na.omit(data), family = fam.gamlss[i], n.cyc = 60, trace = FALSE)$aic
}
## Regresión polinómico local grado 2 con los enlaces
for (i in 1:8) {
  m.q[[fam.gamlss[i]]] <-GAIC(  gamlss(GB.Seqs ~ poly(CVTOT,2) , data=na.omit(data), family = fam.gamlss[i],  n.cyc = 60, trace = FALSE))
}
## Smooth cubic splines (suavización via splines cúbicos)
for (i in 1:8) {
  m.s[[fam.gamlss[i]]] <-GAIC(  gamlss(GB.Seqs ~ cs(CVTOT) , data=na.omit(data), family = fam.gamlss[i],  n.cyc = 60, trace = FALSE))
}
unlist(m.l)
unlist(m.q)
unlist(m.s)

# Recorte de datos
ln_y <- log(y, base = exp(1))
ind_y <- (ln_y == -Inf)
# Nuevos datos
x_i <- x[!ind_y]
y_i <- ln_y[!ind_y]
# Removing outliers
out <- boxplot.stats(y_i)$out
out_ind <- which(y_i %in% c(out))
out_ind
boxplot(y_i,
        ylab = "hwy",
        main = "Num appearences"
)
#Fin preprocesado de datos
lower_bound <- quantile(y_i, 0.00);lower_bound
upper_bound <- quantile(y_i, 0.95);upper_bound

outlier_ind <- which(y_i< lower_bound | y_i > upper_bound)
x_i_n <- x_i[-outlier_ind]
y_i_n <- y_i[-outlier_ind]
par(mfrow = c(1,2))
upper <- c(0.95,0.99)

for (up in upper)
{
  lower_bound <- quantile(y_i, 0.00);lower_bound
  upper_bound <- quantile(y_i, up);upper_bound

  outlier_ind <- which(y_i< lower_bound | y_i > upper_bound)
  x_i_n <- x_i[-outlier_ind]
  y_i_n <- y_i[-outlier_ind]
  data_tmp <- data.frame(x = x_i_n, y = y_i_n)
  # Small visualization
  plot(data_tmp[,1],data_tmp[,2],
       main = 'Smoothing splines', ylab = 'ln(appearances)', xlab = 'CV')
  fit <- smooth.spline(data_tmp$x, data_tmp$y)
  lines(fit, col = 'red')
  summary(sp.cis)
  res <- (fit$yin - fit$y)/(1-fit$lev)      # jackknife residuals
  sigma <- sqrt(var(res))                     # estimate sd

  upper <- fit$y + 2.0*sigma*sqrt(fit$lev)   # upper 95% conf. band
  lower <- fit$y - 2.0*sigma*sqrt(fit$lev)   # lower 95% conf. band
  lines(fit$x,upper, lty = 2, col = 'blue')
  lines(fit$x, lower, lty = 2, col = 'blue')
}