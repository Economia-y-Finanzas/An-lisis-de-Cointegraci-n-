# Análisis de Cointegración
Análisis econométrico de la relación de largo plazo entre los rendimientos de bonos corporativos canadienses AAA y BAA mediante pruebas de cointegración y un modelo de corrección de errores (VECM) en R.

---

## Informe PDF

A continuación puedes **ver el informe directamente**:
[Análisis de cointegración para bonos de Canadá](Análisis_de_cointegración_para_bonos_de_Canada.pdf)


---
#### Librerías utilizadas:
El análisis fue realizado en **R** utilizando los siguientes paquetes:
```r
library(normtest)
library(vars)
library(forecast)
library(urca)
library(moments)
library(vioplot)
library(tsDyn)
library(tseries)
```

#### Definición de la base de datos utilizada en el análisis:
```r
bonds <- data[, c("aaa", "bbb")]
```

#### Gráfico de las series:
```r
par(mfrow = c(2,1), mar = c(2.5,5,1,1), las = 0, cex.axis = 1.05, cex.lab = 1.05)
ts.plot(bonds[,'aaa'],    ylab='bonos corporativos AAA')
ts.plot(bonds[,'bbb'], ylab='bonos corporativos BAA')
```
<img width="629" height="389" alt="image" src="https://github.com/user-attachments/assets/8781b475-c592-4ba9-baae-5ed32e8bea4d" />

#### ACF de las series:
```r
par(mfrow = c(1,2), mar = c(4.2,4.2,1,1), las = 0, cex.axis = 1.05, cex.lab = 1.05)
Acf(bonds[,'aaa'],    ylab='AAA', xlab='Lag', lag.max = 40)
Acf(bonds[,'bbb'], ylab='BAA', xlab='Lag', lag.max = 40)
```
<img width="629" height="389" alt="image" src="https://github.com/user-attachments/assets/04ac3f52-8c3d-406b-8bf9-627d1e64aa65" />


#### Gráficos de violin:
```r
par(mfrow=c(2,1), mar = c(2.5,5,1,1), las = 0, cex.axis = 1.05, cex.lab = 1.05)
for (m in 1:2) {
  vioplot(bonds[,m], ylab=colnames(bonds)[m], col=gray(0.7), horizontal = T)
}
```
<img width="629" height="389" alt="image" src="https://github.com/user-attachments/assets/074ddea6-e50b-4e91-a4fb-998348cdb8de" />

#### Determinación del rezago p para el tes de ADF
```r
p_max <- floor((dim(bonds)[1]-1)^(1/3))
p_max
# [1] 8
```

##### Pruebas de raíz unitaria (ADF)

Se aplicó el test de Dickey-Fuller aumentado para verificar la estacionariedad de las series de rendimientos de bonos.

```r
df_aaa_none <- ur.df(bonds[,"aaa"], type = "none", lags = p_max, selectlags = "BIC")
df_bbb_none <- ur.df(bonds[,"bbb"], type = "none", lags = p_max, selectlags = "BIC")
```

```
ADF statistic (AAA) = -0.168
ADF statistic (BAA) = -0.169
Critical value (5%) = -1.95
```

El test de Dickey-Fuller aumentado contrasta las siguientes hipótesis:

$$
H_0: \rho = 1 \quad \text{(la serie tiene raíz unitaria, no es estacionaria)}
$$

$$
H_1: |\rho| < 1 \quad \text{(la serie es estacionaria)}
$$

No se rechaza la hipótesis nula, lo que sugiere que ambas series son no estacionarias en niveles.


#### Selección del número óptimo de rezagos para el VAR

Para determinar el número de rezagos del modelo VAR se utilizan criterios de información.

```r
VARselect(bonds, lag.max = p_max, type = "none")
```

```
AIC(n)  HQ(n)  SC(n)  FPE(n)
  8      3      3      8
```

##### Estimación del modelo VAR

Con base en los criterios de información se selecciona un VAR con \(p = 3\) rezagos.

```r
var <- VAR(bonds, p = 3, type = "none")
```

#### Prueba de cointegración de Johansen

Para determinar la existencia de relaciones de equilibrio de largo plazo entre las series se aplica el test de cointegración de Johansen utilizando los estadísticos Trace y Maximum Eigenvalue.

```r
summary(ca.jo(bonds, type = "trace", ecdet = "trend", K = 3, spec = "transitory"))
summary(ca.jo(bonds, type = "eigen", ecdet = "trend", K = 3, spec = "transitory"))
```
**Trace test**

```
r = 0   |  test = 30.41   |  critical value (5%) = 25.32
r ≤ 1   |  test = 4.14    |  critical value (5%) = 12.25
```

**eigen test**

```
r = 0   |  test = 26.27   |  critical value (5%) = 18.96
r ≤ 1   |  test = 4.14    |  critical value (5%) = 12.25
```

**Hipótesis del test**

Para cada posible rango de cointegración \(r\):

$$
H_0 : \text{El número de relaciones de cointegración es } \le r
$$

$$
H_1 : \text{El número de relaciones de cointegración es } > r
$$


Los resultados de ambos estadísticos indican que se **rechaza la hipótesis nula de ausencia de cointegración (r = 0)**, pero **no se rechaza la hipótesis de a lo sumo una relación cointegrante (r ≤ 1)**.

#### Estimación del modelo VECM

Dado que el test de Johansen indica la existencia de **una relación de cointegración (r = 1)**, se estima un modelo de corrección del error (VECM).

```r
vecm <- ca.jo(bonds, type = "trace", ecdet = "trend", K = 3, spec = "transitory")
vecm.r1 <- cajorls(vecm, r = 1)
summary(vecm.r1$rlm)
```

```
Ecuación ΔAAA : 0.0242   (p = 0.509)
Ecuación ΔBBB : 0.0915   (p = 0.004)
```

#### Estimación del VECM por máxima verosimilitud

El modelo VECM también se estima mediante máxima verosimilitud para obtener los coeficientes de corto plazo, el vector de cointegración y los parámetros de ajuste.

```r
vecm <- VECM(bonds, lag = 2, r = 1, include = "const", LRinclude = "trend", estim = "ML")
summary(vecm)

coefB(vecm)   # vector de cointegración (β)
coefA(vecm)   # parámetros de ajuste (α)
```

**Vector de cointegración (β)**

```
aaa    1.0000
bbb   -0.8794
trend -0.00059
```

La relación de equilibrio de largo plazo puede expresarse como:

$$
AAA_t - 0.879\,BAA_t - 0.0006,t = 0
$$

**Parámetros de ajuste (α)**

```
ΔAAA : 0.0242
ΔBBB : 0.0915
```
#### Relación de cointegración (ECT) del VECM y grafica en el tiempo.
```r
beta <- vecm.r1$beta
X <- cbind(bonds[, c('aaa', 'bbb')], t=seq(1,dim(bonds)[1]))
coint <- X %*% beta
rel_coint <- ts(as.numeric(coint), start = start(bonds), frequency = frequency(bonds))
plot(rel_coint, main="Relacion de cointegracion (ECT) con tendencia", ylab="ECT", xlab="Tiempo")
```
<img width="629" height="389" alt="image" src="https://github.com/user-attachments/assets/9112ad05-9405-442a-bfd7-01407e795327" />







