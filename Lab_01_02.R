#╒═══════════════════════════════════════════════════════════════════════════════╕
#│ Sesión laboratorio 1                                                          │
#│ Econometría Avanzada I Mención en Economía                                    │
#│ Última Revisión: 24 de junio 2021                                             │
#│ Ojo, este archivo tiene carácteres latinos y griegos, por lo que está         │
#│      codificado en UTF-8.                                                     │
#╘═══════════════════════════════════════════════════════════════════════════════╛
#╒═══════════════════════════════════════════════════════════════════════════════╕
#│ En esta primera sesión construiremos un "modelo verdadero" del tipo:          │
#│ General   :   Y = Xβ + ε                                                      │
#│ Específico:   Y = β1 + β2 x + ε                                               │
#│                                                                               │
#│ Los datos se generan de manera aleatoria, por lo que usamos los comandos      │
#│ set.seed() y rnorm()                                                          │
#╘═══════════════════════════════════════════════════════════════════════════════╛

# Primero las variables "de usuario":
n       <- 20      # Número de Observaciones
beta    <- c(1,1)  # Beta verdadero (desconocido para el investigador)
sigma   <- 1       # Verdadera desviación estándar de ε
set.seed(19000000) # La "semilla" fija la creación de números aleatorios
# De aquí en adelante todo es automático
k       <- length(beta)  # Obtiene el número de variables a partir de beta
# Construimos X: la matriz variables independientes
x       <- matrix(rnorm(n*k), n,k) # Damos estructura de matriz n×k a un vector de 
                                   # números aleatorios de longitud x*k
# Por convención la primera columna es la constante
x[,1]   <- 1                       

# Ahora contruimos el "modelo verdadero" (método del oráculo)
epsilon <- rnorm(n, sd=sigma) # Este es el verdadero ε (desconocido)
y       <- x%*%beta + epsilon # Creamos la var.dep a partir del modelo verdadero

# En la práctica el modelo verdadero es desconocido (β y ε), por lo que se necesita
# un criterio para estimarlos.  El que usamos para estimar β es el de mínimos
# cuadrados ordinarios. A este estimador le llamamos b = (X'X)‾¹X'y

# Estimador de mínimos cuadrados
b       <- solve(t(x)%*%x)%*%t(x)%*%y
# Para verificar su resultado puede usar el comando de R lm():
model1  <- lm(y~0+x)
summary(model1)

# Ahora podemos calcular ŷ y e
# Estimados
y_hat   <- x%*%b
e       <- y - y_hat

# A continuación puede comparar los valores verdaderos y los estimados:
# Gráficos
plot(y, type="l", col="blue")       # plot()  crea un nuevo gráfico
lines(y_hat, type="l", col="red")   # lines() escribe sobre el gráfico anterior
plot(epsilon, type="l", col="blue")
lines(e, type="l", col="red")

# Ahora podemos verificar algunas propiedades de las matrices asociadas:
# (note que el computador no siempre encuentra ceros exactos)
# Para la matriz de proyección:
p  <- x%*%solve(t(x)%*%x)%*%t(x)    # Matriz de proyección
p%*%x - x                           # Cero computacional
p%*%p - p                           # Idempotencia

# Para la matriz de aniquilación:
M  <- diag(n) -  x%*%solve(t(x)%*%x)%*%t(x)
M%*%x                               # Aniquilación
M%*%M - M                           # Idempotencia
M%*%epsilon - e                     # Esta propiedad es muy importante
                                    # Conecta lo desconocido (ε) con lo
                                    #  conocido (M y e)
#╒═══════════════════════════════════════════════════════════════════════════════╕
#│ Sesión laboratorio 2                                                          │
#│ Econometría Avanzada I MAEC/DOEC                                              │
#│ Última Revisión: 1 de septiembre                                              │
#╘═══════════════════════════════════════════════════════════════════════════════╛
#╒═══════════════════════════════════════════════════════════════════════════════╕
#│ En esta sesión realizaremos inferencia estadística del modelo estimado en     │
#│    la sesión anterior                                                         │
#│ General   :   Y = Xβ + ε                                                      │
#│ Específico:   Y = β1 + β2 x + ε                                               │
#╘═══════════════════════════════════════════════════════════════════════════════╛

# Primero calcularemos nuestro estimador de la varianza:
SSE    <- sum(e^2)                 # Suma de errores al cuadrado
df     <- n-k                      # Grados de libertad ("degrees of freedom")
sigma2 <- SSE/df                   # Estimador para σ²
varcov <- sigma2*solve(t(x)%*%x)   # Matriz de Covarianza de b: varcov = σ²(X'X)‾¹
ee_b   <- sqrt(diag(varcov))       # Vector de errores estándar de b 

# Test de significancia  Hₒ: β = 0
#   En este caso realizaremos el test para el vector completo
tc     <- b/ee_b                   # El estadístico de prueba t calculado
pval   <- 2*(1-pt(abs(tc), df=df)) # Valor p

R2     <- 1-SSE/sum(y^2)           # Coeficiente de Determinación
R2adj  <- 1- ((n-1)/(n-k))*(1-R2)  # Coeficiente de Determinación adjustado
# Note que hay una diferencia entre el R2 ajustado calculado acá y el 
# que obtuvimos usando lm().  Esto se debe a que en la estimación lm()
# impusimos un modelo sin intercepto, lo que no permite ajustar el R², por
# lo que lm() usa una versión alternativa. Para obtener el R² "correcto" hay
# que corregir la estimación en lm()
R2adj_sc  <- 1- ((n)/(n-k))*(1-R2) # Coeficiente de Determinación adjustado
                                   # sin constante ¿puede ver la diferencia?
                                   # Tarea: piense por qué la diferencia
# Nuevamente miremos el resultado de la regresión usando lm()
summary(model1)
# Observe los siguietes resultados obtenidos en este laboratorio
quantile(e)
cbind(b, ee_b, tc, pval)

# Usted puede construir en R sus propios "outputs" o resultados. Ejemplos:
cat(paste("Error estándar del residuo:", round(sqrt(sigma2),4), 
          "con", df ,"grados de libertad \n"))
# o por ejemplo, si está escribiendo un informe
cat(paste("La constante es", round(b[1],2), "con un error estándar de",
          round(ee_b[1],2)), "\n")
cat(paste("La pendiente es", round(b[2],2), "con un valor p =",
          round(pval[2],3)), "\n")
# El comando cat() escribe textos en la consola.
# El comando paste() crea textos pegando partes
#    usamos \n para crear una nueva línea
# El comando round() es para redondear:
#            round(b[1],2) redondea el primer elemento de b a 2 decimales
# Fin Laboratorio 2
cat("Fin del Laboratorio 2 \n")

