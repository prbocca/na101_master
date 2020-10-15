# -*- coding: utf-8 -*-
"""ar_hw6.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/18rtZgaq0UI96qYguCIX8fQ6apMYKxQ8W

# Práctico 6 - Modelos de grafos

# 1). Ambientes de trabajo

## 1.a) Ambiente COLAB remoto

1.   Abrir en navegador: https://colab.research.google.com/
2.   Abrir el notebook de la tarea:
     File-> Open Notebook -> Github -> https://github.com/prbocca/na101_master -> homework_6_models
3.   Guardar el notebook en su Google Drive:
     File -> Save a Copy in Drive... 
4.   Renombrar el archivo `"cedula ID"_ar_hw6.ipynb`, por ejemplo *33484022_ar_hw6.ipynb*
5.   Al final usted deberá descargar el notebook. Asegurarse que se están guardando las salidas de ejecución en el notebook: File -> Download .ipynb
6.   Luego estos archivos deberán ser enviados a prbocca@fing.edu.uy 

##1.b) Ambiente RSTUDIO local (opcional)

Abrir el .r de la tarea en: https://github.com/prbocca/na101_master/tree/master/homework_6_models

## 1.c) Cargar Librerias

Todas las librerías deben instalarse correctamente, si el proceso se interrumpe o alguna librería da error en la instalación, entonces habrá problemas en el código más adelante. Si tiene este tipo de problemas pruebe hacer `Runtime -> Factory reset runtime`, y volver a intentar.
"""

# cargar librerias
load_libs <- function(libraries = libs, install=TRUE){
  if (install){ # instalar librerias no instaladas
    new.packages <- libs[!(libs %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
  }
  #cargo librerias  
  for (lib in libraries){
    require(lib, character.only=TRUE, quietly = FALSE)
  } 
} 

libs = c("sand","igraph") 

load_libs(libs)

"""## 1.d) Descargar funciones auxiliares"""

#directorio donde se va a trabajar
data_path = "/content/ar/hw6/"

dir.create(data_path, showWarnings = FALSE, recursive = TRUE)
setwd(data_path)
getwd()
list.files()

# cargo funciones auxiliares
source("https://raw.githubusercontent.com/prbocca/na101_master/master/homeworks_common.r")

"""# 2). Modelos clásicos de grafos en `R`

Seguir las Secciones 5.1 y 5.2 del libro [SANDR], ejecutando el código fuente incluido.
"""

#5.1 Introduction
#5.2 Classical Random Graph Models

#The function erdos.renyi.game in igraph can be used to simulate classical
#random graphs of either type. Figure 5.1 shows a realization of a classical random
#graph, based on the choice of N v = 100 vertices and a probability of p = 0.02 of an
#edge between any pair of vertices. Using a circular layout, we can see that the edges
#appear to be scattered between vertex pairs in a fairly uniform manner, as expected.
set.seed(42)
g.er <- erdos.renyi.game(100, 0.02)
plot(g.er, layout=layout.circle, vertex.label=NA)

# Note that random graphs generated in the manner described above need not be connected.
is.connected(g.er)
#2 [1] FALSE

#Although this particular realization is not connected, it does nevertheless have a
# giant component, containing 71 of the 100 vertices. All other components contain
#between one and four vertices only.
table(sapply(decompose.graph(g.er), vcount))
# 3 1 2 3 4 71
#4 15 2 2 1 1
#In general, a classical random graph G will with high probability have a giant com-
#ponent if p = c/N v for some c > 1.

# Under this same parameterization for p, for c > 0, the degree distribution will
# be well-approximated by a Poisson distribution, with mean c, for large N v .
# Indeed, in our simulated random graph, the mean degree is quite close to the
# expected value of (100 − 1) × 0.02 = 1.98.
mean(degree(g.er))

# Furthermore, we see that the degree distribution is quite homogeneous.
hist(degree(g.er), col="lightblue", xlab="Degree", ylab="Frequency", main="")

#Other properties of classical random graphs include that there are relatively few
#vertices on shortest paths between vertex pairs72
average.path.length(g.er)
#[1] 5.276511
diameter(g.er)
#[1] 14
#and that there is low clustering.
transitivity(g.er)
#[1] 0.01639344
# More specifically, under the conditions above, it can be shown that the diameter
# varies like O(log N v ), and the clustering coefficient, like N v −1.

"""# 3). Distribución de grado de Erdös-Renyi, $G(n, p)$ , en `R`. 

Utilizando la función `erdos.renyi.game()` generar grafos aleatorios de distinto orden.

## 3.a) Crecer la cantidad de vértices

Si fijamos $p = 0.02$, al crecer la cantidad de vértices: 
¿qué sucede con el grado promedio?, y 
¿la distribución de grado se acerca a una Normal o a una Poisson? 

Verificar comparando el histograma con la distribución teórica para $n \in \{10, 100, 1000, 10000\}$.
"""

# Gnp: efecto de crecer la red, se aproxima a una normal
par(mfrow=c(3,2))
p=0.02
for (i in 1:5){
  n = 10^i
  g.er <- erdos.renyi.game(n, p.or.m=p, type="gnp")
  #hist(degree(g.er), probability = T, col="lightblue", xlab="Degree", ylab="Probabilities", main="")
  g.er.dd = degree_distribution(g.er)
  plot(g.er.dd, xlab="Degree", ylab="Probabilities", main=paste("n=",n))
  x = seq(from=1,to=length(g.er.dd), by=1)
  lines(x, dnorm(x, mean = n*p, sd = sqrt(n*p*(1-p))), type="l",col="red")
  print(paste("n =",n,", grado promedio =", mean(degree(g.er))))
}

par(mfrow=c(1,1))

"""## 3.b) Crecer la cantidad de vértices y obtener una distribución de Poisson

¿Cómo debemos cambiar el parámetro $p$ del modelo al crecer $n$ para que la  distribución tienda a una Poisson? 

Verificar comparando el histograma con la distribución de Poisson para $n \in \{10, 100, 1000, 10000\}$ y un grado promedio de $1.8$.

Nota: Los parámtros de la distribución teórica fueron vistos en las clases teóricas.
"""

##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a 3.a), los parámtros de la distribución teórica fueron vistos en las clases teóricas.
#
#
#
##################################################################

"""# 4). Otros modelos de grafos en `R`.

Seguir las Secciones 5.3 y 5.4 del libro [SANDR], ejecutando el código fuente incluido. En ellas se explican los modelos de configuración, acoplamiento preferencial, Barabási-Albert, y mundo pequeño.
"""

# 5.3 Generalized Random Graph Models

set.seed(42)

# The igraph function degree.sequence.game can be used to uniformly
# sample random graphs with fixed degree sequence. Suppose, for example, that we
# are interested in graphs of N v = 8 vertices, half of which have degree d = 2, and the
# other half, degree d = 3. Two examples of such graphs, drawn uniformly from the
# collection of all such graphs, are shown in Fig. 5.2.
degs <- c(2,2,2,2,3,3,3,3)
g1 <- degree.sequence.game(degs, method="vl")
g2 <- degree.sequence.game(degs, method="vl")
plot(g1, vertex.label=NA)
plot(g2, vertex.label=NA)

#Note that these two graphs do indeed differ, in that they are not isomorphic. 
graph.isomorphic(g1, g2)
#[1] FALSE

# For a fixed number of vertices N v , the collection of random graphs with fixed
# degree sequence all have the same number of edges N e , due to the relation d  ̄ =
# 2N e /N v , where d  ̄ is the mean degree of the sequence (d (1) , . . . , d (N v ) ).
c(ecount(g1), ecount(g2))
#[1] 10 10


# On the other hand, it is important to keep in mind that all other characteristics
# are free to vary to the extent allowed by the chosen degree sequence. For example,
# we can generate a graph with the same degree sequence as our network of protein–
# protein interactions in yeast.
data(yeast)
degs <- degree(yeast)
fake.yeast <- degree.sequence.game(degs, method=c("vl"))
all(degree(yeast) == degree(fake.yeast))
# [1] TRUE

#But the original network has twice the diameter of the simulated version
diameter(yeast)
# [1] 15
diameter(fake.yeast)
# [1] 8

# and virtually all of the substantial amount of clustering originally present is now gone.
transitivity(yeast)
#[1] 0.4686178
transitivity(fake.yeast)
#[1] 0.0396504

# 5.4 Network Graph Models Based on Mechanisms

set.seed(42)

# 5.4.1 Small-World Models
# Work on modeling of this type received a good deal of its impetus from the seminal
# paper of Watts and Strogatz [146] and the ‘small-world’ network model introduced
# therein. These authors were intrigued by the fact that many networks in the real
# world display high levels of clustering, but small distances between most nodes.
#
# In order to create a network graph with both of these properties, Watts and Strogatz 
# suggested instead beginning with a graph with lattice structure, and then randomly 
# ‘rewiring’ a small percentage of the edges. More specifically, in this model
# we begin with a set of N v vertices, arranged in a periodic fashion, and join each
# vertex to r of its neighbors to each side. Then, for each edge, independently and
# with probability p, one end of that edge will be moved to be incident to another
# vertex, where that new vertex is chosen uniformly, but with attention to avoid the
# construction of loops and multi-edges.
# An example of a small-world network graph of this sort can be generated in
# igraph using the function watts.strogatz.game.

g.ws <- watts.strogatz.game(1, 25, 5, 0.05)
plot(g.ws, layout=layout.circle, vertex.label=NA)

# For the lattice alone, which we generate by setting p = 0, there is a substantial amount of clustering.
g.lat100 <- watts.strogatz.game(1, 100, 5, 0)
plot(g.lat100, layout=layout.circle, vertex.label=NA)
transitivity(g.lat100)
#[1] 0.6666667
#But the distance between vertices is non-trivial.
diameter(g.lat100)
#[1] 10
average.path.length(g.lat100)
#[1] 5.454545

#The effect of rewiring a relatively small number of edges in a random fashion is to
# noticeably reduce the distance between vertices, while still maintaining a similarly
# high level of clustering.
g.ws100 <- watts.strogatz.game(1, 100, 5, 0.05)
plot(g.ws100, layout=layout.circle, vertex.label=NA)
diameter(g.ws100)
#[1] 5
average.path.length(g.ws100)
#[1] 2.793939393
transitivity(g.ws100)
#[1] 0.5121305


# This effect may be achieved even with very small p. To illustrate, we simulate ac-
# cording to a particular Watts-Strogatz small-world network model, with N v = 1, 000
#and r = 10, and re-wiring probability p, as p varies over a broad range.
steps <- seq(-4, -0.5, 0.1)
len <- length(steps)
cl <- numeric(len)
apl <- numeric(len)
ntrials <- 100
for (i in (1:len)) {
  cltemp <- numeric(ntrials)
  apltemp <- numeric(ntrials)
  for (j in (1:ntrials)) {
    g <- watts.strogatz.game(1, 1000, 10, 1^steps[i])
    cltemp[j] <- transitivity(g)
    apltemp[j] <- average.path.length(g)
  }
  cl[i] <- mean(cltemp)
  apl[i] <- mean(apltemp)
}
# The results shown in Fig. 5.4, where approximate expected values for normalized
# versions of average path length and clustering coefficient are plotted, indicate that
# over a substantial range of p the network exhibits small average distance while main-
# taining a high level of clustering.
plot(steps, cl/max(cl), ylim=c(0, 1), lwd=3, type="l", col="blue", xlab=expression(log[10](p)),
  ylab="Clustering and Average Path Length")
lines(steps, apl/max(apl), lwd=3, col="red")

# 5.4.2 Preferential Attachment Models

# Many networks grow or otherwise evolve in time. The World Wide Web and
# scientific citation networks are two obvious examples. Similarly, many biological
# networks may be viewed as evolving as well, over appropriately defined time scales.
# Much energy has been invested in the development of models that mimic network growth.
# In this arena, typically a simple mechanism(s) is specified for how the network
# changes at any given point in time, based on concepts like vertex preference, fitness,
# copying, age, and the like. A celebrated example of such a mechanism is that of
# preferential attachment, designed to embody the principle that ‘the rich get richer.’
# ...

# Using the igraph function barabasi.game, we can simulate a BA random
# graph of, for example, N v = 100 vertices, with m = 1 new edges added for each new vertex.
set.seed(42)
g.ba <- barabasi.game(100, directed=FALSE)
# A visualization of this graph is shown in Fig. 5.5.
plot(g.ba, layout=layout.circle, vertex.label=NA)

# Note that the edges are spread among vertex pairs in a decidedly less uniform man-
# ner than in the classical random graph we saw in Fig. 5.1. And, in fact, there appear
# to be vertices of especially high degree—so-called ‘hub’ vertices.
# Examination of the degree distribution (also shown in Fig. 5.5)
hist(degree(g.ba), col="lightblue", xlab="Degree", ylab="Frequency", main="")
# confirms this suspicion, and indicates, moreover, that the overall distribution is quite
# heterogeneous. Actually, the vast majority of vertices have degree no more than two
# in this graph, while, on the other hand, one vertex has a degree of 11.
summary(degree(g.ba))
#Min. 1st Qu. Median Mean 3rd Qu. Max.
#1.00 1.00 1.00 1.98 2.00 11.00

# On the other hand, network graphs generated according to the BA model will
# share with their classical counterparts the tendency towards relatively few vertices
# on shortest paths between vertex pairs
average.path.length(g.ba)
# [1] 5.81555555555556
diameter(g.ba)
# [1] 12
#and low clustering.
transitivity(g.ba)
# [1] 0

"""# 5). Distribución de grado de Barabási-Albert en `R`.

Verificar que el modelo de Barabasi-Albert tiende a una distribución power-law de parámetro $\alpha = 3$ al crecer $n$. 
Para esto graficar la función de distribución acumulada complementaria (CCDF) en conjunto con la recta teórica de la power-law. Usar los $n \in \{10, 100, 1000, 10000, 100000, 1000000\}$. 
Para cada caso, encontrar el $\hat{\alpha}$ usando la función de ajuste de igraph, llamada `power.law.fit()`. 
Además entender el test de hipótesis que se realiza en el ajuste, e interpretar el p-value del resultado.
"""

# BA: efecto de crecer la red
set.seed(42)
par(mfrow=c(3,2))

for (i in 1:6){
  n = 10^i
  g.ba = barabasi.game(n)
  
  # Degree distribution is the cumulative frequency of nodes with a given degree
  # this, like degree() can be specified as "in", "out", or "all"
  g.ba.dd = degree.distribution(g.ba, cumulative=T,mode="all")
  g.ba.d = degree(g.ba,v=V(g.ba),mode="all")
  
  # Using the power.law.fit() function I can fit a power law to the degree distribution
  # Plot degree distribution histogram, and theoretical line
  g.ba.power = NA
  ##################################################################
  #                       TU CÓDIGO ACÁ                           
  #
  #
  #
  ##################################################################
  print(paste("n =",n,", alpha =", g.ba.power$alpha, ", xmin =", g.ba.power$xmin))
}

"""# 6) Distribución de grado de una red real.

Mark Newman en 2006 creó un grafo de routers en Internet, con 22963 vértices y 48436 aristas. 
Descargar el grafo en formato GML desde:
https://github.com/prbocca/na101_master/raw/master/homework_6_models/internet_routers-22july06.gml.zip

## 6.a) Descargar los datos

Descargar los datos en formato `GML` en el siguiente link: 
https://github.com/prbocca/na101_master/raw/master/homework_6_models/internet_routers-22july06.gml.zip.

Cargar los datos en `R`.
"""

# download data
download.file(url="https://github.com/prbocca/na101_master/raw/master/homework_6_models/internet_routers-22july06.gml.zip", destfile="internet_routers-22july06.gml.zip", mode="wb")
unzip(zipfile="internet_routers-22july06.gml.zip")
list.files()

# cargar datos
g.real = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: usar funciones de igraph
#
#
#
##################################################################
summary(g.real)

"""## 6.b) (Opcional) Visualizar el grafo en Gephi

Visualizar los vértices de mayor grado con mayor tamaño,
y los vértices con mayor PageRank más oscuros.

El resultado debe ser similar al de la siguiente Figura: ![alt text](https://github.com/prbocca/na101_master/raw/master/homework_6_models/routers_gephi.png)

## 6.c) En `R`, graficar la distribución de grado.

También graficar con ejes log-log para una mejor visualización.
"""

# graficar distribución de grado de g.real
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: usar funciones de igraph
#
#
#
##################################################################

"""## 6.d) Graficar la función de distribución acumulada complementaria (CCDF)."""

# graficar CCDF de g.real
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: usar funciones de igraph
#
#
#
##################################################################

"""## 6.e) Austar una distribución Power-law

Ajustar la distribución a una power-law. 

¿Cuál es el $\hat{\alpha}$ y el $d_{min}$ del ajuste? 

¿Cómo interpretar el $p-value$ del test de hipotesis del resultado del ajuste? 

Verificar las conclusiones comparando la CCDF con la recta ajustada. La gráfica debería ser similar a la siguiente: ![alt text](https://github.com/prbocca/na101_master/raw/master/homework_6_models/routers_powerlaw.png)
"""

# graficar la CCDF de g.real, y la distribución power-low ajustada
g.real.power = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: usar funciones de igraph
#
#
#
##################################################################
print(paste("El alpha estimado es =",g.real.power$alpha))