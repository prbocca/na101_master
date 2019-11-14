#!/usr/bin/env python
# coding: utf-8

# # Práctico 8 - Predecir información faltante en redes
# 
# # 1. Ambientes de trabajo
# 
# ## 1.a) Ambiente COLAB remoto
# 
# 1.   Abrir en navegador: https://colab.research.google.com/
# 2.   Abrir el notebook de la tarea:
#      File-> Open Notebook -> Github -> https://github.com/prbocca/na101_master -> homework_8_prediction
# 3.   Guardar el notebook en su Google Drive:
#      File -> Save a Copy in Drive... 
# 4.   Renombrar el archivo `"cedula ID"_ar_hw8.ipynb`, por ejemplo *33484022_ar_hw8.ipynb*
# 5.   Al final usted deberá descargar el notebook. Asegurarse que se están guardando las salidas de ejecución en el notebook: File -> Download .ipynb
# 6.   Luego estos archivos deberán ser enviados a prbocca@fing.edu.uy 
# 
# ##1.b) Ambiente RSTUDIO local (opcional)
# 
# Abrir el .r de la tarea en: https://github.com/prbocca/na101_master/tree/master/homework_8_prediction

# ## 1.c) Cargar Librerias

# In[ ]:


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

libs = c("ROCR", #evaluacion de problemas de clasificacion
         "rgexf", #igraph to gephi
         "vioplot",  "cowplot",
         # "org.Sc.sgd.db", "GOstats", "GO.db", #parecen no estar disponibles en el ambiente COLAB y los instalaremos luego de otra manera
         "foreach", "doMC", #procesamiento en paralelo
         "sand","igraph", "ggplot2", "reshape2") #las basicas

load_libs(libs)

# ## 1.d) Descargar funciones auxiliares

# In[ ]:


#directorio donde se va a trabajar
data_path = "/content/ar/hw8/"

dir.create(data_path, showWarnings = FALSE, recursive = TRUE)
setwd(data_path)
getwd()
list.files()

# cargo funciones auxiliares
source("https://raw.githubusercontent.com/prbocca/na101_master/master/homeworks_common.r")

# # 2. Predicción de enlaces en `R`.
# 
# Seguir las Secciones 7.1 y 7.2 del libro [SANDR], ejecutando el código fuente incluido.
# 

# In[ ]:


# Chapter 7 Network Topology Inference
# 7.1 Introduction
# 7.2 Link Prediction SANDR (esta es el tutorial de SANDR incluido como ejercicio)
#########################################################
# To illustrate the potential of simple scoring methods like this, recall the network
# fblog of French political blogs. The number of nearest common neighbors for
# each pair of vertices in this network, excluding—if incident to each other—the two
# vertices themselves, may be computed in the following manner.
#library(sand)
fblog = upgrade_graph(fblog)  #evita warnings por ser un grafo en version vieja
nv <- vcount(fblog)
ncn <- numeric()
A <- get.adjacency(fblog)
for(i in (1:(nv-1))){
  ni <- neighborhood(fblog, 1, i)
  nj <- neighborhood(fblog, 1, (i+1):nv)
  nbhd.ij <- mapply(intersect, ni, nj, SIMPLIFY=FALSE)
  temp <- unlist(lapply(nbhd.ij, length)) -
    2 * A[i, (i+1):nv]
  ncn <- c(ncn, temp)
}

# In[ ]:


# In Fig. 7.2 we compare the scores s(i, j) for those vertex pairs that are not incident
# to each other (i.e., ‘no edge’), and those that are (i.e., ‘edge’), using so-called violin
# plots.
#library(vioplot)
Avec <- A[lower.tri(A)]
vioplot(ncn[Avec==0], ncn[Avec==1],
        names=c("No Edge", "Edge"))
title(ylab="Number of Common Neighbors")

# In[ ]:


# It is evident from this comparison that there is a decided tendency towards larger
# scores when there is in fact an edge present. Viewing the calculation we have done
# here as a ‘leave-one-out’ cross-validation, and calculating the area under the curve
# (AUC) of the corresponding ROC curve,
#library(ROCR)
pred <- prediction(ncn, Avec)
perf <- performance(pred, "auc")
slot(perf, "y.values")

#agrego el plot de la curva ROC
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc.perf, col="red")
abline(a=0, b= 1)

# $\def\ccalN{{\mathcal N}}$
# 
# # 3. Predicción de aristas en redes reales.
# 
# En este ejercicio veremos que tan bien funcionan los métodos informales de *score* para predecir enlaces en distintas redes. Este ejercicio fue originalmente propuesto por el Prof. Aaron Clauset, en su curso "Network Analysis and Modeling", CSCI 5352, del Santa Fe Institute.
# 
# Como vimos, los métodos informales de predicción de enlaces basados en *score* definen una medida de enlace para todo par de vértices que no son aristas (es decir, definen el *score*, $s(i,j)$, para todo $(i,j) \notin E^{obs}$), y predicen las aristas faltantes como aquellas parejas con mayor `score`.
# Hay muchas funciones de `score` posibles, en este ejercicio probaremos tres:
# *   camino más corto: $s(i,j)= \frac{1}{dist_{G^{obs}}(i,j)}$, siendo $dist_{G^{obs}}(i, j)$ la distancia geodésica observada entre los vértices $i$ y $j$ en el grafo observado;
# 
# *   vecinos en común: $s(i,j)=|\ccalN_i^{obs}\cap \ccalN_j^{obs}|$; y
# 
# *   producto de grados observados: $s(i,j)= d_i d_j$.

# ## 3.a) Descargar y cargar la red de la enfermedad de Malaria.
# 
# Existen muchas representaciones en red de las relaciones entre los genes para la enfermedad de la Malaria. 
# Usaremos la red de: *D. B. Larremore, A. Clauset, and C. O. Buckee, "A network approach to analyzing highly recombinant malaria parasite genes." PLOS Computational Biology 9(10), e1003268 (2013).* 
# 
# Esta red tiene $307$ vértices y $2684$ aristas.
# Los vértices son genes, y dos genes están conectados si ellos comparten un substring cuyo largo es estadísticamente significativo. Se tiene una etiqueta en los vértices que clasifica los genes.
# Visualizar la red en la siguiente Figura (el color de los vértices es la etiqueta de clasificación, y el tamaño representa el grado): ![alt text](https://github.com/prbocca/na101_master/raw/master/homework_8_prediction/g_malaria.png)
# 
# Descargar el grafo desde: https://github.com/prbocca/na101_master/raw/master/homework_8_prediction/g_malaria.graphml.

# In[ ]:


# descargar y cargar datos
g_malaria = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: -
#
#
#
##################################################################
summary(g_malaria)

#plot
set.seed(25)
E(g_malaria)$color = "gray"
plot(g_malaria, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.color=V(g_malaria)$label, vertex.size=3, vertex.label=NA,
     layout=layout_with_drl(g_malaria), main="Malaria")


# ## 3.b) Descargar y cargar la red de las relaciones entre directivos en empresas públicas de Noruega.
# 
# Usaremos los datos del siguiente trabajo: *C. Seierstad and T. Opsahl, "For the few not the many? The effects of affirmative action on presence, prominence, and social capital of women directors in Norway." Scand. J. Manag. 27(1), 44-54 (2011)*.
# 
# Es una red de afiliaciones entre directores que comparten alguna junta directiva en  las compañias públicas noruegas.
# La red tiene una gran componente, con $854$ vértices (directivos) y $2745$ aristas (juntas directivas en común), y corresponde a la realidad de agosto de 2011.
# Se incluyen algunos metadatos, como los nombres de directores y compañías, la ciudad y código postal para empresas, y el género de los directores. Visualizar la red en la siguiente Figura (el color de los vértices es el género (azul-hombre, rojo-mujer), y el tamaño representa el grado): ![alt text](https://github.com/prbocca/na101_master/raw/master/homework_8_prediction/g_gcdirectores.png)
# 
# Descargar el grafo desde: 
# https://github.com/prbocca/na101_master/raw/master/homework_8_prediction/g_gcdirectores.graphml.
# 

# In[ ]:


# descargar y cargar datos
g_gcdirectores = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: -
#
#
#
##################################################################
summary(g_gcdirectores)

#plot
set.seed(25)
E(g_gcdirectores)$color = "gray"
plot(g_gcdirectores, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.color=V(g_gcdirectores)$gender, vertex.size=3, vertex.label=NA,
     layout=layout_with_drl(g_gcdirectores), main="Red de afiliación de directores noruegos: componente gigante")


# ## 3.c) Predecir enlaces
# 
# Para medir la calidad de predicción de cada método usaremos la métrica $AUC$ (*Area Under Curve*). Ver la página de Wikipedia sobre ROC: http://bit.ly/2ehXHrb. Un predictor aleatorio tiene $AUC=0.5$ y un predictor perfecto tiene $AUC=1.0$ por tanto los métodos serán mejores a mayor valor de $AUC$.
# 
# Utilizando la métrica AUC, evaluar la capacidad predictiva del método de *score* de vecinos en común, sobre una subred observada con el $\%80$ de las aristas de la red de Malaria.
# 
# Se disponen de varias funciones auxiliares precargadas, incluyendo `score()` y  `delete_edges_rand()`. 
# Es posible ver las funciones axiliares en: https://raw.githubusercontent.com/prbocca/na101_master/master/homeworks_common.r.
# 
# El siguiente código resuelve el ejercicio y ejemplifica como usar las funciones auxiliares. Se obtiene una muy buena capacidad predictiva, con un $AUC=0.917$ y la curva ROC de la siguiente Figura: ![Curva ROC para el método de score de vecinos en común, sobre una subred con el %80 de las aristas de la red de Malaria.](https://github.com/prbocca/na101_master/raw/master/homework_8_prediction/roc_malaria.png)
# 
# 

# In[ ]:


# computo las aristas reales y las pongo en un vector (triangular inferior, mismo orden que score())
A_malaria <- get.adjacency(g_malaria)
A_malaria_v <- A_malaria[lower.tri(A_malaria)]

#creo una subred observada
set.seed(42)
g_malaria_obs = delete_edges_rand(g_malaria, p=0.2)

# computo las aristas observadas y las pongo en un vector (triangular inferior, mismo orden que score())
A_malaria_obs = get.adjacency(g_malaria_obs)
A_malaria_obs_v <- A_malaria_obs[lower.tri(A_malaria_obs)]
true_possibleedges_no_obs = A_malaria_v[A_malaria_obs_v==0] #vector real sobre aristas no observadas

# computo aristas predichas 
s_cn  = score(g_malaria_obs, type="common neighbors")
pred_possibleedges_no_obs = s_cn$score[A_malaria_obs_v==0] #vector prediccho sobre aristas no observadas

# computo la performance
pred_malaria <- prediction(pred_possibleedges_no_obs, true_possibleedges_no_obs)
perf_malaria <- performance(pred_malaria, "auc")
printf("El método de score saundo vecinos en cumún, tiene un AUC=%s", slot(perf_malaria, "y.values"))
# el valor de AUC es 0.9172444

#plot de la curva ROC
roc_malaria = performance(pred_malaria, measure = "tpr", x.measure = "fpr")
plot(roc_malaria, col="red")
title("ROC para la red de Malaria y %20 de aristas no observadas")
abline(a=0, b=1)

# ## 3.d) Estudio exahustivo de capacidad predictiva
# 
# Realizar un estudio más exahustivo sobre la calidad predictiva del método de score, probando las tres funciones de score definidas anteriormente, variando el porcentaje de aristas observadas entre \%99 y \%20.
# Buscando eliminar casuística, repetir la prueba sorteando distintos grafos observados.
# Realizar este procedimiento para las siguientes redes:
# *   grafo aleatorio de Erdos Renyi con $100$ vértices y una probabilidad de arista de $0.20$;
# 
# *   modelo Barabási-Albert con $100$ vértices y $m=9$ aristias agregadas en cada paso;
# 
# *   red de Malaria; y 
# 
# *   red de directores.
# 
# 
# Es posible utilizar la función auxiliar cargada previamente `summary_predictions(g, n_sample, p, types=c("1/dist","common neighbors", "degree product"), metric="auc", cores=1)`, que recibe un grafo $g$, crea un grafo observado para cada probabilidad de eliminación de arista del vector $p$, y prueba el métoodo de *score* para distintos tipos de medidas y para la métrica $AUC$. Repite cada caso `n_sample` y promedia resultados.
# 
# 

# In[ ]:


# defino los parametros del experimento

# para cada grafo, sorteo un subgrafo observado, calculo todos los tipos de score, predigo aristas faltantes y evaluo la metrica
p_v = c(0.01, 0.05, 0.10, 0.20, 0.50, 0.80) #las probabilidades de eliminacion de arista a evaluar
n_sample = 10 #cantidad de sorteos de g_obs
types = c("1/dist","common neighbors", "degree product")
metric="auc" 
n_cores = 1 #la ejecución en COLAB solo permite 1 CPU

# In[ ]:


#realizo los calculos para el grafo aleatorio

# sorteo un grafo
g_erdosrenyi = erdos.renyi.game(100,0.20)
write.graph(g_erdosrenyi, file="g_erdosrenyi.graphml", format= "graphml") 
plot(g_erdosrenyi, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.color="black", vertex.size=3, vertex.label=NA,
     layout=layout_with_drl(g_erdosrenyi), main="Erdos-Renyi")

# evaluo la capacidad predictiva sobre ese grafo
results_erdosrenyi = summary_predictions(g_erdosrenyi, n_sample, p_v, types, metric, cores=n_cores)
saveRDS(results_erdosrenyi, "results_erdosrenyi.RDS")

# imprimo resultados
ggplot(results_erdosrenyi$results_summary, aes(100*p, auc, colour = score)) + 
  geom_line() + 
  geom_point(data=results_erdosrenyi$results, aes(100*p, auc, colour = score), alpha=0.01) + 
  labs(title = "AUC para distintos scores y porcentaje de aristas no observadas", 
       subtitle = "grafo Erdos Renyi", x="% aristas no observadas", y="AUC")

# In[ ]:


#realizo los calculos para el grafo Barabasi-Albert

##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a anterior
#
#
#
##################################################################



# In[ ]:


#realizo los calculos para el grafo g_malaria

##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a anterior
#
#
#
##################################################################


# In[ ]:


#realizo los calculos para el grafo g_gcdirectores

##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a anterior
#
#
#
##################################################################

# ## 3.e) Interpretar resultados de sección anterior
# 
# En la sección anterior debieron obtenerse resultados similares a los de la siguiente Figura: ![AUC de los métodos de score para predecir enlaces, usando distintos porcentajes de aristas no observadas y distintas redes.](https://github.com/prbocca/na101_master/raw/master/homework_8_prediction/result_auc_p.png)
# 
# Como es de esperar: 
# *    no puede predecirse aristas en los grafos aleatorios;
# *    el *score* del producto de grado solo es útil en el modelo preferencial de Barabási-Albert debido a que en la construcción del grafo son más probables las aristas entre vértices más populares; y
# *    la capadidad predictiva en redes reales es bastante buena, inclusive para altos porcentajes de aristas no observadas.

# # 4. Predicción de atributos de vértices en `R`
# 
# Seguir las Secciones 8.1 y 8.2 del libro [SANDR], ejecutando el código fuente incluido.
# 

# In[ ]:


# Chapter 8 Modeling and Prediction for Processes on Network Graphs
# 8.1 Introduction
# 8.2 Nearest Neighbor Methods
#############################################################

# We illustrate through the problem of protein function prediction. 
# We saw examples of strong assortative mixing of protein function in
# the underlying network of protein–protein interactions. While the gold standard for
# establishing the functionality of proteins is through direct experimental validation
# (or ‘assays’), results like these have been taken to suggest that, given the vast num-
# ber of proteins yet to be annotated, it is natural to approach this problem from the
# perspective of statistical prediction. Approaches to protein function prediction that
# incorporate network-based information have become standard.
set.seed(42)

data(ppi.CC)
ppi.CC = upgrade_graph(ppi.CC)
summary(ppi.CC)
# IGRAPH 006c919 UN-- 134 241 -- 
#   + attr: name (v/c), ICSC (v/n), IPR000198 (v/n), IPR000403 (v/n), IPR001806 (v/n), IPR001849 (v/n),
# | IPR002041 (v/n), IPR003527 (v/n)
# contains a network data object, called ppi.CC, that consists of a network of 241
# interactions among 134 proteins, as well as various vertex attributes.
#
# The vertex attribute ICSC is a binary vector
V(ppi.CC)$ICSC[1:10]

# A visualization of this network is shown in Fig. 8.1.
V(ppi.CC)[ICSC == 1]$color <- "yellow"
V(ppi.CC)[ICSC == 0]$color <- "blue"
plot(ppi.CC, vertex.size=5, vertex.label=NA)


# In[ ]:


# A simple, but often quite effective, method for producing local predictions is the
# nearest-neighbor method. See Hastie, Tibshirani, and Friedman [71, Chap. 2.3.2],
# for example, for general background on nearest-neighbor methods. For networks,
# the nearest-neighbor method centers on the calculation, for a given vertex i ∈ V , of
# the nearest-neighbor average
# i.e., the average of the values of the vertex attribute vector X in the neighborhood
# N i of i. Here |N i | denotes the number of neighbors of i in G. Calculation of these
# averages over all vertices i ∈ V corresponds to a nearest-neighbor smoothing of X
# across the network.

# In the context of protein function prediction, X is a binary vector, with entries
# indicating whether or not each protein is or is not annotated with a function of
# interest (e.g., ICSC). In predicting binary vertex attributes, the nearest-neighbor av-
#   erages (8.1) typically are compared to some threshold. For example, a threshold
# of 0.5 is commonly used, with a nearest-neighbor average greater than this value
# meaning that a majority of neighbors have the characteristic indicated by X = 1,
# resulting in a prediction for X i of 1 as well. Such methods are also known as ‘guilt-
#   by-association’ methods in some fields.

# me quedo con la componente gigante
clu <- clusters(ppi.CC)
ppi.CC.gc <- induced.subgraph(ppi.CC, clu$membership==which.max(clu$csize))
plot(ppi.CC.gc, vertex.size=5, vertex.label=NA)

# we can calculate the nearest-neighbor average
# for each of the proteins in the giant connected component of our network.
nn.ave <- sapply(V(ppi.CC.gc), function(x) mean(V(ppi.CC.gc)[nei(x)]$ICSC))

# We then plot histograms of the resulting values, separated according to the status
# of the vertex defining each neighborhood, i.e., according to the status of the ‘ego’
# vertex, in the terminology of social networks.
# The results, shown in Fig. 8.2, confirm that ICSC can be predicted with fairly
# good accuracy. 
par(mfrow=c(2,1))
hist(nn.ave[V(ppi.CC.gc)$ICSC == 1], col="yellow",
     ylim=c(0, 30), xlab="Proportion Neighbors w/ ICSC",
     main="Egos w/ ICSC")
hist(nn.ave[V(ppi.CC.gc)$ICSC == 0], col="blue",
     ylim=c(0, 30), xlab="Proportion Neighbors w/ ICSC",
     main="Egos w/out ICSC")

# In particular, using a threshold of 0.5 would yield an error rate
# of roughly 25 %.
nn.pred <- as.numeric(nn.ave > 0.5)
mean(as.numeric(nn.pred != V(ppi.CC.gc)$ICSC))
#[1] 0.2598425

# In[ ]:


# La siguiente sección requiere previamente instalar algunos paquetes especiales

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("GOstats", "GO.db", "org.Sc.sgd.db"))

# In[ ]:


# Interestingly, we can push this illustration a bit further by taking advantage of the
# evolving nature of biological databases like GO. In particular, the proteins annotated
# in GO as not having a given biological function include both (1) those that indeed are
# known not to have that function, and (2) those whose status is simply unknown. As
# a result, by comparing against more recent versions of GO, it is sometimes possible
# to identify proteins whose status has changed for a given functional annotation,
# indicating that in the interim it has been discovered to in fact have that function.
# The R package GOstats, a part of the Bioconductor package, may be used to
# manipulate and analyze Gene Ontology, as contained in the database GO.db.
library(GOstats)
library(GO.db)
# And the annotations specific to the organism yeast can be obtained from the
# org.Sc.sgd.db database.
library(org.Sc.sgd.db)
# At the time of this writing, these annotations were last updated in September of
# 2013, roughly six years after the data in ppi.CC were assembled.

# In[ ]:


# We extract those proteins with the function ICSC—now subsumed under the
# term intercellular signaling transduction (ICST), or GO label 003556—and keep
# only those that have been identified from direct experimental assay, as indicated by
# the evidence code ‘IDA’.
x <- as.list(org.Sc.sgdGO2ALLORFS)
current.icst <- x[names(x) == "GO:0035556"]
ev.code <- names(current.icst[[1]])
icst.ida <- current.icst[[1]][ev.code == "IDA"]
# We then separate out the names of those proteins that had ICSC in our original data
orig.icsc <- V(ppi.CC.gc)[ICSC == 1]$name
# and similarly extract the names of those proteins under the new annotations that
# were present in the giant connected component of our original network.
candidates <- intersect(icst.ida, V(ppi.CC.gc)$name)
# Among these candidates, there are four that have been newly discovered to have
# ICSC, with the following names.
new.icsc <- setdiff(candidates, orig.icsc)
new.icsc
# [1] "YDL159W" "YHL007C" "YIL033C" "YLR362W"
# And among these four, we find that two of them would have been correctly predicted
# by comparing the value of their nearest-neighbor averages to a threshold of 0.5.
nn.ave[V(ppi.CC.gc)$name %in% new.icsc]

# # 5. Predicción de atributos de vértices en redes reales.
# 
# Por lo general, los vértices de las redes reales tiene atributos, los cuales pueden ser categóricos, escalares o reales.
# En muchos casos no se conocen estos atributos para algunos vértices, por ejemplo cuando solo se tiene una muestra de la red, cuando explicitamente el atributo es ocultado por el vértice, cuando el atributo aun no existe y se creará a futuro, etc.
# 
# Un método sencillo para predecir los atributos faltantes es el método de los vecinos más cercanos (*NNM - nearest-neighbor method*), donde se estima el atributo faltante en el vértice $i$ como el promedio de los atributos en el vecindario de $i$. La idea de vecindario y promedio puede variar dependiendo del caso particular. Este método suele funcionar bien cuando existe una marcada homofilia entre los vértices participantes respecto a ese atributo.

# ## 5.a) Cargar los datos
# 
# Cargar la red de la enfermedad de Malaria y la red de relaciones entre directivos en empresas públicas de Noruega siguiendo el procedimiento del ejercicio anterior.
# 

# In[ ]:


# Cargar la red de Malaria
g_malaria = read.graph(file="g_malaria.graphml", format= "graphml")
summary(g_malaria)

# Cargar la red de directores
g_gcdirectores = read.graph(file="g_gcdirectores.graphml", format= "graphml")
summary(g_gcdirectores)


# ## 5.b) Predecir atributos en la red de Malaria
# 
# Evaluar la precisión del método NNM, sobre una subred observada con el $\%80$ de los atributos de vértice `label` de la red de Malaria.
# El atributo `label`, disponibles en la red, es un enumerado (variable categórica) y es representado en `igraph` como números. Por tal motivo, el promedio sobre los vecinos del método NNM se define como el enumerado de mayor frecuencia. 
# Para medir la calidad de predicción del método usaremos la exactitud  (*accuracy* - fracción de predicciones correctas).
# 
# Se disponen de varias funciones auxiliares, incluyendo `nnm()` y `delete_vertex_attr_rand()`. 
# Es posible ver las funciones axiliares en: https://github.com/prbocca/na101_master/raw/master/homeworks_common.r.
# 
# El siguiente código resuelve el ejercicio y ejemplifica como usar las funciones auxiliares. Se obtiene una muy buena capacidad predictiva, con una exactitud de $acc=0.712$.
# 

# In[ ]:


#creo una subred observada
set.seed(42)
g_malaria_obs = delete_vertex_attr_rand(g_malaria,"label", p=0.2)
vertex_attr_deleted = is.na(vertex_attr(g_malaria_obs, "label")) #vector booleando con los atributos borrados

#computo vector con los atributos reales
vertex_attr_actual = vertex_attr(g_malaria, "label")[vertex_attr_deleted] #vector con los atributos reales
vertex_attr_actual = factor(vertex_attr_actual, levels=unique(vertex_attr(g_malaria, "label")))

#computo vector con atributos predichos
vertex_attr_pred  = nnm(g_malaria_obs,"label", fun="freq")[vertex_attr_deleted] #vector con los atributos predichos
vertex_attr_pred = factor(vertex_attr_pred, levels=unique(vertex_attr(g_malaria, "label")))

# computo la exactitud
cm = as.matrix(table(pred=vertex_attr_pred, actual=vertex_attr_actual)) # create the confusion matrix
accuracy = sum(diag(cm)) / sum(cm)
printf("La exactitud de NNM para la red e Malaria es %s", accuracy) # [1] 0.7121212

# otra forma corta de hacerlo
# comparo con funcion de prediccion que hace todo
#r = evaluate_predictions_vertexattr(g_malaria, g_malaria_obs, "label", fun="freq", metric="acc")
#r$metrics


# ## 5.c) Predecir atributos en la red de directores
# 
# En la red de relaciones entre directivos se conoce el género de los directivos, como un atributo binario de vértice llamado `gender`.
# Dado que el atributo es binario se puede utilizar tanto la métrica de exactitud como $AUC$ visto anteriormente. Además el método NNM puede realizarse usando el promedio entre los vecinos (y no la frecuencia más alta).
# Repetir la parte anterior para esta red. 
# ¿Qué se puede concluir sobre la homofilia en esta red?
# 

# In[ ]:


##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a anterior
#
#
#
##################################################################


# ## 5.d) Estudio exahustivo de capacidad predictiva
# 
# Realizar un estudio más exahustivo sobre la calidad predictiva del método NNM, variando el porcentaje de atributos `label` observados entre \%99 y \%20 para la red de Malaria.
# Buscando eliminar casuística, repetir la prueba sorteando distintos grafos observados.
# 
# Es posible utilizar la función auxiliar cargada previamente  `summary_predictions_vertexattr(g, attr, n_sample, p, fun="mean", metric="auc", cores=1)`, que recibe un grafo `g`, crea un grafo observado para cada probabilidad de eliminación de arista del vector `p`, y prueba el métoodo NNM promediando con la función `fun` y para la métrica `acc`. Repite cada caso `n_sample` y promedia resultados.
# 
# Deben obtenerse resultados similares a los de la siguiente Figura: ![precisión del método NNM para predecir atributos de vértices, usando distintos porcentajes de atributos no observados](https://github.com/prbocca/na101_master/raw/master/homework_8_prediction/result_acc_p.png). Observar que en este caso, la capadidad predictiva no es muy buena, y existe gran variación en los resultados entre muestras.
# 
# 

# In[ ]:


# defino los parametros del experimento

# para cada grafo, sorteo un subgrafo observado, calculo todos los tipos de score, predigo aristas faltantes y evaluo la metrica
p_v = c(0.01, 0.05, 0.10, 0.20, 0.50, 0.80) #las probabilidades de eliminacion de arista a evaluar
n_sample = 10 #cantidad de sorteos de g_obs
metric="acc" #metric="auc" #metric="f"
n_cores = 1 #la ejecución en COLAB solo permite 1 CPU

# In[ ]:


#realizo los calculos para el grafo g_malaria

results_malaria_vertexattr = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: ver funciones auxiliares
#
#
#
##################################################################
str(results_malaria_vertexattr)

#imprimo resultados
ggplot(results_malaria_vertexattr$results_summary, aes(100*p, acc)) + 
  geom_line() + 
  geom_point(data=results_malaria_vertexattr$results, aes(100*p, acc), alpha=0.20) + 
  labs(title = "Precisión para distintos porcentajes de atributos no observados", 
       subtitle = "grafo Malaria", x="% atributos de vértices no observados", y="precisión")

