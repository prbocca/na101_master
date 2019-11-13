#!/usr/bin/env python
# coding: utf-8

# # Práctico 4 - Medidas en redes sociales reales
# 
# #1. Ambientes de trabajo
# 
# ## 1.a) Ambiente COLAB remoto
# 
# 1.   Abrir en navegador: https://colab.research.google.com/
# 2.   Abrir el notebook de la tarea:
#      File-> Open Notebook -> Github -> https://github.com/prbocca/na101_master -> homework_4_measures
# 3.   Guardar el notebook en su Google Drive:
#      File -> Save a Copy in Drive... 
# 4.   Renombrar el archivo `"cedula ID"_AR_hw1_IR.ipynb`, por ejemplo *33484022_AR_hw1_IR.ipynb*
# 5.   Al final usted deberá descargar el notebook. Asegurarse que se están guardando las salidas de ejecución en el notebook: File -> Download .ipynb
# 6.   Luego estos archivos deberán ser enviados a prbocca@fing.edu.uy 
# 
# ##1.b) Ambiente RSTUDIO local (opcional)
# 
# Abrir el .r de la tarea en: https://github.com/prbocca/na101_master/tree/master/homework_4_measures

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

libs = c("network", "sna","rgexf", 
         "networkD3","twitteR","stringr","dplyr", "ggplot2",
         "R.matlab", "sand","igraph","igraphdata") #procesamiento en paralelo

load_libs(libs)

# ## 1.d) Descargar funciones auxiliares

# In[ ]:


#directorio donde se va a trabajar
data_path = "/content/ar/hw4/"

dir.create(data_path, showWarnings = FALSE, recursive = TRUE)
setwd(data_path)
getwd()
list.files()

# cargo funciones auxiliares
source("https://raw.githubusercontent.com/prbocca/na101_master/master/homeworks_common.r")

# # 2. Características de vértices y aristas
# Características de vértices y aristas de un grafo en `R`. 
# Seguir las secciones 4.1 y 4.2 del libro [SANDR], ejecutando el código fuente incluido. 
# 

# In[ ]:


# 4.1 Introduction
# 4.2 Vertex and Edge Characteristics
# 4.2.1 Vertex Degree
data(karate)
summary(karate)

# The degree distribution for the karate club network
# We can see that there are three distinct groups of vertices, as measured by degree. 
# The two most highly connected vertices correspond to actors 1 and 34 in the network, 
# representing the instructor and administrator a
hist(degree(karate), col="lightblue", xlim=c(0, 50), xlab="Vertex Degree", ylab="Frequency", main="")

# weighted degree distribution
# For weighted networks, a useful generalization of degree is the notion of vertex
# strength, which is obtained simply by summing up the weights of edges incident to a given vertex.
# The previously observed distinction among the three groups of vertices is lost
hist(graph.strength(karate), col="pink", xlab="Vertex Strength", ylab="Frequency", main="")

# In[ ]:


# Degree distributions can exhibit a variety of shapes. For a network of interactions
# among protein pairs in yeast, for next example, the shape is quite different from that of the karate club. In this case, 
# the distribution of degrees associated with the edges among vertices is quite heterogeneous, as can be seen from examination of the histogram.
data(yeast)
ecount(yeast) #11855
vcount(yeast) #2617
d.yeast <- degree(yeast)
# degree distribution
hist(d.yeast,col="blue", xlab="Degree", ylab="Frequency", main="Degree Distribution")

# Given the nature of the decay in this distribution, a log–log scale is more effective in summarizing the degree information.
dd.yeast <- degree.distribution(yeast)
d <- 1:max(d.yeast)-1
ind <- (dd.yeast != 0)
plot(d[ind], dd.yeast[ind], log="xy", col="blue",
     xlab=c("Log-Degree"), ylab=c("Log-Intensity"),
     main="Log-Log Degree Distribution")
# we see that there is a fairly linear decay in the log-frequency as a function of log-degree.

# Beyond the degree distribution itself, it can be interesting to understand the manner 
# in which vertices of different degrees are linked with each other. Useful in assessing 
# this characteristic is the notion of the average degree of the neighbors of a
# given vertex. For example, a plot of average neighbor degree versus vertex degree
a.nn.deg.yeast <- graph.knn(yeast,V(yeast))$knn
plot(d.yeast, a.nn.deg.yeast, log="xy", col="goldenrod", xlab=c("Log Vertex Degree"),
     ylab=c("Log Average Neighbor Degree"))
# suggests that while there is a tendency for vertices of higher degrees to link with similar vertices, 
# vertices of lower degree tend to link with vertices of both lower and higher degrees.

# In[ ]:


# 4.2.2 Vertex Centrality
# An intuitively appealing way of displaying vertex centralities (for networks of small to moderate size) 
# is to use a radial layout, with more central vertices located closer to the center. 
# The function gplot.target, in the package sna, can be used for this purpose.
A <- get.adjacency(karate, sparse=FALSE)
g <- network::as.network.matrix(A)
par(mfrow=c(2, 2))
sna::gplot.target(g, sna::degree(g), main="degree",
                  #circ.lab = FALSE,
                  circ.col="skyblue",
                  usearrows = FALSE,
                  vertex.col=c("blue", rep("red", 32), "yellow"),
                  edge.col="darkgray")
sna::gplot.target(g, sna::closeness(g), main="closeness",
                  #circ.lab = FALSE, 
                  circ.col="skyblue",
                  usearrows = FALSE,
                  vertex.col=c("blue", rep("red", 32), "yellow"),
                  edge.col="darkgray")
sna::gplot.target(g, sna::betweenness(g), main="betweenness",
                  #circ.lab = FALSE, 
                  circ.col="skyblue",
                  usearrows = FALSE,
                  vertex.col=c("blue", rep("red", 32), "yellow"),
                  edge.col="darkgray")
sna::gplot.target(g, sna::evcent(g), main="evcent",
                  #circ.lab = FALSE, 
                  circ.col="skyblue",
                  usearrows = FALSE,
                  vertex.col=c("blue", rep("red", 32), "yellow"),
                  edge.col="darkgray")

# Extensions of these centrality measures from undirected to directed graphs are
# straightforward. However, in the latter case, there are in addition other useful
# options. For example, the HITS algorithm, based on the concept of ‘hubs and authorities’
# Applying these measures to the AIDS blog network indicates, as seen in Fig. 4.5,
# that only six of the 146 blogs in this network play the role of a hub, while the vast
# majority of the vertices (including some of the hubs) play the role of an authority.
l <- layout.kamada.kawai(aidsblog)
plot(aidsblog, layout=l, main="Hubs", vertex.label="", 
     vertex.size=10 * sqrt(hub.score(aidsblog)$vector))
plot(aidsblog, layout=l, main="Authorities", vertex.label="", 
     vertex.size=10*sqrt(authority.score(aidsblog)$vector))


# In[ ]:


# 4.2.3 Characterizing Edges
# Using edge betweenness with the karate network and examining, for instance,
# the edges with the three largest betweenness values
eb <- edge.betweenness(karate)
E(karate)[order(eb, decreasing=T)[1:3]]
# + 3/78 edges from 4b458a1 (vertex names):
#   [1] Actor 20--John A   Mr Hi   --Actor 20 Mr Hi   --Actor 32
# we are led to note that actor 20 plays a key role from this perspective in facilitating
# the direct flow of information between the head instructor (Mr Hi, vertex 1) and the
# administrator (John A, vertex 34).

# However, many other vertex centrality measures do not extend as easily. One
# way around this problem is to apply vertex centrality measures to the vertices in the
# line graph of a network graph G. The line graph of G, say G' = (V', E'), is obtained
# essentially by changing vertices of G to edges, and edges, to vertices, which in
# igraph may be accomplished using line.graph.
summary(karate)
karate_dual = line.graph(karate)
summary(karate_dual)
A_dual <- get.adjacency(karate_dual, sparse=FALSE)
g_karate_dual <- network::as.network.matrix(A_dual)
sna::closeness(g_karate_dual)
sna::betweenness(g_karate_dual)
sna::evcent(g_karate_dual)

# # 3. Medidas en redes sociales reales.
# 
# Calcular e interpretar medidas de centralidad de vértice en redes sociales reales (Linkedin/Facebook/Twiiter, etc.).
# 

# ## 3.a) Datos de Twitter de `isDotR`
# 
# Utilizando los datos de Twitter publicados por Enmanuel Santana (http://apuntes-r.blogspot.com.uy/2015/07/un-ejemplo-de-pagerank.html),
# cargarlos como un grafo en `R`. Los datos son un pequeño conjunto de la red de ego de la cuenta `isDotR` de Twitter.
# Los datos también pueden descargarse en: 
# https://raw.githubusercontent.com/prbocca/na101_master/master/homework_4_measures/Twitter_network_R.csv
# 
# 
# Utilizando el paquete `igraph` calcular las métricas de centralidad: `degree(), closeness(), betweenness(), evcent(), page.rank()`, 
# y comparar los tres usuarios con mayor centralidad para cada una de ellas.
# El resultado debería ser el de la siguiente Tabla:
# 
# |metric|1|2|3|
# |---|---|---|---|
# | `degree()` |revodavid | cjbayesian | freakonometrics |
# | `closeness()` |statbandit | cjbayesian | freakonometrics | 
# | `betweenness()` |revodavid | hadleywickham | cjbayesian |
# | `evcent()` | freakonometrics | revodavid | statsinthewild |
# | `page.rank()` | isdotr | hadleywickham | revodavid |
# 

# In[ ]:


# download data
download.file(url="https://raw.githubusercontent.com/prbocca/na101_master/master/homework_4_measures/Twitter_network_R.csv", destfile="Twitter_network_R.csv", mode="wb")
list.files()

# load data
datos <- read.csv("Twitter_network_R.csv")

# aprendemos de los datos
head(datos)
unique(c(as.character(datos[,2]),as.character(datos[,1]))) # los usuarios unicos son 34

# In[ ]:


# crea objeto de grafo 
grafo = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: crear grafo desde el data.frame "datos"
#
#
#
##################################################################
summary(grafo)

# grafica de grafo
plot(grafo,layout=layout.fruchterman.reingold, vertex.size=8,
     vertex.label.dist=0.4, vertex.color="red", edge.arrow.size=0.5)


# In[ ]:


# calcular centralidades
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: calcular centralidades sobre "grafo"
#
#
#
##################################################################

# grafica interacciones usando valor page.rank para tamaños en grafico
pr = page.rank(grafo)$vector
V(grafo)$label.cex =  0.6 + pr*3 # tamaño de letra segun page.rank
plot(grafo,layout=layout.fruchterman.reingold, vertex.size=pr*100,
     vertex.label.dist=0.4, vertex.color="red", edge.arrow.size=0.3)


# ## 3.b) (opcional) Centralidad en tu red social
# 
# Lamentablemente desde 2015, las principales redes sociales han cerrado sus APIs para acceder a los datos de redes (amigos, etc). Solo algunas tienen un acceso limitado (gratuito o con suscripción).
# El sitio SociLab (http://socilab.com/) realizaba un análisis básico de la red Linkedin del usuario. De forma excepcional (y por razones históricas) este sitio tuvo acceso a esta API hasta el 2018.
# Actualmente estan todas cerradas, solo existen muchas ofertas de servicios online, que utilizan tu cuenta de usuario para extraer a información (muy parcial) de las redes sociales. Ejemplos son: https://mentionmapp.com/, https://socioviz.net, https://netlytic.org/, etc.
# 
# En este ejercicio utilizaremos Netlytic con datos de twitter. La interfaz no es intuitiva, pero es potente. Debe: 
# i) crearse una cuenta, 
# ii) realizar un nuevo dataset vinculando tu cuenta de Twitter, escribiendo un nombre al dataset y las palabras de búsqueda (ej. `#Uruguay`); el resultado lleva unos minutos, y se puede acceder en la sección `mi dataset`, 
# iii) una de las opciones de análisis es basado en redes, en donde puede visualizar la red y exportarla.
# 
# Se pide realizar los pasos anteriores para generar un dataset sobre algún tema de interés, y exportarlo para visualizarlo con `igraph` o `Gephi`.
# 

# ## 4. Twitter
# 
# Twitter permite acceder parcialmente a datos de la red utilizando una cuenta de desarrollador gratuita. El 30/08/2018 a las 11.30am se descargaron los $5000$ tweets más recientes sobre `#Uruguay`. Los datos pueden descargarse en: 
# https://raw.githubusercontent.com/prbocca/na101_master/master/homework_4_measures/tweets_df.csv

# ## 4.a) Cargar los datos

# In[ ]:


# download and liad data
tweets_df = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a otras descargas
#
#
#
##################################################################
head(tweets_df)

# 
# ## 4.b) Tweets más populares
# 
# Ver los 5 tweets más populares (con más retweets) que no sean spam. Nota: el spam puede ser eliminado viendo los mensajes con más retweets (por ejemplo con más de 15000 retweets).
# 

# In[ ]:


# select top retweeted tweets
tweets_df = tweets_df[order(tweets_df$retweetCount, decreasing = TRUE),]
tweets_df$spam = (tweets_df$retweetCount >= 15000)  #supongo spam lo que tiene mas de 15000 retweets
tweets_df_filter = tweets_df[!tweets_df$spam,] #filtro el spam

# Distribucion de retweets
# mostrar el histograma de la cantidad de tweets  que tienen la misma cantidad de retweets. 
# También puede mostrarse como tabla.
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: usar tweets_df_filter
#
#
#
##################################################################


# Mostrar los 5 tweets mas retweeteados
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: usar tweets_df_filter
#
#
#
##################################################################


# ## 4.c)
# 
# Utilizando el paquete `igraph` crear la red de quién hace retweet de quién. Utilizando los tweets con al menos dos retweets (y que no son spam). Debe obtener una red similar a la siguiente Figura.
# ![](https://github.com/prbocca/na101_master/raw/master/homework_4_measures/rt_top_graph.png)
# 

# In[ ]:


# crear la red de retweets y guardarla en rt_top_graph
rt_top_graph = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: crear un data.frame de quien retweetea a quien, yy luego crear grafo desde el data.frame 
#
#
#
##################################################################

#plot
glay = layout_with_fr(rt_top_graph)
ver_labs = igraph::get.vertex.attribute(rt_top_graph, "name", index=V(rt_top_graph))
par(bg="white", mar=c(1,1,1,1))
plot(rt_top_graph, layout=glay, vertex.color="black", vertex.size=10, vertex.label=ver_labs,
     vertex.label.family="sans", vertex.shape="none", vertex.label.color=hsv(h=.165, s=.28, v=.08, alpha=1),
     vertex.label.cex=0.85, edge.arrow.size=0.8, edge.arrow.width=0.5, edge.width=3,
     edge.color=hsv(h=0, s=1, v=1, alpha=1))
title("5000 last tweets on #Uruguay : Who retweets whom",
      cex.main=1, col.main="black") # add title

# ## 4.d) Centralidad
# 
# Calcular las métricas de centralidad: `betweenness(), hub_score(), authority_score()` y comparar los cinco usuarios con mayor centralidad para cada una de ellas.  El resultado debería ser el de la siguiente Tabla
# 
# |metric|1|2|3|4|5|
# |---|---|---|---|---|---|
# | `betweenness()` | emekavoces | lsm_en_uy | jgamorin | boleadorcharrua | the_view_of_leo |
# | `hub_score()` | jgamorin | emekavoces | ciudadanos_mvd | pablolarraz10 | the_view_of_leo |
# | `authority_score()` | magelameinjua | dibonomiguel | elojochurrinche | arriolaernesto2 | pablolarraz10 | 
# 
# 

# In[ ]:


##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: calcular centralidades usando el grafo rt_top_graph 
#
#
#
##################################################################


# ## 4.e) (opcional)
# 
# En lugar de los datos sugeridos, puede repetirse el proceso utilizando datos propios. Esto puede hacerse mediante el proceso descripto en el ejercicio anterior o utilizando el paquete `twitteR` de `R`. Para hacerlo utilizando `twitteR` es necesario: 
# 1.   crear una cuenta de desarrollo en https://developer.twitter.com/, 
# 2.   crear una aplicación de Twitter para obtener credenciales de acceso al API (`consumer key, consumer secret, access token, access token secret`), 
# 3.   usar el API desde `R`.
# Puede por ejemplo descargarse los 5000 tweets más recientes de `#Uruguay` (o del tópico que se desee) y repetir las partes anteriores.
# 
# 
# 
# 

# In[ ]:


# Step 0: You need to use OAuth from Twitter
# 1. Log into the Twitter Developers section Apply for a account in: https://developer.twitter.com/.
# 2. Get Consumer Key & Consumer Secret, you have to create an app in Twitter via https://apps.twitter.com/app/new.
# 3. Then you’ll be taken to a page containing Consumer Key & Consumer Secret.

# collect recent 5000 tweets of #Uruguay
# you need to use your own key, which can be obtain from tweeter
api_key       = "paste here"
api_secret    = "paste here"
access_token  = "paste here"
access_token_secret = "paste here"
setup_twitter_oauth(api_key,api_secret,access_token,access_token_secret)

my_tweets = searchTwitter("#Uruguay", n=5000, lang=NULL,since=NULL, until=NULL,locale=NULL, 
                          geocode=NULL, sinceID=NULL, maxID=NULL,resultType=NULL, retryOnRateLimit=120)
my_tweets_df = twListToDF(my_tweets) # convert tweets to a data frame
write.csv(my_tweets_df, "my_tweets_df.csv")

