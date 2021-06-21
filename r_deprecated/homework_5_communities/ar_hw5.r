
# coding: utf-8

# # Práctico 5 - Cohesión y detección de comunidades en redes reales
# 
# #1. Ambientes de trabajo
# 
# ## 1.a) Ambiente COLAB remoto
# 
# 1.   Abrir en navegador: https://colab.research.google.com/
# 2.   Abrir el notebook de la tarea:
#      File-> Open Notebook -> Github -> https://github.com/prbocca/na101_master -> homework_5_communities
# 3.   Guardar el notebook en su Google Drive:
#      File -> Save a Copy in Drive... 
# 4.   Renombrar el archivo `"cedula ID"_ar_hw5.ipynb`, por ejemplo *33484022_ar_hw5.ipynb*
# 5.   Al final usted deberá descargar el notebook. Asegurarse que se están guardando las salidas de ejecución en el notebook: File -> Download .ipynb
# 6.   Luego estos archivos deberán ser enviados a prbocca@fing.edu.uy 
# 
# ##1.b) Ambiente RSTUDIO local (opcional)
# 
# Abrir el .r de la tarea en: https://github.com/prbocca/na101_master/tree/master/homework_5_communities

# ## 1.c) Cargar Librerias
# 
# Todas las librerías debe instalarse correctamente, si el proceso se interrumpe o alguna librería da error en la instalación, entonces habrá problemas en el código más adelante. Si tiene este tipo de problemas pruebe hacer `Runtime -> Factory reset runtime`, y volver a intentar.
# 

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

libs = c("network", "sna","rgexf","ape", 
         "R.matlab", "sand","igraph","igraphdata") 

load_libs(libs)


# ## 1.d) Descargar funciones auxiliares

# In[ ]:

#directorio donde se va a trabajar
data_path = "/content/ar/hw5/"

dir.create(data_path, showWarnings = FALSE, recursive = TRUE)
setwd(data_path)
getwd()
list.files()

# cargo funciones auxiliares
source("https://raw.githubusercontent.com/prbocca/na101_master/master/homeworks_common.r")


# # 2. Cohesión en grafos en `R`
# 
# Seguir la Secciones 4.3 del libro [SANDR], ejecutando el código fuente incluido.

# In[ ]:

#4.3 Characterizing Network Cohesion
data(karate)
table(sapply(cliques(karate), length))
# 1 2 3 4 5
# 34 78 45 11 2
# there are 34 nodes (cliques of size one) and 78 edges (cliques of size two), followed by 45 triangles (cliques of size three)
#
# the largest cliques are of size five, of which there are only two. 
# These two both involve four actors in common, including Mr Hi, the head instructor.
cliques(karate)[sapply(cliques(karate), length) == 5]
# [[1]]
# + 5/34 vertices, named, from 4b458a1:
#   [1] Mr Hi    Actor 2  Actor 3  Actor 4  Actor 14
# 
# [[2]]
# + 5/34 vertices, named, from 4b458a1:
#   [1] Mr Hi   Actor 2 Actor 3 Actor 4 Actor 8
#
table(sapply(maximal.cliques(karate), length))
# 2  3  4  5 
# 11 21  2  2 
# In the karate network, the two largest cliques (formally called maximum cliques) are maximal, while, for example, the same can
# be said of only two of the 11 cliques of size four.


# In[ ]:

# representacion en k-core
# esta es una forma relativamente común para visualizar grafos muy grandes (>10k nodos)
cores <- graph.coreness(karate)
A <- get.adjacency(karate, sparse=FALSE)
library(network)
g <- network::as.network.matrix(A)
sna::gplot.target(g, cores, 
                  #circ.lab = FALSE,
                  circ.col="skyblue", usearrows = FALSE,
                  vertex.col=cores, edge.col="darkgray")
#detach("package:sna")
#detach("package:network")


# In[ ]:

# 4.3.2 Density and Related Notions of Relative Frequency
# we see that the sub-graphs corresponding to each of the instructor and the administrator, in union with
# their immediate respective neighborhoods—i.e., the ego-centric networks around
# vertices 1 and 34—are noticeably more dense than the overall network.
ego.instr <- induced.subgraph(karate, neighborhood(karate, 1, 1)[[1]])
ego.admin <- induced.subgraph(karate, neighborhood(karate, 1, 34)[[1]])
graph.density(karate) #0.1390374
graph.density(ego.instr) #0.25
graph.density(ego.admin) #0.2091503
#
# global clustering, summarizing the relative frequency with which connected triples close to form triangles
transitivity(karate) # 0.2556818
# In the case of the instructor and administrator of the karate network, for
# example, we see that their local clustering is only 50–60 % that of the clustering for the network as a whole.
transitivity(karate, "local", vids=c(1,34)) #0.1500000 0.1102941
#
# reciprocity is defined as the total number of reciprocated edges divided by the total number of edges.
# In the AIDS blog network, the reciprocity is quite low by either definition.
data(aidsblog)
reciprocity(aidsblog, mode="default") #[1] 0.03278689
reciprocity(aidsblog, mode="ratio") #[1] 0.01666667


# In[ ]:

# 4.3.3 Connectivity, Cuts, and Flows
data(yeast)
is.connected(yeast)
#
#A census of all connected components within this graph, however, shows that there clearly is a giant component.
# This single component contains 2375/2617 ≈ 90 %
comps <- decompose.graph(yeast)
table(sapply(comps, vcount))
# 2    3    4    5    6    7 2375 
# 63   13    5    6    1    3    1
#
# often attention would be restricted to this giant component alone in carrying out further analysis and modeling.
# small world property:
yeast.gc <- decompose.graph(yeast)[[1]]
average.path.length(yeast.gc) #5.09597
diameter(yeast.gc) #15
# At the same time, the clustering in this network is relatively large
transitivity(yeast.gc)
#
# The vertex (edge) connectivity of G is the largest integer such that G is k-vertex- (k-edge-) connected.
# In the case of the giant component of the yeast network, the vertex and edge connectivity are both equal to one.
vertex.connectivity(yeast.gc) #1
edge.connectivity(yeast.gc) #1
# If the removal of a particular set of vertices (edges) in a graph disconnects the
# graph, that set is called a vertex-cut (edge-cut). A single vertex that disconnects
# the graph is called a cut vertex, or sometimes an articulation point.
# In the giant component of the yeast network, almost 15 % of the vertices are cut vertices.
yeast.cut.vertices <- articulation.points(yeast.gc)
length(yeast.cut.vertices) #[1] 350

# Note that the distinction between strong and weak connectivity can be severe for some digraphs. 
# For example, the AIDS blog network is weakly connected
is.connected(aidsblog, mode=c("weak")) #2 [1] TRUE
# but not strongly connected.
is.connected(aidsblog, mode=c("strong")) #[1] FALSE
# And while there does exist a strongly connected component within the graph, there is only one and it has only four vertices.
aidsblog.scc <- clusters(aidsblog, mode=c("strong"))
table(aidsblog.scc$csize)
# 1   4 
# 142   1 


# # 3. Particionar grafos en `R`. 
# 
# Seguir la secciones 4.4 al 4.6 del libro [SANDR], ejecutando el código fuente incluido.
# 
# 

# In[ ]:

#4.4 Graph Partitioning
#4.4.1 Hierarchical Clustering 

kc <- fastgreedy.community(karate)
length(kc) #[1] 3
sizes(kc) 
membership(kc)
# Community sizes
# 1  2  3 
# 18 11  5 
plot(kc, karate)
#library(ape)
dendPlot(kc, mode="phylo")


# In[ ]:

# 4.4.2 Spectral Partitioning
k.lap <- graph.laplacian(karate)
eig.anal <- eigen(k.lap)
plot(eig.anal$values, col="blue", ylab="Eigenvalues of Graph Laplacian")
f.vec <- eig.anal$vectors[, 33]
faction <- get.vertex.attribute(karate, "Faction")
f.colors <- as.character(length(faction))
f.colors[faction == 1] <- "red"
f.colors[faction == 2] <- "cyan"
plot(f.vec, pch=16, xlab="Actor Number",
     ylab="Fiedler Vector Entry", col=f.colors)
abline(0, 0, lwd=2, col="lightgray")


# In[ ]:

# 4.4.3 Validation of Graph Partitioning
func.class <- get.vertex.attribute(yeast.gc, "Class")
table(func.class)
yc <- fastgreedy.community(yeast.gc)
c.m <- membership(yc)
table(c.m, func.class, useNA=c("no"))


# In[ ]:

# 4.5 Assortativity and Mixing
# The assortativity coefficient with categories
assortativity.nominal(yeast, (V(yeast)$Class=="P")+1, directed=FALSE) #[1] 0.4965229
# assortativity mixing with continuos attributes (Pearson correlation coefficient)    

assortativity.degree(yeast) #0.4610797


# # 4. Particionar la red de blogs políticos de EE.UU.
# 
# Estudiaremos la red de blogs políticos de EE.UU, 
# con el objetivo de particionarla en las dos comunidades políticas existenes: liberales (demócratas) y conservadores (republicanos).
# Los datos son de la elección política de EE.UU. en 2004, 
# fueron recolectados por L. Adamic and N. Glance en 2005, 
# y pueden obtenerse de la colección de Mark Newman en: http://www-personal.umich.edu/~mejn/netdata/. 
# 
# En este ejercicio usaremos una versión no dirigida del grafo dirigido original,
# donde las aristas corresponden a *hyperlinks* entre blogs. 
# La red tiene $N_v=1490$ blogs (vértices), 
# y se conoce la afiliación política de cada blogger (y por tanto de sus blogs), representada por un vector binario (0 es liberal, y 1 es conservador) para cada vértice.
# 

# $\newcommand   \ind   [1] {{\mathbb I \left\{#1\right\}  } }$
# $\def\bbA{{\mathbf A}}$
# $\def\bbB{{\mathbf B}}$
# 
# ## 4.a) Descargar y cargar los datos.
# 
# Descargar los datos en formato `Matlab` de la blogósfera política en el archivo `political_blogs.mat` de la página del curso, o en el siguiente link: 
# https://github.com/prbocca/na101_master/raw/master/homework_5_communities/political_blogs.mat.
# 
# Cargar los datos en `R`, que incluyen la matriz de adyacencia $\bbA \in \{0,1\}^{1490\times 1490}$, 
# y el vector ${\tt nodes}\in\{0,1\}^{1490}$ que contiene la afiliación política de cada blog.

# In[ ]:

# download data
download.file(url="https://github.com/prbocca/na101_master/raw/master/homework_5_communities/political_blogs.mat", destfile="political_blogs.mat", mode="wb")
list.files()

# cargar datos
A = NA
nodes = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: cargar la matriz de adyacencia 'A', y el vector con datos de los nodos 'nodes'
#
#
#
##################################################################
str(A)
str(nodes)


# Crear un grafo no dirigido a partir de los datos (matriz de adyacencia y afiliación). El resultaod debe ser similar al de la siguiente Figura: ![alt text](https://github.com/prbocca/na101_master/raw/master/homework_5_communities/political_blogs.png)
# 

# In[ ]:

# La matriz original tiene algunos loops, no son de nuestro interés y los eliminamos
diag(A) = 0

# cargar el grafo no dirigido, g, con los datos de afinidad en un atributo "nodes" de los vértices
g = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: -
#
#
#
##################################################################
summary(g)

#plot
set.seed(1)
plot(g, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.color=nodes, vertex.size=3, vertex.label=NA,
     layout=layout_components(g))


# ## 4.b) Matriz de modularidad
# 
# Calcular el grado $d_v$ de todos los vértices $v\in V$ y el total de aristas $N_e$. Entonces calcular la matriz de modularidad $\bbB$ de la red. Comparar el resultado con el de la función `modularity_matrix()`.
# 

# In[ ]:

#numero de aristas: 16715
printf("El numero de aristas es: %s", length(E(g)))

# calcular B
B = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: - calcular B a partir del grado
#
#
#
##################################################################
str(B)

# usando igraph
B2 = modularity_matrix(g, membership =rep(1,vcount(g)))
str(B2)
identical(B,B2)


# ## 4.c) Particionado goloso
# 
# Realizar el particionado goloso rápido (clustering jerárquico aglomerativo) en varias comunidades (`fastgreedy.community()`).
# 
# 

# In[ ]:

g.fg = fastgreedy.community(g, merges = T)

printf("Se encontraron %s comunidades", length(g.fg))
printf("La mayoría de las comunidades son muy pequeñas:")
sizes(g.fg)


# Dado que el metodo de particionado es jerárquico es posible definir la cantidad de comunidades que se quiere. Lamentablemente esto no siempre funciona, y este es el caso: no podemos particionar en dos usando este método. Ver el siguiente código.

# In[ ]:

cutat(g.fg, no=2) # no me deja partirlo en dos


# Por tanto, voy a predecir la afiliación utilizando solo las dos comunidades más importantes. El resto de los vértices los dejo sin clasificar.
# 
# Para esto me creo la matriz de confusión entre la afiliación real y la predicha. Conocer más de la matriz de confusión en: https://es.wikipedia.org/wiki/Matriz_de_confusi%C3%B3n

# In[ ]:

# me creo la matriz de confusión entre la afiliación real (nodes), y la predicción (usando solo las dos comunidades mas importantes)
g.fg.mem = membership(g.fg)
g.fg.t = as.data.frame(table(nodes, g.fg.mem), stringsAsFactors = F)
g.fg.t = g.fg.t[order(g.fg.t$Freq, decreasing = T),][1:2,] #obtiene el mejor mapeo de las dos comunidades mas importantes

if (g.fg.t[1,1]==0){
  comm_of_0 = g.fg.t[1,2]
  comm_of_1 = g.fg.t[2,2]
} else{
  comm_of_0 = g.fg.t[2,2]
  comm_of_1 = g.fg.t[1,2]
}
printf("El mejor mapeo para la afiliación 0 es la comunidad %s, y para la afiliación 1 es la comunidad %s", comm_of_0, comm_of_1)

printf("La matriz de confusión es:")
g.fg.mem.bin = ifelse(g.fg.mem %in% c(comm_of_0, comm_of_1), g.fg.mem, "other")
cm.fg = as.matrix(table(nodes, g.fg.mem.bin))
cm.fg = cm.fg[,c(comm_of_0,comm_of_1, "other")] #reordeno las columnas 
cm.fg

printf("La exactitud del método es: %s", sum(diag(cm.fg))/sum(cm.fg)) #porcentaje de bien detectados 0.7583893
# También se pudo calcular con sum(g.fg.t$Freq)/vcount(g)


# Puedo volver a graficar la red, agregando la información de cuando no realizo una buena predicción (vértices cuadrados), con el objetivo de buscar algún patrón en el error.

# In[ ]:

# plot
set.seed(1)
shapes = rep("square",vcount(g)) #dibujo con cuadrado los que estan sin clasificiar o mal clasificados
shapes[(nodes==0)&(g.fg.mem==comm_of_0)] = "circle" # dibujo con circulos los bien clasificados como 0
shapes[(nodes==1)&(g.fg.mem==comm_of_1)] = "circle" # dibujo con circulos los bien clasificados como 1
plot(g, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.color=nodes, vertex.size=3, vertex.label=NA,
     vertex.shape = shapes,
     layout=layout_components(g))

# se ve que principalmente los cuadrados estan afuera de la componente gigante 


# ## 4.d) Particionado espectral
# 
# Realizar el particionado espectral (maximización espectral de modularidad) en varias comunidades (`leading.eigenvector.community()`).

# In[ ]:

g.part = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a anterior
#
#
#
##################################################################

printf("Se encontraron %s comunidades", length(g.part))
printf("La mayoría de las comunidades son muy pequeñas:")
sizes(g.part)


# Dado que el metodo de particionado es jerárquico es posible definir la cantidad de comunidades que se quiere. ¿Es posible generar solo dos comunidades con este método? 

# In[ ]:

##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a anterior
#
#
#
##################################################################


# Entonces, voy a predecir la afiliación utilizando solo las dos comunidades más importantes. El resto de los vértices los dejo sin clasificar.
# 
# Crear la matriz de confusión entre la afiliación real y la predicha. Y calcular la exactitud del método.

# In[ ]:


# calcular la matriz de confusión
cm.part = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a anterior
#
#
#
##################################################################
cm.part

printf("La exactitud del método es: %s", sum(diag(cm.part))/sum(cm.part)) #porcentaje de bien detectados 0.772483221
# También se pudo calcular con sum(g.part.t$Freq)/vcount(g)


# Graficar la red, agregando la información de cuando no realizo una buena predicción (vértices cuadrados), con el objetivo de buscar algún patrón en el error.

# In[ ]:

# plot
set.seed(1)
shapes = rep("square",vcount(g)) #dibujo con cuadrado los que estan sin clasificiar o mal clasificados
shapes[(nodes==0)&(g.part.mem==comm_of_0)] = "circle" # dibujo con circulos los bien clasificados como 0
shapes[(nodes==1)&(g.part.mem==comm_of_1)] = "circle" # dibujo con circulos los bien clasificados como 1
plot(g, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.color=nodes, vertex.size=3, vertex.label=NA,
     vertex.shape = shapes,
     layout=layout_components(g))

# nuevamente se ve que los cuadrados estan principalmente afuera de la componente gigante


# ## 4.e) Particionado espectral con matrices
# 
# ¿Cómo detectar solo dos comunidades con los resultados o funciones anteriores? Si no es posible, entonces implementar el algoritmo espectral de maximización de modularidad visto en teórico. Además, crear la matriz de confusión entre la afiliación real y la predicha. Y calcular la exactitud del método.

# In[ ]:

# crear el dataframe g.mymod.t, con misma estructura de los casos anteriores,
# pero hacerlo calculando la modularidad aplicando el método de teórico
g.mymod.t = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: calcular la modularidad aplicando el método de teórico
#
#
#
##################################################################
g.mymod.t

# crear la matriz de confusión
cm.mymod = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: similar a los anteriores
#
#
#
##################################################################
cm.mymod

printf("La exactitud del método es: %s", sum(diag(cm.mymod))/sum(cm.mymod))  #porcentaje de bien detectados 0.8711409
# También se pudo calcular con sum(g.mymod.t$Freq)/vcount(g)


# Graficar la red, agregando la información de cuando no realizo una buena predicción (vértices cuadrados), con el objetivo de buscar algún patrón en el error.

# In[ ]:

# plot
set.seed(1)
shapes = rep("square",vcount(g)) #dibujo con cuadrado los que estan sin clasificiar o mal clasificados
shapes[(nodes==0)&(g.mymod.mem==comm_of_0)] = "circle" # dibujo con circulos los bien clasificados como 0
shapes[(nodes==1)&(g.mymod.mem==comm_of_1)] = "circle" # dibujo con circulos los bien clasificados como 1
plot(g, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.color=nodes, vertex.size=3, vertex.label=NA,
     vertex.shape = shapes,
     layout=layout_components(g))

#mejora visiblemente respecto a los anteriores


# ## 4.f) Predicción de afinidad política.
# 
# Si usaramos los métodos anteriores para predecir la afinidad política. ¿Cuál sería el mejor método de acuerdo a la exactitud?

# In[ ]:

printf("La exactitud del método Goloso es: %s", sum(diag(cm.fg))/sum(cm.fg)) #porcentaje de bien detectados 0.7583893
printf("La exactitud del método Particionado Espectral es: %s", sum(diag(cm.part))/sum(cm.part)) #porcentaje de bien detectados 0.772483221
printf("La exactitud del método Particionado Espectral con matrices es: %s", sum(diag(cm.mymod))/sum(cm.mymod))  #porcentaje de bien detectados 0.8711409
printf("El mejor método es el desarrollado por nosotros, el de maximizacion de modularidad: %s", sum(diag(cm.mymod))/sum(cm.mymod))

