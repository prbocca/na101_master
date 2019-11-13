#!/usr/bin/env python
# coding: utf-8

# #1. Ambientes de trabajo
# 
# Lo primero es crear el ambiente de trabajo para todos los laboratorios. Hay dos opciones: realizarlo remotamente en COLAB de Google, o hacerlo localmente usando el programa RSTUDIO.
# 
# Resumidamente, RSTUDIO es más costoso al principio pero al final se gana tiempo (principalmente si se va a realizar mucho cómputo). Por otro lado, COLAB permite rápidamente ponerse a trabajar. Supondremos que se realiza el práctico en COLAB, pudiendo cambiar más adelante.

# ##1.a) Ambiente COLAB remoto
# 
# 1.   Abrir en navegador: https://colab.research.google.com/
# 2.   Abrir el notebook de la tarea:
#      File-> Open Notebook -> Github -> https://github.com/prbocca/na101_master -> homework_1_graphs
# 3.   Guardar el notebook en su Google Drive:
#      File -> Save a Copy in Drive... 
# 4.   Renombrar el archivo `"cedula ID"_ar_hw1.ipynb`, por ejemplo *33484022_ar_hw1.ipynb*
# 5.   Al final usted deberá descargar el notebook. Asegurarse que se están guardando las salidas de ejecución en el notebook: File -> Download .ipynb
# 6.   Luego estos archivos deberán ser enviados a prbocca@fing.edu.uy 
# 

# ##1.b) Ambiente RSTUDIO local (opcional)
# 
# Es posible instalar un ambiente de trabajo (IDE) local y los paquetes básicos necesarios para todos los prácticos.
# Este proceso puede demorar y se recomienda hacerlo fuera del práctico.
# 
# Este ambiente de trabajo es opcional y recomendado para realizar mucho cómputo. En lo que resta del curso se asumirá que se utiliza la plataforma en linea: https://colab.research.google.com/, en donde se puede ejectuar y resolver todos los prácticos del curso.
# 
# 
# 1.   Instalar Rstudio. Descargar desde: https://www.rstudio.com/products/rstudio/
# 2.   Instalar Paquetes básicos. En una consola de R:
# ```
# #  cargar librerias
# load_libs <- function(libraries = libs, install=TRUE){
#   if (install){ # instalar librerias no instaladas
#     new.packages <- libs[!(libs %in% installed.packages()[,"Package"])]
#     if(length(new.packages)) install.packages(new.packages)
#   }
#   #cargo librerias  
#   for (lib in libraries){
#     require(lib, character.only=TRUE, quietly = FALSE)
#   } 
# } 
# libs = c("igraph", "igraphdata", "GO.db", "GOstats", "ROCR", "ape", "car", "eigenmodel", "ergm", "fdrtool", "ggplot2", "huge", "kernlab", "lattice", "mixer", "network", "networkDynamic", "networkTomography", "ngspatial", "org.Sc.sgd.db", "sna", "vioplot") 
# load_libs(libs)
# ```
# 
# 3.    Abrir el .r de la tarea en: https://github.com/prbocca/na101_master/tree/master/homework_1_graphs

# ## 1.c) Cargar Librerias
# 
# En esta primera ejecución se asignan los recursos. Además instalar las librerías lleva tiempo. Es de esperar que esto demore unos minutos.
# 
# Lamentablemente, el entorno de Colab se pierde luego de un tiempo de inactividad. Esto incluye que se deben volver a cargar las librerías y a descargar los datos (se borran todos los datos temporales).

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

libs = c("R.matlab", "sand","igraph","igraphdata") #procesamiento en paralelo

load_libs(libs)

# ## 1.d) Descargar funciones auxiliares
# 
# Es necesario repetir esta tarea siempre que se inicia una ejecución.

# In[ ]:


#directorio donde se va a trabajar
data_path = "/content/ar/hw1/"

dir.create(data_path, showWarnings = FALSE, recursive = TRUE)
setwd(data_path)
getwd()
list.files()

# cargo funciones auxiliares
source("https://raw.githubusercontent.com/prbocca/na101_master/master/homeworks_common.r")

# # 2. Tutorial del lenguaje *R*
# 
# Es necesario conocer el lenguaje *R* para ganar fluidez en el práctico. Quien no conoce el lenguaje debe realizar algún tutorial rápido.
# , por ejemplo: `igraph_demo("crashR")` (que lamentablemente no anda en COLAB).
# 

# In[ ]:


# lamentablemente el paquete tcltk de COLAB no permite interfaz interactiva, por lo que este tutorial solo anda localemente (en el ambiente RSTUDIO).
library("tcltk")
if (interactive()){
  igraph_demo("crashR")
  }else{
    print("Lo siento en su ambiente no funciona este tutorial.")
  }


# # 3. Tutorial de manipular grafos de redes
# 
# Si aun no lo tiene, descargar el libro de práctico [SANDR]: Kolaczyk, E.D. and Csárdi, G. "Statistical Analysis of Network Data with R". Use R!, Springer New York, 2014. ISBN 9781493909834. Disponible
# en http://timbo.org.uy/.
# Ir a la descarga directamente siguiendo el link (Noviembre, 2019): https://link-springer-com.proxy.timbo.org.uy/book/10.1007%2F978-1-4939-0983-4
# 
# 
# Seguir el texto del Capítulo 2 del libro [SANDR], ejecutando el código fuente incluido. Utilizar el manual de referencia del paquete igraph para extender el conocimiento: http://igraph.org/r/doc/igraph.pdf
# 
# **TODO poner codigo fuente**

# # 4. Explorar una red de comunicaciones: emails de la empresa Enron
# 
# La red de emails de la empresa Enron se publicó durante una investigación federal de EE.UU (si usted esta interesado en leer más sobre el escándalo de Enron, ver http://en.wikipedia.org/wiki/Enron_scandal). 
# 
# En éste ejercicio haremos cálculos sencillos sobre el grafo de esta red usando *R*.
# 
# La red completa está disponible online en http://www.cs.cmu.edu/~enron/ y contiene $517.431$ correos extraidos de los directorios de correo de $150$ usuarios. Aqui utilizaremos un grafo más pequeño, preparado por C. E. Priebe et al., que consiste en $34.427$ correos enviados entre $184$ correos de empleados de Enron (ver http://cis.jhu.edu/~parky/Enron/enron.html por más detalles). Los correos son del periodo 13/11/1998 al 21/06/2002 ($44$ meses).
# 

# ## 4.a) Descargar los datos y crear matriz de adyacencia
# 
# Descargar el subconjunto de datos de `Enron.zip` de github.com

# In[ ]:


# load data
download.file(url="https://github.com/prbocca/na101_master/blob/master/homework_1_graphs/Enron.zip?raw=true", destfile="Enron.zip", mode="wb")
unzip(zipfile="Enron.zip")

list.files()

# $\newcommand   \ind   [1] {{\mathbb I \left\{#1\right\}  } }$
# $\def\bbA{{\mathbf A}}$
# El .zip contiene dos archivos creados en `Matlab`.
# Para cargar archivos de `Matlab` en `R`, utilizar el paquete `R.matlab`.
# 
# En el archivo `Y.mat` se encuentra el arreglo de tres dimensiones $\underline{\mathbf Y} \in \mathbb{R}^{184 \times 184 \times 44}$, que contiene la secuencia de matrices de adyacencia 
# que representa la red de quién le envió correo a quién para cada uno de los $44$ meses.
# Notar que el grafo de la red es dirigido y ponderado, 
# dado que el elemento $[\underline{\mathbf Y}]_{ijt}$ indica la cantidad de correos enviados por el empleado $i$ al empleado $j$ durante el mes $t$.
# 
# Además, desde el otro archivo, cargar el arreglo `employees` que contiene el nombre, la cuenta y cargo en la empresa para cada uno de los $184$ empleados (vértices) en $G$.
# 

# In[ ]:


# cargo la matriz de correos enviados
Y = readMat("Y.mat")$Y 
str(Y)

# cargo los datos de cada empleado
employees = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: usar readMat() y algunas funciones auxiliares para transformar el texto
#
#
#
##################################################################
str(employees)

# Construir la matriz binaria de adyacencia ${\mathbf A}\in\{0,1\}^{184\times184}$ que agrega todos los correos a lo largo del tiempo
# (es decir, sumar $\underline{\mathbf Y}$ en la tercera dimension, y definir $[{\mathbf A}]_{ij}=\ind{\sum_{t=1}^{44}[\underline{\mathbf Y}]_{ijt}>0}$). 
# 

# In[ ]:


# matriz de adyacencia con peso (cantidad de correos enviados en todos los meses)
A_weight = rowSums(Y, na.rm = FALSE, dims = 2)
str(A_weight)

# matriz de adyacencia binaria
A = ifelse(A_weight>0, 1, 0)
str(A)

# ## 4.b) Construir un grafo a partir de la matriz de adyacencia
# 
# Construir el grafo (con `igraph`) a partir de la matriz de adyacencia $\bbA$ y agregar los atributos de vértices del arreglo `employees`.

# In[ ]:


# podría construise usando la matriz de adyacencia binaria
#   g = graph.adjacency(A)

# pero optamos por hacerlo usando la matriz de adyacencia con pesos, así nos quedan en el grafo
g = graph.adjacency(A_weight, weighted = TRUE)
V(g)$id = employees

summary(g)

# Más adelante veremos como graficar este grafo. Por ahora, la siguiente es una visualización sencilla.

# In[ ]:


plot(g, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.size=3, vertex.label=NA,
     layout=layout_nicely(g))

# Ahora nos familiarizarnos con operaciones básicas sobre grafos.
# 
# Para cada uno de las siguientes partes, realizar los cálculos usando la matriz de adyacencia $\bbA$ (con álegebra de matrices), y usando el grafo (con métodos de `igraph`).
# 
# 

# ## 4.c) Cantidad de aristas (arcos)
# 
# Calcular la cantidad de aristas (arcos) en la red.
# Es decir, la cantidad de parejas únicas ordenadas $(u,v)\in E$ donde $u,v\in V$.

# In[ ]:


# numero de aristas dirigidas: 3007
cant_aristadirigida_matriz = sum(A)
printf("MATRIZ. las aristas dirigidas son: %s", cant_aristadirigida_matriz)

# con grafos
cant_aristadirigida_grafo = length(E(g))
printf("GRAFO. las aristas dirigidas también las calculamos con: %s o con %s", cant_aristadirigida_grafo, ecount(g))

# ## 4.d) Cantidad de aristas no dirigidas
# 
# Calcular la cantidad de aristas no dirigidas en la red (es decir, la cantidad de parejas únicas no ordenadas $(u,v)\in E$
# donde $u,v\in V$. Esto significa que si $(u,v)\in E$ o $(v,u)\in E$, entonces se debe contar el par como una única arista.
# 

# In[ ]:


# con matrices
#numero de aristas no dirigidas: 2097
A.und = ifelse((A + t(A)) > 0, 1, 0)
A.und.triang = A.und
A.und.triang[lower.tri(A.und.triang)] = 0
cant_aristanodirigida_matriz = sum(A.und.triang)
printf("MATRIZ. las aristas no dirigidas son: %s", sum(cant_aristanodirigida_matriz))

#con grafos
cant_aristanodirigida_grafo = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: convertir el grafo a no dirigído                         
#
#
##################################################################
printf("GRAFO. las aristas no dirigidas también las calculamos con: %s", cant_aristanodirigida_grafo)


# ## 4.e) Cantidad de arcos mutuos
# 
# Calcular la cantidad de arcos mutuos en la red (es decir, la cantidad de parejas $(u,v)$ donde $\{(u,v),(v,u)\}\subseteq E$
# y $u,v\in V$). Esto significa que si existen ambas $(u,v)\in E$ y $(v,u)\in E$, debe contarse el par como mutual.

# In[ ]:


# con matrices
cant_aristamutua_matriz = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: usar una matriz triangular luego de transformar la matriz de adyacencia
#
#
##################################################################
printf("MATRIZ. las aristas mutuas son: %s", cant_aristamutua_matriz)


#con grafos
cant_aristamutua_grafo = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: convertir el grafo a no dirigído                         
#
#
##################################################################
printf("GRAFO. las aristas mutuas también las calculamos con: %s", cant_aristamutua_grafo)


# ## 4.f) Empleados que no recibieron ningún correo
# 
# Calcular la cantidad de vértices con $d_v^{\text{in}}=0$, y listar los nombres de los empleados correspondientes.
# 

# In[ ]:


#find people not geting any email: 3
# [1] "Vince Kaminski (j..kaminski) Manager Risk Management Head" "Mary Fischer (mary.fischer) Employee "                    
# [3] "xxx Smith (m..smith) xxx " 
idx.notgetting = colSums(A)==0
cant_notgetting_matriz = sum(idx.notgetting)
printf("MATRIZ. los empleados que no recibieron correo son: %s", cant_notgetting_matriz)
employees[idx.notgetting]

#con grafos
cant_notgetting_grafo = sum(degree(g, mode = "in")==0)
printf("GRAFO. los empleados que no recibieron correo son: %s", cant_notgetting_grafo)
V(g)$id[degree(g, mode = "in")==0]

# ## 4.g) Empleados que no enviaron ningún correo
# 
# Calcular la cantidad de vértices con $d_v^{\text{out}}=0$, y listar los nombres de los empleados correspondientes. 
# 

# In[ ]:


# con matrices
cant_notsending_matriz = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: Similar a quienes no recibieron
#
#
##################################################################
printf("MATRIZ. los empleados que no enviaron correo son: %s", cant_notsending_matriz)


# con grafos
cant_notsending_grafo = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: Similar a quienes no recibieron
#
#
##################################################################
printf("GRAFO. los empleados que no enviaron correo son: %s", cant_notsending_grafo) 


# ## 4.h) Cantidad de empleados contactados por al menos otros 30 empleados
# 
# Calcular la cantidad de empleados que fueron contactados por $30$ o más empleados.

# In[ ]:


# con matrices
cant_indegree30_matriz = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: Similar a quienes no recibieron
#
#
##################################################################
printf("MATRIZ. los contactados por al menos otros 30 empleados son: %s", cant_indegree30_matriz)


# con grafos
cant_indegree30_grafo = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: Similar a quienes no recibieron
#
#
##################################################################
printf("GRAFO. los contactados por al menos otros 30 empleados son: %s", cant_indegree30_grafo)



# ## 4.i) Cantidad de empleados que contactaron a al menos otros 30 empleados
# 
# Calcular la cantidad de empleados que contactaron a $30$ o más empleados.

# In[ ]:


#con matrices
cant_outdegree30_matriz = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: Similar a quienes no recibieron
#
#
##################################################################
printf("MATRIZ. los que contactaron a al menos otros 30 empleados son: %s", cant_outdegree30_matriz)

#con grafos
cant_outdegree30_grafo = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: Similar a quienes no recibieron
#
#
##################################################################
printf("GRAFO. los que contactaron a al menos otros 30 empleados son: %s", cant_outdegree30_grafo)



# ## 4.j) Cantidad de triángulos dirigidos
# 
# Calcular la cantidad de triángulos dirigidos en $G$. 
# Recordar que $G$ es dirigido y la orientación de los triángulos importa, es decir, $i\to j \to k \to i$ es lo mismo que $k \to i \to j \to k $, pero diferente a $i\to k \to j \to i$.
# 

# In[ ]:


#numero de triangulos dirigidos: 6819
printf("MATRIZ. Cantidad de triángulos dirigidos son: %s", sum(diag(A %*% A %*% A/3)))

# con matrices
cant_triangnodirigido_matriz = NA
##################################################################
#
#                       TU CÓDIGO ACÁ                           
# Tip: similar a anterior
#
#
##################################################################
printf("MATRIZ. Cantidad de triángulos no dirigidos son: %s", cant_triangnodirigido_matriz)


# con grafos
cant_triangnodirigido_grafo = NA
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: Buscar una función en igraph que ayude
#
#
##################################################################
printf("GRAFO. Cantidad de triángulos no dirigidos son: %s", cant_triangnodirigido_grafo)


# ## 4.k) ¿Es el grafo conexo?
# 
# Solo usando grafos, determinar si el grafo dirigido es débilmente y fuértemente conexo. 
# En caso de no serlo obtener el subgrafo con la componente conexa más grande (llamada componente gigante).

# In[ ]:


printf("GRAFO. Es débilmente conexo: %s", is.connected(g, mode = "weak"))
printf("GRAFO. Es fuertemente conexo: %s", is.connected(g, mode = "strong"))

g_components = components(g, mode = c("weak", "strong"))
str(g_components)

c_gigante = which.max(g_components$csize) #la componente mas grande
vids_gigante =  which(g_components$membership == c_gigante)
g_gigante = induced_subgraph(g, vids_gigante)
summary(g_gigante)
V(g)$id[g_components$membership != c_gigante] 
# sobran los vertices:
#'Vince Kaminski (j..kaminski) Manager Risk Management Head'
# 'Mary Fischer (mary.fischer) Employee '

plot(g_gigante, edge.color="gray", edge.width=1, edge.lty=1, edge.arrow.size=.5, 
     vertex.size=3, vertex.label=NA,
     #vertex.color=nodes, vertex.shape = shapes,
     layout=layout_nicely(g_gigante))


# ## 4.l) Guardar y cargar grafos
# 
# Usando el paquete `igraph` guardar el grafo $G$ en formato PAJEK.
# Luego cargar una copia del grafo en $G_\text{copia}$.
# ¿Porqué los dos grafos no son identicos?

# In[ ]:


#grabar els grafo G
filename_g = "enron_g.pjk"
##################################################################
#                       TU CÓDIGO ACÁ                           
# Tip: Buscar una función en igraph que ayude
#
#
##################################################################


# cargar el grafo
g_copia = read.graph(file=filename_g, format= "pajek")
summary(g_copia)

#son iguales?
identical(g, g_copia)
