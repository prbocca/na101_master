

# Funciones auxiliares
#########################################################

qpaste = function(...) paste(..., sep="")
printf <- function(...) invisible(print(sprintf(...)))


# Converts the given igraph object to GEXF format and saves it at the given filepath location
#     g: input igraph object to be converted to gexf format
#     filepath: file location where the output gexf file should be saved
saveAsGEXF = function(g, filepath="converted_graph.gexf")
{
  require(igraph)
  require(rgexf)
  
  # gexf nodes require two column data frame (id, label)
  # check if the input vertices has label already present
  # if not, just have the ids themselves as the label
  if(is.null(V(g)$label))
    V(g)$label <- as.character(V(g))
  
  # similarily if edges does not have weight, add default 1 weight
  if(is.null(E(g)$weight))
    E(g)$weight <- rep.int(1, ecount(g))
  
  nodes <- data.frame(cbind(V(g), V(g)$label))
  edges <- t(Vectorize(ends, vectorize.args='es')(g, 1:ecount(g)))
  
  # combine all node attributes into a matrix (and take care of & for xml)
  vAttrNames <- setdiff(vertex_attr_names(g), "label") 
  nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&",vertex_attr(g, attr))))
  
  # combine all edge attributes into a matrix (and take care of & for xml)
  eAttrNames <- setdiff(edge_attr_names(g), "weight") 
  edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&",edge_attr(g, attr))))
  
  # combine all graph attributes into a meta-data
  graphAtt <- sapply(graph_attr_names(g), function(attr) sub("&", "&",graph_attr(g, attr)))
  
  # generate the gexf object
  output <- write.gexf(nodes, edges, 
                       edgesWeight=E(g)$weight,
                       edgesAtt = edgesAtt,
                       nodesAtt = nodesAtt,
                       defaultedgetype = ifelse(is_directed(g), "directed", "undirected"),
                       meta=c(list(creator="Gopalakrishna Palem", description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt))
  
  print(output, filepath, replace=T)
}



score = function(g, type="-dist", sort = FALSE, scale=TRUE){
  
  printf("%s init scoring %s",Sys.time(),type)
  nv <- vcount(g)
  result = data.frame()
  
  if (type=="-dist"){
    for(i in (1:(nv-1))){
      temp = distances(g, v = i, to = (i+1):nv)
      result = rbind(result, data.frame(score=-temp[1,], 
                                        edge=paste(i,"--",(i+1):nv, sep = ""), 
                                        from=rep(i,length(temp[1,])),
                                        to=(i+1):nv))
    }
  }
  
  
  if (type=="1/dist"){
    for(i in (1:(nv-1))){
      temp = distances(g, v = i, to = (i+1):nv)
      result = rbind(result, data.frame(score=1/temp[1,], 
                                        edge=paste(i,"--",(i+1):nv, sep = ""), 
                                        from=rep(i,length(temp[1,])),
                                        to=(i+1):nv))
    }
  }
  
  
  # #common neighbors
  # # ineficiente
  # if (type=="common neighbors2"){
  #   for(i in (1:(nv-1))){
  #     ni <- unlist(neighborhood(g, order=1, nodes=i))
  #     for (j in (i+1):nv){
  #       nj <- unlist(neighborhood(g, order=1, nodes=j))
  #       ni_j = setdiff(ni,c(i,j))
  #       nj_i = setdiff(nj,c(i,j))
  #       nbhd.ij <- length(intersect(ni_j, nj_i))
  #       result = rbind(result, data.frame(score=nbhd.ij, 
  #                                         edge=paste(i,"--",j, sep = ""), 
  #                                         from=i,
  #                                         to=j))
  #     }
  #   }
  # }
  #common neighbors
  if (type=="common neighbors"){
    A <- get.adjacency(g)
    for(i in (1:(nv-1))){
      ni <- neighborhood(g, order=1, nodes=i)
      j_v = (i+1):nv
      nj <- neighborhood(g, 1, nodes=j_v)
      nbhd.ij <- mapply(intersect, ni, nj, SIMPLIFY=FALSE)
      temp <- unlist(lapply(nbhd.ij, length)) - 2 * A[i, j_v]
      result = rbind(result, data.frame(score=temp,
                                        edge=paste(i,"--",j_v, sep = ""),
                                        from=rep(i,length(j_v)),
                                        to=j_v))
      
    }
  }

  # #normalized common neighbors
  # #ineficiente
  # if (type=="normalized common neighbors2"){
  #   for(i in (1:(nv-1))){
  #     ni <- unlist(neighborhood(g, order=1, nodes=i))
  #     for (j in (i+1):nv){
  #       nj <- unlist(neighborhood(g, order=1, nodes=j))
  #       ni_j = setdiff(ni,c(i,j))
  #       nj_i = setdiff(nj,c(i,j))
  #       nbhd.ij <- length(intersect(ni_j, nj_i))
  #       nbhd.ij.union <- length(union(ni_j, nj_i))
  #       result = rbind(result, data.frame(score=ifelse(nbhd.ij.union==0,-Inf,nbhd.ij/nbhd.ij.union), 
  #                                         edge=paste(i,"--",j, sep = ""), 
  #                                         from=i,
  #                                         to=j))
  #     }
  #   }
  # }
  if (type=="normalized common neighbors"){
    A <- get.adjacency(g)
    for(i in (1:(nv-1))){
      ni <- neighborhood(g, order=1, nodes=i)
      j_v = (i+1):nv
      nj <- neighborhood(g, 1, nodes=j_v)
      nbhd.ij <- mapply(intersect, ni, nj, SIMPLIFY=FALSE)
      nbhd.ij.union <- mapply(union, ni, nj, SIMPLIFY=FALSE)
      temp <- unlist(lapply(nbhd.ij, length)) - 2 * A[i, j_v]
      temp.union <- unlist(lapply(nbhd.ij.union, length)) - 2
      result = rbind(result, data.frame(score=temp/temp.union,
                                        edge=paste(i,"--",j_v, sep = ""),
                                        from=rep(i,length(j_v)),
                                        to=j_v))
      
    }
  }
  
  
  #degree product
  if (type=="degree product"){
    for(i in (1:(nv-1))){
      di = degree(g,v= i)
      dj = degree(g, v= (i+1):nv)
      result = rbind(result, data.frame(score=di*dj, 
                                        edge=paste(i,"--",(i+1):nv, sep = ""), 
                                        from=rep(i,length(dj)),
                                        to=(i+1):nv))
    }
  }
  
  if (scale){ #normalizo 0-1
    result$score = (result$score - min(result$score))/max(result$score - min(result$score))
  }
  if (sort){ #ordeno decreciente
    result = result[order(result$score, decreasing = TRUE),] 
  }
  return(result)
}


delete_edges_rand = function(g, p=0.3){
  edges_to_delete = E(g)[ifelse(runif(ecount(g))<p,TRUE,FALSE)]
  g_result = delete_edges(g, edges_to_delete)
  return(g_result)
}

last = function(x, ...){
  tail(x, n=1, ...)
}

linelight <- function(x,y, lty='dashed', col='lightgray', ...) {
  # highlight a point with lines running to the axes.
  left = par('usr')[1]
  bot = par('usr')[3]
  segments(left,y, x,y, lty=lty, col=col, ...)
  segments(x,bot,  x,y, lty=lty, col=col, ...)
}


#https://brenocon.com/blog/2009/04/binary-classification-evaluation-in-r-via-rocr/
binary_eval <- function(pred,labels, cutoff='naive', repar=TRUE, ...) {
  # Various binary classification evaluation plots and metrics
  library(ROCR)
  # plot(performance(prediction(pred,y),'acc'))
  rocr_pred = prediction(pred,labels)
  acc = performance(rocr_pred,'acc')
  f1 = performance(rocr_pred,'f')
  auc = performance(rocr_pred,'auc')@y.values[[1]]
  roc = performance(rocr_pred,'rec','spec')
  bac = if (rocr_pred@n.pos[[1]] != rocr_pred@n.neg[[1]])
    sapply(1:length(roc@x.values[[1]]), function(i)
      mean(c(roc@x.values[[1]][i], roc@y.values[[1]][i])))
  else
    rep(-1,length(pred))
  # sensspec = performance(rocr_pred,'rec','spec')
  pr_curve = performance(rocr_pred,'prec','rec')
  rp_curve = performance(rocr_pred,'rec','prec')
  
  printf("AUC = %.3f\n", auc)
  
  if (cutoff=='naive') {
    if (all(pred>=0) & all(pred<=1)) {
      printf("Predictions seem to be probabilities, so ")
      cutoff = 0.5
    } else if (any(pred<0) & any(pred>0)) {
      printf("Predictions seem to be real-valued scores, so ")
      cutoff = 0
    } else {
      warning("cant tell what naive cutoff should be")
      cutoff = NULL
    }
    printf("using naive cutoff %s:\n", cutoff)
  } else if (class(cutoff)=='character') {
    printf("Using %s-best cutoff ", cutoff)
    if (cutoff=='bac') {
      perf = NULL
      perf_y = bac
    } else {
      perf = performance(rocr_pred, cutoff, ...)
      perf_y = perf@y.values[[1]]
    }
    cutoff_ind = which.max(perf_y)
    cutoff = if (cutoff=='prbe') perf@x.values[[1]][1] else rocr_pred@cutoffs[[1]][cutoff_ind]
    printf("%f\n", cutoff)
  } else {
    printf("For cutoff %s:\n", cutoff)
  }
  cutoff_ind = last(which(rocr_pred@cutoffs[[1]] >= cutoff))
  
  if (repar) par(mfrow=c(2,2))
  
  pp = function(perf)  {
    if (length(cutoff_ind)>0 && is.finite(cutoff_ind)) {
      x=perf@x.values[[1]][cutoff_ind]
      y=perf@y.values[[1]][cutoff_ind]
      points(x,y, col='blue')
      linelight(x,y, col='lightblue')
    }
  }
  plot(acc); pp(acc)
  plot(f1); pp(f1)
  plot(roc); pp(roc)
  abline(a=1,b=-1,lty='dashed',col='gray')
  legend('bottomleft',legend=sprintf("AUC = %.3f",auc))
  plot(rp_curve); pp(rp_curve)
  pp = function(ind,...) points(rp_curve@x.values[[1]][ind], rp_curve@y.values[[1]][ind], ...)
  best_f1 = which.max(f1@y.values[[1]])
  pp(best_f1, pch=2,col='green')
  f05 = performance(rocr_pred,'f',beta=0.5)
  best_f05 = which.max(f05@y.values[[1]])
  pp(best_f05,pch=2,col='green')
  f2 = performance(rocr_pred,'f',beta=2)
  best_f2 = which.max(f2@y.values[[1]])
  pp(best_f2,pch=2,col='green')
  
  prbe = performance(rocr_pred,'prbe')@y.values[[1]]
  linelight(prbe,prbe,col='lightgray')
  
  # printf("Acc = %.3f\n", mean((pred >= cutoff) == (labels > 0)))
  printf("Acc %.3f, ", acc@y.values[[1]][cutoff_ind])
  
  printf("  F %.3f, Prec %.3f, Rec %.3f, Spec %.3f",
         f1@y.values[[1]][cutoff_ind],
         pr_curve@y.values[[1]][cutoff_ind],
         pr_curve@x.values[[1]][cutoff_ind],
         roc@x.values[[1]][cutoff_ind])
  # printf(" Prec = %.3f\n", pr_curve@y.values[[1]][cutoff_ind])
  # printf("  Rec = %.3f\n", pr_curve@x.values[[1]][cutoff_ind])
  # printf(" Spec = %.3f\n", roc@x.values[[1]][cutoff_ind])
  
  if (bac[1] != -1)
    printf(", BalAcc %.3f", mean(bac))
  printf("\n")
  
  
  invisible(rocr_pred)
  
  if (repar) par(mfrow=c(1,1))
}



evaluate_predictions = function(g, g_obs, types=c("1/dist", 
                                               "common neighbors", 
                                               "normalized common neighbors",
                                               "degree product"), metric="auc", plot_roc= FALSE){
  
  require(ROCR)
  
  A_g <- get.adjacency(g)
  A_g_v <- A_g[lower.tri(A_g)]
  A_obs = get.adjacency(g_obs)
  A_obs_v <- A_obs[lower.tri(A_obs)]
  true_possibleedges_no_obs = A_g_v[A_obs_v==0]
  true_possibleedges_no_obs = factor(true_possibleedges_no_obs, levels = c(0,1)) #debe tener ceros y unos
  
  result=data.frame(true = true_possibleedges_no_obs)
  metrics = numeric()
  for (t in types){
    s   = score(g_obs, type=t, sort = FALSE, scale=FALSE)
    pred_possibleedges_no_obs = s$score[A_obs_v==0]
    result[,t] = pred_possibleedges_no_obs
    
    pred <- prediction(pred_possibleedges_no_obs, true_possibleedges_no_obs)
    
    if (metric=="auc"){
      #rocr_pred = prediction(pred,labels)
      #auc = performance(rocr_pred,'auc')@y.values[[1]]
      metrics[t] = tryCatch({
        performance(pred, "auc")@y.values[[1]]
      }, error = function(e) {
        NaN
      })
      # if (length(unique(pred_possibleedges_no_obs))>3){
      #   metrics[t] = performance(pred, "auc")@y.values[[1]]
      # }else{
      #   metrics[t] = NaN
      # }

    } else if (metric=="f"){
      cutoff_ind = tail(which(pred@cutoffs[[1]] >= 0.5), n=1)
      perf = performance(pred,'f')
      metrics[t] = perf@y.values[[1]][cutoff_ind]
    }

  }
  result[,c("edge","from","to")] = s[,c("edge","from","to")]
  
  
  if (plot_roc){
    pred <- prediction(result[,types[1]], true_possibleedges_no_obs)
    roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
    plot(roc.perf, col="red")
    abline(a=0, b= 1)
    for (t in types[2:length(types)]){
      pred <- prediction(result[,t], true_possibleedges_no_obs)
      roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
      lines(roc.perf@x.values[[1]], roc.perf@y.values[[1]])
    }
  }

  return(list(result=result, metrics=metrics))
}

summary_predictions = function(g,
                               n_sample=10, 
                               p_v=seq(from=0.01,to=0.5,by=0.05), 
                               types=c("1/dist","common neighbors", "normalized common neighbors", "degree product"),
                               metric="f",
                               cores=1){
  
  require(reshape2)
  
  if (cores==1){
    results = data.frame()
    #undebug(evaluate_predictions)
    for (p in p_v){
      for (n in 1:n_sample){
        printf("%s %s %s",Sys.time(), p, n)
        g_obs = delete_edges_rand(g,p=p)
        r = evaluate_predictions(g, g_obs, types=types, metric=metric)
        r$metric["p"] = p
        r$metric["n"] = n
        results = rbind(results, r$metric)
      }
    }
    names(results) = c(types,"p","n")
    
  }else{
    require(foreach)
    require(doMC)
    registerDoMC(cores=cores) 
    results = foreach(n = 1:n_sample, .combine=rbind) %dopar% {
      results_in = data.frame()
      for (p in p_v){
        printf("%s %s %s",Sys.time(), p, n)
        g_obs = delete_edges_rand(g,p=p)
        r = evaluate_predictions(g, g_obs, types=types, metric=metric)
        r$metric["p"] = p
        r$metric["n"] = n
        results_in = rbind(results_in, r$metric)
      }
      names(results_in) = c(types,c("p","n"))
      return(results_in)
    }
  }
  
  results_summary = aggregate(. ~ p, results[,setdiff(names(results),"n")], FUN="mean")
  results_summary_long = melt(results_summary, id.vars=c("p"))
  names(results_summary_long) = c("p","score",metric)
  
  results_long = melt(results, id.vars=c("p","n"))
  names(results_long) = c("p","n","score",metric)

  return(list(results=results_long, results_summary=results_summary_long))
}




