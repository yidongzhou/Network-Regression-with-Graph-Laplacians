library(Matrix)
library(matrixcalc)
library(igraph)
library(shapes)
library(gdata)
#library(OpenMx)
#library(CVXR)
#library(ggraph)
#library(rtweet)
library(osqp)
#library(abind)
#library(network)
#library(MASS)
library(GGally)
library(sna)
#####
######

igraph_from_lapl <-function(lapl, wght=TRUE){
  if (wght==FALSE){
    adj <- - as.matrix(lapl)
    diag(adj) <-0
    adj[adj<0]<-0
    grph <- graph_from_adjacency_matrix(as.matrix(adj))
    return ( as.undirected(grph))
  }
  if (wght==TRUE){
    adj <- - as.matrix(lapl)
    diag(adj) <-0 #create adjacency
    adj[adj<0]<-0
    grph <- graph.adjacency(adj, mode="undirected", weighted=TRUE) #igraph from adjacency
    E(grph)$width <- E(grph)$weight 
    return ( grph)
  }
}
#Description: convert laplacian to igraph
######


######
is.laplacian <- function(L, tol = 1e-08) {
  L <- as.matrix(L)
  if (!matrixcalc::is.square.matrix(L)) {
    print("argument L is not a square matrix")
    return(FALSE)
  }

  if (!matrixcalc::is.symmetric.matrix(L)) {
    print("argument L is not a symmetric matrix")
    return(FALSE)
  }
  
  if (!is.numeric(L)) {
    print("argument L is not a numeric matrix")
    return(FALSE)
  }
  upper_diag <- L[upper.tri(L)]
  row_sums <- rowSums(L)
  if (any(row_sums < -tol)) {
    # print("Rows do no sum to zero")
    return(FALSE)
  }
  if (any(row_sums > tol)) {
    # print("Rows do no sum to zero")
    return(FALSE)
  }
  if (any(upper_diag > tol)) {
    # print("Off diagonal elements are not all non-positive")
    return(FALSE)
  }
  return(TRUE)
}
#Description: tests for graph Laplacian objects
#######


#######
ggraph_plot <- function(GL, label_nodes=0, quantile_edges=0, delete_edges=0, tol=-10^-10, colour_node='slateblue', colour_line='springgreen')#delete nodes less then degree of label nodes,
  #set anything belew quantile edges to this quantile, delete edges less then delete edges quantile
{GL[GL> tol]<-0
graph<-igraph_from_lapl(GL, TRUE)
igraph::V(graph)$node_label <- unname(ifelse(strength(graph)[igraph::V(graph)] >label_nodes, names(igraph::V(graph)), "")) 
igraph::V(graph)$node_size <- unname(ifelse(strength(graph)[igraph::V(graph)] > label_nodes, strength(graph), 0))

graph<-delete_edges(graph,E(graph)[E(graph)$weight[E(graph)]< quantile(E(graph)$weight,delete_edges)]) 
E(graph)$test<-E(graph)$weight
E(graph)$test[E(graph)$test<quantile(E(graph)$test,quantile_edges)]<-quantile(E(graph)$test,quantile_edges)

ggraph(graph, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(edge_width=0.125, aes(alpha=(E(graph)$test)^1)) +
  geom_node_label(aes(label=node_label, size=node_size),
                  label.size=0, fill="#ffffff66", segment.colour=colour_line,
                  color=colour_node, repel=TRUE,  fontface="bold") +
  coord_fixed() +
  scale_size_area(trans="sqrt") +
  labs(title="", subtitle="") +
  theme_graph() +
  theme(legend.position="none")}

#Description :plot laplcians
#######





########


proj_osqp <- function(x0, P, A, l, u, diag_index, nodiag_index) {
  #problem specific values
  dim_x <- dim(x0)[1]
  x0_vec <-
    gdata::upperTriangle(as.matrix(x0) / sqrt(2), diag = TRUE, byrow = TRUE)
  x0_vec[diag_index] <- x0_vec[diag_index] / sqrt(2)
  q <- -2 * x0_vec
  
  # print("QP set up!")
  # print(date())
  out <-
    osqp::solve_osqp(P, q, A, l, u, pars = osqp::osqpSettings(eps_abs = 1E-09, eps_rel =
                                                    1E-09, verbose = FALSE))
  theta <- out$x
  theta[nodiag_index] <- theta[nodiag_index] / sqrt(2)
  approx <- matrix(0, dim_x, dim_x)
  gdata::upperTriangle(approx, diag = TRUE, byrow = TRUE) <- theta
  approx <- Matrix::forceSymmetric(approx)
  approx
}



proj_sparse <- function(x0){
  if (is.laplacian(as.matrix(x0)) == TRUE) {
    return(x0)
  }
  M_dim <-dim(x0)[1]
  if (is.laplacian(as.matrix(x0)) == FALSE && M_dim==1000){
    P.save <- readRDS("Projection constraints/P_save_1000.rds")
    A.save <- readRDS("Projection constraints/A_save_1000.rds")
    l.save <- readRDS("Projection constraints/l_save_1000.rds")
    u.save <- readRDS("Projection constraints/u_save_1000.rds")
    diag_index.save <- readRDS("Projection constraints/diag_index_save_1000.rds")
    nodiag_index.save <- readRDS( "Projection constraints/nodiag_index_save_1000.rds")
    return(proj_osqp(x0,P.save,A.save,l.save,u.save,diag_index.save,nodiag_index.save))
  }
  else{
    dim_x<-dim(x0)[1]
    p<-dim_x
    ### This part only needs doing once for dimension M x M
    #sets up the constraints
    #indexes for diagonal and off diagonal
    index_mat <- matrix(0, dim_x, dim_x)
    gdata::upperTriangle(index_mat, diag = TRUE, byrow = TRUE) <-
      seq(1, dim_x * (dim_x + 1) / 2, 1)
    diag_index <- diag(index_mat)
    upper_index <-
      setdiff(seq(1, dim_x * (dim_x + 1) / 2, 1), diag(index_mat))
    nodiag_index <- index_mat + t(index_mat)
    diag(nodiag_index) <- 0
    #####################################
    P<-Matrix::.sparseDiagonal(p+p*(p-1)/2)
    u<-rep(0,p+p*(p-1)/2)
    l<-u
    A<-P
    sq<-1/sqrt(2)
    #print("setting up constraints - only needs doing once - can it be sped up?")
    l[nodiag_index[1:p,]]<- -Inf
    #this loop is slow - how to speed up
    for (i in 1:p){
      A[diag_index[i],nodiag_index[i,]]<- sq
      if ((i/100)==trunc(i/100)) {
        print(i)
      }
    }
    return(proj_osqp(x0,P,A,l,u,diag_index,nodiag_index))
  }
}

#Description: Projection of a square matrix into the graph Laplacian space
#########

########

estSS<-function (S, weights = 1){
  M <- dim(S)[3]
  k <- dim(S)[1]
  H <- defh(k)
  if (length(weights) == 1) {
    weights <- rep(1, times = M)
  }
  Q <- array(0, c(k + 1, k, M))
  for (j in 1:M) {
    Q[, , j] <- t(H) %*% (shapes::rootmat(S[, , j]))
  }
  ans <- fgpa.rot(Q, tol1=1e-05, tol2=1e-05, proc.output = FALSE, 
                  reflect = TRUE)   
  out<-t(H %*% ans$mshape) %*% (H %*% ans$mshape) #changed order of transpose so rotation has affect
  (out+t(out))/2
}
#Description: Calculating the Procrustes mean of an array of graph Laplacians (From Ian Dryden)
#########


#######
Mean_GL <- function(Array,
                    euc = TRUE,
                    sqrt = FALSE,
                    proc = FALSE, project=FALSE)
{
  if (is.list(Array) == TRUE) {
    rownms <- rownames(Array[[1]])
    AArray <-
      array(0, dim = c(dim(Array[[1]])[1], dim(Array[[1]])[2], length(Array)))
    for (i in 1:length(Array))
      AArray[, , i] <- as.matrix(Array[[i]])
    Array <- AArray
    rownames(Array) <-rownms
  }
  n <- dim(Array)[3]
  k <- dim(Array)[2]
  rownms <- rownames(Array[, , 1])
  if (euc == TRUE) {
    a <- Array[,,1]
    for ( i in 2:n)
      a<- a+Array[,,i]
    a<- a/n
    h <- Array[, , 1]
    h[] <- a
    rownames(h) <- rownms
    colnames(h) <- rownms
    return(h)
  }
  if (sqrt == TRUE) {
    #square root euclidean mean
    a <- array(rep(0, k * k * n), dim = c(k, k, n))
    for (i in 1:n)
      a[, , i] <-  shapes::rootmat(Array[, , i])
    b <- a[,,1]
    for (i in 2:n)
      b<- b+a[,,i]
    b<-b/n
    b <- b %*% t(b)
    h <- Array[, , 1]
    h[] <- b
    rownames(b) <- rownms
    colnames(b) <- rownms
    if (project==FALSE){
      return(b)
    }
    if(project==TRUE){
      b <-proj_sparse(b)
      rownames(b) <- rownms
      colnames(b) <- rownms
      return(b)
    }
  }
  if (proc == TRUE) {
    #size-and-shape mean
    c <-estSS(Array)
    h <- Array[, , 1]
    h[] <- c
    rownames(c) <- rownms
    colnames(c) <- rownms
    if (project==FALSE){
      return(c)
    }
    if(project==TRUE){
      c <- proj_sparse(c)
      rownames(c) <- rownms
      colnames(c) <- rownms
      return(c)
    }
  }
}
#Description: Calculating the mean of an array of graph Laplacians, option for different metrics
#and if the mean should be projected into the graph Laplacian space
########

#######
G_2 <- function(mat){
  B<- (mat+t(mat))/2
  H <- shapes::rootmat(t(B)%*%B)
  return((B+H)/2)
}
#Description: function to apply before projection for square root
#######


#######
Interpolation <-
  function(GL1,
           GL2,
           c,
           euc = TRUE,
           sqrt = FALSE,
           proc = FALSE,
           proj_opt = TRUE){
    if ((euc == TRUE) &&
        (sqrt == FALSE) && (proc == FALSE) && (proj_opt == TRUE))
    {
      diff_euc <- GL2 - GL1
      un_proj <- GL1 + c * diff_euc
      int <- proj_sparse(un_proj)
      dist_ofproj <- norm(int - un_proj, type = 'f')
      return(list(
        int = as.matrix(int) ,
        un_proj = un_proj ,
        diff = diff_euc,
        dist = dist_ofproj
      ))
    }
    
    if ((euc == FALSE) &&
        (sqrt == TRUE) && (proc == FALSE) && (proj_opt == TRUE))
    {
      sqrt1 <- shapes::rootmat(GL1)
      sqrt2 <- shapes::rootmat(GL2)
      diff_sqrt <- sqrt2 - sqrt1
      un_proj <- G_2(sqrt1 + c * diff_sqrt)%*%t(G_2(sqrt1 + c * diff_sqrt))
      int <- proj_sparse(un_proj)
      dist_ofunproj <- norm(shapes::rootmat(int) - shapes::rootmat(un_proj), type = 'f')
      return(list(
        int = as.matrix(int) ,
        un_proj = un_proj
        ,
        complete_un_proj = sqrt1 + c * diff_sqrt,
        diff = ((diff_sqrt) %*% t(diff_sqrt)),
        diff_sq = diff_sqrt,
        dist = dist_ofunproj
      ))
    }
    
    if ((euc == FALSE) &&
        (sqrt == FALSE) && (proc == TRUE) && (proj_opt == TRUE))
    {   m <- dim(GL1)[1]
    shape1 <- shapes::rootmat(GL1)
    shape2 <- shapes::rootmat(GL2)
    proc12 <- procOPA(shape1, shape2, scale = FALSE, reflect = TRUE)$Bhat
    sqrt_int <- procOPA(shape1, c*proc12+(1-c)*shape1, scale = FALSE, reflect = TRUE)$Bhat
    diff_proc <- proc12 - shape1
    un_proj <-
      ( t(sqrt_int)%*%sqrt_int)
    int <- proj_sparse(un_proj)
    dist_ofunproj <-
      procdist(
        shapes::rootmat(un_proj),
        shapes::rootmat(int),
        type = 'sizeandshape',
        reflect = TRUE
      )
    return (list(
      int = as.matrix(int) ,
      un_proj = un_proj,
      diff = (((diff_proc)) %*% t((diff_proc))) ,
      dist = dist_ofunproj
    ))
    }
    
    
    if ((euc == FALSE) &&
        (sqrt == FALSE) && (proc == old_way) && (proj_opt == TRUE))
    {
      m <- dim(GL1)[1]
      shape1 <- shapes::rootmat(GL1)
      shape2 <- shapes::rootmat(GL2)
      proc12 <- procOPA(shape1, shape2, scale = FALSE, reflect = TRUE)$Bhat
      sqrt_int <- procOPA(shape1, c*proc12+(1-c)*shape1, scale = FALSE, reflect = TRUE)$Bhat
      diff_proc <- proc12 - shape1
      un_proj <-
        ( sqrt_int%*%t(sqrt_int)) #this is what I used to do, new way has transpose times not transpose
      int <- proj_sparse(un_proj)
      dist_ofunproj <-
        procdist(
          shapes::rootmat(un_proj),
          shapes::rootmat(int),
          type = 'sizeandshape',
          reflect = TRUE
        )
      return (list(
        int = as.matrix(int) ,
        un_proj = un_proj,
        diff = (((diff_proc)) %*% t((diff_proc))) ,
        dist = dist_ofunproj
      ))
    }
    else
      return(0)
  }

########



########
estSS3<-function (S) 
{
  M <- dim(S)[3]
  k <- dim(S)[1]
  Q <- array(0, c(k, k, M))
  for (j in 1:M) {
    Q[, , j] <-  (shapes::rootmat(S[, , j]))
  }
  ans <- fgpa.rot(Q, tol1=1e-05, tol2=1e-05, proc.output = TRUE, 
                  reflect = TRUE)   
  out<-ans$mshape %*% t(ans$mshape)
  ans2<-list(mshape=0,tan=0)
  ans2$mshape<-(out+t(out))/2
  ans2$tan<-matrix(0,k*k,M)
  for (i in 1:M){
    ans2$tan[,i]<-c(t(ans$r.s.r[,,i]-ans$mshape))
  }
  ans2
}
#Description : procrustes tangent coordinate
########


########

PCA_GL <- function(gl_list,
                   euc = TRUE,
                   sqrt = FALSE,
                   proc = FALSE)
{
  amount_words <- dim(gl_list[[1]])[1]
  n <- length(gl_list)
  
  if ((euc == TRUE) && (sqrt == FALSE) && (proc == FALSE))
  {
    eucmean <- (Mean_GL(
      gl_list,
      euc = TRUE,
      sqrt = FALSE,
      proc = FALSE
    ))
    vectorized_euc <- matrix(0, amount_words * amount_words , n)
    for (i in 1:n)
    {
      vectorized_euc[, i] <-as.vector(gl_list[[i]])
    }
    pca <- prcomp(t(vectorized_euc), center = TRUE)
    PC1_euc <- matrix(pca$rotation[, 1], amount_words, amount_words)
    return(list(
      pca = pca,
      PC1 = PC1_euc,
      eigenvl_S = pca$sdev ^ 2
    ))
  }
  if ((euc == FALSE) && (sqrt == TRUE) && (proc == FALSE))
  {
    sqrtmean <- (Mean_GL(
      gl_list,
      euc = FALSE,
      sqrt = TRUE,
      proc = FALSE
    ))
    sqrted_sqrtmean <- shapes::rootmat(sqrtmean)
    sqrted_gl <- gl_list
    for (i in 1:n) {
      sqrted_gl[[i]] <- shapes::rootmat(gl_list[[i]])
    }
    vectorized_sqrt <- matrix(0, amount_words * amount_words , n)
    for (i in 1:n)
    {
      for (j in 1:amount_words)
        vectorized_sqrt[((j - 1) * amount_words + 1):(j * amount_words), i] <-
          sqrted_gl[[i]][, j] - sqrted_sqrtmean[, j]
    }
    pca <- prcomp(t(vectorized_sqrt), center = FALSE)
    PC1_sqrt <- matrix(pca$rotation[, 1], amount_words, amount_words)
    return(list(
      pca = pca,
      PC1 = PC1_sqrt,
      eigenvl_S = pca$sdev ^ 2
    ))
  }
  if ((euc == FALSE) && (sqrt == FALSE) && (proc == TRUE))
  {
    proc_array <- array(0, dim = c(amount_words, amount_words, n))
    for (i in 1:n)
      proc_array[, , i] <- as.matrix(gl_list[[i]])
    proc_tan <- estSS3(proc_array)$tan
    vectorized_proc <- matrix(0, amount_words * amount_words , n)
    for (i in 1:n)
      vectorized_proc[, i] <- as.vector(proc_tan[, i])
    pca <- prcomp(t(vectorized_proc), center = FALSE)
    PC1_proc <- matrix(pca$rotation[, 1], amount_words, amount_words)
    return(list(
      pca = pca,
      PC1 = PC1_proc,
      eigenvl_S = pca$sdev ^ 2
    ))
  }
}  
#Description: PCA

#########
########


distance_between_gl <- function(GL_1,
                                GL_2,
                                euc = TRUE,
                                sqrt = FALSE,
                                proc = FALSE) {
  if ((euc == TRUE) && (sqrt == FALSE) && (proc == FALSE)) {
    return(norm(GL_1 - GL_2, type = 'f'))
  }
  
  if ((euc == FALSE) && (sqrt == TRUE) && (proc == FALSE)) {
    return(norm(shapes::rootmat(GL_1) - shapes::rootmat(GL_2), type = 'f'))
  }
  
  if ((euc == FALSE) && (sqrt == FALSE) && (proc == TRUE)) {
    return(procdist(
      shapes::rootmat(GL_1),
      shapes::rootmat(GL_2),
      type = 'sizeandshape',
      reflect = TRUE
    ))
  }
  
  else
    print('Error, one and only one metric must be TRUE')
}
#Description: Calculate distances between graph Laplacians for the Eulcidean (alpha= 0.5, 1) 
#and Procrustes metric (alpha=0.5)
#######




#######
dist_sq_means <-
  function(GL1_array,
           GL2_array,
           euc = TRUE,
           sqrt = FALSE,
           proc = FALSE) {
    if ((euc == TRUE) && (sqrt == FALSE) && (proc == FALSE)) {
      mean1 <- Mean_GL(GL1_array,
                       euc = TRUE,
                       sqrt = FALSE,
                       proc = FALSE)
      mean2 <- Mean_GL(GL2_array,
                       euc = TRUE,
                       sqrt = FALSE,
                       proc = FALSE)
      return(distance_between_gl(
        mean1,
        mean2,
        euc = TRUE,
        sqrt = FALSE,
        proc = FALSE
      ) ^ 2)
    }
    
    if ((euc == FALSE) && (sqrt == TRUE) && (proc == FALSE)) {
      mean1 <- Mean_GL(GL1_array,
                       euc = FALSE,
                       sqrt = TRUE,
                       proc = FALSE)
      mean2 <- Mean_GL(GL2_array,
                       euc = FALSE,
                       sqrt = TRUE,
                       proc = FALSE)
      return(distance_between_gl(
        mean1,
        mean2,
        euc = FALSE,
        sqrt = TRUE,
        proc = FALSE
      ) ^ 2)
    }
    
    if ((euc == FALSE) && (sqrt == FALSE) && (proc == TRUE)) {
      mean1 <- Mean_GL(GL1_array,
                       euc = FALSE,
                       sqrt = FALSE,
                       proc = TRUE)
      mean2 <- Mean_GL(GL2_array,
                       euc = FALSE,
                       sqrt = FALSE,
                       proc = TRUE)
      return(distance_between_gl(
        mean1,
        mean2,
        euc = FALSE,
        sqrt = FALSE,
        proc = TRUE
      ) ^ 2)
    }
    
    else
      print('Error, one and only one metric must be TRUE')
  }
#Description: Test statistic for certain metrics

