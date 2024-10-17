



rawEnrichmentAnalysis <- function(expr, signatures, genes, parallel.sz = 4, parallel.type = 'SOCK') {

  # Reduce the expression dataset to contain only the required genes
  shared.genes <- intersect(rownames(expr), genes)
  print(paste("Num. of genes:", length(shared.genes)))
  expr <- expr[shared.genes, ]
  if (dim(expr)[1] < 5000) {
    print(paste("ERROR: not enough genes"))
    return(-1)
  }

  # Transform the expression to rank
  expr <- apply(expr, 2, rank)

  # Run ssGSEA analysis for the ranked gene expression dataset
  if(packageVersion("GSVA") >= "1.50.0") {
    exprs <- GSVA::ssgseaParam(expr, signatures,
                         normalize = FALSE)
    scores <- GSVA::gsva(exprs)

  } else {
    exprs <- GSVA::ssgseaParam(expr, signatures,
                         normalize = FALSE)
    scores <- GSVA::gsva(exprs)
  }


  scores = scores - apply(scores,1,min)

  # Combine signatures for same cell types
  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]
  agg <- aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]

  return(scores)
}

transformScores <- function(scores, fit.vals, scale=TRUE) {
  rows <- rownames(scores)[rownames(scores) %in% rownames(fit.vals)]
  tscores <- scores[rows, ]
  minX <- apply(tscores, 1, min)
  A <- rownames(tscores)
  tscores <- (as.matrix(tscores) - minX)/5000
  tscores[tscores < 0] <- 0
  if (scale==FALSE) {
    fit.vals[A,3] = 1
  }

  tscores <- (tscores^fit.vals[A,2])/(fit.vals[A,3]*2)

  return(tscores)
}

spillOver <- function(transformedScores, K, alpha = 0.5) {
  K <- K * alpha
  diag(K) <- 1
  rows <- rownames(transformedScores)[rownames(transformedScores) %in%
                                        rownames(K)]
  scores <- apply(transformedScores[rows, ], 2, function(x) pracma::lsqlincon(K[rows,rows],
                                                                              x, lb = 0))

  scores[scores<0]=0
  rownames(scores) <- rows
  return(scores)
}


microenvironmentScores <- function(adjustedScores) {
  ImmuneScore = apply(adjustedScores[c('B-cells','CD4+ T-cells','CD8+ T-cells','DC','Eosinophils','Macrophages','Monocytes','Mast cells','Neutrophils','NK cells'),],2,sum)/1.5
  StromaScore = apply(adjustedScores[c('Adipocytes','Endothelial cells','Fibroblasts'),],2,sum)/2
  MicroenvironmentScore = ImmuneScore+StromaScore
  adjustedScores = rbind(adjustedScores,ImmuneScore,StromaScore,MicroenvironmentScore)
}

xCellSignifcanceBetaDist = function(scores,beta_params=NULL,rnaseq=T) {
  if (is.null(beta_params)) {
    if (rnaseq==T) {
      beta_params = xCell.data$spill$beta_params
    } else {
      beta_params = xCell.data$spill.array$beta_params
    }
  }
  scores = scores[rownames(scores) %in% colnames(xCell.data$spill$beta_params[[1]]),]
  pvals = matrix(0,nrow(scores),ncol(scores))
  rownames(pvals) = rownames(scores)
  eps = 1e-3

  for (i in 1:nrow(scores)) {
    ct = rownames(scores)[i]
    beta_dist = c()
    for (bp in beta_params) {
      if (sum(bp[,i]==0)) {
        bd = matrix(eps,1,100000)
      } else {
        bd = stats::rbeta(100000,bp[1,ct],bp[2,ct])
        bd = ((1+eps)*(bp[3,ct]))*bd
      }
      beta_dist = c(beta_dist,bd)
    }
    pvals[i,] = 1-mapply(scores[i,],FUN=function(x) mean(x>beta_dist))
  }


  pvals
}

xCellSignifcanceRandomMatrix = function(scores,expr,spill,alpha=0.5,nperm=250) {

  shuff_expr = mapply(seq(1:nperm),FUN=function(x) sample(nrow(expr),nrow(expr)))

  rownames(shuff_expr) = sample(rownames(expr))
  shuff_xcell = xCellAnalysis(shuff_expr,spill=spill,alpha=alpha)

  shuff_xcell = shuff_xcell[rownames(scores),]

  pvals = matrix(0,nrow(scores),ncol(scores))
  beta_dist = matrix(0,nrow(scores),100000)
  eps = 1e-3
  for (i in 1:nrow(scores)) {
    x = shuff_xcell[i,]
    if (stats::sd(x)<eps) {
      beta_dist[i,] = rep(eps,100000)
    } else {
      x = x+eps
      beta_params=MASS::fitdistr(x/((1+2*eps)*(max(x)))+eps,"beta",list(shape1=1,shape2=1),lower=eps)
      beta_dist[i,] = stats::rbeta(100000,beta_params$estimate[1],beta_params$estimate[2])
      beta_dist[i,] = ((1+2*eps)*(max(x)))*beta_dist[i,]
    }
    #sm.density.compare(c(shuff_xcell[i,],beta_dist),factor(c(rep(1,100),rep(2,100000))),xlim=c(0,max(beta_dist)));title(rownames(scores)[i])
    pvals[i,] = 1-unlist(lapply(scores[i,],FUN=function(x) mean(x>beta_dist[i,])))
  }
  rownames(pvals) = rownames(scores)
  colnames(pvals) = colnames(scores)
  rownames(shuff_xcell) = rownames(scores)
  rownames(beta_dist) = rownames(scores)


  list(pvals=pvals,shuff_xcell=shuff_xcell,shuff_expr=shuff_expr,beta_dist=beta_dist)
}

xCellAnalysis <- function(expr, signatures=NULL, genes=NULL, spill=NULL, rnaseq=TRUE, scale=TRUE,
                          alpha = 0.5, parallel.sz = 4, parallel.type = 'SOCK',
                          cell.types.use = NULL) {
  if (is.null(signatures))
    signatures = xCell.data$signatures
  if (is.null(genes))
    genes = xCell.data$genes
  if (is.null(spill)) {
    if (rnaseq==TRUE) {
      spill = xCell.data$spill
    } else {
      spill = xCell.data$spill.array
    }
  }

  # Caulcate average ssGSEA scores for cell types

  if (!is.null(cell.types.use)) {
    A = intersect(cell.types.use,rownames(spill$K))
    if (length(A)<length(cell.types.use)) {
      return ('ERROR - not all cell types listed are available')
    }
  }

  scores <- rawEnrichmentAnalysis(expr,signatures,genes, parallel.sz = parallel.sz, parallel.type = 'SOCK')

  # Transform scores from raw to percentages
  scores.transformed <- transformScores(scores, spill$fv, scale)

  # Adjust scores using the spill over compensation matrix

  if (is.null(cell.types.use)) {
    scores.adjusted <- spillOver(scores.transformed, spill$K, alpha )
    scores.adjusted = microenvironmentScores(scores.adjusted)
  } else {
    scores.adjusted <- spillOver(scores.transformed[cell.types.use,], spill$K, alpha )
  }
  return(scores.adjusted)
}


CoreAlg <- function(X, y){


  ########################
  ## X is the data set
  ## y is labels for each row in X
  ########################



  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}

    #if(i==1){nus <- 0.997}
    #if(i==2){nus <- 0.998}
    #if(i==3){nus <- 0.999}

    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  #Execute In a parallel way the SVM
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <-  mclapply(1:svn_itor, res, mc.cores=svn_itor) #lapply(1:svn_itor, res) #

  #Initiate two variables with 0
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  ##############################
  ## Here CIBERSORT starts    #
  ##############################

  t <- 1
  while(t <= svn_itor) {

    #Get the weights with a matrix multiplications between two vectors. I should get just one number (?)
    #This is done multiplying the coefficients (?) and ???

    #The support vectors
    #are the points of my dataset that lie closely to the plane that separates categories
    #The problem now is that I don't have any category (discrete variable, e.g., "sport", "cinema") but I ave continuous variable
    mySupportVectors <- out[[t]]$SV

    #My coefficients
    myCoefficients <- out[[t]]$coefs

    weights = t(myCoefficients) %*% mySupportVectors

    #set up weight/relevance on each
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)

    #This multiplies the reference profile for the correspondent weigth
    u <- sweep(X,MARGIN=2,w,'*')

    #This does the row sums
    k <- apply(u, 1, sum)

    #Don't know
    nusvm[t] <- sqrt((mean((k - y)^2))) #pitagora theorem
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  print(mn)
  model <- out[[mn]]

  #get and normalize coefficients

  #############################################
  ## THIS IS THE SECRET OF CIBERSORT
  #############################################

  q <- t(model$coefs) %*% model$SV

  #############################################
  #############################################

  q[which(q<0)]<-0

  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}


doPerm <- function(perm, X, Y){


  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #print(itor)

    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){


  #read in data
  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  #Y <- read.table(mixture_file, header=T, sep="\t", row.names=1,check.names=F)
  Y <- mixture_file
  X <- data.matrix(X)
  Y <- data.matrix(Y)

  #order

  ###################################
  ## This is needed to make the two tables consistent in gene
  ###################################

  X <- X[order(rownames(X)),,drop=F]
  Y <- Y[order(rownames(Y)),,drop=F]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file

  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,,drop=F]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,,drop=F]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  # stefano write X and Y
  Y_norm <- apply(Y, 2, function(mc) (mc - mean(mc)) / sd(mc)  )

  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mix <- dim(Y)[2]
  pval <- 9999

  #iterate through mix
  while(itor <= mix){

    ##################################
    ## Analyze the first mixed sample
    ##################################

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1

  }


  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1, drop=F]
  obj <- obj[-1,, drop=F]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")

  list(proportions=obj, mix = Y_norm, signatures = X)
}






#' @title exp2cell
#' @description Function `exp2cell` use gene expression profiles to quantify cell abundance matrix.
#' `exp2cell` provides three methods for estimating the relative infiltration abundance of different cell types in the tumor microenvironment (TME),
#' which including xCell, ssGSEA estimated method proposed by Şenbabaoğlu et al. and CIBERSORT.
#' @param exp The gene expression data set. A matrix with row names as symbols and columns as samples.
#' Gene expression profiles were used to quantify cell abundance matrix.
#' @param method Method must be one of "xCell", "ssGSEA" and "CIBERSORT".
#' @param QN Quantile normalization of input mixture (default = TRUE)
#' @param perm No. permutations; set to >=100 to calculate p-values (default = 100)
#' @importFrom  GSVA ssgseaParam
#' @importFrom  GSVA gsva
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom e1071 svm
#' @importFrom grDevices colorRampPalette
#' @importFrom parallel mclapply
#' @importFrom MASS fitdistr
#' @importFrom pracma lsqlincon
#' @importFrom stats cor
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom utils read.table
#' @importFrom utils packageVersion
#' @importFrom stats aggregate
#' @return Cell abundance matrix.
#' @references 1. Aaron, M, Newman, et al. Robust enumeration of cell subsets from tissue expression profiles.[J]. Nature Methods, 2015.
#' 2. Aran D , Hu Z , Butte A J . xCell: digitally portraying the tissue cellular heterogeneity landscape[J]. Genome Biology, 2017, 18(1):220.
#' 3. Şenbabaoğlu, Yasin, Gejman R S , Winer A G , et al. Tumor immune microenvironment characterization in clear cell renal cell carcinoma identifies prognostic and immunotherapeutically relevant messenger RNA signatures[J]. Genome biology, 2016, 17(1).
#' @export
#' @examples
#' #get breast cancer gene expression profile.
#' exp.example<-GetExampleData("exp.example")
#'
#' #perform the exp2cell method. Method must be one of "xCell","ssGSEA" and "CIBERSORT".
#' cellmatrix<-exp2cell(exp=exp.example,method="ssGSEA") #cell abundance matrix
exp2cell <- function(exp,method="xCell",perm=100,QN=TRUE) {
  if (method=="xCell") {
    ## requireNamespace("xCell")|| stop("package xCell is required,please install package xCell")
    cellmatrix<-xCellAnalysis(exp)
    cellmatrix<-cellmatrix[1:64,]
  }else if (method=="ssGSEA") {
    ## requireNamespace("GSVA")|| stop("package GSVA is required,please install package GSVA")
    exp <- exp[apply(exp, 1, var) > 0, ]
    expr<-ssgseaParam(as.matrix(exp),immunelist,minSize=2)
    cellmatrix<-gsva(expr)
  }else if (method=="CIBERSORT"){
    lm22path<-system.file('extdata', 'LM22.txt', package = 'SMDIC')
    cellmatrix_pre<-CIBERSORT(lm22path,exp,perm = perm,QN = QN)
    cellmatrix_T<-cellmatrix_pre$proportions
    cellmatrix_T<-cellmatrix_T[,1:22]
    cellmatrix<-t(cellmatrix_T)
  } else {
    print("method must be one of xCell, ssGSEA, or CIBERSORT")
  }
  return(cellmatrix)
}
