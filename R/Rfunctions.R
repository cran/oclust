

# roxygen2::roxygenise(roclets = "rd")

#'
#' Mixture of Beta Functions
#'
#'MixBetaDens generates the pdf and cdf of a mixture of beta functions, and calculates the area under the graph between two points.
#'
#' @param n The number of observations in the dataset
#' @param p The dimension
#' @param x A vector of x values to evaluate. Default value is seq(0, 15, by=0.01)
#' @param a Lower bound for area evaluation. Default value is 0
#' @param b Upper bound for area evaluation. Default value is 1
#' @param n_g Vector describing the number of observations in each cluster
#' @param var A list of cluster variance matrices
#'
#' @details
#' The domain for this function is not [0,1] as is typical with a beta function. The domain encompasses the shifted log-likelihoods generated in \code{\link{oclust}}.
#'
#' @return
#' MixBetaDens returns a list with
#' \item{ pdf }{The probability density at each x value}
#' \item{ cdf }{The cumulative density at each x value}
#' \item{ area }{The area under the pdf graph between a and b}
#'
#' @importFrom graphics axis box hist par plot text
#' @importFrom stats rbeta dbeta pbeta dist hclust cutree var runif
#' @importFrom utils flush.console menu
#' @importFrom entropy KL.plugin
#' @importFrom MASS ginv
#' @importFrom mclust plot.mclustBIC mclust1Dplot mclust2Dplot coordProj Mclust emControl priorControl
#' @importFrom mixture gpcm
#' @importFrom mvtnorm rmvnorm
#' @importFrom progress progress_bar
#'
#' @export
MixBetaDens<-function(n,p,x =seq(0, 15, by=0.01),a=0,b=1,n_g=n_g,var=var){
  is.int <- function(no){abs(no-round(no))<1e-15}
  if(!is.numeric(n) ||!is.numeric(p)||!is.numeric(x)||!is.numeric(a)||!is.numeric(b)||!is.numeric(n_g)||sum(unlist(lapply(var,is.numeric))==FALSE)) stop("all arguments must be numeric")
  if(!is.int(n)||!is.int(p)|| n<1 |p<1) stop("n and p must be positive integers")
  if(!is.list(var)) stop("var must be a list")
  if(!is.vector(n_g)) stop("n_g must be a vector")
  if(a>b) stop("a must be less than or equal to b")
  if(length(n_g)!=length(var)) stop("There should be the same number of clusters in n_g and var")


  n_g<-c(n_g)
  pi_g<-n_g/n

  #find determinants of variances
  log_det_var<-pi_g
  for (i in 1:length(pi_g)){
    log_det_var[i]<-determinant(var[[i]])$modulus
    if (log_det_var[i]==-Inf | is.na(log_det_var[i])){
      log_det_var[i]=0
      pi_g[i]=0
    }
  }
  #n_g must be greater than p+1
  for (i in 1:length(pi_g)){
    if (n_g[i]<=(p+1)){pi_g[i]=0
    }
  }
  #generate mixture pdf, cdf, area
  h<-length(x)
  pdf<-rep(0,h)
  for (i in which(pi_g!=0)){
    pdf<-pdf+pi_g[i]*2/(n_g[i]-1)^2*(n_g[i])*dbeta(2*(x+log(pi_g[i])-p/2*log(2*pi)-1/2*log_det_var[i])/(n_g[i]-1)^2*(n_g[i]),p/2,(n_g[i]-p-1)/2)}
  cdf<-rep(0,h)
  for (i in which(pi_g!=0)){
    cdf<-cdf+pi_g[i]*pbeta(2*(x+log(pi_g[i])-p/2*log(2*pi)-1/2*log_det_var[i])/(n_g[i]-1)^2*(n_g[i]),p/2,(n_g[i]-p-1)/2)}
  area<-0
  for (i in which(pi_g!=0)){
    area<-area+
      pi_g[i]*pbeta(2*(b+log(pi_g[i])-p/2*log(2*pi)-1/2*log_det_var[i])/(n_g[i]-1)^2*(n_g[i]),p/2,(n_g[i]-p-1)/2)-
      pi_g[i]*pbeta(2*(a+log(pi_g[i])-p/2*log(2*pi)-1/2*log_det_var[i])/(n_g[i]-1)^2*(n_g[i]),p/2,(n_g[i]-p-1)/2)}

  list(pdf=pdf,cdf=cdf,area=area)
}

#'
#' Simulate Outliers
#'
#'simOuts generates uniform outliers in each dimension in (min- 2.range, max+ 2.range)
#'
#' @param data The data in data frame form
#' @param alpha The proportion of outliers to add in terms of the original data size
#' @param seed Set the seed for reproducibility


#' @return
#' simOuts returns a data frame with the generated outliers appended to the original data.
#'
#'@importFrom stats runif
#' @export

simOuts<-function(data,alpha,seed=123){
  if(!is.numeric(alpha) || !is.numeric(seed)) stop("alpha and seed must be numeric")
  if(alpha<0) stop("alpha must be positive")
  nsim<-round(alpha*nrow(data))

  maxs<-apply(data,2,max)
  mins<-apply(data,2,min)
  means<-apply(data,2,mean)
  bounds<-apply(rbind(maxs-means,means-mins),2,max)
  bounds<-rbind(means-2*bounds,means+2*bounds)

  set.seed(seed)
  data<-rbind(data,apply(bounds,2,function(x) runif(nsim,x[1],x[2])))
  return(data)
}

.dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
}
#'
#' Minimum Mahalanobis Distance
#'
#'minMD calculates the Mahalanobis distance to each cluster and returns the Mahalanobis distance to the closest cluster.
#'
#' @param X A matrix or data frame of the data.
#' @param sigs A list of cluster variance matrices
#' @param mus A list of cluster mean vectors
#'
#' @details
#' This function is used to help identify initial gross outliers.
#'
#' @return
#' minMD returns a vector of length n corresponding to the minimum MD for each point.

#'
#' @export

minMD<-function(X,sigs,mus){
  G=length(mus)
  p=ncol(X)
  X=as.matrix(X)
  if(!is.list(mus) || !is.list(sigs) || (length(mus)!=length(sigs))) stop("sigs and mus must be lists of the same length")
  if(sum(unlist(lapply(sigs,dim))!=p) ||sum(unlist(lapply(sigs,is.matrix))==FALSE)) stop("sigs must be pxp matrices")
  if(sum(unlist(lapply(mus,length))!=p)||sum(unlist(lapply(mus,is.vector))==FALSE)) stop("mus must be vectors of length p")
  if(G==1){invsigs=MASS::ginv(sigs[[1]]); MDs=apply(X,1,function(y) t(y-mus[[1]])%*%invsigs%*%(y-mus[[1]])); minMDs=MDs
  }else{
    invsigs<-lapply(sigs,MASS::ginv); MDs<-apply(X,1,function(y) sapply(1:G,function(i) t(y-mus[[i]])%*%invsigs[[i]]%*%(y-mus[[i]])))
    minMDs<-apply(MDs,2,min)}

  return(minMDs)
}

#'
#' Find Initial Gross Outliers
#'
#'findGrossOuts uses DBSCAN to find areas of high density. Mahalanobis distance to the closest area of high density is calculated for each point. With no elbow specified, the Mahalonis distances are plotted. If the elbow is specified, the indices of the gross outliers are returned.
#'
#' @param X A data matrix
#' @param minPts The minimum number of points in each region of high density. Default is 10
#' @param xlim A vector of form c(xmin,xmax) to specify the domain of the plot. Default is NULL, which sets xmax to 10\% of the data size.
#' @param elbow An integer specifying the location of the elbow in the plot of Mahalanobis distances. Default is NULL, which returns the plot. If elbow is specified, no plot is produced and the gross outliers are returned.
#'
#'
#' @details
#' The function plots Mahalanobis distance to the closest centre in decreasing order or returns the indices of the gross outliers. The elbow location of the plot provides a good indication as to where the gross outliers end. Running the function first without an elbow specified will plot the Mahalonobis distances. Running it again with the elbow specified will return the outliers. It is recommended to choose the elbow conservatively. If the MDs decrease smoothly, there are no gross outliers. Set elbow=1.
#'
#' @return
#' findGrossOuts returns a vector with the indices of the gross outliers. One fewer point is returned than the elbow specified.
#'
#' @importFrom graphics plot
#' @importFrom dbscan kNNdist dbscan
#'
#' @export

findGrossOuts<-function(X,minPts=10,xlim=NULL,elbow=NULL){
  if(!is.null(xlim) && length(xlim)!=2) stop("Upper and lower xlimits must be specified, or leave blank for defaults")
  if(!is.numeric(minPts)) stop("minPts must be numeric")
  if(!is.numeric(xlim) &&!is.null(xlim)) stop("xlim must be numeric or NULL")
  if(!is.numeric(elbow) &&!is.null(elbow)) stop("elbow must be numeric or NULL")
  is.int <- function(no){abs(no-round(no))<1e-15}
  if(!is.null(elbow)){if((!is.int(elbow)||elbow<1)) stop("elbow must be a positive integer or NULL")
  }
  if(!is.int(minPts) || minPts<1) stop("minPts must be a positive integer")
  n=nrow(X)
  if(is.null(xlim)){xlim=c(1,0.1*n)}
  xlim=c(floor(xlim[1]),ceiling(xlim[2]))

  kdist<-dbscan::kNNdist(X,(minPts-1))
  kdist<-sort(kdist,decreasing = F)
  n=nrow(X)

  distx<-sapply(1:n,function(x) .dist2d(c(x,kdist[x]),c(1,kdist[1]),c(n,kdist[n])))
  e<-which.max(distx)

  dbres<-dbscan::dbscan(X,eps=e,minPts=minPts)
  Gd=max(dbres$cluster)
  means=list()
  for(i in 1:Gd){
    means[[i]]=colMeans(X[which(dbres$cluster==i),])
  }

  vars=list()
  for(i in 1:Gd){
    vars[[i]]=var(X[which(dbres$cluster==i),])
  }

  MDs<-minMD(X,vars,means)
  sortedMD<-sort(MDs,decreasing = T)
  if(is.null(elbow)){
  plot(x=(xlim[1]:xlim[2]),y=sortedMD[xlim[1]:xlim[2]],type='l',ylab="MD",xlab="Index")}
  else{
    if(elbow==1){grossOuts=NULL}
    else {grossOuts=order(MDs,decreasing = T)[1:(elbow-1)]}
    return(grossOuts)
  }
}

.pval<-function(p,s,n_gs,likes,like0,nLeft,B){
  likes<-likes[which(!is.na(likes))]
  diff<-likes-like0
  orddif<-sort(diff)
  N=length(likes)


  G<-length(n_gs)

  cdf<-MixBetaDens(n=N,p=p,x=orddif,n_g=n_gs,var=s)$cdf

  empcdf<-(1:N)/N
  diffcdf<-cdf-empcdf

  T0<-abs(min(diffcdf))+max(diffcdf)

  Tstat=NULL
  for (k in 1:B){
    set.seed(k)
    theosamp<-NULL
    for (i in 1:G){
      theosamp<-c(theosamp,rbeta(n_gs[i],p/2,(n_gs[i]-p-1)/2)*(n_gs[i]-1)^2/(2*n_gs[i])-
                    log(n_gs[i]/N)+p/2*log(2*pi)+1/2*log(det(s[[i]])))}
    ordsamp<-sort(theosamp)
    cdf<-MixBetaDens(n=N,p=p,x=ordsamp,n_g=n_gs,var=s)$cdf

    empcdf<-(1:N)/N
    diffcdf<-cdf-empcdf

    Tstat<-c(Tstat,abs(min(diffcdf))+max(diffcdf))
  }

  pval=(sum(Tstat>T0)+1)/(B+1)

  return(pval)
}

.DivLikes<-function(p,s,n_gs,likes,like0,nLeft){
  likes<-likes[which(!is.na(likes))]
  shift_likes=likes-like0

  #generate bins for frequencies
  break_min<-min(floor(min(shift_likes)),0)
  break_max<-max(ceiling(max(shift_likes)),100)
  br=seq(break_min,break_max,by=1)

  #calculate frequencies of likelihoods
  freq<-hist(shift_likes,breaks=br,plot=F)$counts/nLeft

  #calculate frequencies from beta mixture
  freq_null<-rep(0,(length(br)-1))
  for (j in 1:(length(br)-1)){
    freq_null[j]<-MixBetaDens(nLeft,p,x=seq(0,1,by=1),a=br[j],b=br[j+1],n_g=n_gs,var=s)$area
  }

  #Find K-L Divergence#
  KL<-entropy::KL.plugin(freq[which(freq_null!=0)],freq_null[which(freq_null!=0)])

  return(KL)
}

.oLikes<-function(X,maxO,G,B,grossOuts,kuiper=kuiper,pval=1,modelNames="VVV",mc.cores=1,nmax,prnt=FALSE){
  P=pval
  n=nrow(X)
  p=ncol(X)
  nO=length(grossOuts)
  ninitO=nO

  if(prnt){
    total <- maxO+1
    pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = total,clear=FALSE)
    pb$tick(nO)
  }

  # We reset the number of cores to use
  if((n-nO)<mc.cores) mc.cores = n-nO

  # We switch to the right number of cores + a warning if necessary
  max_nb_of_cores = parallel::detectCores()
  if(mc.cores>max_nb_of_cores){
    warning("The argument mc.cores is greater than its maximum.\nmc.cores was set to ", max_nb_of_cores)
    mc.cores = max_nb_of_cores
  }
  #
  # Parallel estimations
  #

  if (mc.cores > 1 && Sys.info()[['sysname']] == 'Windows') {
    cl <- parallel::makeCluster(mc.cores)
    .likes.j.wrapper <- function(i) {
      .likes.j(i, x = newX, G = G, z = z, modelNames = modelNames, nmax = nmax)
    }
    # Load necessary packages on each worker
    parallel::clusterEvalQ(cl, {
      library(oclust)
    })

    #export variables for later
    parallel::clusterExport(cl, varlist = c("G", "modelNames", "nmax"
    ), envir = environment())
  }

  indsLeft<-1:nrow(X)
  initcand<-grossOuts
  outliers<-initcand
  if(!is.null(outliers)) indsLeft<-indsLeft[-outliers]

  #remove the outliers from the dataset
  newX<-X
  if(!is.null(outliers)) newX<-X[-initcand,]


  i=1
  pval<-rep(NA,(maxO-nO+1))
  names(pval)=(nO):maxO
  maxpval=-Inf

  KL<-rep(NA,(maxO-nO+1))
  names(KL)=(nO):maxO
  minKL=Inf

  while(nO<=maxO && maxpval<P){
    dist_mat <- dist(newX, method = 'euclidean')
    hclust_X <- hclust(dist_mat, method = 'ward.D2')
    class<-cutree(hclust_X, k = G)
    #turn into z matrix
    z<-matrix(0,nrow(newX),G)
    for (k in 1:G){
      z[class==k,k]=1
    }
    if (min(table(class))<p+1) z=NULL
    mixtry=NULL
    modfitted=FALSE
    attempt=0
    while(modfitted==FALSE && attempt<5){
      attempt=attempt+1
      if (attempt==1){try(mixtry<-mixture::gpcm(newX,G=G,start=z,mnames = modelNames,nmax=nmax),silent=TRUE)}
      if (attempt==2){try(mixtry<-mixture::gpcm(newX,G=G,start=z,mnames = modelNames,veo=TRUE,nmax=nmax),silent=TRUE)}
      if (attempt==3){try(mixtry<-mixture::gpcm(newX,G=G,start=2,mnames = modelNames,veo=TRUE,nmax=nmax),silent=TRUE)}
      if (attempt==4){try(mixtry<-mclust::Mclust(data=newX,G=G,modelNames = modelNames,control=emControl(itmax=nmax),verbose = FALSE),silent=TRUE)}
      if (attempt==5){try(mixtry<-mclust::Mclust(data=newX,G=G,modelNames = modelNames,prior=priorControl(),control=emControl(itmax=nmax),verbose = FALSE),silent=TRUE)}
      modfitted=!is.null(mixtry)
    }
    if (is.null(mixtry)) {
      message("\nCould not fit model with ", nO, " outliers. Ending OCLUST early.")
      break
    }

    if(attempt<=3){
      like0=mixtry$best_model$loglik
      sigs=mixtry$best_model$model_obj[[1]]$sigs
      n_gs=as.vector(table(mixtry$map))
      pi_gs=mixtry$best_model$model_obj[[1]]$pi_gs
      z=mixtry$z
    }else{
      like0=mixtry$loglik
      sigs=lapply(1:G,function(g)mixtry$parameters$variance$sigma[,,g])
      n_gs=as.vector(table(mixtry$classification))
      pi_gs=matrix(mixtry$parameters$pro,1,G)
      z=mixtry$z
      }


    # compute the subset log-likelihoods
    if(mc.cores==1){
      liks<-unlist(lapply(X=1:nrow(newX), FUN=.likes.j, x=newX, G=G,z=z,modelNames=modelNames,nmax=nmax))
    } else if(Sys.info()[['sysname']] == 'Windows'){
      parallel::clusterExport(cl, varlist = c( "newX", "z"
      ), envir = environment())
      liks = NULL
      try(liks <- unlist(parallel::parLapply(cl, X = 1:nrow(newX), fun = .likes.j.wrapper)))
      if(is.null(liks)) {parallel::stopCluster(cl);
        stop("Unknown error in the parallel computing. Try mc.cores=1 to detect the problem.")}
    } else {
      liks= NULL
      try(liks <- unlist(parallel::mclapply(X=1:nrow(newX), .likes.j, x=newX, G=G,z=z,modelNames=modelNames,nmax=nmax,mc.cores=mc.cores)))
      if(is.null(liks)) stop("Unknown error in the parallel computing. Try mc.cores=1 to detect the problem.")
    }
    out<-which.max(liks)
    ind<-indsLeft[out]
    outliers<-c(outliers,ind)
    indsLeft<-indsLeft[-out]
    KL[i]<-.DivLikes(p=p,s=sigs,n_gs=n_gs,likes=liks,like0=like0,nLeft=nrow(newX))
    if (kuiper){pval[i]<-.pval(p,sigs,n_gs,liks,like0,nLeft =nrow(newX),B=B)}
    newX<-newX[-out,]
    if(!is.null(KL[i])&&!is.na(KL[i])&&KL[i]<minKL){finalModel=mixtry; minKL=KL[i]; numO=nO}
    if(kuiper&&!is.null(pval[i])&&!is.na(pval[i])&&pval[i]>P){finalModel=mixtry; maxpval=pval[i]; numO=nO}
    if(kuiper&&!is.null(pval[i])&&!is.na(pval[i])&&pval[i]>maxpval){maxpval=pval[i]}
    if (prnt) pb$tick(1)
    nO=nO+1
    i=i+1
  }
  if(mc.cores>1 && Sys.info()[['sysname']] == 'Windows') parallel::stopCluster(cl)

  allCand=outliers
  if(numO>0){outliers=outliers[1:numO]}else{outliers=NULL}

  pval=pval[1:(nO-ninitO)]
  KL=KL[1:(nO-ninitO)]

  return(list(numO=numO,outliers=outliers,KL=KL,pval=pval,allCand=allCand,finalModel=finalModel))
}


#' The OCLUST Algorithm
#'
#'oclust is a trimming method in model-based clustering. It iterates over possible values for the number of outliers and returns the model parameters for the best model as determined by the minimum KL divergence. If kuiper=TRUE, oclust calculates an approximate p-value using the Kuiper test and stops the algorithm if the p-value exceeds the specified threhold.

#' @param X A matrix or data frame with n rows of observations and p columns
#' @param maxO An upper bound for the number of outliers
#' @param G The number of clusters
#' @param grossOuts The indices of the initial outliers to remove. Default is NULL.
#' @param modelNames The model to fit using the gpcm function in the mixture package. Default is "VVV" (unconstrained). If modelNames=NULL, all models are fitted for each subset at each iteration. The BIC chooses the best model for each subset.
#' @param mc.cores Number of cores to use if running in parallel. Default is 1
#' @param nmax Maximum number of iterations for each EM algorithm. Decreasing nmax may speed up the algorithm but lose precision in finding the log-likelihoods.
#' @param kuiper A logical specifying whether to use the Kuiper test (Kuiper, 1960) to stop the algorithm when p-value exceeds the specified threshold. Default is FALSE.
#' @param pval The p-value for the Kuiper test. Default is 0.05.
#' @param B Number of samples to calculate the approximate p-value. Default is 100.
#' @param verb A logical specifying whether to print the current iteration number. Default is FALSE
#' @param scale A logical specifying whether to centre and scale the data. Default is TRUE
#'
#' @details
#' Gross outlier indices can be found with the \code{\link{findGrossOuts}} function.
#'
#' N. H. Kuiper, Tests concerning random points on a circle, in: Nederl. Akad. Wetensch. Proc. Ser. A, Vol. 63, 1960, pp. 38–47.
#'
#' @return
#' oclust returns a list of class oclust with
#' \item{data}{A list containing the raw and scaled data}
#' \item{numO}{The predicted number of outliers}
#' \item{outliers}{The most likely outliers in the optimal solution in order of likelihood}
#' \item{class}{The classification for the optimal solution}
#' \item{model}{The model selected for the optimal solution}
#' \item{G}{The number of clusters}
#' \item{pi.g}{The group proportions for the optimal solution}
#' \item{mu}{The cluster means for the optimal solution}
#' \item{sigma}{The cluster variances for the optimal solution}
#' \item{KL}{The KL divergence for each iteration, with the first value being for the initial dataset with the gross outliers removed}
#' \item{allCand}{All outlier candidates in order of likelihood}

#' @examples
#' # simulate 4D dataset
#' library(mvtnorm)
#' set.seed(123)
#' data <- rbind(rmvnorm(250, rep(-3, 4), diag(4)),
#'               rmvnorm(250, rep(3, 4), diag(4)))
#' # add outliers
#' noisy <- simOuts(data = data, alpha = 0.02, seed = 123)
#'
#' # Find gross outliers
#' findGrossOuts(X = noisy, minPts = 10)
#'
#' # Elbow between 5 and 10. Specify limits of graph
#' findGrossOuts(X = noisy, minPts = 10, xlim = c(5, 10))
#'
#' # Elbow at 9
#' gross <- findGrossOuts(X = noisy, minPts = 10, elbow = 9)
#'
#' # run algorithm
#' if (interactive()) {
#' # This example takes a few minutes to run
#'   result <- oclust(X = noisy, maxO = 15, G = 2, grossOuts = gross,
#'                   modelNames = "EEE", mc.cores = 1, nmax = 50,
#'                   kuiper = FALSE, verb = TRUE, scale = TRUE)
#' }
#'@export
#'
oclust<-function(X,maxO,G,grossOuts=NULL,modelNames="VVV",mc.cores=1,nmax=1000,kuiper=FALSE, pval=0.05, B=100,verb=FALSE,scale=TRUE){
  origDat<-as.matrix(X)
  X<-origDat
  if(scale==TRUE){X=scale(X)} else{X=scale(X,center = F,scale = F)}
  p=ncol(X)
  n=nrow(X)
  options(warn=1)
  is.int <- function(no){abs(no-round(no))<1e-15}
  if(!is.numeric(G)||!is.numeric(mc.cores)||!is.numeric(nmax)||!is.numeric(maxO)||!is.numeric(B)) stop("G, B, nmax, maxO and mc.cores must numeric");
  if(!is.int(G)||!is.int(mc.cores)||G<1||mc.cores<1||!is.int(nmax)||nmax<1||!is.int(B)||B<1) stop("G, B, nmax and mc.cores must be positive integers");
  if(!is.int(maxO)||maxO<0||maxO<length(grossOuts)) stop("maxO must be a positive integer and larger than the number of gross outliers specified");
  if(!is.null(grossOuts)){
    if(!is.numeric(grossOuts)) stop("Outlier indices must be numeric")
    if(sum(!is.int(grossOuts))||sum(grossOuts<1)) stop("Outlier indices must be positive integers");
    if(max(grossOuts)>n) stop("Outlier indices cannot be larger than the number of observations");
  }
  if(p<2) stop("oclust is for use on multivariate datasets only")
  if(any(is.na(X))) stop("oclust does not currently accommodate missing data")
  if(!is.null(modelNames) && sum(!(modelNames %in% c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","EEV","VVE","VEV","EVV","VVV")))) stop("modelNames must be NULL or at least one of 'EII','VII','EEI','VEI','EVI','VVI','EEE','VEE','EVE','EEV','VVE','VEV','EVV','VVV'")
  if(!is.logical(verb)||!is.logical(scale)) stop("verb and scale must be logical: TRUE or FALSE");
  if(kuiper==TRUE && !(0<pval && pval<=1)) stop("pval must be between 0 and 1");
  if(n-maxO<G*(p+1)) stop("Specified maxO is too large")
  if(maxO/n>0.5){
    warningmsg<-("Maximum outliers are more than 50%. Consider reducing maxO")
    warning(warningmsg)}
  m<-.oLikes(X=X,maxO=maxO,G=G,B=B,grossOuts=grossOuts,kuiper=kuiper,pval=pval,modelNames = modelNames,mc.cores=mc.cores,nmax=nmax,prnt=verb)

  numO<-as.numeric(m$numO)
  outliers=m$outliers

  pi.g<-m$finalModel$best_model$model_obj[[1]]$pi_gs
  mu<-m$finalModel$best_model$model_obj[[1]]$mus
  if (is.null(outliers)){class<-m$finalModel$map}else{class<-rep(0,n)
  class[-outliers]<-m$finalModel$map}
  sigma<-m$finalModel$best_model$model_obj[[1]]$sigs
  model<-m$finalModel$best_model$cov_type

  if (!(model %in% modelNames)) print(paste0("Requested model could not be estimated. Best model is ",model))
  params<-list(data=list(rawData=origDat,scaledData=X,center=attr(X, 'scaled:center'),scale=attr(X, 'scaled:scale')),numO=numO,outliers=outliers,class=class,model=model,G=G,pi.g=pi.g,mu=mu,sigma=sigma,KL=m$KL,pval=m$pval,allCand=m$allCand)

  class(params)='oclust'
  return (params)
}

.likes.j<-function(j,x,G,z,modelNames,nmax){
  w=NULL
  p=ncol(x)
  if (min(colSums(z))<p+1) z=NULL
  modfitted=FALSE
  attempt=0
  while(modfitted==FALSE && attempt<5){
    attempt=attempt+1
    if (attempt==1){try(w<-mixture::gpcm(x[-j,],G=G,mnames=modelNames,start=z[-j,],nmax = nmax)$best_model$loglik,silent=TRUE)}
    if (attempt==2){try(w<-mixture::gpcm(x[-j,],G=G,mnames=modelNames,start=z[-j,],veo=TRUE,nmax = nmax)$best_model$loglik,silent=TRUE)}
    if (attempt==3){try(w<-mixture::gpcm(x[-j,],G=G,mnames=modelNames,start=2,veo=TRUE,nmax = nmax)$best_model$loglik,silent=TRUE)}
    if (attempt==4){try(w<-mclust::Mclust(data=x[-j,],G=G,modelNames = modelNames,control=emControl(itmax=nmax),verbose = FALSE)$loglik,silent=TRUE)}
    if (attempt==5){try(w<-mclust::Mclust(data=x[-j,],G=G,modelNames = modelNames,prior=priorControl(),control=emControl(itmax=nmax),verbose = FALSE)$loglik,silent=TRUE)}
    modfitted=!is.null(w)
  }
  if(!modfitted) w<-NA
  return(w)
}

#' Print oclust
#'
#' Prints list of available components for \sQuote{oclust} class objects.
#'
#'
#' @method print oclust
#'
#' @param x An \sQuote{oclust} class object obtained by using \code{\link{oclust}}
#' @param ... additional print parameters
#' @export
print.oclust <- function(x,...)
{
  txt <- paste0("\'", class(x)[1], "\' object: ")
  txt <- paste0(txt, "(", "G=",x$G,", numO=",x$numO, ", model=", x$model,")")
  cat(txt)
  cat("\n")
  cat("\nAvailable components:\n")
  print(names(x))
  # str(x, max.level = 1, give.attr = FALSE, strict.width = "wrap")
  invisible(x)
}

#' Summarizes key results for  \sQuote{oclust} class objects.
#'
#'
#' @method summary oclust
#'
#' @param object An \sQuote{oclust} class object obtained by using \code{\link{oclust}}
#' @param ... additional summary arguments
#'
#' @export
summary.oclust <- function(object,...)
{
  # collect info
  G  <- object$G
  numO<-object$numO
  pi.g <- as.numeric(object$pi.g)
  names(pi.g) <- seq(G)
  mean <- object$mu
  names(mean)=seq(G)
  sigma <- object$sigma
  names(sigma)=seq(G)
  model=object$model
  title <- paste("Results of the OCLUST Algorithm")
  #
  obj <- list(title = title,G = G, numO=numO,
              outliers=object$outliers,classification=object$class,
              model=model,pi.g=pi.g, mean = mean, variance = sigma
  )
  class(obj) <- "summary.oclust"
  return(obj)
}
#' Prints the summary of key results for  \sQuote{oclust} class objects.
#'
#'
#' @method print summary.oclust
#'
#' @param x An \sQuote{oclust} class object obtained by using \code{\link{oclust}}
#' @param digits number of digits to print
#' @param ... additional print arguments
#'
#'
#' @export
print.summary.oclust <- function(x, digits = getOption("digits"),...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  message(txt)
  message("\n")
  message(x$title)
  message("\n")
  message(txt)

  message("\n")

  message(paste0("oclust model with ",
                 x$G, ifelse(x$G > 1, " components", " component"),
                 ", covariance structure ",
                 x$model,
                 " and ",x$numO," outliers", ":"))

  message("\n")
  #

  #
  message("\nClustering table:")
  print(table(factor(x$classification,
                     levels = { l <- c(0,seq_len(x$G))
                     })),
        digits = digits)
  #
  message("\nClassification:\n")
  print(x$classification, digits = digits)

  message("\nMixing probabilities:\n")
  print(x$pi.g, digits = digits)
  message("\nMeans:\n")
  print(x$mean, digits = digits)
  message("\nVariances:\n")
  print(x$variance, digits = digits)

  #
  invisible(x)
}

#' Plots results of the \sQuote{oclust} algorithm.
#'
#'
#' @method plot oclust
#'
#' @param x An \sQuote{oclust} class object obtained by using \code{\link{oclust}}
#' @param what A string specifying the type of graph. The options are:
#' "classification"    a plot of the classifications for the optimal solution.
#' For data with p>2, if more than two "dimens" are specified, a pairs plot is produced. If two "dimens" are specified,  a coordinate projection plot is produced with the specified "dimens".
#' Ellipses corresponding to covariances of mixture components are also drawn if "addEllipses = TRUE".
#' "KL"    a plot of Kullback-Leibler divergence for each number of outliers.
#' "pval"    a plot of approximate p-value for each number of outliers.
#' @param dimens a vector specifying the dimensions of the coordinate projections
#' @param xlab,ylab optional argument specifying axis labels for the classsification plot
#' @param ylim optional limits of the y axis of the BIC and KL plots
#' @param addEllipses logical indicating whether to include ellipses corresponding
#' to the covariances of the mixture components
#' @param ... other graphical parameters
#'
#' @export
plot.oclust <- function(x,
                        what = c("classification", "KL","pval"),
                        dimens = NULL, xlab = NULL, ylab = NULL, ylim = NULL,
                        addEllipses = TRUE,
                        ...)
{
  main=F
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "oclust"))
    stop("object not of class \"oclust\"")

  data <- object$data$scaledData
  p <- ncol(data)
  sig<-array(unlist(object$sigma),dim = c(p,p,object$G))
  mu<-matrix(unlist(object$mu),p,object$G,byrow =F )
  if(p == 1)
    colnames(data) <- deparse(x$call$data)
  if(is.null(dimens))
    dimens <- seq(p)
  else
    dimens <- dimens[dimens <= p]
  d <- length(dimens)
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)

  what <- match.arg(what, several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  plot.oclust.KL <- function(...)
    mclust::plot.mclustBIC(matrix(object$KL,dimnames = list(names(object$KL),"EEE")),xlab="Number of Outliers",ylab="KL",legendArgs = list(x = "bottom", ncol = 2, cex = 1, inset = -500),ylim = ylim, ...)

  plot.oclust.pval <- function(...)
  {if(sum(is.na(object$pval))==length(object$pval)) stop("pval can only be plotted if kuiper=TRUE");
    mclust::plot.mclustBIC(matrix(object$pval,dimnames = list(names(object$pval),"EEE")),xlab="Number of Outliers",ylab="pval",legendArgs = list(x = "bottom", ncol = 2, cex = 1, inset = -500),ylim = ylim, ...)
  }

  plot.oclust.classification <- function(...)
  {
    if(p == 1)
    { mclust::mclust1Dplot(data = data,
                           # parameters = object$parameters,
                           what = "classification",
                           classification = object$class,
                           z = NULL,
                           xlab = if(is.null(xlab)) colnames(data)[dimens] else xlab,
                           main = main, ...)
    }
    if(p == 2)
    { mclust::mclust2Dplot(data = data, what = "classification",
                           classification = object$class,
                           parameters = if(addEllipses) list(mean=mu,variance=list(sigma=sig)) else NULL,
                           xlab = if(is.null(xlab)) colnames(data)[1] else xlab,
                           ylab = if(is.null(ylab)) colnames(data)[2] else ylab,
                           main = main, ...)
    }
    if(p > 2)
    {
      if(d == 2)
      {
        mclust::coordProj(data = data, what = "classification",
                          parameters = list(mean=mu,variance=list(sigma=sig)),
                          classification = object$class,
                          addEllipses = addEllipses,
                          dimens = dimens,
                          main = main, ...)
      } else
      {
        par(mfrow = c(d, d),
            mar = rep(c(0.3,0.3/2),each=2),
            oma = c(4, 4, 4, 4))
        on.exit(par(oldpar))
        for(i in seq(d))
        { for(j in seq(d))
        { if(i == j)
        { plot(data[,c(j,i)],type="n",xlab="",ylab="",axes=FALSE)
          text(mean(par("usr")[1:2]),
               mean(par("usr")[3:4]),
               labels = colnames(data[,dimens])[i],
               cex=1.5, adj=0.5)
          box()
        }
          else
          { mclust::coordProj(data = data,
                              what = "classification",
                              parameters = list(mean=mu,variance=list(sigma=sig)),
                              classification = object$class,
                              addEllipses = addEllipses,
                              dimens = dimens[c(j,i)],
                              main = FALSE,
                              xaxt = "n", yaxt = "n", ...)
          }
          if(i == 1 && (!(j%%2))) axis(3)
          if(i == d && (j%%2))    axis(1)
          if(j == 1 && (!(i%%2))) axis(2)
          if(j == d && (i%%2))    axis(4)
        }
        }
      }
    }
  }

  if(interactive() & length(what) > 1)
  { title <- "OCLUST plots:"
  # present menu waiting user choice
  choice <- menu(what, graphics = FALSE, title = title)
  while(choice != 0)
  { if(what[choice] == "classification") plot.oclust.classification(...)
    if(what[choice] == "KL")    plot.oclust.KL(...)
    if(what[choice] == "pval")    plot.oclust.pval(...)
    # re-present menu waiting user choice
    choice <- menu(what, graphics = FALSE, title = title)
  }
  }
  else
  { if(any(what == "classification")) plot.oclust.classification(...)
    if(any(what == "KL"))    plot.oclust.KL(...)
    if(any(what == "pval"))    plot.oclust.pval(...)
  }

  invisible()
}



#' @importFrom("graphics", "axis", "box", "hist", "par", "plot", "text")
#' @importFrom("stats", "dbeta", "pbeta")
#' @importFrom("utils", "flush.console", "menu")
