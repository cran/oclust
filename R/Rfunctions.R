MixBetaDens<-function(n,p,x =seq(0, 15, by=0.01),a=0,b=1,n_g=n_g,var=var){
  n_g<-c(n_g)
  pi_g<-n_g/n
  #turn variance into array
  if (length(dim(var))<3){
    var<-array(var,c(dim(var),1))
  }
  #find determinants of variances
  log_det_var<-pi_g
  for (i in 1:length(pi_g)){
    log_det_var[i]<-determinant(var[,,i])$modulus
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



.DivLikes<-function(func_likes){
  n=func_likes$n
  p=func_likes$p
  s=func_likes$s
  shift_likes=func_likes$shift_likes
  classification=func_likes$class
  G=func_likes$G
  o=func_likes$o

  n_g=matrix(0,G,(o+1))
  for (j in 1:(o+1)){
    for (i in 1:G){
      n_g[i,j]<-length(which(classification[,j]==i))
    }}

  #generate bins for frequencies
  break_min<-floor(min(shift_likes,na.rm=T))
  break_max<-ceiling(max(shift_likes,na.rm=T))
  br=seq(break_min,break_max,by=1)

  #calculate frequencies of likelihoods
  freq<-matrix(0,(break_max-break_min),(o+1))
  for (i in 1:(o+1)){
    freq[,i]<-hist(shift_likes[which(!is.na(shift_likes[,i])),i],breaks=br,plot=F)$counts/(n+1-i)
  }
  #calculate frequencies from beta mixture
  freq_null<-matrix(0,(break_max-break_min),(o+1))
  for (i in 1:(o+1)){
    for (j in 1:(break_max-break_min)){
      freq_null[j,i]<-MixBetaDens((n-i+1),p,x=seq(0,1,by=1),a=br[j],b=br[j+1],n_g=n_g[,i],var=s[,,i,])$area
    }
  }

  #Find K-L Divergence#
  KL<-rep(0,(o+1))
  for (i in 1:(o+1)){
    KL[i]<-entropy::KL.plugin(freq[which(freq_null[,i]!=0),i],freq_null[which(freq_null[,i]!=0),i])}

  list(KL=KL)
}


.oLikes<-function(x,o,G,modelNames=NULL,prior=NULL,mc.cores=1,prnt=FALSE){
  n<-nrow(x)
  p<-ncol(x)
  #find least likely points#
  #Make distribution of likelihoods#
  like0<-rep(0,(o+2))
  bic<-rep(NA,(o+1))
  liks<-matrix(0,n,(o+1))
  index<-rep(0,(o+1))
  trim_data<-array(NA,c(n,p,(o+1)))
  pi_g<-matrix(NA,G,(o+1))
  classification<-matrix(NA,n,(o+1))
  mu<-array(NA,c(p,G,(o+1)))
  var<-array(0,c(p,p,(o+1),G))
  x1<-as.matrix(x)

  # We reset the number of cores to use
  if(n<mc.cores) mc.cores = n

  # We switch to the right number of cores + a warning if necessary
  max_nb_of_cores = parallel::detectCores()
  if(mc.cores>max_nb_of_cores){
    warning("The argument mc.cores is greater than its maximum.\nmc.cores was set to ", max_nb_of_cores)
    mc.cores = max_nb_of_cores
  }
  #
  # Parallel estimations
  #

  if(mc.cores>1 && Sys.info()[['sysname']] == 'Windows'){
    cl = parallel::makeCluster(mc.cores)
    loadMyPackages = function(x){
      # we load the package
      library(oclust)
    }
    par.setup = parallel::parLapply(cl, 1:length(cl), loadMyPackages)
  }

  for (i in 1:(o+1)){
    trim_data[,,i]<-x1
    clust_i<-mclust::Mclust(x1[which(!is.na(x1[,1])),],G,modelNames = modelNames,prior = prior,verbose=F)
    if (is.null(clust_i)){clust_i<-mclust::Mclust(x1[which(!is.na(x1[,1])),],G,prior=prior,verbose=F)}
    bic[i]<-clust_i$bic
    like0[i]<-clust_i$loglik
    pi_g[,i]<-clust_i$parameters$pro
    classification[which(!is.na(x1[,1])),i]<-clust_i$classification
    var[,,i,]<-clust_i$parameters$variance$sigma
    mu[,,i]<-clust_i$parameters$mean
    likes<-rep(NA,nrow(x1))
    if(prnt){
      message("\n","o=",i-1,"\n")
      flush.console()
    }
    # compute the subset log-likelihoods
    if(mc.cores==1){
      likes<-unlist(lapply(1:n, .likes.j, x1=x1, G=G,modelNames=modelNames,prior=prior))
    } else if(Sys.info()[['sysname']] == 'Windows'){
    likes = NULL
    try(likes <- unlist(parallel::parLapply(cl, 1:n, .likes.j, x1=x1, G=G,modelNames=modelNames,prior=prior)))
    if(is.null(likes)) {parallel::stopCluster(cl);
      stop("Unknown error in the parallel computing. Try mc.cores=1 to detect the problem.")}
    } else {
      likes= NULL
      try(likes <- unlist(parallel::mclapply(1:n, .likes.j, x1=x1, G=G,modelNames=modelNames,prior=prior, mc.cores=mc.cores)))
      if(is.null(likes)) stop("Unknown error in the parallel computing. Try mc.cores=1 to detect the problem.")
    }
    liks[,i]<-likes
    max_like<-which(likes==max(likes,na.rm = T))[1]
    index[i]<-max_like
    x1[max_like,]<-rep(NA,p)
  }

  if(mc.cores>1 && Sys.info()[['sysname']] == 'Windows') parallel::stopCluster(cl)

  like0[o+2]<-max(likes,na.rm=T)
  index<-index[1:o]
  #Calculate the variances#
  s<-array(0,c(p,p,(o+1),G))
  for (j in 1:(o+1)){
    for (i in 1:G){
      s[,,j,i]<-var(trim_data[which(classification[,j]==i),,j])
    }}

  #Shift the likelihoods#
  shift_likes<-liks
  for (i in 1:(o+1)){
    shift_likes[,i]<-liks[,i]-like0[i]
  }

  list(index=index,max_like=like0,bic=bic,shift_likes=shift_likes,
       trim_data=trim_data,pi_g=pi_g,mu=mu,s=s,var=var,
       class=classification,n=n,p=p,o=o,G=G)
}

oclust<-function(x,o,G,modelNames=NULL,prior=NULL,mc.cores=1,keepAllRes=F,verb=F){
  is.int <- function(no){abs(no-round(no))<1e-15}
  if(!is.int(o)||!is.int(G)||!is.int(mc.cores)||o<1||G<1||mc.cores<1) stop("o, G and mc.cores must be positive integers");
  x<-as.matrix(x)
  p=ncol(x)
  if(p<2) stop("oclust is for use on multivariate datasets only")
  if(any(is.na(x))) stop("oclust does not currently accommodate missing data")
  m<-.oLikes(x=x,o=o,G=G,modelNames = modelNames,prior=prior,mc.cores=mc.cores,prnt=verb)
  div<-.DivLikes(m)
  KL<-div$KL
  numO<-which.min(KL)-1
  BIC<-m$bic
  bic<-BIC[numO+1]

  pi.g<-m$pi_g[,(numO+1)]
  mu<-m$mu[,,(numO+1)]
  class<-m$class[,(numO+1)]
  class[is.na(class)]<-0
  outs<-m$index[1:(numO+1)]
  sigma<-m$var[,,(numO+1),]

  params<-list(data=as.matrix(x),numO=numO,G=m$G,outs=outs,class=class,pi.g=pi.g,mu=mu,sigma=sigma,KL=KL,BIC=BIC,bic=bic)

  if(keepAllRes){
    all_results = m
    all_results$bic=NULL
    all_results$G=NULL
    params$all_results = all_results
  }
  class(params)='oclust'
  return (params)
}

.likes.j<-function(j,x1,G,modelNames,prior){
    if (!is.na(x1[j,1])){
      k<-as.matrix(x1[-j,]);
      w<-mclust::Mclust(k[which(!is.na(k[,1])),],G,modelNames = modelNames,prior = prior,verbose=F)$loglik;
      if(is.null(w))w<-mclust::Mclust(k[which(!is.na(k[,1])),],G,prior=prior,verbose=F)$loglik}else{w=NA}
    return(w)
}

print.oclust <- function(x,...)
{
  txt <- paste0("\'", class(x)[1], "\' object: ")
  txt <- paste0(txt, "(", "G=",x$G,", numO=",x$numO, ")")
  cat(txt)
  cat("\n")
  cat("\nAvailable components:\n")
  print(names(x))
  # str(x, max.level = 1, give.attr = FALSE, strict.width = "wrap")
  invisible(x)
}

summary.oclust <- function(object,...)
{
  # collect info
  G  <- object$G
  numO<-object$numO
  pi.g <- object$pi.g
  names(pi.g) <- seq(G)
  mean <- object$mu
  colnames(mean)=seq(G)
  sigma <- object$sigma
  title <- paste("Results of the OCLUST Algorithm")
  #
  obj <- list(title = title,G = G, numO=numO,
              outs=object$outs,classification=object$class,
              pi.g=pi.g, mean = mean, variance = sigma,
              classification = object$class,BIC=object$bic)
  class(obj) <- "summary.oclust"
  return(obj)
}

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
                   " and ",x$numO," outliers", ":"))

  message("\n")
  #

  #
  tab <- data.frame("BIC" = x$BIC, "numO" = x$numO,
                    row.names = "", check.names = FALSE)
  print(tab, digits = digits,...)
  #
  message("\nClustering table:")
  print(table(factor(x$classification,
                     levels = { l <- c(0,seq_len(x$G))
                     })),
        digits = digits)
  #

    message("\nMixing probabilities:\n")
    print(x$pi.g, digits = digits)
    message("\nMeans:\n")
    print(x$mean, digits = digits)
    message("\nVariances:\n")
    print(x$variance, digits = digits)

    message("\nClassification:\n")
    print(x$classification, digits = digits)

    #
    invisible(x)
}

plot.oclust <- function(x,
                        what = c("BIC", "classification", "KL"),
                        dimens = NULL, xlab = NULL, ylab = NULL, ylim = NULL,
                        addEllipses = TRUE,
                        ...)
{
  main=F
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "oclust"))
    stop("object not of class \"oclust\"")

  data <- object$data
  p <- ncol(data)
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

  plot.oclust.bic <- function(...)
  mclust::plot.mclustBIC(matrix(object$BIC,dimnames = list(as.character(0:(length(object$BIC)-1)),"EEE")),xlab="Number of Outliers",legendArgs = list(x = "bottom", ncol = 2, cex = 1, inset = -500000000),ylim = ylim, ...)

  plot.oclust.kl <- function(...)
  mclust::plot.mclustBIC(matrix(object$KL,dimnames = list(as.character(0:(length(object$BIC)-1)),"EEE")),xlab="Number of Outliers",ylab="KL",legendArgs = list(x = "bottom", ncol = 2, cex = 1, inset = -500),ylim = ylim, ...)

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
                   parameters = if(addEllipses) list(mean=object$mu,variance=list(sigma=object$sigma)) else NULL,
                   xlab = if(is.null(xlab)) colnames(data)[1] else xlab,
                   ylab = if(is.null(ylab)) colnames(data)[2] else ylab,
                   main = main, ...)
    }
    if(p > 2)
    {
      if(d == 2)
      {
        mclust::coordProj(data = data, what = "classification",
                  parameters = list(mean=object$mu,variance=list(sigma=object$sigma)),
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
                      parameters = list(mean=object$mu,variance=list(sigma=object$sigma)),
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
  { if(what[choice] == "BIC")            plot.oclust.bic(...)
    if(what[choice] == "classification") plot.oclust.classification(...)
    if(what[choice] == "KL")    plot.oclust.kl(...)
    # re-present menu waiting user choice
    choice <- menu(what, graphics = FALSE, title = title)
  }
  }
  else
  { if(any(what == "BIC"))            plot.oclust.bic(...)
    if(any(what == "classification")) plot.oclust.classification(...)
    if(any(what == "KL"))    plot.oclust.kl(...)
  }

  invisible()
}
