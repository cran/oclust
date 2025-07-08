library(mvtnorm)
library(testthat)
data<-rbind(rmvnorm(250,rep(-3,4),diag(4)),
          rmvnorm(250,rep(3,4),diag(4)))
dat<-simOuts(data=data,alpha=0.02,seed=123)
gross<-findGrossOuts(dat,10,elbow=9)

#full oclust without Kuiper test
testfull<-oclust(X=dat,maxO = 10,G=2,grossOuts = gross,modelNames="EEE",mc.cores=2,nmax=20,verb=TRUE,scale=TRUE)
print(testfull)
summary(testfull)
plot(testfull,what="KL")
plot(testfull,what="classification")

#oclust with Kuiper
testkuiper<-oclust(X=dat,maxO = 10,G=2,grossOuts = gross,modelNames="EEE",mc.cores=2,nmax=20,kuiper=TRUE,pval=0.05,B=100,verb=TRUE,scale=TRUE)
print(testkuiper)
summary(testkuiper)
plot(testkuiper,what="KL")
plot(testkuiper,what="pval")
plot(testkuiper,what="classification")

test_that("error for non-int elbow", {
  expect_error(findGrossOuts(dat,10,NULL,.5))
})

test_that("error for zero elbow", {
  expect_error(findGrossOuts(dat,10,NULL,0))
})

test_that("error for only one xlim", {
  expect_error(findGrossOuts(dat,10,xlim=2))
})

test_that("error non-int", {
  expect_error(findGrossOuts(dat,10.9))
})

test_that("error non-pos", {
  expect_error(findGrossOuts(dat,-1))
})

test_that("error non-int", {
  expect_error(findGrossOuts(dat,minPts=NULL))
})

test_that("error non-numeric", {
  expect_error(findGrossOuts(dat,minPts=2,elbow="string"))
})

test_that("error for non-int G", {
  expect_error(oclust(X=dat,maxO=10,G=3.5,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=20,verb=TRUE,scale=TRUE))
})
test_that("error for non-int maxO", {
  expect_error(oclust(X=dat,maxO=10.2,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=20,verb=TRUE,scale=TRUE))
})
test_that("error for invalid modelNames", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = NULL,modelNames="VVX",mc.cores = 1, nmax=20,verb=TRUE,scale=TRUE))
})

test_that("error non-int nmax", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=22.2,verb=TRUE,scale=TRUE))
})

test_that("error neg nmax", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=-2,verb=TRUE,scale=TRUE))
})

test_that("error non-numeric", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax="string",verb=TRUE,scale=TRUE))
})

test_that("error nonbool verb", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=2,verb="string",scale=TRUE))
})

test_that("error nonbool scale", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=2,verb=TRUE,scale="string"))
})

test_that("error pval>1", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=2,kuiper=TRUE,pval=1.01,B=100,verb=TRUE,scale=TRUE))
})

test_that("error nonbool scale", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=2,verb=TRUE,scale="string"))
})

test_that("data is multivariate", {
  expect_error(oclust(X=rnorm(100,0,1),maxO=10,G=3,grossOuts = NULL,modelNames="VVV",mc.cores = 1, nmax=2,verb=TRUE,scale=TRUE))
})

test_that("error more gross outliers than maxO", {
  expect_error(oclust(X=dat,maxO=4,G=3,grossOuts = gross,modelNames="VVV",mc.cores = 1, nmax=2,verb=TRUE,scale=TRUE))
})

test_that("error grossOuts non-numeric", {
  expect_error(oclust(X=dat,maxO=10,G=3,grossOuts = c("one","two","seventeen"),modelNames="VVV",mc.cores = 1, nmax=2,verb=TRUE,scale=TRUE))
})

test_that("error mus not vector minMD", {
  expect_error(minMD(X=dat,sigs=list(diag(13),diag(13)),mus=list(diag(13),diag(13))))
})

test_that("error lists not same length minMD", {
  expect_error(minMD(X=dat,sigs=list(diag(13),diag(13)),mus=list(rep(1,13),rep(2,13),rep(1,13))))
})

test_that("error sigs not matrices minMD", {
  expect_error(minMD(X=dat,sigs=list(rep(1,13),rep(1,13)),mus=list(rep(1,13),rep(2,13))))
})

test_that("error sigs not pxp minMD", {
  expect_error(minMD(X=dat,sigs=list(matrix(1,13,12),matrix(1,13,12)),mus=list(rep(1,13),rep(2,13))))
})

test_that("error not list minMD", {
  expect_error(minMD(X=dat,sigs=array(1,dim=c(13,13,2)),mus=list(rep(1,13),rep(2,13))))
})

test_that("error not int MixBetaDens", {
  expect_error(MixBetaDens(n=101.5,p=2,x=seq(0,15,0.01),a=1,b=2,n_g=c(0.1,0.9),var=list(diag(2),diag(2))))
})

test_that("error not int MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=2.1,x=seq(0,15,0.01),a=1,b=2,n_g=c(0.1,0.9),var=list(diag(2),diag(2))))
})

test_that("error not numeric MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=2,x=NULL,a=1,b=2,n_g=c(0.1,0.9),var=list(diag(2),diag(2))))
})

test_that("error not numeric MixBetaDens", {
  expect_error(MixBetaDens(n=NULL,p=2,x=seq(0,15,0.01),a=1,b=2,n_g=c(0.1,0.9),var=list(diag(2),diag(2))))
})

test_that("error not numeric MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=TRUE,x=seq(0,15,0.01),a=1,b=2,n_g=c(0.1,0.9),var=list(diag(2),diag(2))))
})

test_that("error not numeric MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=2,x=seq(0,15,0.01),a=NULL,b=2,n_g=c(0.1,0.9),var=list(diag(2),diag(2))))
})

test_that("error not numeric MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=2,x=seq(0,15,0.01),a=1,b=NULL,n_g=c(0.1,0.9),var=list(diag(2),diag(2))))
})

test_that("error not numeric MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=2,x=seq(0,15,0.01),a=1,b=2,n_g=TRUE,var=list(diag(2),diag(2))))
})

test_that("error not numeric MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=2,x=seq(0,15,0.01),a=1,b=2,n_g=c(0.1,0.9),var=list(NULL,diag(2))))
})

test_that("a>b MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=2,x=seq(0,15,0.01),a=1,b=0,n_g=c(0.1,0.9),var=list(diag(2),diag(2))))
})

test_that("n_g and var same length MixBetaDens", {
  expect_error(MixBetaDens(n=100,p=2,x=seq(0,15,0.01),a=1,b=2,n_g=c(1),var=list(diag(2),diag(2))))
})

test_that("alpha non-neg simOuts", {
  expect_error(simOuts(wine[,-1],-0.1,123))
})

test_that("alpha non-numeric simOuts", {
  expect_error(simOuts(wine[,-1],NULL,123))
})

test_that("seed non-numeric simOuts", {
  expect_error(simOuts(wine[,-1],0.1,"string"))
})


test_that("error for pval plot when no pval", {
  expect_error(plot(testfull,what="pval"))
})
