fit_alk10 <- function(alk.models, lengths = 1:200, minAge = 1, maxAge = 4) {

  d <- data_frame(LngtCm = 1:200)
  cc <- 1:nrow(d)

  p = matrix(1, nrow = nrow(d), ncol = maxAge)

  punc <- function(k, a) {
    p[k, a] * prod(1 - p[k, 1:(a - 1)])

  }

  W = getOption("warn")
  options(warn = -1)  ## disable warnings temporarily

  for (i in 1:(maxAge - 1)) {
    p[, i] = 1 - predict(alk.models[[i]],
                         newdata = d,
                         type = "response",
                         newdata.guaranteed = TRUE)

  }
  options(warn = W) ## restore warnings

  for (a in 2:maxAge) {
    ## unconditional prob of age_i = Pi_i * Prod(1-Pi_j , j=1..i-1) for i>1. [Rindorf,Lewy p.2]
    p[, a] = sapply(cc, punc, a = a)
  }

  colnames(p) <- minAge:maxAge
  p <-
    p %>%
    as_tibble() %>%
    mutate(length = d$LngtCm) %>%
    gather(age, p, -length, convert = TRUE)

  return(p)
}







## @title Find all variable names in formula which is also in hydro data
## @param formula A formula object
## @param x a DATRASraw object
## @return a vector of variable names
xtraVars <- function(formula,x){
  intersect(all.vars(formula),names(x[[2]]))
}

##' @title Fit part of the continuation ratio model using GAMs (helper function)
##' @param a Condition on being at least this age.
##' @param ages Vector with consecutive ages to consider, last age is plus group.
##' @param AL Age-length data from a DATRASraw object
##' @param model Model formula (string or formula), or a vector of strings specifying the formula for each age group.
##' @param gamma Multiplier for AIC score (see ?gam)
##' @param autoChooseK Automatic choice of the max. dimension for the basis used to represent the smooth term for spatial ALK. See ?s in the mgcv-package.
##' @param useBIC Use Bayesian Information Criterion for smoothness selection instead of AIC.
##' @param varCof Use varying coefficients model for spatial effect.
##' @param maxK Maximum k to use. Only applies if autoChooseK is TRUE.
##' @param verbose Print model summary?
##' @param ... Extra parameters to gam()
##' @return Object of class '"gam"'
fitALKoneXX <-function(a,
                       ages,
                       AL,
                       model,
                       gamma,
                       autoChooseK=FALSE,
                       useBIC=FALSE,
                       varCof=FALSE,
                       maxK=100,
                       verbose=FALSE,...){

  if(length(model)>1){
    idx=which(ages==a);
    f <- as.formula(model[idx])
  } else {
    f <- as.formula(model)
  }
  require(mgcv,quietly=TRUE)

  myd=subset(AL,Age>=a)
  myd$cra=as.factor(myd$Age>a)

  ## automatic choice of k - overrides model arguments
  if(autoChooseK) {
    uniqueCovs = length( unique( paste( myd$lon, myd$lat)))
    if(!varCof){
      k = min(maxK, uniqueCovs-1);
      f = as.formula(paste("cra~LngtCm+s(lon,lat,k=",k,",bs='ts')"))
      if(uniqueCovs<10) f=as.formula("cra~LngtCm")
    } else {
      k = min(maxK, uniqueCovs/2-1);
      f = as.formula(paste("cra~s(lon,lat,by=LngtCm,k=",k,",bs='ts')+s(lon,lat,k=",k,",bs='ts')"))
      if(uniqueCovs<10) f=as.formula("cra~LngtCm")
    }
    if(useBIC) gamma = log( sum( myd$NoAtALK)) / 2
  }

  #m <- tryCatch.W.E( gam(f,data=myd,family="binomial",weights=NoAtALK,gamma=gamma,...) )$value
  m <- gam(f,data=myd,family="binomial",weights=NoAtALK,gamma=gamma,...)
  if(class(m)[2]=="error") {
    print(m)
    stop("Error occured for age ",
         a,
         "\n",
         "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n");}
  if(verbose) print(summary(m))
  return(m)
}
##' @title Fit a continuation-ratio logit model for age given length and possibly other covariates.
##' @param x a DATRASraw object.
##' @param minAge minimum age group to consider
##' @param maxAge maximum age group to consider
##' @param mc.cores The number of cores to use, i.e. how many processes will be spawned (at most)
##' @param model Model formula(string) for ALK, or a vector of strings specifying the formula for each logit (i.e. number of age groups minus one).
##' @param method Use default formula: 1=Normal distributed age-at-length (1st order), 2=Second order polynomial, 3=Spatial ALK, first order length-effect.
##' @param autoChooseK Automatic choice of the max. dimension for the basis used to represent the smooth term for spatial ALK. See ?s in the mgcv-package.
##' @param useBIC Use Bayesian Information Criterion for smoothness selection instead of AIC.
##' @param varCof Use varying coefficients model for spatial effect.
##' @param maxK Maximum k to use. Only applies if autoChooseK is TRUE.
##' @param gamma Multiplier for AIC score (see ?gam)
##' @param verbose Print details about the fitting process.
##' @param ... Optional extra arguments to gam()
##' @return An object of class 'ALKmodel'
fitALKXX<-function(x,
                   minAge,
                   maxAge,
                   mc.cores=1,
                   method=1,
                   model= c( "cra~LngtCm", "cra~poly(LngtCm,2)", "cra~LngtCm+s(lon,lat,bs='ts')" )[method],
                   autoChooseK=FALSE,
                   useBIC=FALSE,
                   varCof=FALSE,
                   maxK=100,
                   gamma=1.4,
                   verbose=FALSE,...){

  #checkSpectrum(x);
  if( (minAge+1)>maxAge ) stop("Invalid age selection.");
  ages=minAge:maxAge;
  nAges=length(ages);
  ## remove last selected age group, as it is automatically added as plus group later
  lastAge=ages[nAges];
  ages=ages[-nAges];

  extraVars=unlist(  lapply( lapply(model,as.formula),xtraVars,x=x) );
  if(length(extraVars)==0) extraVars=NULL;

  x <-
    x %>%
    dplyr::rename(haul.id = id,
                  Age = age,
                  LngtCm = length,
                  NoAtALK = n)
  ##merge covariates from hydro data to age data by haul id.
  #x[[1]]=merge(x[[1]],x[[2]][c("lon","lat","haul.id",extraVars)],by="haul.id",all.x=TRUE,sort=FALSE,suffixes=c("",".y"))
  #x[[1]]=subset(x[[1]],!is.na(Year) & !is.na(Age) & !is.na(LngtCm))

  mylapply<-function(...){
    hasmc=(mc.cores>1 && require(multicore,quietly=TRUE));
    if(!hasmc) return(lapply(...)) else return(mclapply(...,mc.cores=mc.cores))
  }
  if(verbose) cat("Fitting model...");
  #models = mylapply(ages,fitALKone,ages=ages,AL=x[[1]],model=model,gamma=gamma,autoChooseK=autoChooseK,useBIC=useBIC,varCof=varCof,maxK=maxK,verbose=verbose,...);
  models = mylapply(ages,fitALKoneXX,ages=ages,AL=x,model=model,gamma=gamma,autoChooseK=autoChooseK,useBIC=useBIC,varCof=varCof,maxK=maxK,verbose=verbose,...);
  #class(models)<-"ALKmodel";
  #attr(models,"data")<-x;
  #attr(models,"ALKformula")<-model;
  #attr(models,"ages")<-ages;
  return(models)
}



