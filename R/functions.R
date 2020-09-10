#' Test of independence in presence of heteroscedastic treatments
#'
#' NPtest_indept performs the test of independence between the response
#' variable and a single covariate when there is potentially heteroscedastic
#' treatment effects present (see
#' Wang, Tolos and Wang (2010)).
#'
#' @param dat A data frame with three columns named X, trt, and Y, where
#'        X is the covarite, trt is the treatment level, and Y is the response
#'        variable.
#' @param k   An odd integer to specify the number of nearest neighbors
#'          to be used in augmentation. Generally recommend to use 3, 5, or 7.
#'
#' @return A list containing the following variables:
#'
#'        Asys_var: the asymptotic variance for the test statistics
#'
#'        Tstat:  the test statistic
#'
#'        pvalue: the p-value of the test under H0: independence between X and Y.
#'
#' @export
#'
#' @examples
#' n=64;  X=runif(n); trt=gl(2, n/2)
#' e=rnorm(n, 0, 0.1)
#' Y=ifelse(trt==1, 4*(X-0.5)^2+e, 2*X+e)
#' dat=data.frame(X, Y, trt)
#' NPtest_indept(dat, k=7)
#'
#' @references
#' Haiyan Wang, Siti Tolos, and Suojin Wang (2010). A Distribution Free
#'  Nonparametric Test to Detect Dependence Between a Response Variable and
#'  Covariate in Presence of Heteroscedastic Treatment Effects.
#'  The Canadian Journal of Statistics. 38(3), 408433. Doi:10.1002/cjs.10068
#'
NPtest_indept= function(dat, k=7)
{
  X=dat$X;  trt=dat$trt; Y1=dat$Y
  if (stats::var(Y1)> 1e-6 ){
    ranksuse=unlist(tapply(X, trt, rank) )
    alltrt=rbind(Y1, X, ranksuse )
    n=unlist(tapply(rep(1, nrow(dat)), trt, sum))
    N=sum(n); a=length(n)
    for (i1 in 1:a){
      locationi1=trt_position(i1,n);
      orderwant=order(alltrt[2,locationi1[1]:locationi1[2]])+locationi1[1]-1;
      alltrt[,locationi1[1]:locationi1[2] ]=  alltrt[,orderwant]
    }

    psudodat=makepseudo(N,n, k, a, alltrt)
    psudo=psudodat$psudo; index=psudodat$index
    cellmean<-apply(psudo, c(1,2), mean)
    colmean=apply(psudo, 2, mean)
    sig<- stats::cov(t(cellmean))
    # diagonal part gives the \hat\sigma_{1,i}^2
    # and off-diagonal part gives \hat\sigma_{1,i_1, i_2}

    sigXij<-apply(psudo, c(1, 2), stats::var)
    # get a axN matrix with \hat\sigma_i^2(X_{ij}) =sigXij[i, j]

    meanrk=apply(psudo,1, mean)
    MSTphi=k*sum( (cellmean-matrix(rep(meanrk,N),ncol=N))^2)/((N-1)*a)
    MSE=sum((psudo-array(rep(cellmean, k), c(a, N, k)) )^2)/(N*a*(k-1))
    Tss=(sqrt(N) * (MSTphi-MSE))

    ##****************************************
    #Calculate estimate of variance for test statistics
    # count is a matrix; first three columns give the value of i1, j2, i;
    #       the last column gives the number of times X_{ij_2} is used in
    #       construction of windows for all covariate values in group i_1

    count<-matrix(-1, a^2*N, 4)
    whereini=0
    for( i1 in 1:a){
      if (i1==1) lower=1 else lower=sum(n[1:(i1-1)])+1
      upper=sum(n[1:i1])
      for (j2 in 1:N){
        for (i in 1:a){
          whereini=whereini+1 ;  whereisXij2= mapindex(j2, n)
          counti1j2i=sum(index[i,lower:upper, ]==((whereisXij2[1]==i) *whereisXij2[2]) )
          count[whereini, ]=c(i1, j2, i, counti1j2i)
        }}}

    subcount= count[count[,1]!=count[,3], ] # entries in count that i1 \ne i
    prodcount1.1=tapply(subcount[,4],list(subcount[,2], subcount[,3]), sum)
    #above line is sum_{i1} I(i1 \ne i) of count at i1, j2, i
    #prodcount1.1 is \sum_{i_1, i_1 \ne i}^a \frac{n_{i_1}}{n_i} d_{i_1i}(X_{ij})

    tau3=0
    for (i in 1:a){
      starti=trt_position(i,n)[1]-1
      for (jp in starti+(2:n[i])){
        for (j in (max(1, (jp-k+1)): (jp-1)) ){
          Bijjp= (prodcount1.1[j,i]/k+1 ) * (k-jp+j)*(jp-j<=k-1)
          tau3=tau3+(Bijjp^2+Bijjp-2*(jp-j<=((k-1)/2) ))*(jp-j<=k-1)*sigXij[i,j] *sigXij[i,jp] *(j!=jp)
        } } }
    tau3=tau3*4/(N*a^2*(k-1)^2)

    #########
    tauAsys=tau3
    pvalue.sim=1-stats::pnorm(Tss/sqrt(tauAsys))
  } else {Tss=0; tauAsys=1e+8;  pvalue.sim=1}
  list(Asys_var=tauAsys,Tstat=Tss,pvalue=pvalue.sim)
}


#' Nearest neighbor augmentatation based on ranks
#'
#' The function makepseudo performs the nearest neighbor augmentation
#' based on the rank of covariate values according to the scheme discribed
#' on page 410-411 of Wang, Tolos and Wang (2010)
#'
#' @param N total number of covariate values.
#' @param n vector of sample sizes from all treatments
#' @param k number of nearest neighbors
#' @param a number of treatment levels in the data
#' @param alltrt a matrix of dimension 3x\code{N}, whose first two rows are
#'        Y and X, and the third row gives the rank of X values within the
#'        same treatment level.
#'
#' @return A list containing the following:
#'         psudo: a 3-d array of the dimension (k, a, N) that stores the
#'                augmented observations based on k-nearest neighbor rule
#'                in Wang, Tolos and Wang (2010).
#'         index: a 3-d array of the dimension (k, a, N) that stores the
#'                index of which observation was used for augmentation.
#'
#' @export
#'
#' @examples
#'  a=2; n=c(7,9); N=sum(n);  X=runif(N);
#'  trt=c(rep(1,n[1]), rep(2, n[2])); e=rnorm(N, 0, 0.1)
#'  Y=ifelse(trt==1, 4*(X-0.5)^2+e, 2*X+e)
#'  ranksuse=unlist(tapply(X, trt, rank) )
#'  alltrt=rbind(Y, X, ranksuse )
#'  aug=makepseudo(N,n, k=3, a, alltrt)
#'
#' @references
#' Haiyan Wang, Siti Tolos, and Suojin Wang (2010). A Distribution Free
#'  Nonparametric Test to Detect Dependence Between a Response Variable and
#'  Covariate in Presence of Heteroscedastic Treatment Effects.
#'  The Canadian Journal of Statistics. 38(3), 408433. Doi:10.1002/cjs.10068
#'
makepseudo=function(N,n, k, a, alltrt){
  psudo<-array(0, c(a, sum(n), k))
  index<-array(0, c(a,sum(n), k))

  ### Augment observations for each cell
  for (i in 1:a){
    for (j in 1:N){
      if (i==1){
        if ( j<= n[1] ) {
          newtrt<-alltrt[,1:n[1]]
          total<-ncol(newtrt)
          jj<-j
        }

        if  (j>=n[1]+1) {
          newtrt<-cbind(alltrt[,1:n[1]], alltrt[, j])
          total<-jj<- ncol(newtrt)
        }
      }

      if (i>1) {
        sumni=sum(n[1:i]);  sumniminus1=sum(n[1:(i-1)])
        if ((j<=sumni)& (j>=sumniminus1+1) ) {
          newtrt<-alltrt[,(sumniminus1+1): sumni]
          total<- ncol(newtrt)
          jj<- j-sumniminus1
        } else {
          newtrt<-cbind(alltrt[,(sumniminus1+1): sumni], alltrt[,j] )
          total<-jj<-ncol(newtrt)
        }
      }

      newtrt[3, ]<-rank(newtrt[2, ])
      flag<-((jj==total)& (jj>n[i])& c(rep(TRUE, total-1), FALSE)  )  | (jj<=n[i])
      if ((jj==total) & (jj>n[i]) ) {
        newtrt[3, -jj]<- rank(newtrt[2, -jj])
        total<-total-1
      }
      target<-newtrt[3, jj ];  trunc_target=trunc(target)
      newtrt<-newtrt[, flag];

      if (trunc_target <= ((k-1)/2) )
      {ordk=order(newtrt[3, ])[1:k]
      }
      if (trunc_target > (total- ((k-1)/2)))
      {ordk= order(total-newtrt[3, ])[1:k]
      }
      if ((trunc_target  <=(total-(k-1)/2 ) ) & (trunc_target >((k-1)/2) )   )
      {ordk=order((abs(newtrt[3,]-trunc_target) ))[1:k]
      }
      psudo[i,j, ]<-newtrt[1, ordk]
      index[i,j, ]<-seq(1, total)[ ordk]

    }  #end of j
  }   #end of i

  list(psudo=psudo, index=index)
}

#' Index in one vector mapped to treatment and observation index
#'
#' Function mapindex() maps the 1-d index r=1,...,N to 2-d
#' index i=1, ...a, j=1, ..., ni. Generally the covariate values
#' from all treatments are stored together in one vector and
#' r=1,...,N enumerates the values. For any integer between 1 and N,
#' mapindex tells which treatment the rth value belongs to, and which
#' observation in the identified treatment.
#'
#' @param r an integer between 1 and sum(n).
#' @param n a vector of the sample sizes.
#'
#' @return the 2-d index, where the first gives which treatment
#'         the value belongs to and the second gives which observation
#'         in that treatment.
#'
#' @export
#'
#' @examples
#' r=5; n=c(7, 8); mapindex(r, n)
#' r=7; n=c(7, 8); mapindex(r, n)
#' r=9; n=c(7, 8); mapindex(r, n)

mapindex=function(r, n){
  sumn=cumsum(n);  rem=r-sumn
  imap=sum(rem>0)+1
  if (imap<2) jmap=r else jmap=rem[imap-1]
  c(imap, jmap)
}


#' Starting and ending position in a vector
#'
#' Function trt_position() gives the starting and ending
#' index of covariate values in the i1th group if all the
#' covariate values from all treatment groups are together in a vector.
#' E.g., covariate values in group 1 start from 1st value to the
#' n1 th value; those in group 2 start from n1+1 and end at (n1+n2)th
#' value.
#' This function is for retrieving the position of an observations
#' when the covariate values from all treatments are stored together
#' in one vector.
#'
#' @param i1 an integer between 1 and length(n).
#' @param n  the vector of sample sizes.
#'
#' @export
#'
#' @examples
#' i = 2; n=c(7, 8); trt_position(i, n)
#'
trt_position=function(i1, n){
  if (i1==1) lower=1 else lower=sum(n[1:(i1-1)])+1
  upper=sum(n[1:i1])
  c(lower, upper)
}




#' Image structural similarity measure PSSIM based on hypothesis test
#'
#' PSSIM_snow computes image structural similarity PSSIM of Wang, Maldonado and Silwal (2011) using
#'             parallel programming.
#'
#' @param A a grayscale image stored as a matrix.
#' @param A1 grayscale image stored as a matix. Same dimension as A.
#' @param nprocess number of cores (workers) to use for parallel computation.
#'         Note:
#'          In personal computer, nprocess =detectCores() is good to use.
#'          On cluster machine, nprocess need to be specified to a number that is
#'          no more than its number of cores (for courtesy)
#' @param b  Number of columns in each block. Suggest to use default value 64.
#' @param a  Number of rows in each block. Suggest to use default value 2.
#' @param vs Block shift size. Suggest to use default value 32.
#' @param wavecoeff logical of whether the input matrices are wavelet coefficients.
#'            Currently, wavelet version is not implemented.
#'            This parameter is a placeholder for future implementation.
#' @param cs dividing factor to split index.
#' @param dyn logical, whether dynamic scheduling should be used.
#'
#' @return: Image structural similarity based on PSSIM. The value is in [0,1]
#'          with values close to 0 meaning the two images are different
#'          and values close to 1 meaning the two iamges are similar.
#'
#' @export
#'
#' @examples
#'   A=miniimagematrix$A
#'   B=miniimagematrix$B
#'   # see it with image(A, axes=FALSE, col  = gray((0:255)/256) )
#'   PSSIM_snow(A, B, nprocess=2)
#'
#' @references
#' Haiyan Wang, Diego Maldonado, and Sharad Silwal  (2011).   A Nonparametric-Test-Based
#' Structural Similarity Measure for Digital Images.  Computational Statistics and Data Analysis.  55: 2925-2936. Doi:10.1016/j.csda.2011.04.021


PSSIM_snow=function(A, A1, nprocess=min(8, parallel::detectCores()),
                    b=64,  a=2, vs=32, wavecoeff=FALSE,cs=2,dyn=FALSE){
  requireNamespace('parallel')
  # Define a function to compute blockwise p-values with
  # 'snow' parallel programming
  paraPvaluesnow=function(ichunks){
    A=as.matrix(A); A1=as.matrix(A1)
    i=1;   p=numeric()
    contrast=numeric(); luminance=numeric()
    for (i in ichunks) {
      j=1
      p1=numeric()
      while ((j-1)*vs+b<=nrow(A) ) {
        blockA=A[(j-1)*vs+1:b, i+(0:(a-1))*ncol(A)/a  ]
        blockA1=A1[(j-1)*vs+1:b, i+(0:(a-1))*ncol(A)/a  ]
        trtindex=kronecker(seq(a), rep(1,b) )
        dat=data.frame(Y=c(blockA- blockA1), X=c(blockA), trt=trtindex)
        dat2=data.frame(Y=c(blockA- blockA1), X=c(blockA1), trt=trtindex)
        # Consider both A-A1 vs A and A-A1 vs A1 because if A=A1+e,
        # then A-A1 is indept of e but A-A1 is not indept of A.
        pv1=NPtest_indept(dat, k=7)$pvalue
        pv2=NPtest_indept(dat2, k=7)$pvalue
        p1=c(p1, max(pv1, pv2) )
        # Pick the one with bigger pvalue since the test of A-A1 vs A would
        # tend to reject H0 due to depenendence of e with A.
        j=j+1
      }
      p=rbind(p, p1)
    }
    p
  }


  cls=parallel::makeCluster(nprocess)
  chunk = ncol(A)-(a-1)*ncol(A)/a
  ichunks=parallel::splitIndices(chunk,chunk/cs)
  parallel::clusterExport(cls,c('A','A1','b','a','vs','wavecoeff','cs'),envir=environment())
  parallel::clusterExport(cls,c('NPtest_indept','makepseudo','mapindex','trt_position'))

  if (!dyn) {
    out <- parallel::clusterApply(cls,ichunks,paraPvaluesnow)
  } else {
    out <- parallel::clusterApplyLB(cls,ichunks,paraPvaluesnow)
  }
  p=matrix(unlist(out),byrow=TRUE)
  luminance=(2*A*A1+0.001  )/(  A^2 + A1^2 +0.001  )
  result=ifelse(wavecoeff==TRUE, mean(p>0.01, na.rm=TRUE),
                mean(p>0.01, na.rm=TRUE)* mean(luminance) )
  parallel::stopCluster(cls)
  result
}

##
