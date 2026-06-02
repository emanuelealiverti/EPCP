########## R function: ZOSull ##########

# For creation of O'Sullivan-type Z matrices.
#
# Author: Matt Wand <Matt.Wand@uts.edu.au>
# Source: Course materials for 37457 Advanced Bayesian Methods,
#         University of Technology Sydney (UTS).
# Last changed: 15 JAN 2015

#' @export
ZOSull <- function(x,range.x,intKnots,drv=0)
{
   # Check legality of range.x.
  
   if (!missing(range.x))
   {
      if (length(range.x)!=2) stop("range.x must be of length 2.")
      if (range.x[1]>range.x[2]) stop("range.x[1] exceeds range.x[1].")
      if (range.x[1]>min(x)) stop("range.x[1] must be <= than min(x).")
      if (range.x[2]<max(x)) stop("range.x[2] must be >= than max(x).")

   }
   
   if (drv>2) stop("splines not smooth enough for more than 2 derivatives")

   library(splines)

   # Set defaults for `range.x' and `intKnots'

   if (missing(range.x))
      range.x <- c(1.05*min(x)-0.05*max(x),1.05*max(x)-0.05*min(x))
   
   if (missing(intKnots))
   {
      numIntKnots <- min(length(unique(x)),35)
      intKnots <- quantile(unique(x),seq(0,1,length=
                  (numIntKnots+2))[-c(1,(numIntKnots+2))])
   }
   numIntKnots <- length(intKnots) 

   # Obtain the penalty matrix.

   allKnots <- c(rep(range.x[1],4),intKnots,rep(range.x[2],4)) 
   K <- length(intKnots) ; L <- 3*(K+8)
   xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+ 
              rep(allKnots,each=3)[-c(1,2,L)])/2
   wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
   Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                     outer.ok=TRUE)$design  
   Omega     <- crossprod(Bdd*wts,Bdd)     

   # Use the spectral decomposition of Omega to obtain Z.

   svdOmega <- svd(Omega) 
   indsZ <- 1:(numIntKnots+2)
   UZ <- svdOmega$u[,indsZ] 
   LZ <- t(t(UZ)/sqrt(svdOmega$d[indsZ]))

   # Perform stability check.   

   indsX <- (numIntKnots+3):(numIntKnots+4)
   UX <- svdOmega$u[,indsX]   
   L <- cbind(UX,LZ)
   stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))          
   if (sum(stabCheck^2) > 1.0001*(numIntKnots+2))
       print("WARNING: NUMERICAL INSTABILITY ARISING\\
              FROM SPECTRAL DECOMPOSITION")

   # Obtain B and post-multiply by LZ matrix to get Z.

   B <- spline.des(allKnots,x,derivs=rep(drv,length(x)),
                     outer.ok=TRUE)$design  
   
   Z <- B%*%LZ

   # Add the `range.x' and 'intKnots' as attributes
   # of the return object.

   attr(Z,"range.x") <- range.x
   attr(Z,"intKnots") <- intKnots

   # Return Z matrix with 2 attributes.

   return(Z)
}

########## End of ZOSull ##########

