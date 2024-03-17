#' calculates Coancestries and FST from demographic parameters
#' 
#' @usage  get.efst(Nindsperpop,npops,migmat,mutrate,ngen,theta=NULL)
#'
#' @param Nindsperpop vector of population sizes. If only a scalar, 
#' all populations assumed to have the same sizes
#' @param npops number of populations
#' @param migmat migration matrix
#' @param mutrate mutation rate
#' @param ngen number of generations
#' @param theta initial coancestries, assumed to be 0 if not specified
#'
#' @value Coan a matrix of Coancestries
#' @value EFst a matrix of FST, or mean coancestries relative to 
#' the mean off-diagonal elements
#' 
#' @references 
#'
#' \href{https://academic.oup.com/genetics/article/206/4/2085/6072590}{Weir, BS and Goudet J. 2017} A Unified Characterization 
#' of Population Structure and Relatedness. Genetics (2017) 206:2085 
#'
#' \href{https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010871}{Goudet J. and Weir, B.S.} An allele-sharing, moment-based estimator 
#' of global, population-specific and population-pair  FST under a 
#' general model of population structure.  PLOS Genetics (2023) 19:e1010871
#'
#' @export
#########################################

get.efst<-function(Nindsperpop,npops,migmat,mutrate,ngen,theta0=NULL){
na<-2
if (length(Nindsperpop)==1) Nindsperpop=rep(Nindsperpop,npops)

b <- (1-mutrate*(na-1)/na)^2
#  b<-(1-mut)^2
#ID <- diag(npops)
# calculate theta, the coancestry coefficient
# and calculate expected betas, including the expected population-specific Fst
# initialize both at 0, and only save the data at the end of the number of generations
J<-matrix(rep(1,npops^2),nrow=npops)

if(is.null(theta0)) theta <- matrix(numeric(npops^2), nrow = npops) else theta<-theta0
for (i in 1:ngen){
  theta<-b*(migmat%*%((J-diag(1/2/Nindsperpop))*theta+diag(1/2/Nindsperpop))%*%t(migmat))
}
#theta

Mb<-mean(mat2vec(theta))
return(list(Coan=theta,EFst=(theta-Mb)/(1-Mb)))
}

