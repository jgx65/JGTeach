#######################
#' @title Estimates FUNI as ratio of averages
#'
#' @description Estimates FUNI as ratio of averages on a \code{\link[gaston]{bed.matrix}} object
#'
#' @usage get.funiw(dos)
#'
#' @param dos a \code{\link[gaston]{bed.matrix}} object
#'
#' @details
#'
#'
#' @return a list with components
#' \itemize{
#' \item het the sum of loci expected heterozygosities
#' \item Funi ratio of averages \eqn{F_{uni}} estimates
#' }
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references
#'
#' \href{https://doi.org/10.1038/s41437-021-00471-4}{Zhang QS, Goudet J and Weir BS 2022}
#'  Rank-invariant Estimation of Inbreeding coefficients. Heredity
#'
#' @details
#'
#'
#'
#' @examples{\dontrun{}}
#'
#' @export

get.funiw<-function(dos){

  if(!class(dos)[[1]]=="bed.matrix") stop("Argument must be of class bed.matrix. Exiting")
  p<-dos@p
  het<-2*p*(1-p)
  num<-apply(as.matrix(dos),1,function(x) sum(x^2-(1+2*p)*x+2*p^2))
  den<-sum(het)
  return(list(het=den,Funi=num/den))
}

#######################
#' @title Estimates FUNI (FIII) as average of ratios
#'
#' @description Estimates FUNI as average of ratios on a \code{\link[gaston]{bed.matrix}} object
#'
#' @usage get.funiu(dos)
#'
#' @param dos a \code{\link[gaston]{bed.matrix}} object
#'
#' @details
#'
#'
#' @return a list with components
#' \itemize{
#' \item nloc the number of typed loci for each individual
#' \item Funi average of ratios \eqn{F_{uni}} estimates
#' }
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#'
#' @references{
#'
#' \href{https://doi.org/10.1038/s41437-021-00471-4}{Zhang QS, Goudet J and Weir BS 2022}
#'  Rank-invariant Estimation of Inbreeding coefficients. Heredity
#'
#' Original definitions of FUNI and FIII are found in:
#'
#'      Yang J, Benyamin B, McEvoy BP, Gordon S, Henders AK, Nyholt DR,
#' Madden PA, Heath AC, Martin NG, Montgomery GW, Goddard ME,
#' Visscher PM. 2010.  Common SNPs explain a large proportion of the
#' heritability for human height.  Nat Genet. 42(7):565-9. Epub 2010
#' Jun 20.
#'
#' Yang, J., Lee, S. H., Goddard, M. E. & Visscher, P. M.  GCTA: a
#' tool for genome-wide complex trait analysis.  American journal of
#' human genetics 88, 76-82 (2011).
#'
#' }
#'
#' @details
#'
#' This is the original definition of FIII or FUNI given in Yang et al. (2010, 2011)
#'
#'
#' @export
get.funiu<-function(dos){
  if(!class(dos)[[1]]=="bed.matrix") stop("Argument must be of class bed.matrix. Exiting")
  p<-dos@p
  pol<-which(dos@snps$maf>0.0)
  x<-as.matrix(dos[,pol])
  p<-p[pol]
  het<-2*p*(1-p)
  nl<-dim(x)[2]

  if(sum(is.na(x))>0){
    na <- matrix(rep(1,prod(dim(x))),ncol=ncol(x))
    ina<-which(is.na(x))
    na[ina]<-0
    nls<-rowSums(na)
  }
    else nls<-rep(nl,dim(x)[1])
    tmp1<-sweep(x,2,(1+2*p),FUN="*")
    tmp2<-sweep(-tmp1,2,2*p^2,FUN="+")
    res<-rowMeans(sweep(x^2+tmp2,2,het,FUN="/"),na.rm=TRUE)
    return(list(nloc=nls,Funi=res))
}

