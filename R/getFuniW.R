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
#' \item het the sum of loci heterozygosities
#' \item res a ni x 2 matrix, the first column being the numerator
#'  of \eq{F_{uni}} and the second being \eq{F_{uni}}
#' }
#'
#' @export

get.funiw<-function(dos){

  p<-dos@p
  het<-dos@snps$hz
  num<-apply(as.matrix(dos),1,function(x) sum(x^2-(1+2*p)*x+2*p^2))
  den<-sum(het)
  return(list(het=den,res=cbind(num,num/den)))
}
