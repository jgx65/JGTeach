#' Estimates the standard kinship matrix
#'
#' Estimates the standard kinship matrix as
#' an average of ratios
#'
#' @usage kC0(dos,inb=FALSE,matching=FALSE)
#' @param dos either a dosage matrix or a matching matrix
#' @param inb whether to report inbreeding coefficient
#'  or self coancestries on the main diagonal
#' @param matching TRUE if dos is a matching matrix (nixni),
#'  FALSE if it is a dosage (nixnl) matrix
#'
#' @return the standard, ratio of averages, kinship matrix
#'
#' @seealso \code{\link[hierfstat]{beta.dosage}} for allele sharing kinships;
#' \code{\link[gaston]{GRM}} for standard kinship as an average of ratio (default)
#'
#' @export
#'
Kc0<-function(dos,inb=FALSE,matching=FALSE){
  if(!matching){
    tmp<-hiersftat::beta.dosage(dos,inb=inb,Mb=TRUE)
    Mij<-with(tmp,betas*(1-MB)+MB)
    }
  else Mij<-dos
  MT<-mean(Mij)
  n<-nrow(MT)
  CM<-diag(n)-1/n
  CM%*%(1/(1-MT)*(Mij-MT))%*%CM
}

