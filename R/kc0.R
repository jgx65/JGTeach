#' Estimates the standard kinship matrix
#'
#' Estimates the standard kinship Kc0 matrix as
#' a ratio of averages
#'
#' @usage KC0(dos,inb=FALSE,matching=FALSE)
#' @param dos either a dosage matrix or a matching matrix
#' @param inb whether to report inbreeding coefficient
#'  or self coancestries on the main diagonal
#' @param matching TRUE if dos is a matching matrix (nixni),
#'  FALSE (default) if it is a dosage (nixnl) matrix
#'
#' @return the standard, ratio of averages, kinship matrix
#'
#' @seealso \code{\link[hierfstat]{matching} for allele sharing,
#' \code{\link[hierfstat]{beta.dosage}} for allele sharing kinships;
#' \code{\link[gaston]{GRM}} for standard kinship as an average of ratio (default)
#'
#' @export

Kc0<-function(dos,inb=FALSE,matching=FALSE){
  if(!matching){
    tmp<-hierfstat::matching(dos)
    }
  else Mij<-dos
  MT<-mean(Mij)
  n<-nrow(Mij)
  CM<-diag(n)-1/n
  CM%*%(1/(1-MT)*(Mij-MT))%*%CM
}

