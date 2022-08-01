#' @title convert Identity By Descent gametes to dosage data
#'
#' @description Given gametes, the output of \code{\link{drop.along.ped}} on a pedigree and
#' and n.founders, the number of founders in the pedigree,
#' return a matrix ni x nl of dosages after having mapped
#' the 2n.founders alleles at each locus to 0 or 1.
#' The gametes of the founders are given in g.founders or
#' frequencies at each locus are drawn from a \code{\link{beta}}
#' distribution with parameters b1 and b2.
#'
#' @usage ibd2dos(gametes,n.founders,g.founders=NULL,b1=1,b2=10)
#'
#' @param gametes a nixnlx2 array of gametes created by \code{\link{drop.along.ped}}
#' @param n.founders the number of founders in the pedigree
#' @param g.founders a (2*nfounders)xnl array with the bi-allelic gametes of the founders (default to NULL)
#' @param b1 the 1st parameter of the \code{\link{beta}} distribution
#' @param b2 the 2nd parameter of the \code{\link{beta}} distribution
#'
#' @return a dosage matrix with ni rows and nl columns
#'
#' @details The gametes of  the founders can be entered as an (2*nfounders) x nl array. This insures
#' any haplotype structure in the founders can be kept.  If g.founders is NULL, allelic frequencies
#' in the founders are drawn from a beta distribution with parameters b1 (default=1)
#' and b2 (default=10).
#'
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#'
#' @examples
#' \dontrun{
#' ped1<-buildped.2sexes(founders.m=10,founders.f=10,fert=1,death.rate=0.5,breed.prop=0.2,n.tstep=5)
#' ibd.ped1<-drop.along.ped(ped1,nloc=5000,maplength=20)
#' dos.ped1<-ibd2dos(ibd.ped1,n.founders=20)
#' }
#'
#' @export
#'
ibd2dos<-function(gametes,n.founders,g.founders=NULL,b1=1,b2=10){
  nf2<-2*n.founders
  nl<-dim(gametes)[2]
  tmp<-dim(gametes)[3]
  if ((is.na(tmp)) | (tmp!=2))
    stop("Error, gametes must be an array nixnlx2 \n as produced by the drop.along.ped function. Exiting")
  if((n.founders <1) | (n.founders >dim(gametes)[1]))
    stop("Error, n.founders must be the number of founders in the pedigree \n used to generate gametes. Exiting")
  if((!is.null(g.founders)) & (dim(g.founders)[1]!=nf2) & (dim(g.founders)[2]!=nl))
    stop("Error, dim of g.founders must match be 2n.founders x nl. Exiting")

  if (is.null(g.founders)) nA<-round(rbeta(nl,b1,b2)*nf2,digits=0)

  for (il in 1:nl){
    if (!is.null(g.founders))
       alt<-which(g.founders[,il]==1)
    else
       alt<-sample(nf2,size=nA[il])
    x<-gametes[,il,] %in% alt
    gametes[,il,][x]<-1
    gametes[,il,][!x]<-0
  }

  return(gametes[,,1]+gametes[,,2])
}
