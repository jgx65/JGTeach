#' @title convert Identity By Descent gametes to dosage data
#'
#' @description Given gametes, the output of \code{\link{drop.along.ped}} on a pedigree and
#' and n.founders, the number of founders in the pedigree,
#' return a matrix of dosages after having mapped at random
#' the 2n.founders alleles at each locus to 0 or 1.
#' The probability of a 1 (the founders frequency of 1) is drawn from a
#' \code{\link{beta}} distribution with parameters b1 and b2.
#'
#' @usage ibd2dos(gametes,n.founders,b1=1,b2=10)
#'
#' @param gametes a nixnlx2 array of gametes created by \code{\link{drop.along.ped}}
#' @param n.founders the number of founders in the pedigree
#' @param b1 the 1st parameter of the \code{\link{beta}} distribution
#' @param b2 the second parameter of the \code{\link{beta}} distribution
#'
#' @value a dosage matrix with ni rows and nl columns
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#'
#' @examples
#' \dontrun{
#' ped1<-buildped.2sexes(founders.m=10,founders.f=10,fert=1,death.rate=0.5,breed.prop=0.2,n.tstep=5)
#' ibd.ped1<-drop.along.ped(ped1,nloc=5000,maplength=20)
#' dos.ped1<-ibd2dos(ibd.ped1,n.founders=20)
#' kAS.ped1<-hierfstat::beta.dosage(dos.ped1)
#' kc0.ped1<-Kc0(dos.ped1)
#' kgold.ped1<-kinship.gold(ibd.ped1)
#' ARM.ped1<-pedARM(ped1[,2],ped1[,3])
#'
#' round(cor(cbind(hierfstat::mat2vec(ARM.ped1)/2,
#' hierfstat::mat2vec(kgold.ped1),
#' hierfstat::mat2vec(kAS.ped1),
#' hierfstat::mat2vec(kc0.ped1))),digits=3)
#' }
#'
#' @export
#'
ibd2dos<-function(gametes,n.founders,b1=1,b2=10){
  tmp<-dim(gametes)[3]
  if ((is.na(tmp)) | (tmp!=2))
    stop("Error, gametes must be an array nixnlx2 \n as produced by the drop.along.ped function. Exiting")
  if((n.founders <1) | (n.founders >dim(gametes)[1]))
    stop("Error, n.founders must be the number of founders in the pedigree \n used to generate gametes. Exiting")
  nf2<-2*n.founders
  for (il in 1:dim(gametes)[[2]]){
    nA<-rbeta(1,b1,b2)*nf2
    alt<-sample(nf2,size=round(nA,digits=0))
    x<-gametes[,il,] %in% alt
    gametes[,il,][x]<-1
    gametes[,il,][!x]<-0
  }

  return(gametes[,,1]+gametes[,,2])
}
