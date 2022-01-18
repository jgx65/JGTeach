###########################################################################
#' @title drops alleles along pedigree
#'
#' @description Given a pedigree, drops alleles along the pedigree using the given recombination map.
#' Could be seeded with genotypes for the founder individuals, otherwise, uses gene dropping.
#'
#'
#'
#' @usage drop.along.ped(ped=ped,founders.genotypes=NULL,nloc=10000,maplength=20)
#'
#' @param ped a pedigree with at least three columns, the first for the individual id,
#' the second for the dam id and the third for the sire id.
#' dams and sires of founders are entered as NAs
#' founders/unknown parents must appear first in the pedigree
#' @param founders.genotypes (default NULL) if given, an array [nf*2,nloc] containing the (0 or 1)
#' of the nf founders at nloc loci. if only one parent is missing,
#' the missing parent should be considered as a founder
#' @param nloc number of equally spaced loci
#' @param maplength map length in Morgans
#'
#'
#'
#' @return a [ni,nl,2] array of genotypes
#'
#' @details
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references
#'
#' @examples
#' \dontrun{}
#' @export

drop.along.ped<-function(ped=ped,founders.genotypes=NULL,nloc=10000,maplength=20,ndigits=2){
  #assumes a pedigree (3 columns, ind, dam,sire)
  #if gene dropping, no founder genotypes needed.
  #else a data frame dat with founders genotypes
  #encoded as gametes (0,1)
  #make it a c(nfounders*2,nl) array of gametes
  #dam and sire of founders are entered as NA, same for unknown parents
  #maplength in Morgans
  drawgam<-function(dat){
    gs<-sample(2,1) # from which gamete to start
    nbxo<-rpois(1,maplength) # nb cross overs
    nbseg<-nbxo+1 # nb segments
    g<-(1:nbseg+gs)%%2 +1 # which chrom to read
    posxo<-sort(sample(nloc-1,nbxo))+.1 #position cross over
    posxof<-c(floor(posxo),nloc) #end segments
    posxos<-c(1,ceiling(posxo)) #start segments
    res<-unlist(lapply(1:nbseg,function(x) dat[posxos[x]:posxof[x],g[x]]))
    return(res)
  }

  founders<-which(is.na(ped[,2]) & is.na(ped[,3]))
  nfounders<-length(founders)
    noffs<-dim(ped)[1]-nfounders

  if(is.null(founders.genotypes)){
    ufounders<-unique(founders)
  }
  else {
    nloc<-dim(founders.genotypes)[2]
  }

  genos<-array(dim=c((nfounders+noffs),nloc,2))

  if(is.null(founders.genotypes)){
  if (dim(founders.genotypes)[1]<nfounders*2) stop("not enough fouders gametes (at least 2 x nfounders). Exiting")
    genos[1:nfounders,,1]<-1:nfounders
  genos[1:nfounders,,2]<-(1:nfounders)+nfounders
  }
  else {
  genos[1:nfounders,,1]<-founders.genotypes[1:nfounders,]
  genos[1:nfounders,,2]<-founders.genotypes[nfounders+1:nfounders,]
  }
  seqoffs<-(nfounders+1):(nfounders+noffs)
  dam<-match(ped[,2],ped[,1])
  sire<-match(ped[,3],ped[,1])
  mumid<-NA
  dadid<-NA
  for (io in seqoffs){
    if (!is.na(dam[io]) & !is.na(sire[io])){ #parents known
      mumid<-dam[io]#match(ped[io,2],ped[,1])
      dadid<-sire[io]#match(ped[io,3],ped[,1])
      mum<-genos[mumid,,]
      dad<-genos[dadid,,]
      genos[io,,1]<-drawgam(mum)
      genos[io,,2]<-drawgam(dad)
    }
  }
  genos
}
######################################
###########################################################################
#' @title Estimates gold kinship
#'
#' @description Estimates gold kinship from the genotypes generated with \code{\link{drop.along.ped}}
#' assuming the founders are unrelated
#'
#' @usage kinship.gold(dat)
#'
#' @param dat an array (niXnlX2) produced by \code{\link{drop.along.ped}} containing
#' which alleles each individual has inherited from the founders
#'
#' @return 'gold' kinship matrix, including self kinship on the diagonal
#'
#' @details
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references
#'
#' @examples
#' \dontrun{}
#' @export

kinship.gold<-function(dat){
  #kinship from gene dropping on the pedigree
  n<-dim(dat)[1]
  nl<-dim(dat)[2]
  k<-matrix(numeric(n^2),ncol=n)
  for (i in 1:n){
    a1<-dat[i,,1]
    a2<-dat[i,,2]
    for (j in 1:i){
      a3<-dat[j,,1]
      a4<-dat[j,,2]
      k[i,j]<-(sum(a1==a3)+
                 sum(a1==a4)+
                 sum(a2==a3)+
                 sum(a2==a4))/4/nl
      k[j,i]<-k[i,j]
    }}
  k
}
##########################################
kinship.goldm<-function(dat,ncores=1){
#  if (ncores >1) library(parallel)

  psi<-function(i){
    psis<-numeric(n)
    a1<-dat[i,,1]
    a2<-dat[i,,2]
    for (j in 1:i){
      a3<-dat[j,,1]
      a4<-dat[j,,2]
      psis[j]<-(sum(a1==a3)+
                  sum(a1==a4)+
                  sum(a2==a3)+
                  sum(a2==a4))/4/nl
    }
    psis
  }
  #kinship from gene dropping on the pedigree
  n<-dim(dat)[1]
  nl<-dim(dat)[2]
  k<-matrix(numeric(n^2),ncol=n)
  if (ncores==1)
    res<-lapply(1:n,function(x) psi(x))
  else
    res<-parallel::mclapply(1:n,function(x) psi(x),mc.cores=ncores)
  for (i in 1:n)
    k[,i]<-res[[i]]
  for (i in 2:n)
    for (j in 1:(i-1))
      k[i,j]<-k[j,i]
  k
}
library(compiler)
kin.goldm<-cmpfun(kinship.goldm)
