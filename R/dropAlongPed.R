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

drop.along.ped<-function(ped=ped,founders.genotypes=NULL,nloc=10000,maplength=20){
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
     genos[1:nfounders,,1]<-1:nfounders
     genos[1:nfounders,,2]<-(1:nfounders)+nfounders
  }
  else {
   if(dim(founders.genotypes)[1]<nfounders*2)
     stop("not enough founders gametes (at least 2 x n founders). Exiting")
   else{
    genos[1:nfounders,,1]<-founders.genotypes[1:nfounders,]
    genos[1:nfounders,,2]<-founders.genotypes[nfounders+1:nfounders,]
    }
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
#' @title Estimates realized (gold) kinship
#'
#' @description Estimates realized kinship from the genotypes generated with \code{\link{drop.along.ped}}
#' assuming the founders are unrelated
#'
#' @usage kinship.gold(dat)
#'
#' @param dat an array (niXnlX2) produced by \code{\link{drop.along.ped}} containing
#' which alleles each individual has inherited from the founders
#'
#' @return 'gold' kinship matrix,  including self kinship on the diagonal
#'
#' @details 'gold' kinship is the realized proportion of shared alleles between individuals,
#' using the founders as reference (each founder is assigned two unique alleles at each position).
#' the realized inbreeding coefficient of an individual \eqn{F_i} is obtained from its self-kinship \eqn{k_{ii}} as
#' \eqn{F_i=2 k_{ii}-1}
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

kin.goldm<-compiler::cmpfun(kinship.goldm)


############################################################
###########################################################################
#' @title Calculates the 9 Jaccard's Delta coefficients
#'
#' @description Calculates the 9 Jaccard's Delta coefficients from an array of genotypes
#' generated with \code{\link{drop.along.ped}}
#'
#'
#'
#'
#' @usage jaccard(dat)
#'
#' @param dat an nixnlx2 array of genotypes ( ni individuals, nl loci) generated with \code{\link{drop.along.ped}}
#'
#'
#'
#' @return a list of 9 matrices, the lower half of which contain
#' the corresponding 9 condensed Jaccard's coefficients
#'
#' @details given two individuals, X with alleles a and b,
#' and Y with alleles c and d,
#' \itemize{
#' \item{\eqn{\Delta_1}}{ a==b==c==d}
#' \item{\eqn{\Delta_2}}{ a==b, c==d, a!=c}
#' \item{\eqn{\Delta_3}}{ a==b==c, a!=d}
#' \item{\eqn{\Delta_4}}{ a==b, a=!c, a!=d, c!=d}
#' \item{\eqn{\Delta_5}}{ a==c==d, a!=b}
#' \item{\eqn{\Delta_6}}{ c==d, c!=a, c!=b, a!=b}
#' \item{\eqn{\Delta_7}}{ a==c, b==d, a!=b}
#' \item{\eqn{\Delta_8}}{ a==c, a!=b,d, b!=d}
#' \item{\eqn{\Delta_9}}{ a!=b, a!=c, a!=d, b!=c, b!=d, c!=d}
#' }
#' The sum of all 9 coefficients for any pair of individuals is 1
#'
#' Coancestry \eqn{\theta} is
#' \deqn{\theta=\Delta_1+(\Delta_3+\Delta_5+\Delta_7)/2
#' +\Delta_8/4}
#'
#' Inbreeding FX for individual X is:
#' \deqn{F_X=\Delta_1+\Delta_2+\Delta_3+\Delta_4}
#'
#' while inbreeding FY for individual Y is:
#' \deqn{F_Y=\Delta_1+\Delta_2+\Delta_5+\Delta_6}
#'
#' ks (kappas) defined for non inbred individuals only
#' in which case
#'\itemize{
#' \item{\eqn{k_0}}{ \eqn{\Delta_9}}
#' \item{\eqn{k_1}}{ \eqn{\Delta_8}}
#' \item{\eqn{k_2}}{ \eqn{\Delta_7}}
#'}
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references
#'
#' @examples
#' \dontrun{}
#' @export

jaccard<-function(dat){
  #Jaccard's Deltas from gene dropping on the pedigree
  n<-dim(dat)[1]
  nl<-dim(dat)[2]
  ds<-list(length=9);
  for (k in 1:9) ds[[k]]<-matrix(numeric(n^2),ncol=n)
  for (i in 2:n){
    a1<-dat[i,,1] #a
    a2<-dat[i,,2] #b
    for (j in 1:(i-1)){
      a3<-dat[j,,1] #c
      a4<-dat[j,,2] #d
      ds[[1]][i,j]<-sum((a1==a2) & (a3==a4) & (a1==a3))/nl
	  ds[[1]][j,i]<-ds[[1]][i,j]
      ds[[2]][i,j]<-sum((a1==a2) & (a3==a4) & (a1!=a3))/nl
	  ds[[2]][j,i]<-ds[[2]][i,j]
      ds[[3]][i,j]<-(sum((a1==a2) & (a1==a3) & (a3!=a4))+
                       sum((a1==a2) & (a1==a4) & (a3!=a4)))/nl
      ds[[4]][i,j]<-sum((a1==a2) & (a1!=a3) & (a3!=a4) & (a1!=a4))/nl
      ds[[5]][i,j]<-(sum((a3==a4) & (a3==a1) & (a1!=a2))+
                       sum((a3==a4) & (a3==a2) & (a1!=a2)))/nl
      ds[[6]][i,j]<-sum((a3==a4) & (a3!=a1) & (a3!=a2) & (a1!=a2))/nl
	  ds[[3]][j,i]<-ds[[5]][i,j]
	  ds[[5]][j,i]<-ds[[3]][i,j]
	  ds[[4]][j,i]<-ds[[6]][i,j]
	  ds[[6]][j,i]<-ds[[4]][i,j]
      ds[[7]][i,j]<-(sum((a1==a3) & (a2==a4) & (a1!=a2))+
                       sum((a1==a4) & (a2==a3) & (a1!=a2)))/nl
	  ds[[7]][j,i]<-ds[[7]][i,j]
      ds[[8]][i,j]<-(sum((a1==a3) & (a1!=a2) & (a3!=a4) & (a2!=a4))+
                       sum((a1==a4) & (a1!=a2) & (a3!=a4) & (a2!=a3))+
                       sum((a2==a3) & (a1!=a2) & (a3!=a4) & (a1!=a4))+
                       sum((a2==a4) & (a1!=a2) & (a3!=a4) & (a1!=a3)))/nl
	ds[[8]][j,i]<-ds[[8]][i,j]				   
      ds[[9]][i,j]<-sum((a1!=a2) & (a3!=a4) & (a1!=a3) & (a1!=a4) & (a2!=a3) & (a2!=a4))/nl
	  ds[[9]][j,i]<-ds[[9]][i,j]
    }
      ds[[1]][i,i]<-sum(a1==a2)/nl
      ds[[7]][i,i]<-sum(a1!=a2)/nl	  
	
	}
  ds
}
############################################################
