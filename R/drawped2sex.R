###########################################################################
#' @title Builds a two-sex pedigree
#'
#' @description Builds a two-sex pedigree from the passed arguments
#'
#' @usage buildped.2sexes(founders.m=50,founders.f=50,fert=2.0,death.rate=0.2,breed.prop=0.5,male.prop=0.5,n.tstep=10)
#'
#' @param founders.m numbers of male founders
#' @param founders.f number of female founders
#' @param fert female fertility rate per time interval
#' @param death.rate mortality rate per time interval
#' @param breed.prop proportion of breeding males
#' @param male.prop proportion of males (sex ratio)
#' @param n.tstep number of time step
#'
#'
#' @return a pedigree data frame with 4 columns [ind id, dam id, sire id, sex] and as many rows as individuals
#'
#' @details
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references
#'
#' @examples
#' \dontrun{
#'
#' #build a pedigree
#' library(hierfstat)
#' test<-buildped.2sexes(founders.m=10,founders.f=10,fert=1,death.rate=0.5,breed.prop=0.2,n.tstep=5)
#' #range of kinship coefficients
#' range(mat2vec(grm2kinship(pedARM(test[,2],test[,3]))))
#' #range of inbreeding coefficients
#' range(diag(grm2kinship(pedARM(test[,2],test[,3]))))
#' #plot the pedigree
#' library(kinship2)
#' with(test,plot(pedigree(ind,sire,dam,sex)))
#' }
#' @export

#######################################
buildped.2sexes<-function(founders.m=50,founders.f=50,fert=2.0,death.rate=0.2,breed.prop=0.5,male.prop=0.5,n.tstep=10){
  draw.parents<-function(dams,sires,fert,prop.breed,sanc){

    nbabs<-rpois(length(dams),fert)
    allbabs<-sum(nbabs)
    thisgen.ped<-data.frame(ind=integer(allbabs),
                            dam=integer(allbabs),
                            sire=integer(allbabs),
                            sex=integer(allbabs))
    acc.bab<-0
    for (id in 1:length(dams)){
      # monogamy, males sampled at random
      #      mum<-sample(1:dams,size=1) #need to modify 1:dams to have real mum id, same for dad
      dad<-sample(length(sires),size=1)
      if (nbabs[id]>0)
        thisgen.ped[(acc.bab+1):(acc.bab+nbabs[id]),1:3]<-cbind(sanc+(acc.bab+1):(acc.bab+nbabs[id]),
                                                                rep(dams[id],nbabs[id]),rep(sires[dad],nbabs[id]))
      acc.bab<-acc.bab+nbabs[id]

      thisgen.ped$sex<-sample(1:2,size=allbabs,replace=TRUE,prob=c(male.prop,1-male.prop))
    }
    return(thisgen.ped)
  }

  nf<- founders.m+founders.f
  sanc<-nf #sum of ancestors
  ni<-nf #number of inds / time step
  parents<-1:nf
  alive<-rep(TRUE,nf)
  ped<-data.frame(ind=parents,dam=rep(NA,nf),sire=rep(NA,nf),sex=rep(1:2,c(founders.m,founders.f)))
  dams<-founders.f
  sires<-founders.m


  for (igen in 1:n.tstep){
    if (breed.prop <1.0){
    n.breed.males<-round(breed.prop*sum(alive & ped$sex==1),digits=0)
    if (n.breed.males==0) n.breed.males<-1
    repro.males<-sample(which(alive & ped$sex==1),size=n.breed.males)
    }
    else repro.males<-which(alive & ped$sex==1)
    repro.females<-which(alive & ped$sex==2)

    if(length(repro.males)==0 | length(repro.females)==0) stop("No mates available. exiting")

    thisgen.ped<-draw.parents(ped$ind[alive & ped$sex==2],ped$ind[repro.males],fert,prop.breed,sanc)
    #put mortality here
    dead<-sample(c(FALSE,TRUE),size=sum(alive),replace=TRUE,prob=c(death.rate,1-death.rate))
    alive[which(alive)]<-dead
    #    browser()
    print(sum(alive))

    ped<-rbind(ped,thisgen.ped)
    alive<-c(alive,rep(TRUE,dim(thisgen.ped)[1]))
    ni<-c(ni,dim(thisgen.ped)[1])
    sanc<-sum(ni)

  }

  return(ped)
}
#######################################


