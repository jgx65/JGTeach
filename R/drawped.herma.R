###########################################################################
#' @title Builds a hermaphrodite pedigree
#'
#' @description Builds a hermaphrodite pedigree from the passed arguments
#'
#' @usage buildped.herma(founders=20,fert=0.5,death.rate=0.5,self.prop=0.5,n.tstep=10)
#'
#' @param founders numbers of male founders
#' @param fert hermaphrodite fertility rate per time interval
#' @param death.rate mortality rate per time interval
#' @param self.prop proportion of selfing
#' @param n.tstep number of time step
#'
#'
#' @return a pedigree data frame with 3 columns [ind id, parent1 id, parent2 id] and as many rows as individuals
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
#' test<-buildped.herma(founders=20,fert=0.5,death.rate=0.5,self.prop=0.5,n.tstep=5)
#' #range of kinship coefficients
#' range(mat2vec(grm2kinship(pedARM(test[,2],test[,3]))))
#' #range of inbreeding coefficients
#' range(diag(grm2kinship(pedARM(test[,2],test[,3]))))
#' }
#' @export

#######################################
buildped.herma<-function(founders=20,fert=0.5,death.rate=0.5,self.prop=0.5,n.tstep=10){
  draw.parents<-function(parents,fert,self.prop,sanc){
    np<-length(parents)
    nbabs<-rpois(np,fert)
    allbabs<-sum(nbabs)
    thisgen.ped<-data.frame(ind=integer(allbabs),
                            par1=integer(allbabs),
                            par2=integer(allbabs))
    is.self<-sample(c(TRUE,FALSE),replace=TRUE,size=np,prob=c(self.prop,1-self.prop))
    acc.bab<-0
    for (ii in 1:np){
      # hermaphrodites with selfing prop self.prop
      if(is.self[ii]) p2<-ii else p2<-sample((1:np)[-ii],size=1)

      if (nbabs[ii]>0)
        thisgen.ped[(acc.bab+1):(acc.bab+nbabs[ii]),1:3]<-cbind(sanc+(acc.bab+1):(acc.bab+nbabs[ii]),
                                                                rep(parents[ii],nbabs[ii]),rep(parents[p2],nbabs[ii]))
      acc.bab<-acc.bab+nbabs[ii]

    }
    return(thisgen.ped)
  }

  nf<- founders
  sanc<-nf #sum of ancestors
  ni<-nf #number of inds / time step
  parents<-1:nf
  alive<-rep(TRUE,nf)
  ped<-data.frame(ind=parents,par1=rep(NA,nf),par2=rep(NA,nf))

  for (igen in 1:n.tstep){
    thisgen.ped<-draw.parents(ped$ind[alive],fert,self.prop,sanc)
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
