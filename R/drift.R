#' Simulates drift
#'
#' Simulates drift
#'
#' @usage drift (nind=100, ngen=100, p0=0.5, nrep=1000, PlotIt=FALSE)
#'
#' @param nind number of individuals
#' @param ngen number of generations
#' @param p0 initial allel frequency
#' @param nrep number of replicates
#' @param PlotIt whether to plot the allelic trajectories
#'
#' @return a matrix with the frequencies of replicates in column for each generation in rows
#'
#' @example
#' \notrun{
#' res<-drift(ngen=500)
#' par(mfrow=c(3,3))
#' for (i in c(1,2,5,10,20,50,100,200,500))
#' hist(res[i,],breaks=0:100/100,main=paste0('gen: ',i),xlab="p",ylab="")
#' }
#'
#' @export

drift<-function(nind=100, ngen=100, p0=0.5, nrep=1000, PlotIt=FALSE){
  #nind: is the number of individulas
  #ngen: number of generations
  #p0: initial allele frequency
  #nrep: number of replicates
  #PlotIt: whether
  freq<-matrix(numeric(ngen*nrep), ncol=nrep)
  freq[1,]<-p0
  for (j in 2:ngen){
    next.gene.draws = rbinom(nrep, 2*nind, freq[j-1,])
    freq[j,]<-next.gene.draws/2/nind
  }
  if (PlotIt){
    plot(1:ngen, freq[,1], ylim=c(0,1), type="l", lty=1,
         xlab="Generations", ylab="p",main=paste("p0=",p0," , popsize = ",nind))
    for (i in 2:nrep){
      lines(1:ngen, freq[,i], lty=i)
    }
  }
  return(freq)
}
