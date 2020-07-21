#' Expected coancestries and betas from a two population model
#'
#' Quantifies expected coancestries and betas (coancestries relative to
#' coancestries between population, see
#' \href{http://www.genetics.org/content/206/4/2085}{Weir and Goudet
#' (2017) Genetics}) in a two populations model.
#'
#'
#'
#' @usage thet.bet.2pops((mu=10^-8,m1=10^-4,m2=10^-4,n1=1000,n2=1000,ngen=100,plotit=TRUE))
#'
#' @param mu mutation rate
#' @param m1 immigration into pop1
#' @param m2 immigration into pop2
#' @param n1 size of pop1
#' @param n2 size of pop2
#' @param ngen number of generations
#' @param plotit
#'
#' @return a list with elements
#' \itemize{
#' \item Th a matrix of coancestries \eqn{theta} for pop1, pop2 and pop1-2 (columns) and generations (rows)
#' \item Be a matrix of \eqn{beta}s for pop1, pop2 and average (columns) and generations (rows)
#' }
#'
#'  @reference{
#'  \href{http://www.genetics.org/content/206/4/2085}{Weir and Goudet
#' (2017) Genetics 206:2085}) A Unified Characterization
#' of Population Structure and Relatedness.
#'  }
#'
#'  @example
#'  \notrun{
#'  # reproduces fig 1 in Weir and Goudet (2017)
#'  #top row
#'  thet.bet.2pops(mu=0.0,m2=0.0,m1=0.0,n1=10000,n2=100,ngen=10000)
#'  #middle row
#'  thet.bet.2pops(mu=0.001,m2=0.0,m1=0.0,n1=10000,n2=100,ngen=10000)
#'  #bottom row
#'  thet.bet.2pops(mu=0.001,m2=0.0,m1=0.01,n1=10000,n2=100,ngen=10000)
#'  }
#'
#'
#' @export

thet.bet.2pops<-function(mu=10^-8,m1=10^-4,m2=10^-4,n1=10000,n2=100,ngen=100,plotit=TRUE){
  a<-(1-mu)^2
  b1<-(1-m1)^2
  b2<-(1-m2)^2
  c1<-1/2/n1
  c2<-1/2/n2

  Xv<-a*c(b1*c1+m1^2*c2,m2^2*c1+b2*c2,(1-m1)*m2*c1+m1*(1-m2)*c2)
  Xm<-a*matrix(c(
    b1*(1-c1),m1^2*(1-c2),2*m1*(1-m1),
    m2^2*(1-c1),b2*(1-c2),2*m2*(1-m2),
    (1-m1)*m2*(1-c1),m1*(1-m2)*(1-c2),(1-m1)*(1-m2)+m1*m2
  ),nrow=3,byrow=TRUE)
  y<-matrix(numeric(ngen*3),ncol=3)
  betas<-matrix(numeric(ngen*3),ncol=3)
  y[1,]<-0
  betas[1,]<-NA
  for (igen in 2:ngen){
    y[igen,]<-Xv+Xm%*%y[(igen-1),]
  }
  betas[,1]<-(y[,1]-y[,3])/(1-y[,3])
  betas[,2]<-(y[,2]-y[,3])/(1-y[,3])
  betas[,3]<-(betas[,1]+betas[,2])/2
  if(plotit){
    par(mfrow=c(1,2))
    plot(1:ngen,y[,1],ylim=range(y),xlab="Generations",
         ylab=expression(theta),col="red",type="l",lwd=2)
    lines(1:ngen,y[,2],col="blue",lwd=2)
    lines(1:ngen,y[,3],col="orange",lwd=2)
    legend("bottomright",c(expression(theta^{11}),expression(theta^{22}),
                           expression(theta^{12})),col = c("red","blue","orange"),lwd=2)
    title(paste("mu: ",mu,"; m1: ",m1,"; m2: ",m2,"\n",
                "N1: ",n1,"; N2: ",n2,sep="",collapse=""))
    plot(1:ngen,betas[,1],ylim=range(betas,na.rm=TRUE),xlab="Generations",
         ylab=expression(beta),type="l",col="red",lwd=2)
    lines(1:ngen,betas[,2],col="blue",lwd=2)
    lines(1:ngen,betas[,3],col="black",lwd=2)
    legend("bottomright",c(expression(beta^{11}),expression(beta^{22}),
                           expression(beta^{W})),col = c("red","blue","black"),lwd=2)
    par(mfrow=c(1,1))
  }
  return(list(Th=y,Be=betas))
}

