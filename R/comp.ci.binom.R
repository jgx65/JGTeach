#' Compares Confidence intervals of the binomial distribution
#'
#' computes 95\% confidence intervals of the probability of
#' success using the normal approximation, bootstrapping
#' and an exact confidence interval
#'
#' @usage comp.ci.binom(p=0.5,n=1000,nboot=1000,digits=4)
#'
#' @param p the binomial probabilty of success (frequency)
#' @param n How many draws
#' @param nboot  number of bootstrap samples
#' @param digits how many digits to print out
#'
#' @return a list with components
#' \itemize{
#' \item p.hat
#' \item norm.ci Wald confidence interval
#' \item boot.ci
#' \item exact.ci
#' }
#'
#' @seealso the corresponding \href{https://jgx65.shinyapps.io/BinomNorm/}{shiny app}
#'
#' @examples {
#' comp.ci.binom()
#' comp.ci.binom(n=10,p=0.1)
#' }
#'
#' @export
comp.ci.binom<-function(p=0.5,n=1000,nboot=1000,digits=4){
  #illustrates slides 48 and followings
  x<-rbinom(1,n,p)
  p.hat<-x/n
  sd.p<-(p.hat*(1-p.hat)/n)^0.5
  norm.ci<-round(c(p.hat-1.96*sd.p,p.hat+1.96*sd.p),digits=digits) #wald ci
  exact.ci<-round(binom.test(x,n)$conf.int[1:2],digits=digits) #exact from binomial
  bx<-numeric(nboot)
  ox<-numeric(n);ox[1:x]<-1
  for (i in 1:nboot) bx[i]<-sum(sample(ox,replace=TRUE))
  boot.ci<-round(quantile(bx,c(0.025,0.975))/n,digits=digits) #empirical bootstrap ci
  return(list(p.hat=p.hat,norm.ci=norm.ci,boot.ci=boot.ci,exact.ci=exact.ci))
}


