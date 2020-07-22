#' Creates a phenotype from randomly chosen SNPs
#'
#'
#'
#' @usage make.traits(bed,n.causal=1000,h2=0.5,minp=0.01,sel=0.0,pop=NULL,popid=NULL)
#'
#' @param bed a bed object
#' @param n.causal number of causal loci
#' @param h2 intended trait heritability
#' @param minp minimum maf of a causal allele
#' @param sel if !=0 effect size will be scaled by \eqn{1/(2*p*(1-p))^(sel/2)}
#' @param pop variance of population effet, relative to VP
#' @param popid if NULL, popid taken from bed@ped$famid
#'
#' @return A list with components
#'
#' \itemize{
#' \item h2 he intended heritability
#' \item h2hat the realized heritability
#' \item trait a data frame of two columns, BV being the individual breeding values and pheno their phenotypes
#' \item causal a data frame with three columns, chr the chromosome where the causl loci are located, pos their position
#' and efs their effect size
#'
#' }
#'
#'
#' @export
make.traits<-function(bed,n.causal=1000,h2=0.5,minp=0.01,sel=0.0,pop=NULL,popid=NULL){
  #If sel is !=0, effect size will be scaled by
  #pop is var of the population effect, relative to VP.
  # if popid is NULL, pop taken from bed@ped$famid
  causal<-sample(which(bed@snps$maf>minp),size=n.causal)
  ef.size<-stats::rnorm(n.causal,sd=1/n.causal^.5)
  if (sel!=0.0){
    caus.p<-bed@p[causal]
    w.f.eff.siz<-(2*caus.p*(1-caus.p))^(sel/2)
    ef.size<-ef.size*w.f.eff.siz #effect size inversely prop to st.dev of allele freq
  }
  gval<-as.numeric(gaston::as.matrix(bed[,causal])%*%ef.size)
  va<-var(gval)
  phenog<-gval+stats::rnorm(nrow(bed),sd=(va*(1/h2-1))^.5)
  if(!is.null(pop)){
    if (is.null(popid)) popid<-bed@ped$famid
    popid<-factor(popid)
    npop<-nlevels(popid)
    pop.efsize<-stats::rnorm(npop,sd=sd(phenog)*pop^0.5)
    phenog<-phenog+pop.efsize[popid]
  }
  phenog<-phenog-mean(phenog)
  h2hat<-va/var(phenog) #h2
  return(list(h2=h2,h2hat=h2hat,
              trait=data.frame(BV=gval-mean(gval),pheno=phenog),
              causal=data.frame(chr=bed@snps$chr[causal],id=bed@snps$id[causal],efs=ef.size)))
}
######################################
