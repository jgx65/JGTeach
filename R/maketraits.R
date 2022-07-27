#' Creates a phenotype from randomly chosen SNPs
#'
#' @description
#'
#' Generates a phenotype from a random subset of SNPs with the required heritability.
#' possibly including
#' a (random) population specific effect.
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
#' \item trait a data frame of two columns,
#' gval being the individual genetic values and
#' pheno their phenotypes
#' \item causal a data frame with three columns,
#' chr the chromosome where the causal loci are located;
#' id the SNP id; and efs their effect size
#'
#' }
#'
#' @details
#'
#' Traits are built using additive effects only. Genetic
#' values are the sum of the products of dosages at causal loci by their effect sizes.
#'
#' Additive variance VA is measured as the variance of genetic values. A normal
#' deviate with variance \eqn{VA * (1/h^2-1)}, where  \eqn{h^2} is the intended trait
#' heritability, is  added to the genetic values to make the phenotypes. If pop!=NULL,
#' a random deviate with variance proportional to the phenotypic variance VP (i.e.
#' \eqn{Vpop=pop*VP}) is added to the individual phenotypes in each population.
#'
#'
#' If sel is !=0, effect size will be scaled by  \eqn{1/(2*p*(1-p))^(sel/2)}, where p is the
#' sample frequency of the derived allele.
#' pop is the variance of the population effect, relative to VP.
#' if popid is NULL, popid taken from bed@ped$famid
#' @export

########################3
make.traits<-function(bed=bed,n.causal=1000,h2=0.5,minp=0.01,sel=0.0,pop=NULL,popid=NULL){

    causal<-sample(which(bed@snps$maf>minp),size=n.causal)
  ef.size<-stats::rnorm(n.causal,sd=1/n.causal^.5)
  if (sel!=0.0){
    caus.p<-bed@p[causal]
    w.f.eff.siz<-(2*caus.p*(1-caus.p))^(sel/2)
    ef.size<-ef.size*w.f.eff.siz
  }
  gval<-as.numeric(gaston::as.matrix(bed[,causal])%*%ef.size)
  va<-var(gval)
  phenog<-gval+stats::rnorm(nrow(bed),sd=(va*(1/h2-1))^.5)
  if(!is.null(pop)){
    if (is.null(popid)) popid<-bed@ped$famid
    popid<-factor(popid)
    npop<-nlevels(popid)
    pop.efsize<-stats::rnorm(npop,sd=sd(phenog)*pop^0.5)
    for (i in 1:npop){
      x<-levels(popid)[i]
      phenog[popid==x]<-phenog[popid==x]+pop.efsize[i]
    }
  }
  phenog<-phenog-mean(phenog)
  h2hat<-va/var(phenog) #h2
  return(list(h2=h2,h2hat=h2hat,
              trait=data.frame(gval=gval-mean(gval),pheno=phenog),
              causal=data.frame(chr=bed@snps$chr[causal],id=bed@snps$id[causal],efs=ef.size)))
}
######################################
