% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maketraits.R
\name{make.traits}
\alias{make.traits}
\title{Creates a phenotype from randomly chosen SNPs}
\usage{
make.traits(bed,n.causal=1000,h2=0.5,minp=0.01,sel=0.0,pop=NULL,popid=NULL)
}
\arguments{
\item{bed}{a bed object}

\item{n.causal}{number of causal loci}

\item{h2}{intended trait heritability}

\item{minp}{minimum maf of a causal allele}

\item{sel}{if !=0 effect size will be scaled by \eqn{1/(2*p*(1-p))^(sel/2)}}

\item{pop}{variance of population effet, relative to VP}

\item{popid}{if NULL, popid taken from bed@ped$famid}
}
\value{
A list with components

\itemize{
\item h2 he intended heritability
\item h2hat the realized heritability
\item trait a data frame of two columns,
BV being the individual breeding values and
pheno their phenotypes
\item causal a data frame with three columns,
chr the chromosome where the causal loci are located,
id the SNP id and efs their effect size

}
}
\description{
If sel is !=0, effect size will be scaled by
pop is var of the population effect, relative to VP.
if popid is NULL, pop taken from bed@ped$famid
}
