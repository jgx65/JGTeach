###########################################################################
#' @title Estimates allele sharing (matching) among diploid individuals
#'
#' @description Estimates from genomic data the allele sharing, or matching allelic proportion, between pairs of diploid individuals, 
#'   
#'
#' @usage get.matchingT(bed,nloc=100000,nb.cores=1)
#'
#' @param bed a bed matrix class object from the gaston package
#' @param nloc  the number of loci to process simultaneously
#' @param nb.cores if running on a multiple cores system, the number of cores to use. Default to 1 for
#' compatibility across all systems, but 4 in multi core system works well
#' 
#' @return a ni X ni matrix of allele sharing between pairs of individuals 
#'
#' @details  the matrix returned contains the proportion of loci where alleles are matching:
#' 1 if the 2 individuals are homozygous for the same allele, 0 if the two individuals are homozygous for different alleles, and 0.5 if at least one of the two individuals is a heterozygote.
#'  If M is the object where the result is stored, a genomic relationship matrix can be obtained 
#'  as kinship2GRM(beta.dosage(M, MATCHING=TRUE))
#'
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references
#'
#'
#' @examples
#' \dontrun{}
#' @export

get.matchingT<-function(bed,nloc=100000,nb.cores=1){
matching<-function(dos){
    dos <- gaston::as.matrix(dos)
	NAs<-sum(is.na(dos))
    if (NAs > 0) {
        na <- matrix(rep(1, prod(dim(dos))), ncol = ncol(dos))
        ina <- which(is.na(dos))
        na[ina] <- 0
        dos[ina] <- 1
		mNAs<-tcrossprod(na)
        Mij <- 1/2 * (1 + 1/mNAs * tcrossprod(dos -
            1))
    }
    else {
        nl <- dim(dos)[2]
        Mij <- 1/2 * (1 + tcrossprod(dos - 1)/nl)
    }
	if (NAs==0) nl<-nl else nl<-mNAs
    return(list(Mij=Mij,nl=nl))
}

if (class(bed)[[1]] != "bed.matrix") stop("bed must be of class bed.matrix. Exiting")

if (dim(bed)[2]<=nloc){ 
  matching(bed)$Mij 
  }
else {
   nl<-dim(bed)[2]
   ni<-dim(bed)[1]
   a<-split(1:nl,cut(1:nl,nl %/% nloc))

   if (nb.cores==1){
     tmp_all<-lapply(a,function(x) matching(bed[,x]))
   }
   else {
     tmp_all<-parallel::mclapply(a,function(x) matching(bed[,x]),mc.cores=nb.cores)
   }
   MijT<-matrix(numeric(ni^2),ncol=ni)
   nlt<-matrix(numeric(ni^2),ncol=ni)
   for (i in 1:length(tmp_all)) {
     MijT<-MijT+tmp_all[[i]]$Mij*tmp_all[[i]]$nl
     nlt<-nlt+tmp_all[[i]]$nl
}

   1/nlt*MijT

}
}
