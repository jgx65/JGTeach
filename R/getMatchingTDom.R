###########################################################################
#' @title Estimates dominance effects sharing among individuals
#'
#' @description Estimates from genomic data the dominance effects sharing between pairs of individuals, 
#'  the proportion of loci where both individuals are heterozygous 
#'
#' @usage get.matchingT.dom(bed,nloc=100000,nb.cores=1)
#'
#' @param bed a bed matrix class object from the gaston package
#' @param nloc  the number of loci to process simultaneously
#' @param nb.cores if running on a multiple cores system, the number of cores to use. Default to 1 for
#' compatibility across all systems, but 4 in multi core system works well
#' 
#' @return a ni X ni matrix of dominance sharing between pairs of individuals 
#'
#' @details  the matrix returned contains the proportion of loci where both inddividuals of the pair are heterozygous
#'  if MD is the object where the result is stored, a dominance relationship matrix can be obtained 
#'  as kinship2dist(beta.dosage(MD, MATCHING=TRUE))
#'
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references
#'
#'
#' @examples
#' \dontrun{}
#' @export

get.matchingT.dom<-
function(bed,nloc=100000,nb.cores=1){
matchingD<-function(dos){
    dos <- gaston::as.matrix(dos)
        NAs<-sum(is.na(dos))
    if (NAs > 0) {
        na <- matrix(rep(1, prod(dim(dos))), ncol = ncol(dos))
        ina <- which(is.na(dos))
        na[ina] <- 0
        dos[ina] <- 1
                mNAs<-tcrossprod(na)
        Mij <- 1/mNAs*tcrossprod(dos*(2-dos))# 1/2 * (1 + 1/mNAs * tcrossprod(dos -1))
    }
    else {
        nl <- dim(dos)[2]
        Mij <- tcrossprod(dos*(2-dos))/nl #1/2 * (1 + tcrossprod(dos - 1)/nl)
    }
        if (NAs==0) nl<-rep(nl,ni) else nl<-diag(mNAs)
    return(list(Mij=Mij,nl=nl))
}

if (class(bed)[[1]] != "bed.matrix") stop("bed must be of class bed.matrix. Exiting")

if (dim(bed)[2]<=nloc){
  matchingD(bed)$Mij
  }
else {
   nl<-dim(bed)[2]
   ni<-dim(bed)[1]
   a<-split(1:nl,floor(1:nl/nloc))
   if (nb.cores>1)
     tmp_all<-parallel::mclapply(a,function(x) matchingD(bed[,x]),mc.cores=nb.cores)
   else tmp_all<-lapply(a,function(x) matchingD(bed[,x]))

   MijT<-matrix(numeric(ni^2),ncol=ni)
   nlt<-rep(0,ni)
   for (i in 1:length(tmp_all)) {
     MijT<-MijT+tmp_all[[i]]$Mij*tmp_all[[i]]$nl
     nlt<-nlt+tmp_all[[i]]$nl
}

   MijT/nlt

}
}
