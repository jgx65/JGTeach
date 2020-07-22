#' calculates the additive relationship matrix from a pedigree
#'
#' given a vector of sires and dams from a pedigree, calculates
#' the Additive relationship matrix
#'
#' @usage pedARM(s,d,fcoan=NULL)
#'
#' @param s a vector with sire identies
#' @param d a vector wit dam identities
#' @param fcoan  a matrix of coancestries for the founders, if available
#'
#' @return a matrix with the expected relatedness
#' among all airs of individuals and twice self coancestry
#' on the diagonal
#'
#' @export
pedARM<-function(s, d,fcoan=NULL){
  #fcoan is the matrix of relatedness among founders

  if (nargs()==1){
    stop("sire vector and dam vector are required")
  }

  if (length(s) != length(d)){
    stop("size of the sire vector and dam vector are different!")
  }

  n <- length(s)
  N <- n + 1
  A <- matrix(0, ncol=N, nrow=N)

  #all founders to last column
  if (is.null(fcoan)){
    s <- (is.na(s))*(N) + replace(s,is.na(s),0)
    d <- (is.na(d))*N + replace(d,is.na(d),0)
    n0<-1
  }
  else{
    s <- (is.na(s))*(N) + replace(s,is.na(s),0)
    d <- (is.na(d))*N + replace(d,is.na(d),0)
    nf<-dim(fcoan)[1]
    A[1:nf,1:nf]<-fcoan
    n0<-nf+1
  }


  for(i in 1:n){
    if (i>=n0)
      A[i,i] <- 1 + A[s[i], d[i]]/2

    for(j in (i+1):n){
      if (j > n) break
      if (j>=n0){
        A[i,j] <- ( A[i, s[j]] + A[i,d[j]] )/2
        A[j,i] <- A[i,j]
      }
    }
  }
  return(A[1:n, 1:n])

}
