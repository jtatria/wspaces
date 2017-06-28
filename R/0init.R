#' @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp evalCpp
#' @useDynLib wspaces
NULL

#' @importFrom rJava .jpackage
.onLoad <- function( libname, pkgname ) {
    .jpackage( pkgname, lib.loc = libname )
}
