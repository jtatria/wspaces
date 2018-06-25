#
# Copyright (C) 2017 José Tomás Atria <jtatria at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

####################################################################################################
#### Package-wide initialization                                                                ####
####################################################################################################

# @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib wspaces
NULL

JAVA_XMX = 2048

LECTOR_JAR = 'lector.jar'

java_setup <- function( libname, pkgname ) {
    # TODO: remove pending rJava memory management
    options( java.parameters=sprintf( "-Xmx%dm", JAVA_XMX ) )
    rJava::.jpackage( pkgname, lib.loc = libname )
}

# TODO: spawn lector backend to its own package
.onLoad <- function( libname, pkgname ) {
    if( file.exists( file.path( libname, pkgname, 'java', LECTOR_JAR ) ) ) {
        java_setup( libname, pkgname )    
    } else {
        message( "Lector back-end not found. All 'lector_' functions are not available" )
    }
}
