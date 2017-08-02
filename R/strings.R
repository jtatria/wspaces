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

# String manipultion tools

#' Perl-style concatenation operator.
#' 
#' Wrapper for paste( ..., sep='' )
#' 
#' @export
"%.%" <- function( ... ) paste( ..., sep ='' )

#' Normalize strings.
#'
#' This function offers a consistent procedure for string normalization combining several common
#' operations into one call. Term transofrmations across a corpus should always use a consistent
#' normalization function, e.g. currying this function with a fixed set of parameters.
#'
#' @param x     A string (i.e. a character vector of length 1).
#' @param trim  Logical. Remove trailing and leading whitespace.
#' @param lower Logical. Reduce everything to lower case.
#' @param white Logical. Eliminate duplicate whitespace.
#' @param nl    Loigcal. Replace newlines with spaces.
#' @param punct Logicel. Remove punctuation.
#' @param word  Logical. Remove non-word characters.
#'
#' @return a reasonably normalized string.
#'
#' @examples
#' \code{
#'     x <- "some Ugly    dirty\nstring!"
#'     cat( x )
#'     cat( normalizeString( x ) )
#' }
#' \code{
#'     std_str_normalize <- function( x ) {
#'         str_normalize( x, trim=TRUE, lower=TRUE, white=TRUE, nl=TRUE, punct=FALSE )
#'     }
#'     normals <- std_str_normalize( ugly_strings )
#' }
str_normalize <- function( x, trim=TRUE, lower=TRUE, white=TRUE, nl=FALSE, punct=FALSE ) {
    if( !any( trim, lower, white, nl, punct ) ) return( x ) # All flags are false, nothing to do.
    if( lower ) x <- tolower( x )
    if( trim  ) x <- gsub( '^\\s+(.*?)\\s+$', '\\1', x, perl=TRUE )
    if( white ) x <- gsub( '\\s+', ' ', x )
    if( nl    ) x <- gsub( '\\n+', ' ', x )
    if( punct ) x <- gsub( '[[:punct:]]', '', x )
    if( word )  x <- gsub( '\\W+', '', x, perl=TRUE )
    return( x )
}