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

# IO functions to read and write wspaces data to/from disk

# entry sizes by difference to account for vector overhead.
comp_t <- object.size( c( complex( 1 ), complex( 1 ) ) ) - object.size( complex( 1 ) )
real_t <- object.size( c( 1.0, 1.0 ) ) - object.size( 1.0 )


# io -----------------------------------------------------------------------------------------------
#' Load corpus data from the given directory.
#'
#' Loads corpus data from the files found in the given directory. Corpus data includes a lexicon
#' file, frequency counts for all corpus partitions (if any), POS counts for every word (lemma)
#' in the lexicon and a cooccurrence matrix.
#'
#' See details below for the different file formats used. DSV files are used for lexicon,
#' frequencies and POS counts data. Cooccurrence counts are saved as a sparse matrix stored in
#' triplet format as a plain array of (int, int, float) triplets. Each data file's format is
#' documented in its respective loading function read_*.
#'
#' @param dir         Directory to find corpus data files in.
#' @param lexicon     Lexicon file name, default 'lexicon.tsv'
#' @param frequencies Frequencies file name, default 'frequencies.tsv'
#' @param pos_table   POS counts file name, default 'pos_counts.tsv'
#' @param cooccur     Cooccurrence matrix file name, default 'cooccur.bin'
#' @param quiet       Logical indicating whether progress and general stats messages should be
#'                    silenced.
#'
#' @return A list of length four containing the loaded corpus datasets.
#'
#' @seealso \code{\link{read_lexicon}}, \code{\link{read_frequencies}},
#'          \code{\link{read_pos_counts}} and \code{\link{read_cooccur}}.
#'
#' @export
load_corpus <- function(
    dir=getwd(),
    lexicon="lxcn.dsv", frequencies="freq.dsv", pos_counts="posc.dsv", cooccur="cooc.bin",
    quiet=FALSE, attach=FALSE, env=.GlobalEnv
) {
    if( !quiet ) message( sprintf( "Loading wspaces corpus data from %s", dir ) )

    lxcn = read_lexicon( file.path( dir, lexicon ) )
    if( !quiet ) {
        message( sprintf(
            'Lexicon contains %d terms. Total TF: %d; total DF: %d',
            nrow( lxcn ), sum( lxcn$tf ), sum( lxcn$df )
        ) )
    }

    freq = read_frequencies( file.path( dir, frequencies ) )
    if( nrow( lxcn ) != nrow( freq ) ) corpus_io_error( "freq", nrow( lxcn ), nrow( freq ) )
    if( !quiet ) {
        cmrg <- colSums( freq )
        message( sprintf( 'Frequency data recorded for %d corpus splits.', length( cmrg ) ) )
        message( sprintf( 'Frequencies: min: %d; mean: %4.2f; stdv: %4.2f; max: %d',
            min( cmrg ), mean( cmrg ), sd( cmrg ), max( cmrg )
        ) )
    }

    posc = read_pos_counts( file.path( dir, pos_counts ) )
    if( nrow( lxcn ) != nrow( posc ) ) corpus_io_error( "posc", nrow( lxcn ), nrow( posc ) )
    if( !quiet ) {
        message( sprintf(
            'POS counts recorded for %d tags. Tagset: [ %s ].',
            ncol( posc ), paste( colnames( posc ), collapse=' ' )
        ) )
    }

    if( file.exists( file.path( dir, cooccur ) ) ) {
        cooc_file = file.path( dir, cooccur )
        if( !quiet ) {
            message( sprintf( "Loading cooccurrence counts found in %s", cooc_file ) )
        }
        cooc = read_cooccur( file.path( dir, cooccur ), lxcn=lxcn, shrink=FALSE )
        if( nrow( lxcn ) != nrow( cooc ) ) {
            stop( sprintf(
                "Wrong dimensions in global cooccurence file: exp %d, got %d",
                nrow( lxcn ), nrow( cooc )
            ) )
        }
        if( !quiet ) message( sprintf(
            "Cooccurrence counts loaded as %dx%d matrix with %d non-zero entries (%5.4f%% full)",
            nrow( cooc ), ncol( cooc ), Matrix::nnzero( cooc ),
            100 * (
                Matrix::nnzero( cooc ) / ( as.numeric( nrow( cooc ) ) * as.numeric( ncol( cooc ) ) )
            )
        ) )
    }

    corpus <- list(
        "lxcn" = lxcn,
        "freq" = freq,
        "posc" = posc
    )
    if( exists( 'cooc' ) ) corpus$cooc <- cooc;

    if( attach ) list2env( corpus, envir=env )
    return( corpus )
}

#' Wrapper for \code{\link{data.table::fread}} for lexical datasets.
#'
#' This function is used to load all lexical datasets in order to enforce all conventions assumed
#' in the rest of the lexical dataset manipulations in this package. See
#' \code{\link{lexical_dataset}} for details.
#'
#' @param file   Location of the file containing the lexical dataset to load.
#' @param sep    Column separator. Defaults to '@'.
#' @param header Logical. Take variable names from the first row. Defaults to TRUE.
#' @param nas    Value to replace NA's in loaded dataset. Defaults to 0. Set to NA to keep NA's.
#'
#' @return a data.frame constructed from the loaded file, with columns and keys set in the
#'         appropriate way for further lexical analysis.
#'
#' @importFrom data.table fread
read_dataset <- function( file, sep='@', header=TRUE, nas=0 ) {
    d <- data.table::fread( file, sep=sep, header=header, quote='', data.table = FALSE )
    if( !is.na( nas ) ) d[ is.na( d ) ] <- nas
    return( lexical_dataset( d ) )
}

#' Read lexicon data from DSV file.
#'
#' Reads lexicon data from a DSV file. Lexicon data is recorded on disk as a series of
#' records containings each term's string form, its term total frequency and its total document
#' frequency. Additional columns may be used for additional corpus-wide term statistics. Terms
#' themselves are \emph{not} stored as data frame columns, but are instead used as row names.
#'
#' The character vector terms <- rownames( lexicon ) should be considered the canonical list of
#' terms in a given corpus, and all other lexical data sets are guaranteed to have no more elements
#' than its length (though they may be shorter if any filters are in effect)
#'
#' @param file   Location of the lexicon file on disk.
#' @param header Logical indicating whether the given file contains a header as first row. Default
#'               TRUE.
#' @param sep    Character indicating the field separator value. Deafault '@'.
#'
#' @return A data frame containing the list of terms as its row names, a tf column for total term
#'         frequencies and a df column for total document frequency for each term. The returned data
#'         frame may contain additional columns for extra corpus-wide statistics depending on the
#'         corpus backend implementation.
#'
#' @export
read_lexicon <- function( file, header=TRUE, sep='@' ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    d <- read_dataset( file, sep=sep, header=header )
    return( d )
}

#' Read term frequency data from DSV file.
#'
#' Reads term frequency data for all corpus partitions from a DSV file. Frequency data
#' is recorded on disk as a series of records containing each term's string form and an array of
#' frequency counts for all corpus partitions. Terms themselves are \emph{not} stored as data frame
#' columns, but are instead used as \emph{row names}.
#'
#' @param file   Location of the frequencies file on disk.
#' @param header Logical indicating whether the given file contains a header as first row. Default
#'               TRUE.
#' @param sep    Character indicating the field separator value. Deafault '@'.
#' @param lxcn   (optional) A lexicon data frame used for consistency checks.
#'
#' @return A data frame containing the list of terms as its row names, and as many columns as there
#'         are corpus partitions containing each term's frequency counts for the respective corpus
#'         partition.
#'
#' @export
read_frequencies <- function( file, header=TRUE, sep='@', lxcn=NULL ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    d <- read_dataset( file, sep=sep, header=header )
    return( d )
}

#' Read POS count data from DSV file.
#'
#' Reads POS count data from a DSV file. POS count data is recorded on disk as a series
#' of records containings each term's string form, followed by one column for each major POS group
#' in the tagset used by the corpus backend for its NLP pipeline. Terms themselves are
#' \emph{not} stored as data frame columns, but are instead used as \emph{row names}.
#'
#' @param file   Location of the POS count data file on disk.
#' @param header Logical indicating whether the given file contains a header as first row. Default
#'               TRUE.
#' @param sep    Character indicating the field separator value. Deafault '@'.
#' @param drop   Logical. Drop POS columns with 0 occurrences (e.g. PUNCT in most cases). Default
#'               TRUE.
#'
#' @return A data frame containing the list of terms as its row names, and as many columns as there
#'         are POS groups in the corresponding tagset, with each term's counts on each POS group.
#'
#' @export
read_pos_counts <- function( file, header=TRUE, sep='@', drop=TRUE ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    d <- read_dataset( file, sep=sep, header=header )
    d <- if( drop ) d[ , which( colSums( d ) > 0 ) ] else d
    return( d )
}

#' Read cooccurrence matrix from disk.
#'
#' Reads cooccurrence data from a file on disk, stored as a sparse tensor saved in
#' tuple format; i.e. a plain array of ( coord, value ) structs, where coord is itself an array
#' of integer coordinates for each dimension in the represented tensor and value is a floating point
#' number.
#'
#' The default corpus backend implementation uses a term-term context definition and stores
#' cooccurrences in a rank-2 tensor, using long integers for coordinates and double precision
#' floats for values.
#'
#' The specific relationship between values in the returned tensor and the actual
#' distributional patterns for each term depends on the corpus backend cooccurrence
#' counting strategy.
#'
#' The default implementation uses a sliding window with harmonic weights equal to the inverse of
#' the distance between focus and context terms.
#'
#' This is similar to the strategy used by the Glove model for the computation of cooccurrence
#' counts. This format is for now hardcoded, but this is subject to change in future versions,
#' including both the addition of higher ranks for storing partial sum components for each
#' term-context tensor as well as summary statistics that may be useful for additional
#' transformations of the cooccurrence data.
#'
#' For now, the stored data consists of scalar values for the weighted sum of all
#' term-context observations, with no data retained about partial results.
#'
#' @param file   Location of the cooccurrence matrix file on disk.
#' @param lxcn   A corpus lexicon. Optional in most cases, but necessary for shrinking partial sets.
#' @param shrink Logical. Reduce matrix by removing empty rows and columns. Only useful with partial
#'               sets.
#'
#' @return A Matrix::sparseMatrix object containing occurrence counts for terms in contexts.
#'
#' @export
#' @importFrom Matrix rowSums colSums
read_cooccur <- function( file, lxcn=NULL, shrink=is.null( lxcn ) ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    m <- load_spm( path.expand( file ) )
    if( !is.null( lxcn ) ) {
        if( ( exp = nrow( lxcn ) ) != ( obs = nrow( m ) ) ) {
            stop( sprintf( "Wrong dimensions for cooccurence matrix: exp %d, got %d", exp, obs ) )
        }
        rownames( m ) <- rownames( lxcn )[ 1:nrow( m ) ]
        colnames( m ) <- rownames( lxcn )[ 1:ncol( m ) ]
    } else {
        rownames( m ) <- as.character( 1:nrow( m ) )
        colnames( m ) <- as.character( 1:ncol( m ) )
    }
    if( shrink ) m <- m[ Matrix::rowSums( m ) > 0, Matrix::colSums( m ) > 0 ]
    return( m )
}

#' @export
read_vectors <- function(
    filename, n=n, dim=dim, type=c( 'real', 'complex', 'binary' ), combine=FALSE
) {
    if( missing( type ) ) type = 'real'
    size <- as.integer(
        switch( type,
            real    = real_t,
            complex = comp_t,
            binary  = stop( "Binary vectors not supported yet" )
        )
    )

    bytes <- file.info( filename )$size
    if( bytes %% ( 2 * size ) ) stop( "Truncated or corrupt file!" )

    if( !missing( n ) ) {
        dim <- ( bytes / size / 2 / n ) - 1
    } else if( !missing( dim ) ) {
        n <- ( bytes / size / 2 ) / ( dim + 1 )
    } else stop( "Either dim or n must be defined" )

    f <- file( filename, 'rb' )
    d <- array(
        readBin( f,
            switch( type,
                real = 'numeric',
                complex = 'complex',
                binary = stop( "Binary vectors not supported" )
            ),
            n = ( dim + 1 ) * n * 2
        ),
        dim = c( dim + 1, n, 2 ) )
    close( f )

    if( combine ) {
        d <- d[,,1] + d[,,2]
    }

    return( d )
}


# datasets -----------------------------------------------------------------------------------------
#' Enforce internal conventions for lexical datasets.
#'
#' Lexical datasests created by upstream in memory or dumped to files always: 1) are sorted in
#' descending term frequency order, 2) have as many rows as there are terms in a corpus' lexicon,
#' even if the requested dataset contains no observations for some terms (values will be NA) and
#' 3) contain a column indicating the row's term, named as \code{conf()$termId}.
#'
#' @param d A dataset created by an obo object or read from disk.
#'
#' @return d with internal conventions enforced: a 'key' column preserving original sort order and
#'         a \code{conf()$termId} column indicating each's row term, with all data columns
#'         following.
#'
#' @export
#' @importFrom data.table setkeyv setcolorder
# TODO: define standard structure for data and metadata
lexical_dataset <- function( d ) {
    rn <- obo_mkconf()$termId()
    rownames( d ) <- d[[rn]]
    d <- d[, ! names( d ) %in% rn] # FUCK R.
    return( d )
}


#' Join lexical datasets.
#'
#' wspaces lexical datasets maintain row identity in rownames. This function wraps
#' \code{\link{dplyr::left_join}} to join the given datasets on rownames.
#'
#' @param d1 A lexical dataset
#' @param d2 A lexical dataset
#'
#' @return the left join between d1 and d2, on rownames.
#'
#' @export
#' @importFrom dplyr left_join select
lexical_join <- function( d1, d2 ) {
    d1$key_ <- rownames( d1 )
    d2$key_ <- rownames( d2 )
    out <- dplyr::left_join( d1, d2, by='key_' ) # TODO remove dplyr
    rownames( out ) <- out$key_
    out %<>% dplyr::select( -key_ )
    return( out )
}

#' Extract term data from the given lexical dataset.
#'
#' Extracts term data from the given term character vector by matching the given terms against the
#' rownames in the given dataset. All lexical datasets constructed by wspaces keep their term
#' identifier as rownames.
#'
#' @param d A lexical dataset, i.e. a data frame or matrix with terms as row names.
#' @param terms A character vector with terms to get data for.
#'
#' @return A subset of d containing only entries that are in terms.
#'
#' @export
term_data <- function( d, terms ) {
    filter <- term_idx( d, terms )
    return( if( is.vector( d ) ) d[ filter ] else d[ filter, ] )
}

#' Get index vector into the given data for the given terms.
#'
#' Generates an index vector for each entry in the given terms character vector into the elements
#' or rows of the given named vector, matrix or lexical dataset.
#'
#' Use this method to consistently extract data from different sources for a given subset of terms.
#'
#' @param d Data into which t compute index vector. Can be a named vector, named matrix or lexical
#'          dataset.
#' @param terms A character vector with the terms for which to compute indices.
#'
#' @return An index vector with the location of data entries in d for the terms in terms.
#'
#' @export
term_idx <- function( d, terms ) {
    u <- if( is.vector( d ) && !is.null( names( d ) ) ) {
        names( d )
    } else if( is.matrix( d ) || is.data.frame( d ) ) {
        rownames( d )
    } else {
        stop( sprintf(
            "d must be a vector, matrix or data.frame, got %s",
            paste( class( d ), collapse=", " )
        ) )
    }
    return( which( u %in% terms ) )
}
