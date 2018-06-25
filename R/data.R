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
#' file, frequency counts for all corpus partitions (if any), POS class counts for every term
#' in the lexicon and a cooccurrence matrix.
#'
#' See details below for the different file formats used. DSV files are used for lexicon,
#' frequencies and POS counts data. Cooccurrence counts are saved as a sparse matrix stored in
#' triplet format as a plain array of (coordinates,value) tuples. Each data file's format is
#' documented in its respective loading function read_*.
#' 
#' TODO: Create S3 classes and methods.
#'
#' @param dir         Directory in which to find corpus data files. Defaults to the current working
#'                    dir.
#' @param lxcn.file   Lexicon file name. Default 'lxcn.dsv'.
#' @param freq.file   Frequencies file name. Default 'freq.dsv'.
#' @param posc.file   POS counts file name. Default 'posc.dsv'.
#' @param cooc.file   Cooccurrence matrix file name. Default 'cooc.bin'.
#' @param quiet       Logical indicating whether progress and general stats messages should be
#'                    silenced. Defaults to \code{FALSE}.
#' @param attach      Logical indicating whtether the symbols associated to the loaded objects 
#'                    should be attached to the given environment. Defaults to \code{FALSE}.
#' @param env         Environment in which to attach symbols, if \code{attach} is \code{TRUE}.
#'                    Defaults to \code{.GlobalEnv}, ignored if \code{attach} is \code{FALSE}.
#'
#' @return A list with named entries containing the loaded corpus datasets.
#'
#' @seealso \code{\link{read_lexicon}}, \code{\link{read_frequencies}},
#'          \code{\link{read_pos_counts}} and \code{\link{read_cooccur}}.
#'
#' @export
load_corpus <- function(
    dir=getwd(),
    lxcn.file="lxcn.dsv", freq.file="freq.dsv", posc.file="posc.dsv", cooc.file="cooc.bin",
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

#' Load corpus sample data from the given directory.
#' 
#' A corpus sample is defined a subset of documents from a given corpus. All corpus data sets for a 
#' given sample are compiled considering only documents in the subset in the sample.
#' Samples can be loaded autonomously, in which case this function does the same as 
#' \link{load_corpus}, or associated to a set of universe corpus data sets, in which case the 
#' lexicon from the global corpus is needed to ensure proper term indices for sparse cooccurrence 
#' matrices. This is irrelevant if no cooccurrence counts are loaded.
#' 
#' TODO: merge this function with \link{load_corpus}.
#' 
#' @param dir         Directory in which to find corpus data files. Defaults to the current working
#'                    dir.
#' @param lxcn        An optional lexicon for the global corpus from which the sample is taken. 
#'                    Defaults to \code{NULL}.
#' @param lxcn.file   Lexicon file name for global lexicon. Default 'lxcn.dsv'.
#' @param freq.file   Sample frequencies file name. Default 'freq.dsv'.
#' @param posc.file   Sample POS counts file name. Default 'posc.dsv'.
#' @param cooc.file   Sample cooccurrence matrix file name. Default 'cooc.bin'.
#' @param quiet       Logical indicating whether progress and general stats messages should be
#'                    silenced. Defaults to \code{FALSE}.
#' @param attach      Logical indicating whtether the symbols associated to the loaded objects 
#'                    should be attached to the given environment. Defaults to \code{FALSE}.
#' @param env         Environment in which to attach symbols, if \code{attach} is \code{TRUE}.
#'                    Defaults to \code{.GlobalEnv}, ignored if \code{attach} is \code{FALSE}.
#'
#'@return A list with named entries containing the loaded corpus datasets.
#'
#' @export
load_sample <- function(
    dir, lxcn=NULL,
    lxcn.file='lxcn.dsv', freq.file='freq.dsv', posc.file='posc.dsv', cooc.file='cooc.bin',
    quiet=FALSE, attach=FALSE, env=.GlobalEnv
) {
    if( !quiet ) message( sprintf( 'Loading corpus sample data from %s', dir ) )
    
    add.lxcn <- FALSE
    if( is.null( lxcn ) ) {
        lxcn.p <- file.path( dir, lxcn.file )
        if( !file.exists( lxcn.p ) ) stop( sprintf( 'Can\'t read lexicon from %s: file not found', lxcn.p ) )
        message( sprintf( 'Loading sample lexicon from %s', lxcn.p ) )
        lxcn <- read_lexicon( lxcn.p )
        add.lxcn <- TRUE
    }
    
    freq <- read_frequencies( file.path( dir, freq.file ) )
    posc <- read_pos_counts( file.path( dir, posc.file ) )
    cooc <- read_cooccur( file.path( dir, cooc.file ), lxcn=lxcn )
    out <- list( freq=freq, posc=posc, cooc=cooc, dir=dir )
    if( add.lxcn ) out$lxcn <- lxcn
    class( out ) <- 'corpus_sample'
    if( attach ) list2env( out, envir=env )
    return( out )
}

#' Load all samples found in the given base directory.
#' 
#' Files corresponding to each sample should be located in subdirectories within the base directory. 
#' The name of each subdirectory will be used as sample identifier; this function will scan for all 
#' samples in the given base dir and run \link{load_sample} on each.
#' 
#' If no lexicon is given as \code{lxcn}, this function will look for a global lexicon in the file 
#' named as \code{lxcn.file} \emph{in the base directory}. Other sample data sets will be searched 
#' for inside each sample's subdirectory.
#' 
#' If \code{attach} is \code{TRUE}, this function will attach the symbols corresponding 
#' \emph{to each sample}, each of which will be associated to a list containing each sample's 
#' data sets, as per \link{load_sample}.
#' 
#' TODO: revise design for load_* functions.
#' 
#' @param dir         Directory in which to find corpus data files. Defaults to the current working
#'                    dir.
#' @param lxcn        An optional lexicon for the global corpus from which the sample is taken. 
#'                    Defaults to \code{NULL}.
#' @param lxcn.file   Lexicon file name for global lexicon. Default 'lxcn.dsv'.
#' @param freq.file   Sample frequencies file name. Default 'freq.dsv'.
#' @param posc.file   Sample POS counts file name. Default 'posc.dsv'.
#' @param cooc.file   Sample cooccurrence matrix file name. Default 'cooc.bin'.
#' @param quiet       Logical indicating whether progress and general stats messages should be
#'                    silenced. Defaults to \code{FALSE}.
#' @param attach      Logical indicating whtether the symbols associated to the loaded objects 
#'                    should be attached to the given environment. Defaults to \code{FALSE}. See
#'                    details.
#' @param env         Environment in which to attach symbols, if \code{attach} is \code{TRUE}.
#'                    Defaults to \code{.GlobalEnv}, ignored if \code{attach} is \code{FALSE}.
#'                    
#' @return A list containing an element for each sample, each of which will be equal to the value 
#'         of \link{load_sample}.
#' 
#' @seealso \code{link{load_sample}} \code{\link{load_corpus}}.
#' 
#' @export
load_sample_set <- function(
    dir=getwd(), lxcn=NULL,
    lxcn.file='lxcn.dsv', freq.file='freq.dsv', posc.file='posc.dsv', cooc.file='cooc.bin',
    quiet=FALSE, attach=FALSE, env=.GlobalEnv
) {
    if( !quiet ) message( sprintf( "Loading corpus samples from %s", dir ) )
    if( is.null( lxcn ) ) {
        lxcn.p <- file.path( dir, lxcn.file )
        if( !file.exists( lxcn.p ) ) stop( sprintf( 'Can\'t read lexicon from %s: file not found', lxcn.p ) )
        message( sprintf( 'Loading sample lexicon from %s', lxcn.p ) )
        lxcn <- read_lexicon( lxcn.p )
        add.lxcn <- TRUE
    }
    sdirs <- list.dirs( dir, recursive=FALSE )
    samples <- list()
    for( sdir in sdirs ) {
        sname <- gsub( ".*/", "", sdir )
        samples[[sname]] <- load_sample( sdir, lxcn=lxcn, quiet=quiet, attach=FALSE )
    }
    samples$lxcn <- lxcn
    class( samples ) <- 'corpus_sample_set'
    if( attach ) list2env( samples, envir=env )
    return( samples )
}

# low-level io --------------------------------------------------------------------------------

#' Wrapper for \code{\link{data.table::fread}} for lexical datasets.
#'
#' This function is used to load all lexical datasets in order to enforce all conventions assumed
#' in the rest of the lexical dataset manipulations in this package. See
#' \code{\link{lexical_dataset}} for details.
#'
#' TODO: get rid of data.table::fread
#'
#' @param file   Location of the file containing the lexical dataset to load.
#' @param sep    Column separator. Defaults to '@'.
#' @param header Logical. Take variable names from the first row. Defaults to TRUE.
#' @param nas    Value to replace NA's in loaded dataset. Defaults to 0. Set to NA to keep NA's.
#'
#' @return a data.frame constructed from the loaded file, with columns and keys set in the
#'         appropriate way for further lexical analysis.
#'
#' @export
#' @importFrom data.table fread
read_dataset <- function( file, sep='@', header=TRUE, nas=0 ) {
    d <- data.table::fread( file, sep=sep, header=header, quote='', data.table = FALSE )
    if( !is.na( nas ) ) d[ is.na( d ) ] <- nas
    return( lexical_dataset( d ) )
}

#' Read lexicon data from DSV file.
#'
#' Reads lexicon data from a DSV file. Lexicon data is recorded on disk as a series of
#' records containings each term's string form, its total term frequency and its total document
#' frequency.
#' 
#' Additional columns may be used for additional corpus-wide term statistics. Terms
#' themselves are \emph{not} stored as data frame columns, but are instead used as row names.
#'
#' The character vector \code{terms <- rownames( lexicon )} should be considered the canonical list 
#' of terms in a given corpus, and all other lexical data sets from the same corpus are guaranteed 
#' to have no more elements than its length (though they may be shorter if any term filters are in 
#' effect).
#'
#' @param file   Location of the lexicon file on disk.
#' @param header Logical indicating whether the given file contains a header as first row. Default
#'               \code{TRUE}.
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
read_frequencies <- function( file, header=TRUE, sep='@' ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    d <- read_dataset( file, sep=sep, header=header )
    return( d )
}

#' Read POS classes count data from DSV file.
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
#' @param drop   Logical. Drop POS class columns with 0 occurrences (e.g. PUNCT in most cases).
#'               Default TRUE.
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
#' cooccurrences in a rank-2 tensor (i.e. a matrix), using long integers for coordinates and double
#' precision floats for values.
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
        if( ( exp = nrow( lxcn ) ) < ( obs = nrow( m ) ) ) {
            stop( sprintf( "Extra rows in cooccurence matrix: exp %d, got %d", exp, obs ) )
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
    file, n=n, dim=dim, type=c( 'real', 'complex', 'binary' ), combine=FALSE
) {
    if( missing( type ) ) type = 'real'
    size <- as.integer(
        switch( type,
            real    = real_t,
            complex = comp_t,
            binary  = stop( "Binary vectors not supported yet" )
        )
    )

    bytes <- file.info( file )$size
    if( bytes %% ( 2 * size ) ) stop( "Truncated or corrupt file!" )

    if( !missing( n ) ) {
        dim <- ( bytes / size / 2 / n ) - 1
    } else if( !missing( dim ) ) {
        n <- ( bytes / size / 2 ) / ( dim + 1 )
    } else stop( "Either dim or n must be defined" )

    f <- file( file, 'rb' )
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
#' @param d A dataset created by a Lector object or read from disk.
#'
#' @return d with internal conventions enforced: a 'key' column preserving original sort order and
#'         a \code{conf()$termId} column indicating each's row term, with all data columns
#'         following.
#'
#' @export
#' @importFrom data.table setkeyv setcolorder
# TODO: define standard structure for data and metadata
lexical_dataset <- function( d ) {
    rn <- lector_mkconf()$termId()
    rownames( d ) <- d[[rn]]
    d <- d[, ! names( d ) %in% rn] # FUCK R.
    return( d )
}

#' Join lexical datasets.
#'
#' wspaces lexical datasets maintain row identity in rownames. This function wraps joining
#' operations to enforce wspaces rowname conventions.
#'
#' If no key names are given, joins are performed on the rownames of the given data sets. If a key
#' name is given for either d1 or d2, the column named with that key name is used for the
#' corresponding data set.
#'
#' By default, the effective join key (rownames or named vector) will be added as rownames in the
#' results data frame, unless the no.rn parameter is TRUE or its values in the result are not
#' unique, which will result in a warning.
#'
#' The mode parameter indicates the kind of join: 'left' or 'right' will return all rows in d1
#' or d2 with all columns from d1 and d2. 'inner' will return only the intersection of both d1 and
#' d2 with columns from d1 and d2. 'full' will return the union of d1 and d2, with columns from d1
#' and d2. The default is 'left': add data from d2 to all rows in d1.
#' 
#' TODO: get rid of dplyr.
#'
#' @param d1    A lexical dataset.
#' @param d2    A lexical dataset.
#' @param k1    An optional character vector with the name of a key for d1 to use instead of
#'              rownames. Defaults to NULL (join on rownames).
#' @param k2    An optional character vector with the name of a key for d2 to use instead of
#'              rownames. Defaults to NULL (join on rownames).
#' @param mode  One of 'left', 'right', 'inner', or 'full' specifying the kind of join to perform.
#'              Defaults to 'left' (add data from d2 to d1).
#' @param no.rn Logical, don't attempt to add rownmaes to result. Defaults to FALSE.
#'
#' @return A data frame with the results of the join.
#'
#' @export
#' @importFrom dplyr left_join full_join select
lexical_join <- function(
    d1, d2, k1=NULL, k2=NULL, mode=c( 'left','right', 'inner','full' ), no.rn=FALSE
) {
    d1$key_ <- if( is.null( k1 ) ) rownames( d1 ) else d1[[k1]]
    d2$key_ <- if( is.null( k2 ) ) rownames( d2 ) else d1[[k2]]
    mode = match.arg( mode )
    out <- switch( mode,
         # TODO remove dplyr
        left  = dplyr::left_join(  d1, d2, by='key_' ),
        right = dplyr::left_join(  d2, d1, by='key_' ),
        inner = dplyr::inner_join( d1, d2, by='key_' ),
        full  = dplyr::full_join(  d1, d2, by='key_' ),
        stop( "Unsupported join mode requested" )
    )
    if( no.rn || length( unique( out$key_ ) ) != length( out$key_ ) ) {
        if( !no.rn ) warning( 'Can\'t add non-unique key as rownames to result' )
    } else {
        rownames( out ) <- out$key_
    }
    out %<>% dplyr::select( -key_ )
    return( out )
}

#' Lexical sampling.
#'
#' Extracts samples from a TF vector such that the sampled terms account for the given \emph{theta}
#' proportion of the total mass in the TF vector.
#'
#' This function assumes that the TF vector is sorted in descending frequency following Zipf's
#' law such that sampling a relatively small number of entries from the top of the given vector
#' will suffice to cover the requested total mass percentage.
#'
#' This strategy will fail miserably if a) the given TF vector is not sorted in descending order
#' (though it can be sorted internally) or b) the entries in the given vector do not follow Zipf's
#' law (i.e. do not correspond to natural term frequencies). Use accordingly.
#'
#' This function also allows for censoring the given TF vector according to a filtering vector. In
#' this case, the 'univ' parameter controls whether mass coverage should be computed against the
#' mass of the full tf vector, or against the mass of the censored tf vector. Setting 'univ' to
#' TRUE will result in larger samples, as a higher number of elements is needed to account for the
#' same mass proportion.
#'
#' @param tf     A term frequency vector, typically sorted in descending order.
#' @param filter Logical vector of length equal to tf to exclude entries based on prior conditions
#'               (e.g. POS tag). Defaults to TRUE for all entries.
#' @param theta  A numeric value s.t. 0 < theta < 1 indicating the total tf mass that should be
#'               covered by the resulting sample.
#' @param univ   Logical indicating whether the covered mass should be computed against all entries
#'               in tf or only against those that pass the filter. Ignored if filter is NULL,
#'               defaults to FALSE (compute against censored tf).
#' @param sort   Logical indicating whether the given vector should be resorted. Use with caution,
#'               as most lexical sets in wspaces assume a fixed sort order.
#' @param value  Logical. Return term names instead of a logical vector. Requires 'tf' to be named.
#'
#' @return A logical vector \emph{of the same length as tf} indicating whether the corresponding
#'         entry in tf is included in the sample or not.
#'
#' @export
lexical_sample <- function( tf, filter=NULL, theta=.95, univ=FALSE, sort=FALSE, value=FALSE ) {
    if( theta >= 1 || theta <= 0 ) stop( sprintf(
        "Invalid theta %s given: Must be between 0 and 1 (both excluded)"
    ) )
    if( sort ) tf <- tf[ order( tf ) ]
    filter <- if( !is.null( filter ) ) {
        # TODO allow recylcing? int vectors?
        if( length( filter ) != length( tf ) ) stop( 'Wrong dimension for filter vector' )
        else filter
    } else rep( TRUE, length( tf ) )
    cum <- cumsum( ifelse( filter, tf, 0 ) ) / sum( if( univ ) tf else tf[filter] )
    return( ( cum <= theta ) & filter )
}

#' Extract term data from the given lexical dataset.
#'
#' Extracts term data from the given terms character vector from the given named vector, matrix or 
#' data.frame \link{lexical_dataset}.
#' 
#' Use this method to consistently extract data from different datasets for the same lexicon for a 
#' given subset of terms. Calls \link{term_idx} internally.
#'
#' @param d A \link{lexical_dataset}, i.e. a data frame or matrix with lexicon terms as row names.
#' @param terms A character vector with terms to get data for.
#'
#' @return A proper subset of \code{d} containing only entries in terms.
#'
#' @export
term_data <- function( d, terms ) {
    filter <- term_idx( d, terms )
    return( if( is.vector( d ) ) d[ filter ] else d[ filter, ] )
}

#' Get index vector for the given terms into the given lexical dataset.
#'
#' Generates an index vector for the given terms character vector into the elements or rows of the 
#' given named vector, matrix or lexical dataset.
#'
#' Use this method to consistently extract data from different datasets for the same lexicon for a 
#' given subset of terms.
#'
#' @param d Dataset into which to get index vector. Can be a named vector, named matrix or a
#'          \link{lexical_dataset}.
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
            "d must be a named vector, matrix or data.frame, got %s",
            paste( class( d ), collapse=", " )
        ) )
    }
    return( which( u %in% terms ) )
}

