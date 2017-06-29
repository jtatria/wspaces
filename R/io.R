# IO functions to read and write wspaces data to/from disk

# entry sizes by difference to account for vector overhead.
comp_t <- object.size( c( complex( 1 ), complex( 1 ) ) ) - object.size( complex( 1 ) )
real_t <- object.size( c( 1.0, 1.0 ) ) - object.size( 1.0 )

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
    dir = getwd(), lexicon = "lxcn.dsv", frequencies = "freq.dsv",
    pos_counts = "posc.dsv", cooccur = "cooc.bin",
    quiet = FALSE, attach = FALSE, env = .GlobalEnv
) {
    if( !quiet ) message( sprintf( "Loading corpus data from %s", dir ) )

    lxcn = read_lexicon( file.path( dir, lexicon ) )
    if( !quiet ) message( sprintf(
        'Lexicon contains %d terms. Total TF: %d; total DF: %d',
        nrow( lxcn ), sum( lxcn$tf ), sum( lxcn$df )
    ) )

    freq = read_frequencies( file.path( dir, frequencies ) )
    if( nrow( lxcn ) != nrow( freq ) ) corpus_io_error( "freq", nrow( lxcn ), nrow( freq ) )
    if( !quiet ) message( sprintf(
        'Frequency data recorded for %d corpus splits. Mean size: %4.2f; std. dev.: %4.2f',
        ncol( freq ), mean( colSums( freq ) ), sd( colSums( freq ) )
    ) )

    posc = read_pos_counts( file.path( dir, pos_counts ) )
    if( nrow( lxcn ) != nrow( posc ) ) corpus_io_error( "posc", nrow( lxcn ), nrow( posc ) )
    if( !quiet ) {
        message( sprintf(
            'POS counts recorded for %d tags. Tagset: [ %s ]. Mean pos conf.: %4.2f; std. dev.: %4.2f',
            ncol( posc ) - 2, paste( colnames( posc )[ 1:( ncol( posc ) - 2 ) ], collapse = ' ' ),
            mean( posc$POS_conf, na.rm = TRUE ), sd( posc$POS_conf, na.rm = TRUE )
        ) )
        values = summary( posc$POS )
        for( pos in names( values ) ) {
            message( sprintf( "%s:\t%d", pos, values[pos] ) )
        }
    }

    if( file.exists( file.path( dir, cooccur ) ) ) {
        cooc_file = file.path( dir, cooccur )
        if( !quiet ) {
            message( sprintf( "Loading cooccurrence counts found in %s", cooc_file ) )
        }
        cooc = read_cooccur( file.path( dir, cooccur ), lxcn = lxcn )
        if( nrow( lxcn ) != nrow( cooc ) ) corpus_io_error( "cooc", nrow( lxcn ), nrow( cooc ) )
        if( !quiet ) message( sprintf(
            "Cooccurrence counts loaded as %dx%d matrix with %d non-zero entries (%5.4f%% fill)",
            nrow( cooc ), ncol( cooc ), Matrix::nnzero( cooc ),
            100 *( Matrix::nnzero( cooc ) / ( as.numeric( nrow( cooc ) ) * as.numeric( ncol( cooc ) ) ) )
        ) )
    }

    corpus <- list(
        "lxcn" = lxcn,
        "freq" = freq,
        "posc" = posc
    )
    if( exists( 'cooc' ) ) corpus$cooc <- cooc;

    if( attach ) list2env( corpus, envir = env )
    return( corpus )
}

#' Read lexicon data from DSV file.
#'
#' Reads lexicon data from a DSV file. Lexicon data is recorded on disk as a series of
#' records containings each term's string form, its term total frequency and its total document
#' frequency. Additional columns may be used for additional corpus-wide term statistics. Terms
#' themselves are \emph{not} stored as data frame columns, but are instead used as row names.
#' Indices of terms in the vector produced by names( lexicon ) are used across all corpus datasets
#' as term keys.
#'
#' @param file   Location of the lexicon file on disk.
#' @param header Logical indicating whether the given file contains a header as first row. Default
#'               TRUE.
#' @param sep    Character indicating the field separator value. Deafault '@'.
#'
#' @return A data frame containing the list of terms as its row names, a tf column for total term
#' frequencies and a df column for total document frequency for each term. The returned data frame
#' may contain additional columns for extra corpus-wide statistics depending on the corpus source
#' implementation.
#'
#' @export
read_lexicon <- function( file, header = TRUE, sep = '@' ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    d <- read.table( file, quote = '', comment.char = '', sep = sep, header = header, row.names = 1 )
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
#' are corpus partitions containing each term's frequency counts for the respective corpus segment.
#' The returned data frame.
#'
#' @export
read_frequencies <- function( file, header = TRUE, sep = '@', lxcn = NULL ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    d <- read.table( file, quote = '', comment.char = '', sep = sep, header = header, row.names = 1 )
    d[ is.na( d ) ] <- 0
    return( d )
}

#' Read POS count data from DSV file.
#'
#' Reads POS count data from a DSV file. POS count data is recorded on disk as a series
#' of records containings each term's string form, followed by one column for each major POS group
#' in the tagset used by the upstream corpus source for its NLP pipeline. Terms themselves are
#' \emph{not} stored as data frame columns, but are instead used as \emph{row names}.
#'
#' @param file   Location of the POS count data file on disk.
#' @param header Logical indicating whether the given file contains a header as first row. Default
#'               TRUE.
#' @param sep    Character indicating the field separator value. Deafault '@'.
#'
#' @return A data frame containing the list of terms as its row names, and as many columns as there
#' are POS groups in the corresponding tagset, with each term's counts on each POS group.
#'
#' @export
read_pos_counts <- function( file, header = TRUE, sep = '@' ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    d <- read.table( file, header = header, sep = sep, quote = "", comment.char = "", row.names = 1 )
    d[ is.na( d ) ] <- 0
    d$POS <- counts_to_factor( d, drop = TRUE )
    d$POS_conf <- apply( d[,-length( d )], 1, function( x ) max( x ) / sum( x ) )
    return( d )
}

#' Read cooccurrence matrix from disk.
#'
#' Reads cooccurrence data from a file on disk, stored as a sparse tensor saved in
#' tuple format; i.e. a plain array of ( coord, value ) structs, where coord is itself an array
#' of integer coordinates for each dimension in the represented tensor and value is a floating point
#' number.
#'
#' The default corpus source implementation uses a term-term context definition and stores
#' cooccurrences in a rank-2 tensor, using 64 bit integers for coordinates and double precision
#' floats for values.
#'
#' The specific relationship between values in the returned tensor and the actual
#' distributional patterns for each term depends on the upstream corpus source cooccurrence
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
#' term-context observations, with no information about partial results.
#'
#' @param file Location of the cooccurrence matrix file on disk.
#'
#' @return A Matrix::sparseMatrix object containing occurrence counts for terms in contexts.
#'
#' @export
#' @importFrom Matrix rowSums colSums
read_cooccur <- function( file, lxcn = NULL, shrink = !is.null( lxcn ) ) {
    if( !file.exists( file ) ) stop( sprintf( "%s: file not found", file ) )
    m <- load_spm( file )
    if( !is.null( lxcn ) ) {
      if( nrow( m ) > nrow( lxcn ) ) {
        stop(
          paste(
            "Wrong number of dimensions in cooccurrence matrix:", nrow( m ),
            "is larger than the number of terms in the lexicon", nrow( lxcn )
          )
        )
      }
      rownames( m ) <- rownames( lxcn )[ 1:nrow( m ) ]
      colnames( m ) <- rownames( lxcn )[ 1:ncol( m ) ]
    } else {
      rownames( m ) <- as.character( 1:nrow( m ) )
      colnames( m ) <- as.character( 1:ncol( m ) )
    }
    if( shrink ) m <- m[ rowSums( m ) > 0, colSums( m ) > 0 ]
    return( m )
}

#' @export
read_vectors <- function( filename, n = n, dim = dim, type = c( 'real', 'complex', 'binary' ), combine = FALSE ) {
    if( missing( type ) ) type = 'real'

    size <- as.integer(
        switch( type,
            real    = real_t,
            complex = comp_t,
            binary  = stop( "Binary vectors not supported" )
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
