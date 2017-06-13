
#include "SparseMatrix.h"
#include <R.h>
#include <Rinternals.h>

SparseMatrix * sexp_to_spm( SEXP x, SEXP i, SEXP j ) {
    int lx = length( x ), li = length( i ), lj = length( j );
    if( lx != li || lx != lj || li != lj ) {
        error( "i, j, x vectors must be of the same length!" );
    }

    size_t size = sizeof( SparseMatrix )
        + li * sizeof( int )
        + lj * sizeof( int )
        + lx * sizeof( double )
    ;

    SparseMatrix *spm;
    if( ( spm = (SparseMatrix *) malloc( size ) ) == 0 ) {
        fprintf( stderr, "Not enough memory!\n" );
    }

    spm->entry_count = lx;
    spm->row_indices = INTEGER( i );
    spm->col_indices = INTEGER( j );
    spm->entries = REAL( x );

    return spm;
}

SEXP spm_to_sexp( SparseMatrix *spm ) {

    SEXP ivec = PROTECT( allocVector( INTSXP, spm->entry_count ) );
    SEXP jvec = PROTECT( allocVector( INTSXP, spm->entry_count ) );
    SEXP xvec = PROTECT( allocVector( REALSXP, spm->entry_count ) );

    int    *pi = INTEGER( ivec );
    int    *pj = INTEGER( jvec );
    double *px = REAL( xvec );

    for( int i = 0; i < spm->entry_count; i++ ) {
        pi[i] = spm->row_indices[i];
        pj[i] = spm->col_indices[i];
        px[i] = spm->entries[i];
    }

    spm_free( spm );

    SEXP out = PROTECT( allocVector( VECSXP, 3 ) );
    SET_VECTOR_ELT( out, 0, ivec );
    SET_VECTOR_ELT( out, 1, jvec );
    SET_VECTOR_ELT( out, 2, xvec );

    SEXP nms = PROTECT( allocVector( STRSXP, 3 ) );
    const char *names[] = { "i", "j", "x", "" };
    for( int i = 0; i < 3; i++ ) SET_STRING_ELT( nms, i,  mkChar( names[i] ) );
    setAttrib( out, R_NamesSymbol, nms );

    UNPROTECT( 5 );

    return out;
}

SEXP read_spm( SEXP file ) {
    const char *fname = CHAR( asChar( file ) );
    SparseMatrix *spm = spm_read( fname );
    return spm_to_sexp( spm );
}
