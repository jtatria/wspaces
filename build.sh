#!/bin/bash

( cd obo-java; mvn clean install );
ln -s obo-java/target/obo-1.0-SNAPSHOT.jar inst/java/;

Rscript -e 'roxygen2::roxygenize( ".", roclets = c( "rd", "collate", "namespace" ) )'
Rscript -e 'Rcpp::compileAttributes()'
R CMD INSTALL --no-multiarch --with-keep.source .
exit 0;
