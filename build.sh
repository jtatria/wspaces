#!/bin/bash

Rscript -e 'roxygen2::roxygenize( ".", roclets = c( "rd", "collate", "namespace" ) )'
R CMD INSTALL --no-multiarch --with-keep.source .
exit 0;
