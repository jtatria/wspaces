# wspaces

Construction and analysis of semantic space models.

This package has been developed as part of my doctoral dissertation. It provides facilities for the
compilation of distributional information as cooccurrence matrices, and the projection of the
distributional spaces described by them into different representations of semantic spaces:
graph-based representations that we call semantic networks, and real- or complex- valued vector
spaces.

The package has been developed with a rather obsessive focus on speed and generality/modularity.
It includes facilities for accessing and manipulating a Lucene index for the compilation of
distributional information, but this is not required as the package also provides I/O facilities
for reading and writing corpus data from external sources (using its own custom format for now).

Most numerical procedures are coded (and optimized) in C++ using Rcpp, but the computation model
used has been designed to allow for custom functions defined in user code (though this has not been
optimized yet and is much slower than the native functions provided in the package).

WARNING: This package is under heavy development and refactoring, things are broken and will stay
broken for a while while I decide on a distribution model for the java back-end and sort out names,
stabilize an API and finish related administratrivia.

Consider it pre-alpha, and use with care (or better yet: not at all!).

## Java binary blob

As mentioned above, the package uses an external java library to communicate with a Lucene index.
I have not yet sorted out the best way to include (and build!) the necessary java dependencies, so
this external library is now included as a binary blob. Package may not build correctly, as I have
not figured out how to include the java toolchain (with the corresponding dependency hell) in the
R compilation process.

# Author

J. T. Atria (jtatria@nomoi.org)

# License

GPL v3 or later.
