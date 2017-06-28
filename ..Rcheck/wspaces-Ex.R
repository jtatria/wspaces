pkgname <- "wspaces"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('wspaces')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("counts_to_factor")
### * counts_to_factor

flush(stderr()); flush(stdout())

### Name: counts_to_factor
### Title: Build a factor from max/min values in a series of variables
### Aliases: counts_to_factor

### ** Examples

d <- as.data.frame( matrix( rnorm( 300 ), nrow = 100, ncol = 3 ) )
names( d ) <- c( 'red', 'blue', 'green' )
count_to_factor( d )




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
