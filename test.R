RcppParallel::setThreadOptions( numThreads = 1 )
require( wspaces )

K = 100

nd = sample( (9*K):(11*K), 1 )
md = sample( (9*K):(11*K), 1 )
pd = sample( (9*K):(11*K), 1 )

mode = 0

cat( sprintf( "nd: %d\n", nd ) )
cat( sprintf( "md: %d\n", md ) )
cat( sprintf( "pd: %d\n", pd ) )

ma <- matrix( runif( nd*md ), nrow=nd, ncol=md )
mb <- matrix( runif( pd*md ), nrow=pd, ncol=md )

obs0 <- innerp( ma[1,], mb[1,], -1, -1, mode )
exp0 <- ma[1,] %*% mb[1,]
cat( sprintf( "%s: max error: %10.8f\n", 'vector', abs( obs0 - exp0 ) ) )

obs1 <- marg_innerp( ma, ma, mode )
exp1 <- vapply( 1:nrow( ma ), function( x ) ma[x,] %*% ma[x,], 0.0 )
cat( sprintf( "%s: max error: %10.8f\n", 'marginal', max( abs( obs1 - exp1 ) ) ) )

obs2 <- full_innerp( ma, mb, mode )
exp2 <- ma %*% t( mb )
cat( sprintf( "%s: max error: %10.8f\n", 'full', max( abs( obs2 - exp2 ) ) ) )

obs3 <- self_innerp( ma, mode )
exp3 <- ma %*% t( ma )
cat( sprintf( "%s: max error: %10.8f\n", 'self', max( abs( obs3 - exp3 ) ) ) )

