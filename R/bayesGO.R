#' @title bayesGO function
#'
#' @description Fit the bayesGO model
#'
#' @export
#' @import rjags stats label.switching
#' @param pmat A matrix of hypergeometric p-values
#' @param nRow Number of rows to analyze (default: use all rows)
#' @param nCol Number of columns to analyze (default: use all columns)
#' @param Vmax Maximum number of row clusters (default: 0.1 * nRow)
#' @param Kmax Maximum number of column clusters (default: 0.1 * nCol)
#' @param nBurnin MCMC: Number of burn-in iterations (default: 20000)
#' @param nMain MCMC: Number of main iterations (default: 10000)
#' @param thin Thinning (default: 10)
#' @param seed Random seed (default: 12345)
#'
#' @examples
#' data(simul)
#' fit.bayesGO <- bayesGO( simul, nBurnin=10, nMain=10, thin=1 )
#' fit.bayesGO
#' predict( fit.bayesGO, thresGene=0.9, thresTerm=0.9 )
#' plot( fit.bayesGO, thresGene=0.9, thresTerm=0.9 )

bayesGO <- function( pmat, nRow=NA, nCol=NA, Vmax=NA, Kmax=NA, 
  nBurnin=20000, nMain=10000, thin=10, seed=12345 ) {
    
    #library(rjags); pmat=simul; nRow=NA; nCol=NA; Vmax=NA; Kmax=NA; nBurnin=10; nMain=10; thin=10; seed=12345
    
    # gene list
  
    pmat.sub <- pmat
    zmat.sub <- qnorm(pmat.sub)
    zmat.sub[ zmat.sub == Inf ] <- NA
    zmat.sub[ zmat.sub == -Inf ] <- min( zmat.sub[ zmat.sub != -Inf ] )
    
    # sort rows and columns based on numbers of missing cells & subset rows and columns
    
    if ( !is.na(nRow) ) {
      nNArow <- apply( zmat.sub, 1, function(x) length(which( is.na(x) )) )
      zmat.sub <- zmat.sub[ nNArow <= sort(nNArow)[nRow], ]
    }
    if ( !is.na(nCol) ) {
      nNAcol <- apply( zmat.sub, 2, function(x) length(which( is.na(x) )) )
      zmat.sub <- zmat.sub[ , nNAcol <= sort(nNAcol)[nCol] ]
    }
    
    # convert matrices to vectors
    
    rindmat <- matrix( rep( 1:nrow(zmat.sub), ncol(zmat.sub) ), nrow=nrow(zmat.sub) )
    cindmat <- matrix( rep( 1:ncol(zmat.sub), nrow(zmat.sub) ), ncol=ncol(zmat.sub), byrow=TRUE )
    
    zvec.sub <- as.vector(zmat.sub[ !is.na(zmat.sub) ])
    rind.sub <- as.vector(rindmat[ !is.na(zmat.sub) ])
    cind.sub <- as.vector(cindmat[ !is.na(zmat.sub) ])
    
    # calculate K & V
    
    if ( !is.na(Kmax) ) {
      K <- Kmax
    } else {
      K <- round( ncol(zmat.sub) / 10 )
    }
    
    if ( !is.na(Vmax) ) {
      V <- Vmax
    } else {
      V <- round( nrow(zmat.sub) / 10 )
    }
    
    # gene clustering analysis for initialization of M
    
    cormat <- matrix( NA, ncol(zmat.sub), ncol(zmat.sub) )
    
    for ( i in 1:ncol(zmat.sub) ) {
      for ( j in 1:ncol(zmat.sub) ) {
        cor.ij <- cor(na.omit( cbind( zmat.sub[ , i ], zmat.sub[ , j ] ) ))[1,2]
        if ( !is.na(cor.ij) ) {
          cormat[i,j] <- cor.ij
        } else {
          cormat[i,j] <- 0
        }
      }
    }
    
    Minit <- cutree( hclust(as.dist( 1 - cormat )), K )
    
    # term clustering analysis for initialization of L
    
    cormat <- matrix( NA, nrow(zmat.sub), nrow(zmat.sub) )
    
    for ( i in 1:nrow(zmat.sub) ) {
      for ( j in 1:nrow(zmat.sub) ) {
        cor.ij <- cor(na.omit( cbind( zmat.sub[ i, ], zmat.sub[ j, ] ) ))[1,2]
        if ( !is.na(cor.ij) ) {
          cormat[i,j] <- cor.ij
        } else {
          cormat[i,j] <- 0
        }
      }
    }
    
    Linit <- cutree( hclust(as.dist( 1 - cormat )), V )
    
    # generate JAGS data 
    
    jdata <- list(
      N = length(zvec.sub),
      T = nrow(zmat.sub),
      G = ncol(zmat.sub),
      Vmax = V,
      Kmax = K,
      paramEpsilon = c( 1, 1 ),
      paramBeta0 = c( 0.1, 0.1 ),
      paramEta = c( 1, 1 ),
      paramAlpha0 = c( 0.1, 0.1 ),
      paramTheta01 = c( 0.1, 0.1 ),
      paramTheta02 = c( 0.1, 0.1 ),
      paramZeta0 = c( 0, 0.0001 ),
      paramXi0 = c( 2, 1 ),
      paramZeta1 = c( 0, 0.0001 ),
      paramXi1 = c( 2, 1 ),
      paramTau0 = c( 2, 1 ),
      paramTau1 = c( 2, 1 ),
      Y = zvec.sub,
      rind = rind.sub,
      cind = cind.sub
    )
    
    # generate JAGS initial value
    
    code.jags <- system.file( file.path("JAGS","JAGS_model.txt"), package="bayesGO")
    
    message("---------------------\n")
    message("Initializing rjags...\n")
    message("---------------------\n")
    
    jinit <- list(
      M = Minit,
      L = Linit,
      epsilon = 0.5,
      lambda = rep(1,K),
      beta0 = 1,
      eta = 0.5,
      phi = rep(1,K),
      alpha0 = 1,
      mu0 = rep( -1, ncol(zmat.sub) ),
      tau0 = rep( 1, ncol(zmat.sub) ),
      mu1 = rep( -10, ncol(zmat.sub) ),
      tau1 = rep( 0.2, ncol(zmat.sub) )
    )
    
    # initialize rjags
    
    set.seed(seed)
    
    jmodel <- jags.model(
      file=code.jags,
      data=jdata,
      inits=jinit,
      n.chains=2,
      n.adapt=0
    )
    
    # burn-in iteration
    
    message("---------------------\n")
    message("MCMC: Burn-in\n")
    message("---------------------\n")
    
    set.seed(seed)
    
    update( jmodel, nBurnin )
    
    # main MCMC iteration
    
    message("---------------------\n")
    message("MCMC: Main iteration\n")
    message("---------------------\n")
    
    set.seed(seed)
    
    out <- coda.samples( jmodel,
      variable.names=c("E","M","L","lambdaSum","phiSum"),
      n.iter=nMain, thin=thin )
    
    # post-processing: extract objects
    
    chain1.phiSum <- out[[1]][ , grep( "phiSum", colnames(out[[1]]) ) ]
    chain1.lambdaSum <- out[[1]][ , grep( "lambdaSum", colnames(out[[1]]) ) ]
    chain1.M <- out[[1]][ , grep( "M", colnames(out[[1]]) ) ]
    chain1.L <- out[[1]][ , grep( "L", colnames(out[[1]]) ) ]
    chain1.E <- out[[1]][ , grep( "E", colnames(out[[1]]) ) ]
    
    chain2.phiSum <- out[[2]][ , grep( "phiSum", colnames(out[[2]]) ) ]
    chain2.lambdaSum <- out[[2]][ , grep( "lambdaSum", colnames(out[[2]]) ) ]
    chain2.M <- out[[2]][ , grep( "M", colnames(out[[2]]) ) ]
    chain2.L <- out[[2]][ , grep( "L", colnames(out[[2]]) ) ]
    chain2.E <- out[[2]][ , grep( "E", colnames(out[[2]]) ) ]
  
    # post-processing: renumbering labels
    
    chain1.new.M <- chain1.M
    for ( i in 1:nrow(chain1.M) ) {
      label.i <- names(sort( table(chain1.M[i,]), decreasing=TRUE ))
      for ( j in 1:length(label.i) ) {
        chain1.new.M[ i, chain1.M[i,] == label.i[j] ] <- j
      }
    }
    
    chain2.new.M <- chain2.M
    for ( i in 1:nrow(chain2.M) ) {
      label.i <- names(sort( table(chain2.M[i,]), decreasing=TRUE ))
      for ( j in 1:length(label.i) ) {
        chain2.new.M[ i, chain2.M[i,] == label.i[j] ] <- j
      }
    }
    
    chain1.new.L <- chain1.L
    for ( i in 1:nrow(chain1.L) ) {
      label.i <- names(sort( table(chain1.L[i,]), decreasing=TRUE ))
      for ( j in 1:length(label.i) ) {
        chain1.new.L[ i, chain1.L[i,] == label.i[j] ] <- j
      }
    }
    
    chain2.new.L <- chain2.L
    for ( i in 1:nrow(chain2.L) ) {
      label.i <- names(sort( table(chain2.L[i,]), decreasing=TRUE ))
      for ( j in 1:length(label.i) ) {
        chain2.new.L[ i, chain2.L[i,] == label.i[j] ] <- j
      }
    }
    
    # post-processing: ECR algorithm
    
    set.seed(seed)
  
    chain.phiSum <- c( chain1.phiSum, chain2.phiSum )
    nclust.col <- as.numeric(names(table(chain.phiSum))[ which.max(table(chain.phiSum)) ])
    
    chain.lambdaSum <- c( chain1.lambdaSum, chain2.lambdaSum )
    nclust.row <- as.numeric(names(table(chain.lambdaSum))[ which.max(table(chain.lambdaSum)) ])
    
    zpivot.col <- kmeans( t(rbind( chain1.new.M, chain2.new.M )), 
      centers=nclust.col, iter.max=1000, nstart=100 )$cluster
    
    zpivot.row <- kmeans( t(rbind( chain1.new.L, chain2.new.L )), 
      centers=nclust.row, iter.max=1000, nstart=100 )$cluster
    
    perm.col <- ecr( zpivot=zpivot.col, 
      z=rbind( chain1.new.M, chain2.new.M ), K=max( max(chain1.new.M), max(chain2.new.M) ) )
    
    chain1.final.M <- chain1.new.M
    for ( i in 1:nrow(chain1.final.M) ) {
      for ( j in 1:ncol(perm.col[[1]]) ) {
        chain1.final.M[ i, chain1.new.M[i,] == perm.col[[1]][i,j] ] <- j
      }
    }
    
    chain2.final.M <- chain2.new.M
    for ( i in 1:nrow(chain2.final.M) ) {
      for ( j in 1:ncol(perm.col[[1]]) ) {
        chain2.final.M[ i, chain2.new.M[i,] == perm.col[[1]][(nrow(chain1.final.M)+i),j] ] <- j
      }
    }
    
    perm.row <- ecr( zpivot=zpivot.row, 
      z=rbind( chain1.new.L, chain2.new.L ), K=max( max(chain1.new.L), max(chain2.new.L) ) )
    
    chain1.final.L <- chain1.new.L
    for ( i in 1:nrow(chain1.final.L) ) {
      for ( j in 1:ncol(perm.row[[1]]) ) {
        chain1.final.L[ i, chain1.new.L[i,] == perm.row[[1]][i,j] ] <- j
      }
    }
    
    chain2.final.L <- chain2.new.L
    for ( i in 1:nrow(chain2.final.L) ) {
      for ( j in 1:ncol(perm.row[[1]]) ) {
        chain2.final.L[ i, chain2.new.L[i,] == perm.row[[1]][(nrow(chain1.final.L)+i),j] ] <- j
      }
    }
  
    # cluster assignment & frequency
    
    chain.col <- apply( rbind( chain1.final.M, chain2.final.M ), 2, 
      function(x) names(table(x))[ which.max(table(x))] )
    
    chain.row <- apply( rbind( chain1.final.L, chain2.final.L ), 2, 
      function(x) names(table(x))[ which.max(table(x))] )
    
    chain.col.prop <- apply( rbind( chain1.final.M, chain2.final.M ), 2, 
      function(x) length(which( x == names(table(x))[ which.max(table(x)) ] )) ) / 
      nrow(rbind( chain1.final.M, chain2.final.M ))
    
    chain.row.prop <- apply( rbind( chain1.final.L, chain2.final.L ), 2, 
      function(x) length(which( x == names(table(x))[ which.max(table(x)) ] )) ) / 
      nrow(rbind( chain1.final.L, chain2.final.L ))
  
    # enrichment matrix
    
    chain.evec <- apply( rbind( chain1.E, chain2.E ), 2, mean )
    chain.emat <- matrix( NA, nrow(zmat.sub), ncol(zmat.sub) )
    for ( i in 1:length(chain.evec) ) {
      chain.emat[ rind.sub[i], cind.sub[i] ] <- chain.evec[i]
    }
    
    chain.row.sel <- unique(chain.row)
    chain.col.sel <- unique(chain.col)
    
    chain.enrichment <- matrix( NA, length(chain.row.sel), length(chain.col.sel) )
    for ( i in 1:length(chain.row.sel) ) {
      for ( j in 1:length(chain.col.sel) ) {
        chain.enrichment[i,j] <- mean(round( 
          chain.emat[ chain.row == chain.row.sel[i], chain.col == chain.col.sel[j] ] ), 
        na.rm=TRUE )
      }
    }
    
    rownames(chain.enrichment) <- chain.row.sel
    colnames(chain.enrichment) <- chain.col.sel
    
    chain.enrichment <- chain.enrichment[ order(rownames(chain.enrichment)), order(colnames(chain.enrichment)) ]
    
    # export results
  
    methods::new( "BayesGO",
      data = list( pmat=pmat, zmat.sub=zmat.sub, 
        zvec.sub=zvec.sub, rind.sub=rind.sub, cind.sub=cind.sub ),
      param = list( nRow=nRow, nCol=nCol, Vmax=Vmax, Kmax=Kmax, 
        nBurnin=nBurnin, nMain=nMain, thin=thin, seed=seed ),
      init = list( V=V, K=K, Minit=Minit, Linit=Linit ),
      fit = list( jinit=jinit, jmodel=jmodel, out=out ),
      out = list( chain.row=chain.row, chain.col=chain.col, 
        chain.row.prop=chain.row.prop, chain.col.prop=chain.col.prop,
        chain.emat=chain.emat, chain.enrichment=chain.enrichment, 
        chain.lambdaSum=chain.lambdaSum, chain.phiSum=chain.phiSum )
  )
}
