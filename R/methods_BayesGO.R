
#' @import stats pheatmap
#' @exportMethod predict plot

#' @title Summary of bayesGO results
#' 
#' @description Summary of bayesGO results.
#' 
#' @aliases show,BayesGO-method
#' 
#' @param object output of bayesGO function

setMethod(
  f="show",
  signature="BayesGO",
  definition=function( object ) {
    
    emat <- object@out$chain.enrichment
    vDigit <- 100
    
    chain.gene <- object@out$chain.col
    chain.term <- object@out$chain.row
    
    cat( "Summary: Bayesian ontology fingerprint analysis results (class: BayesGO)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Model settings:\n")
    cat( "Number of genes to be analyzed: ", ncol(object@data$zmat.sub), "\n", sep="" )
    cat( "Number of GO terms to be analyzed: ", nrow(object@data$zmat.sub), "\n", sep="" )
    cat( "Maximum possible number of gene clusters: ", object@init$K, "\n", sep="" )
    cat( "Maximum possible number of GO term clusters: ", object@init$V, "\n", sep="" )
    cat( "--------------------------------------------------\n" )
    cat( "Data analysis results:\n")
    #cat( "Median number of gene clusters: ", median(object@out$chain.phiSum), "\n", sep="" )
    #cat( "Median number of GO term clusters: ", median(object@out$chain.lambdaSum), "\n", sep="" )
    cat( "Median number of gene clusters: ", length(unique(chain.gene)), "\n", sep="" )
    cat( "Median number of GO term clusters: ", length(unique(chain.term)), "\n", sep="" )
		cat( "Association between GO terms (rows) and genes (columns):\n" )
		for ( i in 1:nrow(emat) ) { 
			cat( "\t    ", paste( format( round(emat[i,]*vDigit)/vDigit, nsmall=2 ), collapse="\t" ), "\n", sep="" )
		}
    cat( "--------------------------------------------------\n" )
  }
)

#' @title Gene and GO term cluster prediction
#'
#' @description Return gene and GO term clustering results and their association patterns.
#'
#' @aliases predict,BayesGO-method
#'
#' 
#' @param object output of bayesGO function
#' @param thresGene threshold for high confidence assignment to gene clusters, which is between 0 and 1 (default: 0.9)
#' @param thresTerm threshold for high confidence assignment to GO term clusters, which is between 0 and 1

setMethod(
  f="predict",
  signature="BayesGO",
  definition=function( object, thresGene=0, thresTerm=0 ) {

    # extract objects
    
    gene.name <- colnames(object@data$zmat.sub)
    if ( is.null(gene.name) ) {
      gene.name <- 1:ncol(object@data$zmat.sub)
    }
    term.name <- rownames(object@data$zmat.sub)
    if ( is.null(term.name) ) {
      term.name <- 1:nrow(object@data$zmat.sub)
    }
    
    chain.gene <- object@out$chain.col
    chain.term <- object@out$chain.row
    
    chain.gene.prop <- object@out$chain.col.prop
    chain.term.prop <- object@out$chain.row.prop
    
    # gene and GO term clustering results

    mat.gene <- data.frame( gene.name, chain.gene.prop )[ chain.gene.prop > thresGene, ]
    mat.term <- data.frame( term.name, chain.term.prop )[ chain.term.prop > thresTerm, ]
    
    pred.gene <- split( mat.gene, chain.gene[ chain.gene.prop > thresGene ] )
    pred.term <- split( mat.term, chain.term[ chain.term.prop > thresTerm ] )
    
    pred.gene <- lapply( pred.gene, function(x) x[ order( x[,2], decreasing=TRUE ), ] )
    pred.term <- lapply( pred.term, function(x) x[ order( x[,2], decreasing=TRUE ), ] )

    return(list(
      geneClust = pred.gene,
      termClust = pred.term,
      association=object@out$chain.enrichment
    ))
  }
)

#' @title Plotting function
#'
#' @description Plot a heatmap of association between GO terms and genes.
#'
#' @aliases plot,BayesGO-method
#'
#' @param x output of bayesGO function
#' @param y (unused argument)
#' @param thresGene threshold for high confidence assignment to gene clusters, which is between 0 and 1 (default: 0.9)
#' @param thresTerm threshold for high confidence assignment to GO term clusters, which is between 0 and 1

setMethod(
  f="plot",
  signature=c("BayesGO","missing"),
  definition=function( x, y, thresGene=0, thresTerm=0 ) {

    # extract objects
    
    chain.gene <- x@out$chain.col
    chain.term <- x@out$chain.row
    chain.emat <- x@out$chain.emat
    
    chain.gene.uniq <- sort(unique(chain.gene))
    chain.term.uniq <- sort(unique(chain.term))
    
    n.gene.cluster <- length(chain.gene.uniq)
    n.term.cluster <- length(chain.term.uniq)
    
    for ( i in 1:n.gene.cluster ) {
      chain.gene[ chain.gene == chain.gene.uniq[i] ] <- i
    }
    for ( i in 1:n.term.cluster ) {
      chain.term[ chain.term == chain.term.uniq[i] ] <- i
    }
    
    chain.gene.prop <- x@out$chain.col.prop
    chain.term.prop <- x@out$chain.row.prop
    
    # heatmap
    
    colrow <- chain.term
    colcol <- chain.gene
    
    colrow[ chain.term.prop < thresTerm ] <- NA
    colcol[ chain.gene.prop < thresGene ] <- NA

    E2 <- chain.emat
    E2[ is.na(E2) ] <- 0.5
    rownames(E2) <- paste("term",1:nrow(chain.emat),sep="")
    colnames(E2) <- paste("gene",1:ncol(chain.emat),sep="")
    
    colrow2 <- as.data.frame( colrow, stringsAsFactors=FALSE )
    rownames(colrow2) <- paste("term",1:nrow(chain.emat),sep="")
    colnames(colrow2) <- "L   "
    
    colcol2 <- as.data.frame( colcol, stringsAsFactors=FALSE )
    rownames(colcol2) <- paste("gene",1:ncol(chain.emat),sep="")
    colnames(colcol2) <- "M  "
    
    E2 <- E2[ order(colrow2[,1]), order(colcol2[,1]) ]
    
    pheatmap( E2, cluster_rows=FALSE, cluster_cols=FALSE, scale="none", 
      show_rownames = FALSE, show_colnames = FALSE,
      annotation_row = colrow2, annotation_col = colcol2 )
  }
)
