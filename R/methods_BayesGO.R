
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
    
    cat( "Summary: Bayesian ontology fingerprint analysis results (class: BayesGO)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Model settings:\n")
    cat( "Number of genes to be analyzed: ", ncol(object@data$zmat.sub), "\n", sep="" )
    cat( "Number of GO terms to be analyzed: ", nrow(object@data$zmat.sub), "\n", sep="" )
    cat( "Maximum possible number of gene clusters: ", object@init$K, "\n", sep="" )
    cat( "Maximum possible number of GO term clusters: ", object@init$V, "\n", sep="" )
    cat( "--------------------------------------------------\n" )
    cat( "Data analysis results:\n")
    cat( "Median number of gene clusters: ", median(object@out$chain.phiSum), "\n", sep="" )
    cat( "Median number of GO term clusters: ", median(object@out$chain.lambdaSum), "\n", sep="" )
		cat( "Association between GO terms (rows) and genes (columns):\n" )
		for ( i in 1:nrow(emat) ) { 
			cat( "\t    ", paste( round(emat[i,]*vDigit)/vDigit, collapse="\t" ), "\n", sep="" )
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
  definition=function( object, thresGene=0.9, thresTerm=0.9 ) {

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
    
    pred.gene <- split( gene.name[ chain.gene.prop > thresGene ], chain.term[ chain.gene.prop > thresGene ] )
    pred.term <- split( term.name[ chain.term.prop > thresTerm ], chain.term[ chain.term.prop > thresTerm ] )

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
  definition=function( x, y, thresGene=0.9, thresTerm=0.9 ) {

    # extract objects
    
    chain.gene <- x@out$chain.col
    chain.term <- x@out$chain.row
    chain.emat <- x@out$chain.emat
    
    chain.gene.prop <- x@out$chain.col.prop
    chain.term.prop <- x@out$chain.row.prop
    
    # heatmap
    
    colrow <- chain.term
    colcol <- chain.gene
    
    colrow[ chain.term.prop < thresTerm ] <- NA
    colcol[ chain.gene.prop < thresGene ] <- NA

    E2 <- chain.emat
    rownames(E2) <- paste("term",1:nrow(chain.emat),sep="")
    colnames(E2) <- paste("gene",1:ncol(chain.emat),sep="")
    
    colrow2 <- as.data.frame( colrow, stringsAsFactors=FALSE )
    rownames(colrow2) <- paste("term",1:nrow(chain.emat),sep="")
    colnames(colrow2) <- "L   "
    
    colcol2 <- as.data.frame( colcol, stringsAsFactors=FALSE )
    rownames(colcol2) <- paste("gene",1:ncol(chain.emat),sep="")
    colnames(colcol2) <- "M  "
    
    pheatmap( E2, cluster_rows=FALSE, cluster_cols=FALSE, scale="none", 
      show_rownames = FALSE, show_colnames = FALSE,
      annotation_row = colrow2, annotation_col = colcol2 )
  }
)
