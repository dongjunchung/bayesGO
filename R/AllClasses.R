
#' An S4 class to represent bayesGO model fitting results.
#'
#' @slot data data
#' @slot param model parameter settings
#' @slot init model initialization
#' @slot fit rjags fitting results
#' @slot out clustering and association analysis results

setClass( Class="BayesGO",
    representation=representation(
        data="list",
        param="list",
        init="list",
        fit="list",
        out="list"
    )
)
