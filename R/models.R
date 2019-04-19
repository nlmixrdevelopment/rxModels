## Need to have at least 1 package export
##' @useDynLib rxModels, .registration=TRUE
##' @importFrom RxODE RxODE
##' @export
RxODE::RxODE

##' @importFrom RxODE et
##' @export
RxODE::et

##' @importFrom RxODE etRep
##' @export
RxODE::etRep

##' @importFrom RxODE etSeq
##' @export
RxODE::etSeq

##' @importFrom RxODE as.et
##' @export
RxODE::as.et

##' @importFrom RxODE eventTable
##' @export
RxODE::eventTable

##' @importFrom RxODE add.dosing
##' @export
RxODE::add.dosing

##' @importFrom RxODE add.sampling
##' @export
RxODE::add.sampling

##' @importFrom RxODE rxControl
##' @export
RxODE::rxControl

##' @importFrom RxODE rxClean
##' @export
RxODE::rxClean

##' @importFrom RxODE rxUse
##' @export
RxODE::rxUse

##' @importFrom RxODE rxShiny
##' @export
RxODE::rxShiny

##' @importFrom RxODE genShinyApp.template
##' @export
RxODE::genShinyApp.template

##' @importFrom RxODE cvPost
##' @export
RxODE::cvPost

# This is actually from `magrittr` but allows less imports
##' @importFrom RxODE %>%
##' @export
RxODE::`%>%`


##' @importFrom RxODE rxSolve
##' @export
RxODE::rxSolve
