##' Jones 2013 Human PBPK model
##'
##' @examples
##' library(RxODE)
##' Jones2013 %>% solve(et(amt=100,cmt=D,ii=4,until=24) %>% et(seq(0,24,length.out=100))) %>% plot
##'
##' @references 1: Jones H, Rowland-Yeo K. Basic concepts in
##'     physiologically based pharmacokinetic modeling in drug
##'     discovery and development. CPT Pharmacometrics Syst
##'     Pharmacol. 2013 Aug 14;2:e63. doi: 10.1038/psp.2013.41. PubMed
##'     PMID: 23945604; PubMed Central PMCID: PMC3828005.
##'
##' https://github.com/metrumresearchgroup/mrgsolve/blob/9ff7240da48ec133635a6f0db2143f13623f0aff/inst/models/pbpk.cpp
"Jones2013"


##' RxODE one compartment model (solved)
##' @examples
##' oral1cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot
"oral1cmt"

##' RxODE two compartment model (solved)
##' @examples
##' oral2cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot
"oral2cmt"

##' RxODE three compartment model (solved)
##' @examples
##' oral3cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot
"oral3cmt"




##'@importFrom RxODE RxODE
##' @export
RxODE::RxODE

##'@importFrom RxODE et
##' @export
RxODE::et

##'@importFrom RxODE add.sampling
##' @export
RxODE::add.sampling

##'@importFrom RxODE add.dosing
##' @export
RxODE::add.dosing

##'@importFrom RxODE eventTable
##' @export
RxODE::eventTable

##'@importFrom RxODE rxClean
##' @export
RxODE::rxClean

##'@importFrom RxODE rxSolve
##' @export
RxODE::rxSolve


##'@importFrom RxODE rxControl
##' @export
RxODE::rxControl


