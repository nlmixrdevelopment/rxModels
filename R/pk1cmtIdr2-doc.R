##' pk1cmtIdr2 RxODE model
##'
##' @format An \emph{RxODE} model with 22 parameters, 3 ODE states, and 5 calc vars.
##'
##'\emph{Parameters (pk1cmtIdr2$params)}
##'
##' \describe{
##'   \item{popCl}{ (default=1)}
##'   \item{popV}{ (default=20)}
##'   \item{popKa}{ (default=1)}
##'   \item{bsvCl}{ (default=0)}
##'   \item{bsvV}{ (default=0)}
##'   \item{bsvKa}{ (default=0)}
##'   \item{popLagDepot}{ (default=0)}
##'   \item{popLagCentral}{ (default=0)}
##'   \item{popRateCentral}{ (default=0)}
##'   \item{popDurCentral}{ (default=0)}
##'   \item{bsvLagDepot}{ (default=0)}
##'   \item{bsvLagCentral}{ (default=0)}
##'   \item{bsvRateCentral}{ (default=0)}
##'   \item{bsvDurCentral}{ (default=0)}
##'   \item{bsvImax}{ (default=0)}
##'   \item{popImax}{ (default=0.9999)}
##'   \item{popIc50}{ (default=100)}
##'   \item{bsvIc50}{ (default=0)}
##'   \item{popKin}{ (default=9)}
##'   \item{bsvKin}{ (default=0)}
##'   \item{popKout}{ (default=0.3)}
##'   \item{bsvKout}{ (default=0)}
##'}
##'
##' \emph{State pk1cmtIdr2$state}
##'
##' \describe{
##'   \item{R}{ (=1)}
##'   \item{depot}{ (=2)}
##'   \item{central}{ (=3)}
##' }
##'
##' \emph{Calculated Variables pk1cmtIdr2$lhs}
##'
##' \describe{
##'   \item{cp}{}
##'   \item{Imax}{}
##'   \item{ic50}{}
##'   \item{kin}{}
##'   \item{kout}{}
##' }
##'
##' @seealso \code{\link[RxODE]{eventTable}}, \code{\link[RxODE]{et}}, \code{\link[RxODE]{rxSolve}}, \code{\link[RxODE]{RxODE}}
##' 
##' @examples
##' ## Showing the model code
##' summary(pk1cmtIdr2)
##'
"pk1cmtIdr2"
