##' pk3cmtIdr2 RxODE model
##'
##' @format An \emph{RxODE} model with 31 parameters, 3 ODE states, and 5 calc vars.
##'
##'\emph{Parameters (pk3cmtIdr2$params)}
##'
##' \describe{
##'   \item{popCl}{ (default=1)}
##'   \item{popV}{ (default=20)}
##'   \item{popKa}{ (default=1)}
##'   \item{popVp}{ (default=10)}
##'   \item{popQ}{ (default=2)}
##'   \item{popQ2}{ (default=2)}
##'   \item{popVp2}{ (default=100)}
##'   \item{bsvCl}{ (default=0)}
##'   \item{bsvV}{ (default=0)}
##'   \item{bsvKa}{ (default=0)}
##'   \item{bsvVp}{ (default=0)}
##'   \item{bsvQ}{ (default=0)}
##'   \item{bsvQ2}{ (default=0)}
##'   \item{bsvVp2}{ (default=0)}
##'   \item{popLagDepot}{ (default=0)}
##'   \item{popLagCentral}{ (default=0)}
##'   \item{popRateCentral}{ (default=0)}
##'   \item{popDurCentral}{ (default=0)}
##'   \item{bsvLagDepot}{ (default=0)}
##'   \item{bsvLagCentral}{ (default=0)}
##'   \item{bsvRateCentral}{ (default=0)}
##'   \item{bsvDurCentral}{ (default=0)}
##'   \item{pi}{ (default=0)}
##'   \item{bsvImax}{ (default=0.9999)}
##'   \item{popImax}{ (default=100)}
##'   \item{popIc50}{ (default=0)}
##'   \item{bsvIc50}{ (default=9)}
##'   \item{popKin}{ (default=0)}
##'   \item{bsvKin}{ (default=0.3)}
##'   \item{popKout}{ (default=0)}
##'   \item{bsvKout}{ (default=3.14159265358979)}
##'}
##'
##' \emph{State pk3cmtIdr2$state}
##'
##' \describe{
##'   \item{R}{ (=1)}
##'   \item{depot}{ (=2)}
##'   \item{central}{ (=3)}
##' }
##'
##' \emph{Calculated Variables pk3cmtIdr2$lhs}
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
##' summary(pk3cmtIdr2)
##'
"pk3cmtIdr2"
