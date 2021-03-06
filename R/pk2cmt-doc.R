##' pk2cmt RxODE model
##'
##' @format An \emph{RxODE} model with 18 parameters, 2 ODE states, and 1 calc vars.
##'
##'\emph{Parameters (pk2cmt$params)}
##'
##' \describe{
##'   \item{popCl}{ (default=1)}
##'   \item{popV}{ (default=20)}
##'   \item{popKa}{ (default=1)}
##'   \item{popVp}{ (default=10)}
##'   \item{popQ}{ (default=2)}
##'   \item{bsvCl}{ (default=0)}
##'   \item{bsvV}{ (default=0)}
##'   \item{bsvKa}{ (default=0)}
##'   \item{bsvVp}{ (default=0)}
##'   \item{bsvQ}{ (default=0)}
##'   \item{popLagDepot}{ (default=0)}
##'   \item{popLagCentral}{ (default=0)}
##'   \item{popRateCentral}{ (default=0)}
##'   \item{popDurCentral}{ (default=0)}
##'   \item{bsvLagDepot}{ (default=0)}
##'   \item{bsvLagCentral}{ (default=0)}
##'   \item{bsvRateCentral}{ (default=0)}
##'   \item{bsvDurCentral}{ (default=0)}
##'}
##'
##' \emph{State pk2cmt$state}
##'
##' \describe{
##'   \item{depot}{ (=1)}
##'   \item{central}{ (=2)}
##' }
##'
##' \emph{Calculated Variables pk2cmt$lhs}
##'
##' \describe{
##'   \item{cp}{}
##' }
##'
##' \emph{Model Code}
##'
##' @seealso \code{\link[RxODE]{eventTable}}, \code{\link[RxODE]{et}}, \code{\link[RxODE]{rxSolve}}, \code{\link[RxODE]{RxODE}}
##'
##' @examples
##' ## To see the model code:
##' summary(pk2cmt)
##'
##' pk2cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot # Oral
##' pk2cmt %>% solve(et(amt=100,cmt=central,timeUnits="hr")) %>% plot # IV
##'
"pk2cmt"
