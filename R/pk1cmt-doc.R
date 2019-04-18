##' pk1cmt RxODE model
##'
##' @format An \emph{RxODE} model with 14 parameters, 2 ODE states, and 1 calc vars.
##'
##'\emph{Parameters (pk1cmt$params)}
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
##'}
##'
##' \emph{State pk1cmt$state}
##'
##' \describe{
##'   \item{depot}{ (=1)}
##'   \item{central}{ (=2)}
##' }
##'
##' \emph{Calculated Variables pk1cmt$lhs}
##'
##' \describe{
##'   \item{cp}{}
##' }
##'
##' \emph{Model Code}
##'
##' RxODE({
##'     popCl = 1
##'     popV = 20
##'     popKa = 1
##'     bsvCl = 0
##'     bsvV = 0
##'     bsvKa = 0
##'     cl ~ popCl * exp(bsvCl)
##'     v ~ popV * exp(bsvV)
##'     ka ~ popKa * exp(bsvKa)
##'     popLagDepot = 0
##'     popLagCentral = 0
##'     popRateCentral = 0
##'     popDurCentral = 0
##'     bsvLagDepot = 0
##'     bsvLagCentral = 0
##'     bsvRateCentral = 0
##'     bsvDurCentral = 0
##'     rx_ka ~ ka
##'     rx_rate ~ popRateCentral * exp(bsvRateCentral)
##'     rx_dur ~ popDurCentral * exp(bsvDurCentral)
##'     rx_tlag ~ popLagDepot * exp(bsvLagDepot)
##'     rx_tlag2 ~ popLagCentral * exp(bsvLagCentral)
##'     rx_F ~ 1
##'     rx_F2 ~ 1
##'     rx_v ~ v
##'     rx_k ~ cl/v
##'     rx_alpha ~ rx_k
##'     rx_A ~ rx_ka/(rx_ka - rx_alpha)/rx_v
##'     rx_A2 ~ 1/rx_v
##'     rx_beta ~ 0
##'     rx_B ~ 0
##'     rx_B2 ~ 0
##'     rx_gamma ~ 0
##'     rx_C ~ 0
##'     rx_C2 ~ 0
##'     cp = solveLinB(rx__PTR__, t, 0, rx_A, rx_A2, rx_alpha, rx_B,
##'         rx_B2, rx_beta, rx_C, rx_C2, rx_gamma, rx_ka, rx_tlag,
##'         rx_tlag2, rx_F, rx_F2, rx_rate, rx_dur)
##' })
##'
##' @seealso \code{\link[RxODE]{eventTable}}, \code{\link[RxODE]{et}}, \code{\link[RxODE]{rxSolve}}, \code{\link[RxODE]{RxODE}}
##'
##' @examples
##'
##' pk1cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot # Oral
##' pk1cmt %>% solve(et(amt=100,cmt=central,timeUnits="hr")) %>% plot # IV
##'
"pk1cmt"