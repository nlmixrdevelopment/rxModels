##' pk3cmt RxODE model
##'
##' @format An \emph{RxODE} model with 22 parameters, 2 ODE states, and 1 calc vars.
##'
##'\emph{Parameters (pk3cmt$params)}
##'
##' \describe{
##'   \item{popCl}{Population Clearance value (default=1)}
##'   \item{popV}{Population Volume of distribution (default=20)}
##'   \item{popKa}{Population Absorption coefficient (default=1)}
##'   \item{popVp}{Population Peripheral Volume of Distribution #1 (default=10)}
##'   \item{popVp2}{Population Peripheral Volume of Distribution #2 (default=100)}
##'   \item{popQ}{Population Inter-Compartmental Clearance of Compartment #1 (default=2)}
##'   \item{popQ2}{Population Inter-Compartmental Clearance of Compartment #2 (default=2)}
##'   \item{bsvCl}{Between Subject Variability for Clearance (default=0)}
##'   \item{bsvV}{Between Subject Variability for Volume (default=0)}
##'   \item{bsvKa}{Between Subject Variability for Ka (default=0)}
##'   \item{bsvVp}{Between Subject Variability for Peripheral Volume of Distribution #1 (default=0)}
##'   \item{bsvVp2}{Between Subject Variability for Peripheral Volume for Compartment #2 (default=0)}
##'   \item{bsvQ}{Between Subject Variability for Inter-compartmental Clearance #1 (default=0)}
##'   \item{bsvQ2}{Between Subject Variability for Inter-compartmental Clearance #2 (default=0)}
##'   \item{popLagDepot}{Population Value for Lagged doses on Central Compartment (default=0)}
##'   \item{popLagCentral}{Population Value for Lagged doses on Central compartment (default=0)}
##'   \item{popRateCentral}{Population Value for modeled rates (default=0)}
##'   \item{popDurCentral}{Population Value for modeled duration (default=0)}
##'   \item{bsvLagDepot}{Between Subject variability for depot (default=0)}
##'   \item{bsvLagCentral}{Between Subject Variability for Central Lag-time (default=0)}
##'   \item{bsvRateCentral}{Between Subject Variability for Central Rate (default=0)}
##'   \item{bsvDurCentral}{Between Subject Variability for Central Duration (default=0)}
##'}
##'
##' \emph{State pk3cmt$state}
##'
##' \describe{
##'   \item{depot}{Depot Compartment (=1)}
##'   \item{central}{Central compartment (=2)}
##' }
##'
##' \emph{Calculated Variables pk3cmt$lhs}
##'
##' \describe{
##'   \item{cp}{Concentration in the plasma}
##' }
##'
##' \emph{Model Code}
##'
##' RxODE({
##'     popCl = 1
##'     popV = 20
##'     popKa = 1
##'     popVp = 10
##'     popQ = 2
##'     popQ2 = 2
##'     popVp2 = 100
##'     bsvCl = 0
##'     bsvV = 0
##'     bsvKa = 0
##'     bsvVp = 0
##'     bsvQ = 0
##'     bsvQ2 = 0
##'     bsvVp2 = 0
##'     cl ~ popCl * exp(bsvCl)
##'     v ~ popV * exp(bsvV)
##'     ka ~ popKa * exp(bsvKa)
##'     q ~ popQ * exp(bsvQ)
##'     vp ~ popVp * exp(bsvVp)
##'     q2 ~ popQ2 * exp(bsvQ2)
##'     vp2 ~ popVp2 * exp(bsvVp2)
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
##'     rx_k12 ~ q/v
##'     rx_k21 ~ q/vp
##'     rx_beta ~ 0.5 * (rx_k12 + rx_k21 + rx_k - sqrt((rx_k12 +
##'         rx_k21 + rx_k) * (rx_k12 + rx_k21 + rx_k) - 4 * rx_k21 *
##'         rx_k))
##'     rx_alpha ~ rx_k21 * rx_k/rx_beta
##'     rx_A ~ rx_ka/(rx_ka - rx_alpha) * (rx_alpha - rx_k21)/(rx_alpha -
##'         rx_beta)/rx_v
##'     rx_B ~ rx_ka/(rx_ka - rx_beta) * (rx_beta - rx_k21)/(rx_beta -
##'         rx_alpha)/rx_v
##'     rx_A2 ~ (rx_alpha - rx_k21)/(rx_alpha - rx_beta)/rx_v
##'     rx_B2 ~ (rx_beta - rx_k21)/(rx_beta - rx_alpha)/rx_v
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
##' pk3cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot # Oral
##' pk3cmt %>% solve(et(amt=100,cmt=central,timeUnits="hr")) %>% plot # IV
"pk3cmt"
