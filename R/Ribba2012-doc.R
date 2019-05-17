##' Ribba2012 RxODE model
##'
##' A tumor growth inhibition model for low-grade glioma treated with
##' chemotherapy or radiotherapy
##'
##' @format An \emph{RxODE} model with 9 parameters, 4 ODE states, and 1 calc vars.
##'
##'\emph{Parameters (Ribba2012$params)}
##'
##' Typical parameters:
##' \describe{
##'   \item{tk}{Maximum tumor size (default=100) mm}
##'   \item{tkde}{Decay rate of PCV concentration in plasma (default=0.3)}
##'   \item{tkpq}{Transition rate of proliferative tissue to quiescence tissue (default=0.025)}
##'   \item{tkqpp}{Rate constant for damaged quiescent tissue transferring x1to proliferative tissue (default=0.004)}
##'   \item{tlambdap}{Rate of tissue growth for proliferative tissue (default=0.12)}
##'   \item{tgamma}{Damage rate in all tissue (default=1)}
##'   \item{tdeltaqp}{Rate constant for elimination of the damaged quiescent tissue (default=0.01)}
##'   \item{tpt0}{Initial proliferative tissue amount (default=5)}
##'   \item{tq0}{Initial quiescent tissue level (default=40)}
##'}
##'
##' Eta parameters are all zero by default.
##'
##' \emph{State Ribba2012$state}
##'
##' \describe{
##'   \item{c}{Concentration compartment (=1)}
##'   \item{pt}{Proliferative Tissue (=2)}
##'   \item{q}{Nonproliferative or quiescent tissue (=3)}
##'   \item{qp}{DNA-Damaged quiescent tissue (=4)}
##' }
##'
##'
##' \emph{Calculated Variables Ribba2012$lhs}
##'
##' \describe{
##'   \item{pstar}{The total tumor tissue ie p + pt + q}
##' }
##'
##' A tumor growth inhibition model for low-grade glioma treated with chemotherapy or radiotherapy
##'
##' @seealso \code{\link[RxODE]{eventTable}}, \code{\link[RxODE]{et}}, \code{\link[RxODE]{rxSolve}}, \code{\link[RxODE]{RxODE}}
##'
##' @references Ribba B, Kaloshi G, et al.  A. tumor growth inhibition
##'     model for low-grade glioma treated with chemotherapy or
##'     radiotherapy. Clin Cancer Res. 2012 Sep
##'     15;18(18):5071-80. Epub 2012 Jul 3.
##'
##' https://www.ncbi.nlm.nih.gov/pubmed/22761472
##' http://simulx.webpopix.org/videos/simulx-video2/
##'
##' @examples
##' ## Showing the model code
##' summary(Ribba2012)
##'
##' ## Model without doses
##'
##' Ribba2012 %>%
##'   et(time.units="months") %>%
##'   et(0, 250, by=0.5) %>%
##'   rxSolve() %>%
##'   plot(pt, q, qp, pstar)
##'
##' ## Add dose of "1" from 50 to 57.5 months by 1.5
##'
##' Ribba2012 %>%
##'   et(time.units="months") %>%
##'   et(0, 150, by=0.5) %>%
##'   et(amt=1, time=50, until=58, ii=1.5) %>%
##'   rxSolve() %>%
##'   plot(pt, q, qp, pstar)
##'
##' ## Add CVs from paper for individual simulation
##' ## Uses exact formula:
##'
##' lognCv = function(x){log((x/100)^2+1)}
##'
##' library(lotri)
##' ## Now create omega matrix
##' omega <- lotri(eta.pt0 ~ lognCv(94),
##'                eta.q0 ~ lognCv(54),
##'                eta.lambdap ~ lognCv(72),
##'                eta.kqp ~ lognCv(76),
##'                eta.qpp ~ lognCv(97),
##'                eta.deltaqp ~ lognCv(115),
##'                eta.kde ~ lognCv(70))
##' set.seed(3708)
##'
##' ## simulate 3 subjects
##' Ribba2012 %>%
##'   et(time.units="months") %>%
##'   et(0, 150, by=0.5) %>%
##'   et(amt=1, time=50, until=58, ii=1.5) %>%
##'   rxParams(omega=omega, nSub=3) %>%
##'   rxSolve() %>%
##'   plot(pt, q, qp, pstar)
##'
##' ## simulate 3 subjects w/ 3 studies
##' Ribba2012 %>%
##'   et(time.units="months") %>%
##'   et(0, 150, by=0.5) %>%
##'   et(amt=1, time=50, until=58, ii=1.5) %>%
##'   rxParams(omega=omega, nSub=3, nStud=3, dfSub=24) %>%
##'   rxSolve() %>%
##'   plot(pt, q, qp, pstar)
##'
##'
"Ribba2012"


