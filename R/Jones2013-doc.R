##' Jones2013 Human PBPK RxODE model
##'
##' @format An \emph{RxODE} model with 53 parameters, 16 ODE states, and 1 calc vars.
##'
##'\emph{Parameters (Jones2013$params)}
##'
##' \describe{
##'   \item{BW}{ (default=70)}
##'   \item{FVad}{ (default=0.213)}
##'   \item{FVbo}{ (default=0.085629)}
##'   \item{FVbr}{ (default=0.02)}
##'   \item{FVgu}{ (default=0.0171)}
##'   \item{FVhe}{ (default=0.0047)}
##'   \item{FVki}{ (default=0.0044)}
##'   \item{FVli}{ (default=0.021)}
##'   \item{FVlu}{ (default=0.0076)}
##'   \item{FVmu}{ (default=0.4)}
##'   \item{FVsk}{ (default=0.0371)}
##'   \item{FVsp}{ (default=0.0026)}
##'   \item{FVte}{ (default=0.01)}
##'   \item{FVve}{ (default=0.0514)}
##'   \item{FVar}{ (default=0.0257)}
##'   \item{FVpl}{ (default=0.0424)}
##'   \item{FVrb}{ (default=0.0347)}
##'   \item{FVre}{ (default=0.099771)}
##'   \item{FQad}{ (default=0.05)}
##'   \item{FQbo}{ (default=0.05)}
##'   \item{FQbr}{ (default=0.12)}
##'   \item{FQgu}{ (default=0.146462)}
##'   \item{FQhe}{ (default=0.04)}
##'   \item{FQki}{ (default=0.19)}
##'   \item{FQh}{ (default=0.215385)}
##'   \item{FQlu}{ (default=1)}
##'   \item{FQmu}{ (default=0.17)}
##'   \item{FQsk}{ (default=0.05)}
##'   \item{FQsp}{ (default=0.017231)}
##'   \item{FQte}{ (default=0.01076)}
##'   \item{FQre}{ (default=0.103855)}
##'   \item{Kpad}{ (default=0.191)}
##'   \item{Kpbo}{ (default=0.374)}
##'   \item{Kpbr}{ (default=0.606)}
##'   \item{Kpgu}{ (default=0.578)}
##'   \item{Kphe}{ (default=0.583)}
##'   \item{Kpki}{ (default=0.597)}
##'   \item{Kpli}{ (default=0.57)}
##'   \item{Kplu}{ (default=0.62)}
##'   \item{Kpmu}{ (default=0.622)}
##'   \item{Kpsk}{ (default=0.6)}
##'   \item{Kpsp}{ (default=0.591)}
##'   \item{Kpte}{ (default=0.6)}
##'   \item{Kpre}{ (default=0.6)}
##'   \item{fup}{ (default=0.681)}
##'   \item{BP}{ (default=0.98)}
##'   \item{fumic}{ (default=1)}
##'   \item{HLM_CLint}{ (default=8)}
##'   \item{CLrenal}{ (default=0)}
##'   \item{Ka}{ (default=2.18)}
##'   \item{F}{ (default=1)}
##'   \item{CO}{ (default=108.33)}
##'   \item{MPPGL}{ (default=45)}
##'}
##'
##' \emph{State Jones2013$state}
##'
##' \describe{
##'   \item{Aad}{ (=1)}
##'   \item{Abo}{ (=2)}
##'   \item{Abr}{ (=3)}
##'   \item{Agu}{ (=4)}
##'   \item{Ahe}{ (=5)}
##'   \item{Aki}{ (=6)}
##'   \item{Ali}{ (=7)}
##'   \item{Alu}{ (=8)}
##'   \item{Amu}{ (=9)}
##'   \item{Ask}{ (=10)}
##'   \item{Asp}{ (=11)}
##'   \item{Ate}{ (=12)}
##'   \item{Ave}{ (=13)}
##'   \item{Aar}{ (=14)}
##'   \item{Are}{ (=15)}
##'   \item{D}{ (=16)}
##' }
##'
##' \emph{Calculated Variables Jones2013$lhs}
##'
##' \describe{
##'   \item{Cp}{Concentration in plasma}
##' }
##'
##' \emph{Model Code}
##'
##' RxODE({
##'     BW = 70
##'     FVad = 0.213
##'     FVbo = 0.085629
##'     FVbr = 0.02
##'     FVgu = 0.0171
##'     FVhe = 0.0047
##'     FVki = 0.0044
##'     FVli = 0.021
##'     FVlu = 0.0076
##'     FVmu = 0.4
##'     FVsk = 0.0371
##'     FVsp = 0.0026
##'     FVte = 0.01
##'     FVve = 0.0514
##'     FVar = 0.0257
##'     FVpl = 0.0424
##'     FVrb = 0.0347
##'     FVre = 0.099771
##'     FQad = 0.05
##'     FQbo = 0.05
##'     FQbr = 0.12
##'     FQgu = 0.146462
##'     FQhe = 0.04
##'     FQki = 0.19
##'     FQh = 0.215385
##'     FQlu = 1
##'     FQmu = 0.17
##'     FQsk = 0.05
##'     FQsp = 0.017231
##'     FQte = 0.01076
##'     FQre = 0.103855
##'     Kpad = 0.191
##'     Kpbo = 0.374
##'     Kpbr = 0.606
##'     Kpgu = 0.578
##'     Kphe = 0.583
##'     Kpki = 0.597
##'     Kpli = 0.57
##'     Kplu = 0.62
##'     Kpmu = 0.622
##'     Kpsk = 0.6
##'     Kpsp = 0.591
##'     Kpte = 0.6
##'     Kpre = 0.6
##'     fup = 0.681
##'     BP = 0.98
##'     fumic = 1
##'     HLM_CLint = 8
##'     CLrenal = 0
##'     Ka = 2.18
##'     F = 1
##'     CO = 108.33
##'     Vad ~ BW * FVad
##'     Vbo ~ BW * FVbo
##'     Vbr ~ BW * FVbr
##'     Vgu ~ BW * FVgu
##'     Vhe ~ BW * FVhe
##'     Vki ~ BW * FVki
##'     Vli ~ BW * FVli
##'     Vlu ~ BW * FVlu
##'     Vmu ~ BW * FVmu
##'     Vsk ~ BW * FVsk
##'     Vsp ~ BW * FVsp
##'     Vte ~ BW * FVte
##'     Vve ~ BW * FVve
##'     Var ~ BW * FVar
##'     Vpl ~ BW * FVpl
##'     Vrb ~ BW * FVrb
##'     Vre ~ BW * FVre
##'     Vplas_ven ~ Vpl * Vve/(Vve + Var)
##'     Vplas_art ~ Vpl * Var/(Vve + Var)
##'     QC ~ CO/1000 * 60 * 60
##'     Qad ~ QC * FQad
##'     Qbo ~ QC * FQbo
##'     Qbr ~ QC * FQbr
##'     Qgu ~ QC * FQgu
##'     Qhe ~ QC * FQhe
##'     Qki ~ QC * FQki
##'     Qh ~ QC * FQh
##'     Qha ~ Qh - Qgu - Qsp
##'     Qlu ~ QC * FQlu
##'     Qmu ~ QC * FQmu
##'     Qsk ~ QC * FQsk
##'     Qsp ~ QC * FQsp
##'     Qte ~ QC * FQte
##'     Qre ~ QC * FQre
##'     Cadipose ~ Aad/Vad
##'     Cbone ~ Abo/Vbo
##'     Cbrain ~ Abr/Vbr
##'     Cgut ~ Agu/Vgu
##'     Cheart ~ Ahe/Vhe
##'     Ckidney ~ Aki/Vki
##'     Cliver ~ Ali/Vli
##'     Clung ~ Alu/Vlu
##'     Cmuscle ~ Amu/Vmu
##'     Cskin ~ Ask/Vsk
##'     Cspleen ~ Asp/Vsp
##'     Ctestes ~ Ate/Vte
##'     Cvenous ~ Ave/Vve
##'     Carterial ~ Aar/Var
##'     Crest ~ Are/Vre
##'     Absorption ~ Ka * D * F
##'     Cliverfree ~ Cliver * fup
##'     Ckidneyfree ~ Ckidney * fup
##'     MPPGL = 45
##'     CLmet ~ (HLM_CLint/fumic) * MPPGL * Vli * 60/1000
##'     Venous ~ Qad * (Cadipose/Kpad * BP) + Qbo * (Cbone/Kpbo *
##'         BP) + Qbr * (Cbrain/Kpbr * BP) + Qhe * (Cheart/Kphe *
##'         BP) + Qki * (Ckidney/Kpki * BP) + Qh * (Cliver/Kpli *
##'         BP) + Qmu * (Cmuscle/Kpmu * BP) + Qsk * (Cskin/Kpsk *
##'         BP) + Qte * (Ctestes/Kpte * BP) + Qre * (Crest/Kpre *
##'         BP)
##'     d/dt(Aad) ~ Qad * (Carterial - Cadipose/Kpad * BP)
##'     d/dt(Abo) ~ Qbo * (Carterial - Cbone/Kpbo * BP)
##'     d/dt(Abr) ~ Qbr * (Carterial - Cbrain/Kpbr * BP)
##'     d/dt(Agu) ~ Absorption + Qgu * (Carterial - Cgut/Kpgu * BP)
##'     d/dt(Ahe) ~ Qhe * (Carterial - Cheart/Kphe * BP)
##'     d/dt(Aki) ~ Qki * (Carterial - Ckidney/Kpki * BP) - CLrenal *
##'         Ckidneyfree
##'     d/dt(Ali) ~ Qha * Carterial + Qgu * (Cgut/Kpgu * BP) + Qsp *
##'         (Cspleen/Kpsp * BP) - Qh * (Cliver/Kpli * BP) - Cliverfree *
##'         CLmet
##'     d/dt(Alu) ~ Qlu * Cvenous - Qlu * (Clung/Kplu * BP)
##'     d/dt(Amu) ~ Qmu * (Carterial - Cmuscle/Kpmu * BP)
##'     d/dt(Ask) ~ Qsk * (Carterial - Cskin/Kpsk * BP)
##'     d/dt(Asp) ~ Qsp * (Carterial - Cspleen/Kpsp * BP)
##'     d/dt(Ate) ~ Qte * (Carterial - Ctestes/Kpte * BP)
##'     d/dt(Ave) ~ Venous - Qlu * Cvenous
##'     d/dt(Aar) ~ Qlu * (Clung/Kplu * BP) - Qlu * Carterial
##'     d/dt(Are) ~ Qre * (Carterial - Crest/Kpre * BP)
##'     d/dt(D) ~ -Absorption
##'     Cvenous ~ Ave/Vve
##'     Cp = Cvenous/BP
##' })
##' @references Jones H, Rowland-Yeo K. Basic concepts in
##'     physiologically based pharmacokinetic modeling in drug
##'     discovery and development. CPT Pharmacometrics Syst
##'     Pharmacol. 2013 Aug 14;2:e63. doi: 10.1038/psp.2013.41. PubMed
##'     PMID: 23945604; PubMed Central PMCID: PMC3828005.
##'
##' https://github.com/metrumresearchgroup/mrgsolve/blob/9ff7240da48ec133635a6f0db2143f13623f0aff/inst/models/pbpk.cpp
##'
##' @seealso \code{\link[RxODE]{eventTable}}, \code{\link[RxODE]{et}}, \code{\link[RxODE]{rxSolve}}, \code{\link[RxODE]{RxODE}}
##'
##' @examples
##' Jones2013 %>% solve(et(amt=100,cmt=D,ii=4,until=24) %>% et(seq(0,24,length.out=100))) %>% plot
"Jones2013"
