library(RxODE)
message("Jones2013")
Jones2013 <- RxODE({
    ## 1: Jones H, Rowland-Yeo K. Basic concepts in physiologically based
    ## pharmacokinetic modeling in drug discovery and development. CPT Pharmacometrics
    ## Syst Pharmacol. 2013 Aug 14;2:e63. doi: 10.1038/psp.2013.41. PubMed PMID:
    ## 23945604; PubMed Central PMCID: PMC3828005.

    BW = 70	##; BW (kg)

    ## {Fractional tissue volumes}

    FVad = 0.213     ##; adipose
    FVbo = 0.085629  ##; bone
    FVbr = 0.02      ##; brain
    FVgu = 0.0171    ##; gut
    FVhe = 0.0047    ##; heart
    FVki = 0.0044    ##; kidney
    FVli = 0.021     ##; liver
    FVlu = 0.0076    ##; lung
    FVmu = 0.4       ##; muscle
    FVsk = 0.0371    ##; skin
    FVsp = 0.0026    ##; spleen
    FVte = 0.01      ##; testes
    FVve = 0.0514    ##; venous
    FVar = 0.0257    ##; arterial
    FVpl = 0.0424    ##; plasma
    FVrb = 0.0347    ##; erythrocytes
    FVre = 0.099771  ##; rest of body

    ## {Fractional tissue blood flows}

    FQad = 0.05      ##; adipose
    FQbo = 0.05      ##; bone
    FQbr = 0.12      ##; brain
    FQgu = 0.146462  ##; gut
    FQhe = 0.04      ##; heart
    FQki = 0.19      ##; kidney
    FQh  = 0.215385  ##; hepatic (venous side)
    FQlu = 1         ##; lung
    FQmu = 0.17      ##; muscle
    FQsk = 0.05      ##; skin
    FQsp = 0.017231  ##; spleen
    FQte = 0.01076   ##; testes
    FQre = 0.103855  ##; rest of body

    ## {COMPOUND SPECIFIC PARAMETERS}

    ## {Tissue to plasma partition coefficients}

    Kpad = 0.191  ##; adipose
    Kpbo = 0.374  ##; bone
    Kpbr = 0.606  ##; brain
    Kpgu = 0.578  ##; gut
    Kphe = 0.583  ##; heart
    Kpki = 0.597  ##; kidney
    Kpli = 0.570  ##; liver
    Kplu = 0.620  ##; lung
    Kpmu = 0.622  ##; muscle
    Kpsk = 0.600  ##; skin
    Kpsp = 0.591  ##; spleen
    Kpte = 0.600  ##; testes
    Kpre = 0.600  ##; rest of body

    ## {In vitro binding data}

    fup   = 0.681  ##; fraction unbound in plasma
    BP    = 0.98  ##; blood to plasma ratio
    fumic = 1  ##; fraction unbound in microsomes

    ## {Clearances}

    HLM_CLint =  8   ##; HLM CLint apparent (ul/min/mg)
    CLrenal   =  0   ##; CLint renal (L/hr)

    ## {Absorption}

    Ka = 2.18     ##; Ka (hr-1)
    F  = 1.00     ##; fraction absorbed
    CO = 108.33   ##; cardiac output (ml/s)

    ## {Total tissue volumes - L}
    Vad ~ BW*FVad;  ## adipose
    Vbo ~ BW*FVbo;  ## bone
    Vbr ~ BW*FVbr;  ## brain
    Vgu ~ BW*FVgu;  ## gut
    Vhe ~ BW*FVhe;  ## heart
    Vki ~ BW*FVki;  ## kidney
    Vli ~ BW*FVli;  ## liver
    Vlu ~ BW*FVlu;  ## lung
    Vmu ~ BW*FVmu;  ## muscle
    Vsk ~ BW*FVsk;  ## skin
    Vsp ~ BW*FVsp;  ## spleen
    Vte ~ BW*FVte;  ## testes
    Vve ~ BW*FVve;  ## venous blood
    Var ~ BW*FVar;  ## arterial blood
    Vpl ~ BW*FVpl;  ## plasma
    Vrb ~ BW*FVrb;  ## erythrocytes
    Vre ~ BW*FVre;  ## rest of body

    Vplas_ven ~ Vpl*Vve/(Vve + Var) 	;  ## venous plasma
    Vplas_art ~ Vpl*Var/(Vve + Var) 	;  ## arterial plasma

    ## {Total tissue blood flows - L/hr}

    QC  ~ CO/1000*60*60 ;  ## cardiac output (L/hr)
    Qad ~ QC*FQad       ;  ## adipose
    Qbo ~ QC*FQbo       ;  ## bone
    Qbr ~ QC*FQbr       ;  ## brain
    Qgu ~ QC*FQgu       ;  ## gut
    Qhe ~ QC*FQhe       ;  ## heart
    Qki ~ QC*FQki       ;  ## kidney
    Qh  ~ QC*FQh        ;  ## hepatic (venous side)
    Qha ~ Qh - Qgu - Qsp;  ## hepatic artery
    Qlu ~ QC*FQlu       ;  ## lung
    Qmu ~ QC*FQmu       ;  ## muscle
    Qsk ~ QC*FQsk       ;  ## skin
    Qsp ~ QC*FQsp       ;  ## spleen
    Qte ~ QC*FQte       ;  ## testes
    Qre ~ QC*FQre       ;  ## rest of body

    ## D     ~  dose
    ## Aad   ~ adipose
    ## Abo   ~ bone
    ## Abr   ~ brain
    ## Agu   ~ gut
    ## Ahe   ~ heart
    ## Aki   ~ kidney
    ## Ali   ~ liver
    ## Alu   ~ lung
    ## Amu   ~ muscle
    ## Ask   ~ skin
    ## Asp   ~ spleen
    ## Ate   ~ testes
    ## Ave   ~ venous blood
    ## Aar   ~ arterial blood
    ## Are   ~ rest of body
    Cadipose  ~ Aad/Vad;  ## adipose
    Cbone     ~ Abo/Vbo;  ## bone
    Cbrain    ~ Abr/Vbr;  ## brain
    Cgut      ~ Agu/Vgu;  ## gut
    Cheart    ~ Ahe/Vhe;  ## heart
    Ckidney   ~ Aki/Vki;  ## kidney
    Cliver    ~ Ali/Vli;  ## liver
    Clung     ~ Alu/Vlu;  ## lung
    Cmuscle   ~ Amu/Vmu;  ## muscle
    Cskin     ~ Ask/Vsk;  ## skin
    Cspleen   ~ Asp/Vsp;  ## spleen
    Ctestes   ~ Ate/Vte;  ## testes
    Cvenous   ~ Ave/Vve;  ## venous blood
    Carterial ~ Aar/Var;  ## arterial blood
    Crest     ~ Are/Vre;  ## rest of body
    Absorption ~ Ka*D*F;

    ## {Calculation of free concentrations - mg/L}

    Cliverfree ~ Cliver*fup;  ## liver
    Ckidneyfree ~ Ckidney*fup; ## kidney

    ## {Clearance calculations}

    MPPGL = 45; ## mg microsomal protein per g liver
    CLmet ~ (HLM_CLint/fumic)*MPPGL*Vli*60/1000; ## CLint scaled (L/hr)

    Venous ~ Qad*(Cadipose/Kpad*BP) + Qbo*(Cbone/Kpbo*BP)   +
        Qbr*(Cbrain/Kpbr*BP)   + Qhe*(Cheart/Kphe*BP)  + Qki*(Ckidney/Kpki*BP) +
        Qh*(Cliver/Kpli*BP)    + Qmu*(Cmuscle/Kpmu*BP) + Qsk*(Cskin/Kpsk*BP) +
        Qte*(Ctestes/Kpte*BP)  + Qre*(Crest/Kpre*BP);



    d/dt(Aad) ~ Qad*(Carterial - Cadipose/Kpad*BP);    ## adipose
    d/dt(Abo) ~ Qbo*(Carterial - Cbone/Kpbo*BP);       ## bone
    d/dt(Abr) ~ Qbr*(Carterial - Cbrain/Kpbr*BP);      ## brain
    d/dt(Agu) ~ Absorption +
        Qgu*(Carterial - Cgut/Kpgu*BP);        ## gut
    d/dt(Ahe) ~ Qhe*(Carterial - Cheart/Kphe*BP);      ## heart
    d/dt(Aki) ~ Qki*(Carterial - Ckidney/Kpki*BP) -
        CLrenal*Ckidneyfree;                   ## kidney
    d/dt(Ali) ~ Qha*Carterial +
        Qgu*(Cgut/Kpgu*BP) +
        Qsp*(Cspleen/Kpsp*BP) -
        Qh*(Cliver/Kpli*BP) -
        Cliverfree*CLmet;                      ## liver
    d/dt(Alu) ~ Qlu*Cvenous - Qlu*(Clung/Kplu*BP);     ## lung
    d/dt(Amu) ~ Qmu*(Carterial - Cmuscle/Kpmu*BP);     ## muscle
    d/dt(Ask) ~ Qsk*(Carterial - Cskin/Kpsk*BP);       ## skin
    d/dt(Asp) ~ Qsp*(Carterial - Cspleen/Kpsp*BP);     ## spleen
    d/dt(Ate) ~ Qte*(Carterial - Ctestes/Kpte*BP);     ## testes
    d/dt(Ave) ~ Venous - Qlu*Cvenous;                  ## venous blood
    d/dt(Aar) ~ Qlu*(Clung/Kplu*BP) - Qlu*Carterial;   ## arterial blood
    d/dt(Are) ~ Qre*(Carterial - Crest/Kpre*BP);       ## rest of body
    d/dt(D  ) ~ - Absorption;                          ## oral dosing

    Cvenous ~ Ave/Vve
    Cp = Cvenous / BP
})

rxUse(Jones2013);
message("oral1cmt")
oral1cmt <- RxODE({
    popCl <- 1
    popV <- 20
    popKa <- 1
    bsvCl <- 0
    bsvV  <- 0
    bsvKa <- 0
    cl ~ popCl * exp(bsvCl)
    v ~ popV * exp(bsvV)
    ka ~ popKa * exp(bsvKa)
    popLagDepot <- 0
    popLagCentral <- 0
    popRateCentral <- 0
    popDurCentral <- 0
    bsvLagDepot <- 0
    bsvLagCentral <- 0
    bsvRateCentral <- 0
    bsvDurCentral <- 0
    lag(depot) <- popLagDepot * exp(bsvLagDepot)
    lag(central) <- popLagCentral * exp(bsvLagCentral)
    rate(central) <- popRateCentral *  exp(bsvRateCentral)
    dur(central) <- popDurCentral * exp(bsvDurCentral)
    cp <- linCmt()
});
rxUse(oral1cmt);

message("oral2cmt")
oral2cmt <- RxODE({
    popCl <- 1
    popV <- 20
    popKa <- 1
    popVp <- 10
    popQ <- 2
    bsvCl <-0
    bsvV <- 0
    bsvKa <-0
    bsvVp <- 0
    bsvQ <-0
    cl ~ popCl * exp(bsvCl)
    v ~ popV * exp(bsvV)
    ka ~ popKa * exp(bsvKa)
    q ~ popQ * exp(bsvQ)
    vp ~ popVp * exp(bsvVp)
    popLagDepot <- 0
    popLagCentral <- 0
    popRateCentral <- 0
    popDurCentral <- 0
    bsvLagDepot <- 0
    bsvLagCentral <- 0
    bsvRateCentral <- 0
    bsvDurCentral <- 0
    lag(depot) <- popLagDepot * exp(bsvLagDepot)
    lag(central) <- popLagCentral * exp(bsvLagCentral)
    rate(central) <- popRateCentral * exp(bsvRateCentral)
    dur(central) <- popDurCentral * exp(bsvDurCentral)
    cp <- linCmt()
});

rxUse(oral2cmt);

message("oral3cmt")
oral3cmt <- RxODE({
    popCl <- 1
    popV <- 20
    popKa <- 1
    popVp <- 10
    popQ <- 2
    popQ2 <- 2
    popVp2 <- 100
    bsvCl <- 0
    bsvV <- 0
    bsvKa <- 0
    bsvVp <- 0
    bsvQ <- 0
    bsvQ2 <- 0
    bsvVp2 <- 0
    cl ~ popCl * exp(bsvCl)
    v ~ popV * exp(bsvV)
    ka ~ popKa * exp(bsvKa)
    q ~ popQ * exp(bsvQ)
    vp ~ popVp * exp(bsvVp)
    q2 ~ popQ2 * exp(bsvQ2)
    vp2 ~ popVp2 * exp(bsvVp2)
    popLagDepot <- 0
    popLagCentral <- 0
    popRateCentral <- 0
    popDurCentral <- 0
    bsvLagDepot <- 0
    bsvLagCentral <- 0
    bsvRateCentral <- 0
    bsvDurCentral <- 0
    lag(depot)    <- popLagDepot * exp(bsvLagDepot)
    lag(central)  <- popLagCentral * exp(bsvLagCentral)
    rate(central) <- popRateCentral * exp(bsvRateCentral)
    dur(central)  <- popDurCentral * exp(bsvDurCentral)
    cp <- linCmt()
});
rxUse(oral3cmt);

