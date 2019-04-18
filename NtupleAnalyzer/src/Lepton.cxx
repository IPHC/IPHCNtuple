#include "../include/Lepton.h"

Lepton::Lepton()
{
    pt                 = 0.;
    conept                 = 0.;
    //ptCor              = 0.;
    //ptUnc              = 0.;
    eta                = 0.;
    phi                = 0.;
    E                  = 0.;

    p4.SetPtEtaPhiM(0,0,0,0);

    id                 =  0;

    idx                = -1;
    
    isFakeableTTH      = false;
    isLooseTTH		= false;
    isTightTTH         = false;

    lepMVA         = 0.;
    
    passTightCharge    = false;
    cutEventSel        = false;

    charge             =  0;
    
    hasMCMatch		= 0;
    hasChargeMCMatch	= 0;
    hasPhotonMCMatch	= 0;
}

Lepton::~Lepton()
{
}
