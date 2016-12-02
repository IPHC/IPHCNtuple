#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(Tau)

Tau::Tau()
{
}

Tau::~Tau()
{
}

void Tau::read()
{
    _ID = idx;
  
    // general informations
    if( CHECK(ntP->tau_E)                             ) _E      = ntP->tau_E->at(idx);
    if( CHECK(ntP->tau_pt)                            ) _pt     = ntP->tau_pt->at(idx);
    if( CHECK(ntP->tau_pt)                            ) _ptCor  = ntP->tau_pt->at(idx);
    if( CHECK(ntP->tau_pt)                            ) _ptUnc  = ntP->tau_pt->at(idx);
    if( CHECK(ntP->tau_eta)                           ) _eta    = ntP->tau_eta->at(idx);
    if( CHECK(ntP->tau_phi)                           ) _phi    = ntP->tau_phi->at(idx);
    if( CHECK(ntP->tau_m)                             ) _m      = ntP->tau_m->at(idx);
    if( CHECK(ntP->tau_leadingTrackDxy)               ) _dxy  = ntP->tau_leadingTrackDxy->at(idx);
    //if( CHECK(ntP->tau_leadingTrackDxy)               ) _dxy    = ntP->tau_leadingTrackDxy->at(idx);
    if( CHECK(ntP->tau_leadingTrackDz)                ) _dz   = ntP->tau_leadingTrackDz->at(idx);
    //if( CHECK(ntP->tau_leadingTrackDz)                ) _dz     = ntP->tau_leadingTrackDz->at(idx);
    if( CHECK(ntP->tau_charge)                        ) _charge = ntP->tau_charge->at(idx);
    if( CHECK(ntP->tau_id)                            ) _id     = ntP->tau_id->at(idx);

    // selection variables
    
    
    // more variables
    _decayMode = ntP->tau_decayMode->at(idx); 
    _hasLeadChargedHadrCand = ntP->tau_hasLeadChargedHadrCand->at(idx); 
    _leadingTrackPt = ntP->tau_leadingTrackPt->at(idx); 
    _leadingTrackCharge = ntP->tau_leadingTrackCharge->at(idx); 
    
 
    _decayModeFindingOldDMs = ntP->tau_decayModeFindingOldDMs->at(idx);
    _decayModeFindingOldDMs = ntP->tau_decayModeFindingOldDMs->at(idx);
  
    _byLooseCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(idx);
    _byMediumCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(idx);
    _byTightCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(idx);
 
    //_byLooseIsolationMVA3newDMwLT = ntP->tau_byLooseIsolationMVA3newDMwLT->at(idx); //byVLooseIsolationMVArun2v1DBnewDMwLT
    //_byMediumIsolationMVA3newDMwLT = ntP->tau_byMediumIsolationMVA3newDMwLT->at(idx);
    //_byTightIsolationMVA3newDMwLT = ntP->tau_byTightIsolationMVA3newDMwLT->at(idx);

    //_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03  = ntP->tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03->at(idx);
    //_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = ntP->tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03->at(idx);
    //_byTightCombinedIsolationDeltaBetaCorr3HitsdR03  = ntP->tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03->at(idx);


 
    _byLooseIsolationMVArun2v1DBdR03oldDMwLT         = ntP->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
    _byMediumIsolationMVArun2v1DBdR03oldDMwLT        = ntP->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
    _byTightIsolationMVArun2v1DBdR03oldDMwLT         = ntP->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
    _byVTightIsolationMVArun2v1DBdR03oldDMwLT        = ntP->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
 
 
    _byCombinedIsolationDeltaBetaCorrRaw3Hits        = ntP->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(idx);
    _chargedIsoPtSum = ntP->tau_chargedIsoPtSum->at(idx);
    _neutralIsoPtSum = ntP->tau_neutralIsoPtSum->at(idx);
    _puCorrPtSum = ntP->tau_puCorrPtSum->at(idx);
    
   
    _againstMuonLoose3 = ntP->tau_againstMuonLoose3->at(idx);
    _againstMuonTight3 = ntP->tau_againstMuonTight3->at(idx);
 
  
    //AC8X
   // _againstElectronVLoose = ntP->tau_againstElectronVLooseMVA6->at(idx);
   // _againstElectronLoose  = ntP->tau_againstElectronLooseMVA6->at(idx);
   // _againstElectronMedium = ntP->tau_againstElectronMediumMVA6->at(idx);
   // _againstElectronTight = ntP->tau_againstElectronMediumMVA6->at(idx);
    
  
    _pfEssential_jet_pt 	 = ntP->tau_pfEssential_jet_pt->at(idx);
    _pfEssential_jet_eta	 = ntP->tau_pfEssential_jet_eta->at(idx);
    _pfEssential_jet_phi	 = ntP->tau_pfEssential_jet_phi->at(idx);
    _pfEssential_jet_m  	 = ntP->tau_pfEssential_jet_m->at(idx);
    _pfEssential_jetCorr_pt	 = ntP->tau_pfEssential_jetCorr_pt->at(idx);
    _pfEssential_jetCorr_eta	 = ntP->tau_pfEssential_jetCorr_eta->at(idx);
    _pfEssential_jetCorr_phi	 = ntP->tau_pfEssential_jetCorr_phi->at(idx);
    _pfEssential_jetCorr_m	 = ntP->tau_pfEssential_jetCorr_m->at(idx);
    _pfEssential_hasSV  	 = ntP->tau_pfEssential_hasSV->at(idx);
    _pfEssential_sv_x		 = ntP->tau_pfEssential_sv_x->at(idx);
    _pfEssential_sv_y		 = ntP->tau_pfEssential_sv_y->at(idx);
    _pfEssential_sv_z		 = ntP->tau_pfEssential_sv_z->at(idx);
    _pfEssential_flightLengthSig = ntP->tau_pfEssential_flightLengthSig->at(idx);
    _pfEssential_dxy		 = ntP->tau_pfEssential_dxy->at(idx);
    _pfEssential_dxy_error	 = ntP->tau_pfEssential_dxy_error->at(idx);
    _pfEssential_dxy_Sig	 = ntP->tau_pfEssential_dxy_Sig->at(idx);
    
    
}

void Tau::init()
{  
    _fakeType = -1;

    // general informations
    _E        = -100;
    _pt       = -100;
    _ptUnc    = -100;
    _eta      = -100;
    _phi      = -100;
    _m        = -100;
    _charge   =    0;
    _id       =    0;

    // Id
    _isLoose        = 0;
    _isFakeableTTH  = 0; //dummy !!!
    _isTightTTH     = 0; //dummy !!!
    _lepMVA_TTH     = 0.; //dummy !!!

    // variables for Id
    _dxy      = -100;
    _dz       = -100;

    _passTightCharge    = true;
    _cutEventSel        = true;
    _noLostHits         = true;

    // more variables
    _decayMode              = -1;
    _hasLeadChargedHadrCand = false;
    _leadingTrackPt         = -1.;
    _leadingTrackCharge     = -1;

    _decayModeFindingOldDMs = -1.;
    _decayModeFindingOldDMs = -1.;

    _byLooseCombinedIsolationDeltaBetaCorr3Hits  = -1.;
    _byMediumCombinedIsolationDeltaBetaCorr3Hits = -1.;
    _byTightCombinedIsolationDeltaBetaCorr3Hits  = -1.;

    _byLooseIsolationMVA3newDMwLT  = -1.;  
    _byMediumIsolationMVA3newDMwLT = -1.; 
    _byTightIsolationMVA3newDMwLT  = -1.; 

    _byLooseCombinedIsolationDeltaBetaCorr3HitsdR03  = -1.;
    _byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -1.;
    _byTightCombinedIsolationDeltaBetaCorr3HitsdR03  = -1.;

    _byLooseIsolationMVArun2v1DBdR03oldDMwLT  = -1;
    _byMediumIsolationMVArun2v1DBdR03oldDMwLT = -1;
    _byTightIsolationMVArun2v1DBdR03oldDMwLT  = -1;
    _byVTightIsolationMVArun2v1DBdR03oldDMwLT = -1;

    _byCombinedIsolationDeltaBetaCorrRaw3Hits = -1.;
    _chargedIsoPtSum = -1.;
    _neutralIsoPtSum = -1.;
    _puCorrPtSum     = -1.;

    _againstMuonLoose3 = -1.;
    _againstMuonTight3 = -1.;

    _againstElectronVLoose = -1.;
    _againstElectronLoose  = -1.;
    _againstElectronMedium = -1.;
    _againstElectronTight = -1.;

    _pfEssential_jet_pt          = -1.;
    _pfEssential_jet_eta         = -1.;
    _pfEssential_jet_phi         = -1.;
    _pfEssential_jet_m           = -1.;
    _pfEssential_jetCorr_pt      = -1.;
    _pfEssential_jetCorr_eta     = -1.;
    _pfEssential_jetCorr_phi     = -1.;
    _pfEssential_jetCorr_m       = -1.;
    _pfEssential_hasSV           = -1.;
    _pfEssential_sv_x            = -1.;
    _pfEssential_sv_y            = -1.;
    _pfEssential_sv_z            = -1.;
    _pfEssential_flightLengthSig = -1.;
    _pfEssential_dxy             = -1.;
    _pfEssential_dxy_error       = -1.;
    _pfEssential_dxy_Sig         = -1.;
}

bool Tau::sel()
{   

    // selection
    bool pass_pt  = (_pt        > 20.  );
    bool pass_eta = (fabs(_eta) < 2.3  );
    bool pass_dxy = (fabs(_dxy) < 1000 );
    bool pass_dz  = (fabs(_dz)  < 0.2  );

    bool pass_decayModeFindingOldDMs                        = ( ntP->tau_decayModeFindingOldDMs->at(idx)                          > 0.5 );
    bool pass_byLooseIsolationMVArun2v1DBdR03oldDMwLT       = ( ntP->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(idx)         > 0.5 );

    bool pass_muOverlap = 1;
    int nMuon = nt->NtMuon->size();
    for(int im=0;im<nMuon;im++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtMuon->at(im).eta(),nt->NtMuon->at(im).phi());
        if( dr < 0.4 ) pass_muOverlap = 0; //&& nt->NtMuon->at(im).isLoose() ) pass_muOverlap = 0;
    }  

    bool pass_elOverlap = 1;
    int nEl = nt->NtElectron->size();
    for(int iEl=0;iEl<nEl;iEl++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtElectron->at(iEl).eta(),nt->NtElectron->at(iEl).phi());
        if( dr < 0.4 ) pass_elOverlap = 0; //&& nt->NtElectron->at(iEl).isLoose() ) pass_elOverlap = 0;
    }  

    bool isSelectionTau     = ( pass_pt                                         &&
                                pass_eta                                        &&
                                pass_dxy                                        &&
                                pass_dz                                         &&
                                pass_decayModeFindingOldDMs                     &&
                                pass_byLooseIsolationMVArun2v1DBdR03oldDMwLT    &&
                                pass_muOverlap                                  &&
                                pass_elOverlap                                   );

    bool pass_byMediumIsolationMVArun2v1DBdR03oldDMwLT = (ntP->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(idx)         > 0.5 );

    _isTightTTH         = ( isSelectionTau                                  &&
                            pass_byMediumIsolationMVArun2v1DBdR03oldDMwLT   );

    cout<<std::setiosflags(ios::fixed)<<setprecision(5);
    // synchronization printout
    if( false )     std::cout   << "Taus: "
                                << nt->NtEvent->at(0).id()                                      << " "
                                << _pt                                                          << " "
                                << _eta                                                         << " "
                                << _phi                                                         << " "
                                << _E                                                           << " "
                                << _dxy                                                         << " "
                                << _dz                                                          << " "
                                << ntP->tau_decayModeFindingOldDMs->at(idx)                     << " "
                                << ntP->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(idx)    << " "
                                << ntP->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(idx)   << std::endl;


    return isSelectionTau;
}
