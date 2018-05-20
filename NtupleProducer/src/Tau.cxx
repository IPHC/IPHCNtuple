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
   
   if( CHECK(ntP->tau_E)                             ) _E      = ntP->tau_E->at(idx);
   if( CHECK(ntP->tau_pt)                            ) _pt     = ntP->tau_pt->at(idx);
   if( CHECK(ntP->tau_pt)                            ) _ptCor  = ntP->tau_pt->at(idx);
   if( CHECK(ntP->tau_pt)                            ) _ptUnc  = ntP->tau_pt->at(idx);
   if( CHECK(ntP->tau_eta)                           ) _eta    = ntP->tau_eta->at(idx);
   if( CHECK(ntP->tau_phi)                           ) _phi    = ntP->tau_phi->at(idx);
   if( CHECK(ntP->tau_m)                             ) _m      = ntP->tau_m->at(idx);
   if( CHECK(ntP->tau_leadingTrackDxy)               ) _dxy  = ntP->tau_leadingTrackDxy->at(idx);
   if( CHECK(ntP->tau_leadingTrackDz)                ) _dz   = ntP->tau_leadingTrackDz->at(idx);
   if( CHECK(ntP->tau_charge)                        ) _charge = ntP->tau_charge->at(idx);
   if( CHECK(ntP->tau_id)                            ) _id     = ntP->tau_id->at(idx);

   _decayMode = ntP->tau_decayMode->at(idx); 
   _hasLeadChargedHadrCand = ntP->tau_hasLeadChargedHadrCand->at(idx); 
   _leadingTrackPt = ntP->tau_leadingTrackPt->at(idx); 
   _leadingTrackCharge = ntP->tau_leadingTrackCharge->at(idx); 

   _decayModeFinding = ntP->tau_decayModeFinding->at(idx);
//   _decayModeFindingOldDMs = ntP->tau_decayModeFindingOldDMs->at(idx);
   _decayModeFindingNewDMs = ntP->tau_decayModeFindingNewDMs->at(idx);
   _byLooseCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(idx);
   _byMediumCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(idx);
   _byTightCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(idx);
   _againstElectronVLooseMVA6 = ntP->tau_againstElectronVLooseMVA6->at(idx);
   _againstElectronLooseMVA6 = ntP->tau_againstElectronLooseMVA6->at(idx);
   _againstElectronMediumMVA6 = ntP->tau_againstElectronMediumMVA6->at(idx);
   _againstElectronTightMVA6 = ntP->tau_againstElectronTightMVA6->at(idx);

   _byVLooseIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
   _byLooseIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
   _byMediumIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
   _byTightIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
   _byVTightIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(idx);

   _byCombinedIsolationDeltaBetaCorrRaw3Hits = ntP->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(idx);
   _againstMuonLoose3 = ntP->tau_againstMuonLoose3->at(idx);
   _againstMuonTight3 = ntP->tau_againstMuonTight3->at(idx);

   _chargedIsoPtSum = ntP->tau_chargedIsoPtSum->at(idx);
   _neutralIsoPtSum = ntP->tau_neutralIsoPtSum->at(idx);
   _puCorrPtSum = ntP->tau_puCorrPtSum->at(idx);
   
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

   _hasMCMatch = ntP->el_hasMCMatch->at(idx);
   _gen_pt = ntP->el_gen_pt->at(idx);
   _gen_eta = ntP->el_gen_eta->at(idx);
   _gen_phi = ntP->el_gen_phi->at(idx);
   _gen_m = ntP->el_gen_m->at(idx);
   _gen_E = ntP->el_gen_E->at(idx);
   _gen_status = ntP->el_gen_status->at(idx);
   _gen_id = ntP->el_gen_id->at(idx);
   _gen_charge = ntP->el_gen_charge->at(idx);
   _gen_dr = ntP->el_gen_dr->at(idx);
}

void Tau::init()
{  
   _fakeType = -1;
   
   _E        = -100;
   _pt       = -100;
   _ptUnc    = -100;
   _eta      = -100;
   _phi      = -100;
   _m        = -100;
   _charge   =    0;
   _id       =    0;

   _isFakeableTTH  = 0;

   _dxy      = -100;
   _dz       = -100;
   
   _passTightCharge    = true;
   _cutEventSel        = true;
   _noLostHits         = true;

   _decayMode              = -1;
   _hasLeadChargedHadrCand = false;
   _leadingTrackPt         = -1.;
   _leadingTrackCharge     = -1;

   _chargedIsoPtSum = -1.;
   _neutralIsoPtSum = -1.;
   _puCorrPtSum     = -1.;

   _decayModeFinding = -1;
//   _decayModeFindingOldDMs = -1;
   _decayModeFindingNewDMs = -1;
   _byLooseCombinedIsolationDeltaBetaCorr3Hits = -1;
   _byMediumCombinedIsolationDeltaBetaCorr3Hits = -1;
   _byTightCombinedIsolationDeltaBetaCorr3Hits = -1;
   _againstElectronVLooseMVA6 = -1;
   _againstElectronLooseMVA6 = -1;
   _againstElectronMediumMVA6 = -1;
   _againstElectronTightMVA6 = -1;

   _byVLooseIsolationMVArun2v1DBdR03oldDMwLT = -1;
   _byLooseIsolationMVArun2v1DBdR03oldDMwLT = -1;
   _byMediumIsolationMVArun2v1DBdR03oldDMwLT = -1;
   _byTightIsolationMVArun2v1DBdR03oldDMwLT = -1;
   _byVTightIsolationMVArun2v1DBdR03oldDMwLT = -1;

   _byCombinedIsolationDeltaBetaCorrRaw3Hits = -1;
   _againstMuonLoose3 = -1;
   _againstMuonTight3 = -1;
   
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

   _hasMCMatch = false;
   _gen_pt = -100;
   _gen_eta = -100;
   _gen_phi = -100;
   _gen_m = -100;
   _gen_E = -100;
   _gen_status = -100;
   _gen_id = -100;
   _gen_charge = -100;
   _gen_dr = -100;
}

void Tau::sel()
{   
   bool pass_pt  = ( _pt > 20. );
   bool pass_eta = ( fabs(_eta) < 2.3 );
   bool pass_dxy = ( fabs(_dxy) < 1000 );
   bool pass_dz  = ( fabs(_dz) < 0.2 );
   
   bool pass_decayModeFinding = ( _decayModeFinding );
   bool pass_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = ( _byVLooseIsolationMVArun2v1DBdR03oldDMwLT );
   
   bool pass_muOverlap = 1;
   int nMuon = nt->NtMuonLoose->size();
   for(int im=0;im<nMuon;im++)
     {
        float dr = GetDeltaR(_eta,_phi,nt->NtMuonLoose->at(im).eta(),nt->NtMuonLoose->at(im).phi());
        if( dr < 0.3 ) pass_muOverlap = 0;
     }  
   
   bool pass_elOverlap = 1;
   int nEl = nt->NtElectronLoose->size();
   for(int iEl=0;iEl<nEl;iEl++)
     {
        float dr = GetDeltaR(_eta,_phi,nt->NtElectronLoose->at(iEl).eta(),nt->NtElectronLoose->at(iEl).phi());
        if( dr < 0.3 ) pass_elOverlap = 0;
     }  
   
   _isFakeableTTH = ( pass_pt &&
		      pass_eta &&
		      pass_dxy &&
		      pass_dz &&
		      pass_decayModeFinding &&
		      pass_byVLooseIsolationMVArun2v1DBdR03oldDMwLT &&
		      pass_muOverlap &&
		      pass_elOverlap );
   
/*   if( evdebug->at(0) == nt->NtEvent->at(0).id() )
     {
	std::cout << "pass_pt = " << pass_pt << " (" << _pt << ")" << std::endl;
	std::cout << "pass_eta = " << pass_eta << " (" << _eta << ")" << std::endl;
	std::cout << "pass_dxy = " << pass_dxy << " (" << _dxy << ")" << std::endl;
	std::cout << "pass_dz = " << pass_dz << " (" << _dz << ")" << std::endl;
	std::cout << "pass_decayModeFinding = " << pass_decayModeFinding << " (" << _decayModeFinding << ")" << std::endl;
	std::cout << "pass_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = " << pass_byVLooseIsolationMVArun2v1DBdR03oldDMwLT << " (" << _byVLooseIsolationMVArun2v1DBdR03oldDMwLT << ")" << std::endl;
	std::cout << "pass_muOverlap = " << pass_muOverlap << std::endl;
	std::cout << "pass_elOverlap = " << pass_elOverlap << std::endl;
     }*/
}

