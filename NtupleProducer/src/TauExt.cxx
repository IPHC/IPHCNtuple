#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(TauExt)

TauExt::TauExt()
{
}

TauExt::~TauExt()
{
}

void TauExt::read(bool isdata)
{
   ID = idx;
   
   matchedJetId = ntP->tau_matchedJetId->at(idx);
   
   E      = ntP->tau_E->at(idx);
   pt     = ntP->tau_pt->at(idx);
   ptCor  = ntP->tau_pt->at(idx);
   ptUnc  = ntP->tau_pt->at(idx);
   eta    = ntP->tau_eta->at(idx);
   phi    = ntP->tau_phi->at(idx);
   m      = ntP->tau_m->at(idx);
   dxy  = ntP->tau_leadingTrackDxy->at(idx);
   dz   = ntP->tau_leadingTrackDz->at(idx);
   charge = ntP->tau_charge->at(idx);
   id     = ntP->tau_id->at(idx);
   
   decayMode = ntP->tau_decayMode->at(idx); 
   hasLeadChargedHadrCand = ntP->tau_hasLeadChargedHadrCand->at(idx); 
   leadingTrackPt = ntP->tau_leadingTrackPt->at(idx); 
   leadingTrackCharge = ntP->tau_leadingTrackCharge->at(idx); 

   decayModeFinding = ntP->tau_decayModeFinding->at(idx);
   decayModeFindingNewDMs = ntP->tau_decayModeFindingNewDMs->at(idx);
   byLooseCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(idx);
   byMediumCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(idx);
   byTightCombinedIsolationDeltaBetaCorr3Hits = ntP->tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(idx);
   againstElectronVLooseMVA6 = ntP->tau_againstElectronVLooseMVA6->at(idx);
   againstElectronLooseMVA6 = ntP->tau_againstElectronLooseMVA6->at(idx);
   againstElectronMediumMVA6 = ntP->tau_againstElectronMediumMVA6->at(idx);
   againstElectronTightMVA6 = ntP->tau_againstElectronTightMVA6->at(idx);

   byVLooseIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
   byLooseIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
   byMediumIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
   byTightIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(idx);
   byVTightIsolationMVArun2v1DBdR03oldDMwLT = ntP->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(idx);

   byCombinedIsolationDeltaBetaCorrRaw3Hits = ntP->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(idx);
   againstMuonLoose3 = ntP->tau_againstMuonLoose3->at(idx);
   againstMuonTight3 = ntP->tau_againstMuonTight3->at(idx);

   chargedIsoPtSum = ntP->tau_chargedIsoPtSum->at(idx);
   neutralIsoPtSum = ntP->tau_neutralIsoPtSum->at(idx);
   puCorrPtSum = ntP->tau_puCorrPtSum->at(idx);
   
   pfEssential_jet_pt 	 = ntP->tau_pfEssential_jet_pt->at(idx);
   pfEssential_jet_eta	 = ntP->tau_pfEssential_jet_eta->at(idx);
   pfEssential_jet_phi	 = ntP->tau_pfEssential_jet_phi->at(idx);
   pfEssential_jet_m  	 = ntP->tau_pfEssential_jet_m->at(idx);
   pfEssential_jetCorr_pt	 = ntP->tau_pfEssential_jetCorr_pt->at(idx);
   pfEssential_jetCorr_eta	 = ntP->tau_pfEssential_jetCorr_eta->at(idx);
   pfEssential_jetCorr_phi	 = ntP->tau_pfEssential_jetCorr_phi->at(idx);
   pfEssential_jetCorr_m	 = ntP->tau_pfEssential_jetCorr_m->at(idx);
   pfEssential_hasSV  	 = ntP->tau_pfEssential_hasSV->at(idx);
   pfEssential_sv_x		 = ntP->tau_pfEssential_sv_x->at(idx);
   pfEssential_sv_y		 = ntP->tau_pfEssential_sv_y->at(idx);
   pfEssential_sv_z		 = ntP->tau_pfEssential_sv_z->at(idx);
   pfEssential_flightLengthSig = ntP->tau_pfEssential_flightLengthSig->at(idx);
   pfEssential_dxy		 = ntP->tau_pfEssential_dxy->at(idx);
   pfEssential_dxy_error	 = ntP->tau_pfEssential_dxy_error->at(idx);
   pfEssential_dxy_Sig	 = ntP->tau_pfEssential_dxy_Sig->at(idx);        

   if( !isdata ) 
   {
	hasMCMatch = ntP->tau_hasMCMatch->at(idx);
	hasChargeMCMatch = ntP->tau_hasChargeMCMatch->at(idx);
	gen_pt = ntP->tau_gen_pt->at(idx);
	gen_eta = ntP->tau_gen_eta->at(idx);
	gen_phi = ntP->tau_gen_phi->at(idx);
	gen_m = ntP->tau_gen_m->at(idx);
	gen_E = ntP->tau_gen_E->at(idx);
	gen_status = ntP->tau_gen_status->at(idx);
	gen_id = ntP->tau_gen_id->at(idx);
	gen_charge = ntP->tau_gen_charge->at(idx);
	gen_dr = ntP->tau_gen_dr->at(idx);
   }  
}

void TauExt::init()
{  
   matchedJetId = -100;
   
   E        = -100;
   pt       = -100;
   ptUnc    = -100;
   eta      = -100;
   phi      = -100;
   m        = -100;
   charge   =    0;
   id       =    0;
   
   isFakeableTTH  = 0;
   isLooseTTH  = 0;
   isMediumTTH  = 0;
   isTightTTH  = 0;

   dxy      = -100;
   dz       = -100;
   
   passTightCharge    = true;
   cutEventSel        = true;
   noLostHits         = true;

   decayMode              = -1;
   hasLeadChargedHadrCand = false;
   leadingTrackPt         = -1.;
   leadingTrackCharge     = -1;

   chargedIsoPtSum = -1.;
   neutralIsoPtSum = -1.;
   puCorrPtSum     = -1.;

   decayModeFinding = -1;
   decayModeFindingNewDMs = -1;
   byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   againstElectronVLooseMVA6 = 0;
   againstElectronLooseMVA6 = 0;
   againstElectronMediumMVA6 = 0;
   againstElectronTightMVA6 = 0;

   byVLooseIsolationMVArun2v1DBdR03oldDMwLT = 0;
   byLooseIsolationMVArun2v1DBdR03oldDMwLT = 0;
   byMediumIsolationMVArun2v1DBdR03oldDMwLT = 0;
   byTightIsolationMVArun2v1DBdR03oldDMwLT = 0;
   byVTightIsolationMVArun2v1DBdR03oldDMwLT = 0;

   byCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   againstMuonLoose3 = 0;
   againstMuonTight3 = 0;
   
   pfEssential_jet_pt          = -1.;
   pfEssential_jet_eta         = -1.;
   pfEssential_jet_phi         = -1.;
   pfEssential_jet_m           = -1.;
   pfEssential_jetCorr_pt      = -1.;
   pfEssential_jetCorr_eta     = -1.;
   pfEssential_jetCorr_phi     = -1.;
   pfEssential_jetCorr_m       = -1.;
   pfEssential_hasSV           = -1.;
   pfEssential_sv_x            = -1.;
   pfEssential_sv_y            = -1.;
   pfEssential_sv_z            = -1.;
   pfEssential_flightLengthSig = -1.;
   pfEssential_dxy             = -1.;
   pfEssential_dxy_error       = -1.;
   pfEssential_dxy_Sig         = -1.;

   hasMCMatch = false;
   hasChargeMCMatch = false;
   hasPhotonMCMatch = false;
   gen_pt = -100;
   gen_eta = -100;
   gen_phi = -100;
   gen_m = -100;
   gen_E = -100;
   gen_status = -100;
   gen_id = -100;
   gen_charge = -100;
   gen_dr = -100;
}

void TauExt::sel(bool DEBUG,int year)
{
   bool pass_pt  = ( pt > 20. );
   bool pass_eta = ( fabs(eta) < 2.3 );
   bool pass_dxy = ( fabs(dxy) < 1000 );
   bool pass_dz  = ( fabs(dz) < 0.2 );
   
   bool pass_decayModeFinding = ( decayModeFinding );
   bool pass_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = ( byVLooseIsolationMVArun2v1DBdR03oldDMwLT );
   bool pass_byLooseIsolationMVArun2v1DBdR03oldDMwLT = ( byLooseIsolationMVArun2v1DBdR03oldDMwLT );
   bool pass_byMediumIsolationMVArun2v1DBdR03oldDMwLT = ( byMediumIsolationMVArun2v1DBdR03oldDMwLT );
   
   bool pass_muOverlap = 1;
   int nMuon = nt->NtMuonLooseExt->size();
   for(int im=0;im<nMuon;im++)
     {
        float dr = GetDeltaR(eta,phi,nt->NtMuonLooseExt->at(im).eta,nt->NtMuonLooseExt->at(im).phi);
        if( dr < 0.3 ) pass_muOverlap = 0;
     }  
   
   bool pass_elOverlap = 1;
   int nEl = nt->NtElectronLooseExt->size();
   for(int iEl=0;iEl<nEl;iEl++)
     {
        float dr = GetDeltaR(eta,phi,nt->NtElectronLooseExt->at(iEl).eta,nt->NtElectronLooseExt->at(iEl).phi);
        if( dr < 0.3 ) pass_elOverlap = 0;
     }  
   
   isFakeableTTH = ( pass_pt &&
		     pass_eta &&
		     pass_dxy &&
		     pass_dz &&
		     pass_decayModeFinding &&
		     pass_byVLooseIsolationMVArun2v1DBdR03oldDMwLT &&
		     pass_muOverlap &&
		     pass_elOverlap );
   
   isMediumTTH = ( isFakeableTTH &&
		   pass_byMediumIsolationMVArun2v1DBdR03oldDMwLT );
   
   isLooseTTH = ( isFakeableTTH &&
		  pass_byLooseIsolationMVArun2v1DBdR03oldDMwLT );
		  
		  

	if(DEBUG)
	  {
	     std::cout << "------------------------------" << std::endl;
	     std::cout << "Event #" << std::setprecision(12) << ntP->ev_id << std::endl;
	     std::cout << "  Tau #" << ID << std::endl;
	     std::cout << "  pt = " << pt << std::endl;
	     std::cout << "  eta = " << eta << std::endl;
	     std::cout << "  phi = " << phi << std::endl;
	     std::cout << "  pass_pt = " << pass_pt << std::endl;
	     std::cout << "  pass_eta = " << pass_eta << std::endl;
	     std::cout << "  pass_dxy = " << pass_dxy << std::endl;
	     std::cout << "  pass_elOverlap = " << pass_elOverlap << std::endl;
	     std::cout << "  pass_muOverlap = " << pass_muOverlap << std::endl;
	     std::cout << "  isFakeableTTH = " << isFakeableTTH << std::endl;
	     std::cout << "  isLooseTTH = " << isLooseTTH << std::endl;
	     std::cout << "  hasMCMatch = " << hasMCMatch << std::endl;
	  }		  
}

