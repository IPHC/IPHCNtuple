#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(MuonExt)

MuonExt::MuonExt()
{
}

MuonExt::~MuonExt()
{
}

void MuonExt::read(bool isdata)
{
   ID                             = idx;

   matchedJetId                   = ntP->mu_matchedJetId->at(idx);
   
   E                              = ntP->mu_E->at(idx);
   pt                             = ntP->mu_pt->at(idx);
   ptUnc                          = ntP->mu_pt->at(idx);
   eta                            = ntP->mu_eta->at(idx);
   phi                            = ntP->mu_phi->at(idx);
   m                              = ntP->mu_m->at(idx);
   conept                         = ntP->mu_conept->at(idx);
   charge                         = ntP->mu_charge->at(idx);
   id                             = ntP->mu_id->at(idx);
   
   isLoose                        = ntP->mu_isLooseMuon->at(idx);
   isMedium	                    = ntP->mu_isMediumMuon->at(idx);
   isTight	                    = ntP->mu_isTightMuon->at(idx);
   isPFMuon                       = ntP->mu_isPFMuon->at(idx);

   dxy                            = ntP->mu_innerTrack_PV_dxy->at(idx);
   dz                             = ntP->mu_bestTrack_PV_dz->at(idx); //changed
   //dz                             = ntP->mu_innerTrack_PV_dz->at(idx); //switched to bestTrack, as in nanoAOD
   iso                            = ntP->mu_miniIsoTTH->at(idx);
   sip3d                          = ntP->mu_ip3d->at(idx) / ntP->mu_ip3dErr->at(idx);
   bestTrack_pt                   = ntP->mu_bestTrack_pt->at(idx);
   bestTrack_ptError              = ntP->mu_bestTrack_ptError->at(idx);

   tightCharge = (bestTrack_pt) ? (bestTrack_ptError/bestTrack_pt < 0.2) : 0;
   
   lepMVA	                        = ntP->mu_lepMVA->at(idx);
      
   lepMVA_miniRelIsoCharged       = ntP->mu_lepMVA_miniRelIsoCharged->at(idx);
   lepMVA_miniRelIsoNeutral       = ntP->mu_lepMVA_miniRelIsoNeutral->at(idx);
   lepMVA_jetPtRelv2              = ntP->mu_lepMVA_jetPtRelv2->at(idx);
   lepMVA_jetPtRatio              = ntP->mu_lepMVA_jetPtRatio->at(idx);
   jetRelIso                      = ntP->mu_jetRelIso->at(idx);
   lepMVA_jetBTagCSV              = ntP->mu_lepMVA_jetBTagCSV->at(idx);
   lepMVA_jetBTagDeepCSV          = ntP->mu_lepMVA_jetBTagDeepCSV->at(idx);
   lepMVA_jetBTagDeepFlavour      = ntP->mu_lepMVA_jetBTagDeepFlavour->at(idx);
   lepMVA_sip3d                   = ntP->mu_lepMVA_sip3d->at(idx);
   lepMVA_dxy                     = ntP->mu_lepMVA_dxy->at(idx);
   lepMVA_dz                      = ntP->mu_lepMVA_dz->at(idx);
   lepMVA_mvaId                   = ntP->mu_lepMVA_mvaId->at(idx);
   lepMVA_eta                     = ntP->mu_lepMVA_eta->at(idx);
   lepMVA_jetNDauChargedMVASel    = ntP->mu_lepMVA_jetNDauChargedMVASel->at(idx);


   if( !isdata )
   {
      hasMCMatch = ntP->mu_hasMCMatch->at(idx);
      hasChargeMCMatch = ntP->mu_hasChargeMCMatch->at(idx);
      gen_pt = ntP->mu_gen_pt->at(idx);
      gen_eta = ntP->mu_gen_eta->at(idx);
      gen_phi = ntP->mu_gen_phi->at(idx);
      gen_m = ntP->mu_gen_m->at(idx);
      gen_E = ntP->mu_gen_E->at(idx);
      gen_status = ntP->mu_gen_status->at(idx);
      gen_id = ntP->mu_gen_id->at(idx);
      gen_charge = ntP->mu_gen_charge->at(idx);
      gen_dr = ntP->mu_gen_dr->at(idx);
   }
   
   PFRelIso04 = ntP->mu_PFRelIso04->at(idx);
}

void MuonExt::init()
{
   fakeType          = -100.;
   
   matchedJetId      = -100;
   
   E                 = -100.;
   pt                = -100.;
   ptUnc             = -100.;
   eta               = -100.;
   phi               = -100.;
   m                 = -100.;
   conept            = -100.;
   charge            = 0;
   id                = 0;
   
   isLoose           = 0;
   isMedium          = 0;
   isTight           = 0;
   isPFMuon          = 0;
   
   isLooseTTH        = 0;
   isFakeableTTH     = 0;
   isMediumTTH       = 0;
   isTightTTH        = 0;
   
   dxy                = -100.;
   dz                 = -100.;
   iso                = -100.;
   PFRelIso04         = -100.;
   sip3d              = -100.;
   bestTrack_pt       = -100.;
   bestTrack_ptError  = -100.;

   cutEventSel        = 1;
   noLostHits         = 1;
   
   tightCharge = 0;
   
   lepMVA                      = -100.;
   
   lepMVA_miniRelIsoCharged    = -100.;
   lepMVA_miniRelIsoNeutral    = -100.;
   lepMVA_jetPtRelv2           = -100.;
   lepMVA_jetPtRatio           = -100.;
   jetRelIso                   = -100.;
   lepMVA_jetBTagCSV           = -100.;
   lepMVA_jetBTagDeepCSV       = -100.;
   lepMVA_jetBTagDeepFlavour   = -100.;
   lepMVA_sip3d                = -100.;
   lepMVA_dxy                  = -100.;
   lepMVA_dz                   = -100.;
   lepMVA_mvaId                = -100.;
   lepMVA_eta                  = -100.;
   lepMVA_jetNDauChargedMVASel = -100.;   

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

void MuonExt::sel(bool DEBUG,int year)
{
   bool pass_pt      = ( pt > 5 );
   bool pass_eta     = ( fabs(eta) < 2.4 );
   bool pass_dxy     = ( fabs(dxy) < 0.05 );
   bool pass_dz      = ( fabs(dz) < 0.1 );
   bool pass_miniIso = ( iso < 0.4 );
   bool pass_SIP     = ( fabs(sip3d) < 8 );
   bool pass_isLoose = ( isLoose );
   
   float EffArea = getEffArea(eta,year);
   
   //Changed definition, ttH uses "deltaBeta" for Muons and "rhoArea" correction for electrons -- see : https://github.com/peruzzim/cmg-cmssw/blob/heppy_94X_dev_ttH/PhysicsTools/Heppy/python/analyzers/objects/LeptonAnalyzer.py#L331
   float isoR04 = (pt > 0.) ? (ntP->mu_pfIso04_sumChargedHadronPt->at(idx) + std::max( 0.0, double(ntP->mu_pfIso04_sumNeutralHadronEt->at(idx)+ntP->mu_pfIso04_sumPhotonEt->at(idx) - ntP->mu_pfIso04_sumPUPt->at(idx) / 2. )))/pt : -9999;
//   float isoR04 = (pt > 0.) ? (ntP->mu_pfIso04_sumChargedHadronPt->at(idx) + std::max( 0.0, double(ntP->mu_pfIso04_sumNeutralHadronEt->at(idx)+ntP->mu_pfIso04_sumPhotonEt->at(idx) - ntP->ev_rho*EffArea )))/pt : -9999; //previous
   
   isLooseTTH = ( pass_pt      &&
		  pass_eta     &&
		  pass_dxy     &&
		  pass_dz      &&
		  pass_miniIso &&
		  pass_SIP     &&
		  pass_isLoose );
   
   bool pass_fakeable_lepMVA = 1;
   bool pass_lepMVA = ( lepMVA >= 0.85 );
   bool pass_jetSel = 1;
   float btagCut = -777;

   bool pass_clJet = 1;
   
   if( lepMVA < 0.85 )
     {
	float jetpt = 0.9*pt*(1.+jetRelIso);
	btagCut = smoothBFlav(jetpt,20.,45.,year);
	pass_jetSel = (lepMVA_jetBTagDeepFlavour < btagCut);
	if( !(jetRelIso < 0.5 && pass_jetSel) ) pass_fakeable_lepMVA = 0;
     }
   else
     {
	if( year == 2016 )
	  if( lepMVA_jetBTagDeepFlavour > 0.3093 ) pass_clJet = 0;
	if( year == 2017 )
	  if( lepMVA_jetBTagDeepFlavour > 0.3033 ) pass_clJet = 0;
	if( year == 2018 )
	  if( lepMVA_jetBTagDeepFlavour > 0.2770 ) pass_clJet = 0;
     }   
   
   bool pass_conept = ( conept > 10. );
   
   isFakeableTTH = ( isLooseTTH &&
		     pass_fakeable_lepMVA &&
		     pass_clJet &&
		     pass_conept );
   
   isTightTTH = ( isFakeableTTH && 
		  pass_lepMVA && 
		  isMedium );

	if(DEBUG)
	  {
	     std::cout << "------------------------------" << std::endl;
	     std::cout << "Event #" << std::setprecision(12) << ntP->ev_id << std::endl;
	     std::cout << "  muon #" << ID << std::endl << std::endl;
	     std::cout << "  conept = " << conept << std::endl;
	     std::cout << "  pt = " << pt << std::endl;
	     std::cout << "  eta = " << eta << std::endl;
	     std::cout << "  phi = " << phi << std::endl;
	     std::cout << "  E = " << E << std::endl << std::endl;
	     std::cout << "  charge = " << charge << std::endl << std::endl;
	     std::cout << "  dxy = " << dxy << std::endl;
	     std::cout << "  dz = " << dz << std::endl;
	     std::cout << "  iso = " << iso << std::endl;
	     std::cout << "  isoR04 = " << isoR04 << std::endl;
	     std::cout << "  btagCut = " << btagCut << std::endl;
	     std::cout << "  miniIsoCharged = " << lepMVA_miniRelIsoCharged << std::endl;
	     std::cout << "  mu_pfIso04_sumChargedHadronPt = " << ntP->mu_pfIso04_sumChargedHadronPt->at(idx)/pt << std::endl;
	     std::cout << "  miniIsoNeutral = " << lepMVA_miniRelIsoNeutral << std::endl;
	     std::cout << "  sip3d = " << sip3d << std::endl;
	     std::cout << "  lepMVA_jetPtRatio = " << lepMVA_jetPtRatio << std::endl;
	     std::cout << "  jetRelIso = " << jetRelIso << std::endl;
	     std::cout << "  lepMVA = " << lepMVA << std::endl;
	     std::cout << "  lepMVA_jetBTagDeepFlavour = " << lepMVA_jetBTagDeepFlavour << std::endl;
	     std::cout << "  lepMVA_mvaId = " << lepMVA_mvaId << std::endl;
	     std::cout << "  PFRelIso04 = " << PFRelIso04 << std::endl;
	     std::cout << "  isLoose = " << isLoose << std::endl;
	     std::cout << "  isLooseTTH = " << isLooseTTH << std::endl << std::endl;
	     std::cout << "  pass_fakeable_lepMVA = " << pass_fakeable_lepMVA << std::endl;
	     std::cout << "  pass_clJet = " << pass_clJet << std::endl;
	     std::cout << "  pass_conept = " << pass_conept << std::endl;
	     std::cout << "  = isFakeableTTH = " << isFakeableTTH << std::endl << std::endl;
	     std::cout << "  pass_lepMVA = " << pass_lepMVA << std::endl;
	     std::cout << "  isMedium = " << isMedium << std::endl;
	     std::cout << "  isTightTTH = " << isTightTTH << std::endl;
	     std::cout << "  hasMCMatch = " << hasMCMatch << std::endl;
 	}		  
}

float MuonExt::getEffArea(float eta,int year)
{   
   float ea = -1;
   
   if( year == 2016 )
     {	
	if(fabs(eta) < 0.8)      ea = 0.0735;
	else if(fabs(eta) < 1.3) ea = 0.0619;
	else if(fabs(eta) < 2.0) ea = 0.0465;
	else if(fabs(eta) < 2.2) ea = 0.0433;
	else                     ea = 0.0577;
     }
   else if( year == 2017 || year == 2018 )
     {	
	if(fabs(eta) < 0.8)      ea = 0.0566;
	else if(fabs(eta) < 1.3) ea = 0.0562;
	else if(fabs(eta) < 2.0) ea = 0.0363;
	else if(fabs(eta) < 2.2) ea = 0.0119;
	else                     ea = 0.0064;
     }
      
   ea*=16./9.;
   
   return ea;
}

float MuonExt::smoothBFlav(float jetpt,float ptmin,float ptmax,int year)
{   
   float wploose[3] = {0.0614, 0.0521, 0.0494};
   float wpmedium[3] = {0.3093, 0.3033, 0.2770};
   
   float x = std::min(std::max(0.f, jetpt - ptmin)/(ptmax-ptmin), 1.f);
   
   return x*wploose[year-2016] + (1-x)*wpmedium[year-2016]; 
}
