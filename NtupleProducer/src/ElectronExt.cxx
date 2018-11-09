#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(ElectronExt)

ElectronExt::ElectronExt()
{
}

ElectronExt::~ElectronExt()
{
}

void ElectronExt::read(bool isdata)
{
   ID                                = idx;
   
   E 	                            = ntP->el_E->at(idx);
   pt	                            = ntP->el_pt->at(idx);
   ptUnc                             = ntP->el_pt->at(idx);
   eta	                            = ntP->el_eta->at(idx);
   phi	                            = ntP->el_phi->at(idx);
   m 	                            = ntP->el_m->at(idx);
   conept	                    = ntP->el_conept->at(idx);
   charge	                        = ntP->el_charge->at(idx);
   id	                            = ntP->el_id->at(idx);
   
   dxy	                            = ntP->el_gsfTrack_PV_dxy->at(idx);
   dz	                            = ntP->el_gsfTrack_PV_dz->at(idx);
   nlosthits                         = ntP->el_numberOfLostHits->at(idx);
   ip3d	                            = ntP->el_ip3d->at(idx);
   ip3dErr	                        = ntP->el_ip3dErr->at(idx);
   
   isLoose	                        = ntP->el_looseCBId->at(idx);
   isMedium                          = ntP->el_mediumCBId->at(idx);
   passCV	                        = ntP->el_passConversionVeto->at(idx);
   isGsfCtfScPixChargeConsistent     = ntP->el_isGsfCtfScPixChargeConsistent->at(idx);
   isGsfScPixChargeConsistent     = ntP->el_isGsfScPixChargeConsistent->at(idx);

   sip3d                          = (ntP->el_ip3dErr->at(idx)) ? ntP->el_ip3d->at(idx) / ntP->el_ip3dErr->at(idx) : -9999;
   ooEmooP                        = ntP->el_ooEmooP->at(idx);
   
   lepMVA	                        = ntP->el_lepMVA->at(idx);
   mvaIso                      = ntP->el_mvaIso->at(idx);
   mvaNoIso                      = ntP->el_mvaNoIso->at(idx);
   
   lepMVA_miniRelIsoCharged          = ntP->el_lepMVA_miniRelIsoCharged->at(idx);
   lepMVA_miniRelIsoNeutral          = ntP->el_lepMVA_miniRelIsoNeutral->at(idx);
   lepMVA_jetPtRelv2                 = ntP->el_lepMVA_jetPtRelv2->at(idx);
   lepMVA_jetPtRatio                 = ntP->el_lepMVA_jetPtRatio->at(idx);
   lepMVA_jetBTagCSV                 = ntP->el_lepMVA_jetBTagCSV->at(idx);
   lepMVA_jetBTagDeepCSV             = ntP->el_lepMVA_jetBTagDeepCSV->at(idx);
   lepMVA_sip3d                      = ntP->el_lepMVA_sip3d->at(idx);
   lepMVA_dxy                        = ntP->el_lepMVA_dxy->at(idx);
   lepMVA_dz                         = ntP->el_lepMVA_dz->at(idx);
   lepMVA_mvaId                      = ntP->el_lepMVA_mvaId->at(idx);
   lepMVA_eta                        = ntP->el_lepMVA_eta->at(idx);
   lepMVA_jetNDauChargedMVASel       = ntP->el_lepMVA_jetNDauChargedMVASel->at(idx);
   
   miniIso			                = ntP->el_miniIsoTTH->at(idx);
   
   sigmaIetaIeta		                = ntP->el_full5x5_sigmaIetaIeta->at(idx);
   superCluster_eta  	            = ntP->el_superCluster_eta->at(idx);
   hadronicOverEm		            = ntP->el_hadronicOverEm->at(idx);
   deltaEtaSuperClusterTrackAtVtx    = ntP->el_deltaEtaSuperClusterTrackAtVtx->at(idx);
   deltaPhiSuperClusterTrackAtVtx    = ntP->el_deltaPhiSuperClusterTrackAtVtx->at(idx);
   eSuperClusterOverP	            = ntP->el_eSuperClusterOverP->at(idx);
   correctedEcalEnergy	            = ntP->el_correctedEcalEnergy->at(idx);
   ecalEnergy	                    = ntP->el_ecalEnergy->at(idx);
   
   trackMomentumError	            = ntP->el_trackMomentumError->at(idx);
   tightCharge = (isGsfCtfScPixChargeConsistent + isGsfScPixChargeConsistent > 1);

   if( !isdata ) 
   {
	hasMCMatch = ntP->el_hasMCMatch->at(idx);
	gen_pt = ntP->el_gen_pt->at(idx);
	gen_eta = ntP->el_gen_eta->at(idx);
	gen_phi = ntP->el_gen_phi->at(idx);
	gen_m = ntP->el_gen_m->at(idx);
	gen_E = ntP->el_gen_E->at(idx);
	gen_status = ntP->el_gen_status->at(idx);
	gen_id = ntP->el_gen_id->at(idx);
	gen_charge = ntP->el_gen_charge->at(idx);
	gen_dr = ntP->el_gen_dr->at(idx);
   } 
}

void ElectronExt::init()
{   
   E                              = -100.;
   pt                             = -100.;
   ptUnc                          = -100.;
   eta                            = -100.;
   phi                            = -100.;
   m                              = -100.;
   conept                         = -100.;
   charge                         = 0;
   id                             = 0;      
   
   isLooseCBId           	        = 0;
   isMediumCBId          	        = 0;
   isLoose               	        = 0;
   isMedium                  	    = 0;
   isTight                   	    = 0;
   isLooseMVA                	    = 0;
   isTightMVA             	    = 0;
   isLooseTTH            	        = 0;
   isFakeableTTH             	    = 0;
   isMediumTTH                	    = 0;
   isTightTTH                	    = 0;
   
   dxy                            = -100;
   dz                             = -100;
   miniIso                        = -100;
   isoR04                         = -100;
   nlosthits                      = -100;
   sip3d                          = -100;
   ooEmooP                        = -100;
   
   passCV                         = 0;
   isGsfCtfScPixChargeConsistent  = 0;
   isGsfScPixChargeConsistent     = 0;
   passPtEta                      = 0;
   ip3d                           = -100.;
   ip3dErr                        = -100.;
   
   lepMVA                         = -100;

   lepMVA_miniRelIsoCharged       = -100.;
   lepMVA_miniRelIsoNeutral       = -100.;
   lepMVA_jetPtRelv2              = -100.;
   lepMVA_jetPtRatio              = -100.;
   lepMVA_jetBTagCSV              = -100.;
   lepMVA_jetBTagDeepCSV          = -100.;
   lepMVA_sip3d                   = -100.;
   lepMVA_dxy                     = -100.;
   lepMVA_dz                      = -100.;
   lepMVA_mvaId                   = -100.;
   
   lepMVA_eta                     = -100.;
   lepMVA_jetNDauChargedMVASel    = -100.;
   
   passChargeFlip                 = 0;
   hasMatchedConversion           = 0;
   isGsfCtfScPixChargeConsistent  = 0;
   
   sigmaIetaIeta                  = -100.;
   hadronicOverEm                 = -100.;
   correctedEcalEnergy            = -100.;
   ecalEnergy                     = -100.;
   eSuperClusterOverP             = -100.;
   deltaEtaSuperClusterTrackAtVtx = -100.;
   deltaPhiSuperClusterTrackAtVtx = -100.;
   see                            = -100.;
   superCluster_eta               = -100.;
   
   trackMomentumError             = -100;
   tightCharge                    = 0;
   cutEventSel                    = false;
   noLostHits                     = false;
   mvaIso                   = -100;
   mvaNoIso                   = -100;

   hasMCMatch = false;
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

void ElectronExt::sel()
{   
   bool isLoose = false;

   if( pt <= 10. )
     {
	if( fabs(eta) < 0.8 ) isLoose = ( mvaNoIso > -0.13285867293779202 );
	else if( fabs(eta) < 1.479 ) isLoose = ( mvaNoIso > -0.31765300958836074 );
	else isLoose = ( mvaNoIso > -0.0799205914718861 );
     }
   else
     {
	if( fabs(eta) < 0.8 ) isLoose = ( mvaNoIso > -0.856871961305474 );
	else if( fabs(eta) < 1.479 ) isLoose = ( mvaNoIso > -0.8107642141584835 );
	else isLoose = ( mvaNoIso > -0.7179265933023059 );
     }   
   
   float EffArea = getEffArea(superCluster_eta);
   isoR04 = (pt > 0.) ? (ntP->el_pfIso_sumChargedHadronPt->at(idx) + std::max( 0.0, double(ntP->el_pfIso_sumNeutralHadronEt->at(idx)+ntP->el_pfIso_sumPhotonEt->at(idx) - ntP->ev_rho*EffArea )))/pt : -9999;
   
   bool pass_pt       = ( pt > 7 );
   bool pass_eta      = ( fabs(eta) < 2.5 );
   bool pass_dxy      = ( fabs(dxy) < 0.05 );
   bool pass_dz       = ( fabs(dz) < 0.1 );
   bool pass_miniIso  = ( miniIso < 0.4 );
   bool pass_SIP      = ( fabs(sip3d) < 8 );
   bool pass_isLoose  = ( isLoose );
   bool pass_losthits = ( nlosthits < 2 );
   
   bool passMuOverlap = 1;
   int nMuonLoose = nt->NtMuonLooseExt->size();
   for(int im=0;im<nMuonLoose;im++)
     {
        float dr = GetDeltaR(eta,phi,nt->NtMuonLooseExt->at(im).eta,nt->NtMuonLooseExt->at(im).phi);
        if( dr < 0.05 ) passMuOverlap = 0;
     }
   
   isLooseTTH = ( pass_pt          &&
		  pass_eta         &&
		  pass_dxy         &&
		  pass_dz          &&
		  pass_miniIso     &&
		  pass_SIP         &&
		  pass_isLoose     &&
		  pass_losthits    &&
		  passMuOverlap   );

   bool pass_fakeable_lepMVA = 0;
   bool pass_conept = ( conept > 10. );
   bool pass_lepMVA = ( lepMVA >= 0.90 );
   
   if( pass_lepMVA )
     {	
	if( lepMVA_jetBTagDeepCSV < 0.4941 ) pass_fakeable_lepMVA = 1;
     }   
   else
     {	
	if( lepMVA_jetPtRatio > 0.6 && lepMVA_jetBTagDeepCSV < 0.07 && lepMVA_mvaId > 0.5 ) pass_fakeable_lepMVA = 1;
     }   
   
   bool pass_sc = 0;
   float eInvMinusPInv  = (ecalEnergy > 0) ? (1./ecalEnergy-eSuperClusterOverP/ecalEnergy) : 99;
   if( fabs(superCluster_eta) < 1.479 )
     {
	if( sigmaIetaIeta < 0.011 && hadronicOverEm < 0.10 && eInvMinusPInv > -0.04 ) pass_sc = 1;
     }
   else if( fabs(superCluster_eta) < 2.5 )
     {
	if( sigmaIetaIeta < 0.030 && hadronicOverEm < 0.10 && eInvMinusPInv > -0.04 ) pass_sc = 1;
     }
   
   bool pass_nlosthits = (nlosthits == 0);

   isFakeableTTH = ( isLooseTTH &&
		     pass_sc &&
		     pass_nlosthits &&
		     passCV &&
		     pass_fakeable_lepMVA &&
		     pass_conept );
   
   isTightTTH = ( isFakeableTTH &&
		  pass_lepMVA );

   for(int d=0;d<evdebug->size();d++)
     {		       
	double evId = ntP->ev_id;
	if( evId == evdebug->at(d) )
	  {
	     std::cout << "------------------------------" << std::endl;
	     std::cout << "Event #" << evId << std::endl;
	     std::cout << "  electron #" << ID << std::endl;
	     std::cout << "  conept = " << conept << std::endl;
	     std::cout << "  pt = " << pt << std::endl;
	     std::cout << "  isGsfCtfScPixChargeConsistent = " << isGsfCtfScPixChargeConsistent << std::endl;
	     std::cout << "  isGsfScPixChargeConsistent = " << isGsfScPixChargeConsistent << std::endl;
	     std::cout << "  tightCharge = " << tightCharge << std::endl;
	     std::cout << "  isLooseTTH = " << isLooseTTH << std::endl;
	     std::cout << "  pass_sc = " << pass_sc << std::endl;
	     std::cout << "  pass_fakeable_lepMVA = " << pass_fakeable_lepMVA << std::endl;
	     std::cout << "  pass_conept = " << pass_conept << std::endl;
	     std::cout << "  = isFakeableTTH = " << isFakeableTTH << std::endl;
	  }		  
     }		  
}

float ElectronExt::getEffArea(float eta)
{   
   float ea = -1;
   
   if(fabs(eta) < 1.0)        ea = 0.1566;
   else if(fabs(eta) < 1.479) ea = 0.1626;
   else if(fabs(eta) < 2.0)   ea = 0.1073;
   else if(fabs(eta) < 2.2)   ea = 0.0854;
   else if(fabs(eta) < 2.3)   ea = 0.1051;
   else if(fabs(eta) < 2.4)   ea = 0.1204;
   else                       ea = 0.1524;
   
   return ea;
}
