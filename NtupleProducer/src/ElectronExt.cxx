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

   E_preCorr                        = ntP->el_ecalTrkEnergyPreCorr->at(idx);
   E 	                            = ntP->el_ecalTrkEnergyPostCorr->at(idx);
   pt_preCorr                       = ntP->el_pt->at(idx);
   
   float smearCorr = E / E_preCorr;
   float pt_postCorr = pt_preCorr * smearCorr;
   
   pt 	 			    = pt_postCorr;
  
   eta	                            = ntP->el_eta->at(idx);
   phi	                            = ntP->el_phi->at(idx);
   m 	                            = ntP->el_m->at(idx);
   conept	                    = ntP->el_conept->at(idx);
   charge	                    = ntP->el_charge->at(idx);
   id	                            = ntP->el_id->at(idx);
   
   dxy	                            = ntP->el_gsfTrack_PV_dxy->at(idx);
   dz	                            = ntP->el_gsfTrack_PV_dz->at(idx);
   nlosthits                         = ntP->el_numberOfLostHits->at(idx);
   ip3d	                            = ntP->el_ip3d->at(idx);
   ip3dErr	                        = ntP->el_ip3dErr->at(idx);
   
//   isLoose	                        = ntP->el_looseCBId->at(idx);
   isLoose	                        = ntP->el_NoIsoLooseMVAId->at(idx);
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
   lepMVA_jetBTagDeepFlavour         = ntP->el_lepMVA_jetBTagDeepFlavour->at(idx);
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
	hasChargeMCMatch = ntP->el_hasChargeMCMatch->at(idx);
	hasPhotonMCMatch = ntP->el_hasPhotonMCMatch->at(idx);
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
   
   PFRelIso04 = ntP->el_PFRelIso04->at(idx);
}

void ElectronExt::init()
{   
   E                              = -100.;
   pt                             = -100.;
   pt_preCorr                     = -100.;
   E_preCorr			  = -100.;
   
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
   PFRelIso04                     = -100;
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
   lepMVA_jetBTagDeepFlavour      = -100.;
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

void ElectronExt::sel(bool DEBUG,int year)
{   
//   bool isLoose = false;

/*   if( pt <= 10. )
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
     }*/
   
   float EffArea = getEffArea(superCluster_eta,year);
   
   //CHANGED -- was using cone size of 0.3 instead of 0.4 ==> Need to call different variables, and rescale effArea
   //isoR04 = (pt > 0.) ? (ntP->el_chargedHadronIso->at(idx) + std::max( 0.0, double(ntP->el_neutralHadronIso->at(idx)+ntP->el_photonIso->at(idx) - ntP->ev_rho*EffArea )))/pt : -9999;
   
   //isoR04 = (pt > 0.) ? (ntP->el_pfIso_sumChargedHadronPt->at(idx) + std::max( 0.0, double(ntP->el_pfIso_sumNeutralHadronEt->at(idx)+ntP->el_pfIso_sumPhotonEt->at(idx) - ntP->ev_rho*EffArea )))/pt : -9999; //previous
   
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
        if( dr < 0.3 ) passMuOverlap = 0; //CHANGED from 0.05 (agreed in ttH meeting)
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

	if(DEBUG)
	  {
	     std::cout << "------------------------------" << std::endl;
	     std::cout << "Event #" << std::setprecision(12) << ntP->ev_id << std::endl;
	     std::cout << "  electron #" << ID << std::endl << std::endl;
	     std::cout << "  conept = " << conept << std::endl;
	     std::cout << "  pt = " << pt << std::endl;
	     std::cout << "  eta = " << eta << std::endl;
	     std::cout << " phi = " << phi << std::endl << std::endl;
	     std::cout << " E = " << E << std::endl << std::endl;
	     std::cout << " charge = " << charge << std::endl << std::endl;
	     std::cout << " isLoose = " << isLoose << std::endl;
	     std::cout << " passMuOverlap = " << passMuOverlap << std::endl;
	     std::cout << " sip3d = " << sip3d << std::endl;
	     std::cout << " dxy = " << dxy << std::endl;
	     std::cout << " dz = " << dz << std::endl;
	     std::cout << " miniIso = " << miniIso << std::endl;
	     std::cout << " PFRelIso04 = " << PFRelIso04 << std::endl;
	     std::cout << " mvaNoIso = " << mvaNoIso << std::endl;
	     std::cout << " lepMVA = " << lepMVA << std::endl;
	     std::cout << " lepMVA_jetPtRatio = " << lepMVA_jetPtRatio << std::endl;
	     std::cout << " lepMVA_jetBTagDeepCSV = " << lepMVA_jetBTagDeepCSV << std::endl;
	     std::cout << " lepMVA_mvaId = " << lepMVA_mvaId << std::endl;
	     std::cout << " isGsfCtfScPixChargeConsistent = " << isGsfCtfScPixChargeConsistent << std::endl;
	     std::cout << " isGsfScPixChargeConsistent = " << isGsfScPixChargeConsistent << std::endl;
	     std::cout << " tightCharge = " << tightCharge << std::endl;
	     std::cout << " isLooseTTH = " << isLooseTTH << std::endl << std::endl;
	     std::cout << " pass_sc = " << pass_sc << std::endl;
     	     std::cout << " nlosthits = " << nlosthits << std::endl;
     	     std::cout << " passCV = " << passCV << std::endl;
	     std::cout << " pass_fakeable_lepMVA = " << pass_fakeable_lepMVA << std::endl;
	     std::cout << " pass_conept = " << pass_conept << std::endl;
	     std::cout << " isFakeableTTH = " << isFakeableTTH << std::endl<<std::endl;
	     std::cout << " pass_lepMVA = " << pass_lepMVA << std::endl;
	     std::cout << " isTightTTH = " << isTightTTH << std::endl;
	     std::cout << " hasMCMatch = " << hasMCMatch << std::endl;
	  }		  
}


//EffArea values for DR=0.3, rescaled for DR=0.4 (also done in FlatTreePeroducer.cc -- for electrons, different for muons)
float ElectronExt::getEffArea(float eta,int year)
{   
   float ea = -1;

   if( year == 2016 )
     {
	if(fabs(eta) < 1.0)        ea = 0.1752;
	else if(fabs(eta) < 1.479) ea = 0.1862;
	else if(fabs(eta) < 2.0)   ea = 0.1411;
	else if(fabs(eta) < 2.2)   ea = 0.1534;
	else if(fabs(eta) < 2.3)   ea = 0.1903;
	else if(fabs(eta) < 2.4)   ea = 0.2243;
	else                       ea = 0.2687;		
     }
   else if( year == 2017 || year == 2018 )
     {	   
	if(fabs(eta) < 1.0)        ea = 0.1566;
	else if(fabs(eta) < 1.479) ea = 0.1626;
	else if(fabs(eta) < 2.0)   ea = 0.1073;
	else if(fabs(eta) < 2.2)   ea = 0.0854;
	else if(fabs(eta) < 2.3)   ea = 0.1051;
	else if(fabs(eta) < 2.4)   ea = 0.1204;
	else                       ea = 0.1524;
     }   
   
   //Warning: EAs not computed for cone DR=0.4, use the values for DR=0.3 scaled by 16/9 instead
   //NB : dr=0.4 only used for PFRelIso04 ? other variables use dr=0.3 ?
   ea*=16./9.;
   
   return ea;
}
