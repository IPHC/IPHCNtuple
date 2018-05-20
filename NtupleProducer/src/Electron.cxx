#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(Electron)

Electron::Electron()
{
}

Electron::~Electron()
{
}

void Electron::read()
{
   _ID                                = idx;
   
   _E 	                            = ntP->el_E->at(idx);
   _pt	                            = ntP->el_pt->at(idx);
   _ptUnc                             = ntP->el_pt->at(idx);
   _eta	                            = ntP->el_eta->at(idx);
   _phi	                            = ntP->el_phi->at(idx);
   _m 	                            = ntP->el_m->at(idx);
   _conept	                    = ntP->el_conept->at(idx);
   _charge	                        = ntP->el_charge->at(idx);
   _id	                            = ntP->el_id->at(idx);
   
   _dxy	                            = ntP->el_gsfTrack_PV_dxy->at(idx);
   _dz	                            = ntP->el_gsfTrack_PV_dz->at(idx);
   _nlosthits                         = ntP->el_numberOfLostHits->at(idx);
   _ip3d	                            = ntP->el_ip3d->at(idx);
   _ip3dErr	                        = ntP->el_ip3dErr->at(idx);
   
   _isLoose	                        = ntP->el_looseCBId->at(idx);
   _isMedium                          = ntP->el_mediumCBId->at(idx);
   _passCV	                        = ntP->el_passConversionVeto->at(idx);
   _isPCC	                            = ntP->el_isGsfCtfScPixChargeConsistent->at(idx);
   _isGsfScPixChargeConsistent	            = ntP->el_isGsfScPixChargeConsistent->at(idx);

   _sip3d                          = (ntP->el_ip3dErr->at(idx)) ? ntP->el_ip3d->at(idx) / ntP->el_ip3dErr->at(idx) : -9999;
   _ooEmooP                        = ntP->el_ooEmooP->at(idx);
   
   _lepMVA	                        = ntP->el_lepMVA->at(idx);
   _mvaIso                      = ntP->el_mvaIso->at(idx);
   _mvaNoIso                      = ntP->el_mvaNoIso->at(idx);
   
   _lepMVA_miniRelIsoCharged          = ntP->el_lepMVA_miniRelIsoCharged->at(idx);
   _lepMVA_miniRelIsoNeutral          = ntP->el_lepMVA_miniRelIsoNeutral->at(idx);
   _lepMVA_jetPtRelv2                 = ntP->el_lepMVA_jetPtRelv2->at(idx);
   _lepMVA_jetPtRatio                 = ntP->el_lepMVA_jetPtRatio->at(idx);
   _lepMVA_jetBTagCSV                 = ntP->el_lepMVA_jetBTagCSV->at(idx);
   _lepMVA_sip3d                      = ntP->el_lepMVA_sip3d->at(idx);
   _lepMVA_dxy                        = ntP->el_lepMVA_dxy->at(idx);
   _lepMVA_dz                         = ntP->el_lepMVA_dz->at(idx);
   _lepMVA_mvaId                      = ntP->el_lepMVA_mvaId->at(idx);
   _lepMVA_eta                        = ntP->el_lepMVA_eta->at(idx);
   _lepMVA_jetNDauChargedMVASel       = ntP->el_lepMVA_jetNDauChargedMVASel->at(idx);
   
   _miniIso			                = ntP->el_miniIsoTTH->at(idx);
   
   _sigmaIetaIeta		                = ntP->el_full5x5_sigmaIetaIeta->at(idx);
   _superCluster_eta  	            = ntP->el_superCluster_eta->at(idx);
   _hadronicOverEm		            = ntP->el_hadronicOverEm->at(idx);
   _deltaEtaSuperClusterTrackAtVtx    = ntP->el_deltaEtaSuperClusterTrackAtVtx->at(idx);
   _deltaPhiSuperClusterTrackAtVtx    = ntP->el_deltaPhiSuperClusterTrackAtVtx->at(idx);
   _eSuperClusterOverP	            = ntP->el_eSuperClusterOverP->at(idx);
   _correctedEcalEnergy	            = ntP->el_correctedEcalEnergy->at(idx);
   _ecalEnergy	                    = ntP->el_ecalEnergy->at(idx);
   
   _trackMomentumError	            = ntP->el_trackMomentumError->at(idx);
   _tightCharge                       = ntP->el_isGsfCtfScPixChargeConsistent->at(idx) + ntP->el_isGsfScPixChargeConsistent->at(idx);
   
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

void Electron::init()
{   
   _E                              = -100.;
   _pt                             = -100.;
   _ptUnc                          = -100.;
   _eta                            = -100.;
   _phi                            = -100.;
   _m                              = -100.;
   _conept                         = -100.;
   _charge                         = 0;
   _id                             = 0;      
   
   _isLooseCBId           	        = 0;
   _isMediumCBId          	        = 0;
   _isLoose               	        = 0;
   _isMedium                  	    = 0;
   _isTight                   	    = 0;
   _isLooseMVA                	    = 0;
   _isTightMVA             	    = 0;
   _isLooseTTH            	        = 0;
   _isFakeableTTH             	    = 0;
   _isTightTTH                	    = 0;
   
   _dxy                            = -100;
   _dz                             = -100;
   _miniIso                        = -100;
   _isoR04                         = -100;
   _nlosthits                      = -100;
   _sip3d                          = -100;
   _ooEmooP                        = -100;
   
   _passCV                         = 0;
   _isPCC                          = 0;
   _isGsfScPixChargeConsistent     = 0;
   _passPtEta                      = 0;
   _ip3d                           = -100.;
   _ip3dErr                        = -100.;
   
    _lepMVA                         = -100;

   _lepMVA_miniRelIsoCharged       = -100.;
   _lepMVA_miniRelIsoNeutral       = -100.;
   _lepMVA_jetPtRelv2              = -100.;
   _lepMVA_jetPtRatio              = -100.;
   _lepMVA_jetBTagCSV              = -100.;
   _lepMVA_sip3d                   = -100.;
   _lepMVA_dxy                     = -100.;
   _lepMVA_dz                      = -100.;
   _lepMVA_mvaId                   = -100.;
   
   _lepMVA_eta                     = -100.;
   _lepMVA_jetNDauChargedMVASel    = -100.;
   
   _passChargeFlip                 = 0;
   _hasMatchedConversion           = 0;
   _isGsfCtfScPixChargeConsistent  = 0;
   
   _sigmaIetaIeta                  = -100.;
   _hadronicOverEm                 = -100.;
   _correctedEcalEnergy            = -100.;
   _ecalEnergy                     = -100.;
   _eSuperClusterOverP             = -100.;
   _deltaEtaSuperClusterTrackAtVtx = -100.;
   _deltaPhiSuperClusterTrackAtVtx = -100.;
   _see                            = -100.;
   _superCluster_eta               = -100.;
   
   _trackMomentumError             = -100;
   _tightCharge                    = -100;
   _cutEventSel                    = false;
   _noLostHits                     = false;
   _mvaIso                   = -100;
   _mvaNoIso                   = -100;

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

void Electron::sel()
{   
   bool isLoose = false;

   if( _pt <= 10. )
     {	
	if( fabs(_eta) < 0.8 ) isLoose = ( _mvaNoIso > -0.13285867293779202 );
	else if( fabs(_eta) < 1.479 ) isLoose = ( _mvaNoIso > -0.31765300958836074 );
	else isLoose = ( _mvaNoIso > -0.0799205914718861 );
     }
   else
     {
	if( fabs(_eta) < 0.8 ) isLoose = ( _mvaNoIso > -0.856871961305474 );
	else if( fabs(_eta) < 1.479 ) isLoose = ( _mvaNoIso > -0.8107642141584835 );
	else isLoose = ( _mvaNoIso > -0.7179265933023059 );
     }   

   _isoR04 = (_pt > 0.) ? (ntP->el_pfIso_sumChargedHadronPt->at(idx) + std::max( 0.0, ntP->el_pfIso_sumNeutralHadronEt->at(idx)+ntP->el_pfIso_sumPhotonEt->at(idx) - 0.5*ntP->el_pfIso_sumPUPt->at(idx) ))/_pt : -9999;
   
   bool pass_pt       = ( _pt > 7 );
   bool pass_eta      = ( fabs(_eta) < 2.5 );
   bool pass_dxy      = ( fabs(_dxy) < 0.05 );
   bool pass_dz       = ( fabs(_dz) < 0.1 );
   bool pass_miniIso  = ( _miniIso < 0.4 );
   bool pass_SIP      = ( fabs(_sip3d) < 8 );
   bool pass_isLoose  = ( isLoose );
   bool pass_losthits = ( _nlosthits < 2 );
   
   bool _passMuOverlap = 1;
   int nMuonLoose = nt->NtMuonLoose->size();
   for(int im=0;im<nMuonLoose;im++)
     {
        float dr = GetDeltaR(_eta,_phi,nt->NtMuonLoose->at(im).eta(),nt->NtMuonLoose->at(im).phi());
        if( dr < 0.05 ) _passMuOverlap = 0;
     }
   
   bool isLooseTTH = ( pass_pt          &&
		       pass_eta         &&
		       pass_dxy         &&
		       pass_dz          &&
		       pass_miniIso     &&
		       pass_SIP         &&
		       pass_isLoose     &&
		       pass_losthits    &&
		       _passMuOverlap   );
   
   _isLooseTTH = isLooseTTH;

   bool pass_fakeable_lepMVA = 0;
   bool pass_conept = ( _conept > 10. );
   bool pass_lepMVA = ( _lepMVA >= 0.90 );
   
   if( pass_lepMVA )
     if( _lepMVA_jetBTagCSV < 0.4941 ) pass_fakeable_lepMVA = 1;
   else
     if( _lepMVA_jetPtRatio > 0.6 && _lepMVA_jetBTagCSV < 0.07 && _lepMVA_mvaId > 0.5 ) pass_fakeable_lepMVA = 1;
   
   bool pass_sc = 1;
   float eInvMinusPInv  = (_ecalEnergy > 0) ? (1./_ecalEnergy-_eSuperClusterOverP/_ecalEnergy) : 99;
   if( _hadronicOverEm >= (0.10-0.03*(fabs(_superCluster_eta)>1.479)) ) pass_sc = 0;
   if( fabs(_deltaEtaSuperClusterTrackAtVtx) >= (0.01-0.002*(fabs(_superCluster_eta)>1.479)) ) pass_sc = 0;
   if( fabs(_deltaPhiSuperClusterTrackAtVtx) >= (0.04+0.03*(fabs(_superCluster_eta)>1.479)) ) pass_sc = 0;
   if( eInvMinusPInv <= -0.05 ) pass_sc = 0;
   if( eInvMinusPInv >= (0.01-0.005*(fabs(_superCluster_eta)>1.479)) ) pass_sc = 0;
   if( _sigmaIetaIeta >= (0.011+0.019*(fabs(_superCluster_eta)>1.479)) ) pass_sc = 0;
   
   _isFakeableTTH = ( _isLooseTTH &&
		      pass_sc &&
		      pass_fakeable_lepMVA &&
		      pass_conept );

   _isTightTTH = ( _isFakeableTTH &&
		   pass_lepMVA &&
		   _passCV );
}
