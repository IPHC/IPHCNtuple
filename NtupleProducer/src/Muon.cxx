#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(Muon)

Muon::Muon()
{
}

Muon::~Muon()
{
}

void Muon::read()
{
   _ID                             = idx;
   _E                              = ntP->mu_E->at(idx);
   _pt                             = ntP->mu_pt->at(idx);
   _ptUnc                          = ntP->mu_pt->at(idx);
   _eta                            = ntP->mu_eta->at(idx);
   _phi                            = ntP->mu_phi->at(idx);
   _m                              = ntP->mu_m->at(idx);
   _conept                         = ntP->mu_conept->at(idx);
   _charge                         = ntP->mu_charge->at(idx);
   _id                             = ntP->mu_id->at(idx);
   
   _isLoose                        = ntP->mu_isLooseMuon->at(idx);
   _isMedium	                    = ntP->mu_isMediumMuon->at(idx);
   _isTight	                    = ntP->mu_isTightMuon->at(idx);
   _isPFMuon                       = ntP->mu_isPFMuon->at(idx);

   _dxy                            = ntP->mu_innerTrack_PV_dxy->at(idx);
   _dz                             = ntP->mu_innerTrack_PV_dz->at(idx);
   _iso                            = ntP->mu_miniIsoTTH->at(idx);
   _sip3d                          = ntP->mu_ip3d->at(idx) / ntP->mu_ip3dErr->at(idx);
   _bestTrack_pt                   = ntP->mu_bestTrack_pt->at(idx);
   _bestTrack_ptError              = ntP->mu_bestTrack_ptError->at(idx);
   
   _lepMVA	                        = ntP->mu_lepMVA->at(idx);
   
   _lepMVA_miniRelIsoCharged       = ntP->mu_lepMVA_miniRelIsoCharged->at(idx);
   _lepMVA_miniRelIsoNeutral       = ntP->mu_lepMVA_miniRelIsoNeutral->at(idx);
   _lepMVA_jetPtRelv2              = ntP->mu_lepMVA_jetPtRelv2->at(idx);
   _lepMVA_jetPtRatio              = ntP->mu_lepMVA_jetPtRatio->at(idx);
   _lepMVA_jetBTagCSV              = ntP->mu_lepMVA_jetBTagCSV->at(idx);
   _lepMVA_sip3d                   = ntP->mu_lepMVA_sip3d->at(idx);
   _lepMVA_dxy                     = ntP->mu_lepMVA_dxy->at(idx);
   _lepMVA_dz                      = ntP->mu_lepMVA_dz->at(idx);
   _lepMVA_mvaId                   = ntP->mu_lepMVA_mvaId->at(idx);
   _lepMVA_eta                     = ntP->mu_lepMVA_eta->at(idx);
   _lepMVA_jetNDauChargedMVASel    = ntP->mu_lepMVA_jetNDauChargedMVASel->at(idx);
}

void Muon::init()
{
   _fakeType          = -100.;
   
   _E                 = -100.;
   _pt                = -100.;
   _ptUnc             = -100.;
   _eta               = -100.;
   _phi               = -100.;
   _m                 = -100.;
   _conept            = -100.;
   _charge            = 0;
   _id                = 0;
   
   _isLoose           = 0;
   _isMedium          = 0;
   _isTight           = 0;
   _isPFMuon          = 0;
   
   _isLooseTTH        = 0;
   _isFakeableTTH     = 0;
   _isTightTTH        = 0;
   
   _dxy                = -100.;
   _dz                 = -100.;
   _iso                = -100.;
   _isoR04             = -100.;
   _sip3d              = -100.;
   _bestTrack_pt       = -100.;
   _bestTrack_ptError  = -100.;

   _cutEventSel        = 1;
   _noLostHits         = 1;
   
   _lepMVA                      = -100.;
   
   _lepMVA_miniRelIsoCharged    = -100.;
   _lepMVA_miniRelIsoNeutral    = -100.;
   _lepMVA_jetPtRelv2           = -100.;
   _lepMVA_jetPtRatio           = -100.;
   _lepMVA_jetBTagCSV           = -100.;
   _lepMVA_sip3d                = -100.;
   _lepMVA_dxy                  = -100.;
   _lepMVA_dz                   = -100.;
   _lepMVA_mvaId                = -100.;
   _lepMVA_eta                  = -100.;
   _lepMVA_jetNDauChargedMVASel = -100.;   
}

void Muon::sel()
{
   bool pass_pt      = ( _pt > 5 );
   bool pass_eta     = ( fabs(_eta) < 2.4 );
   bool pass_dxy     = ( fabs(_dxy) < 0.05 );
   bool pass_dz      = ( fabs(_dz) < 0.1 );
   bool pass_miniIso = ( _iso < 0.4 );
   bool pass_SIP     = ( fabs(_sip3d) < 8 );
   bool pass_isLoose = ( _isLoose );

   _isoR04 = (_pt > 0.) ? (ntP->mu_pfIso04_sumChargedHadronPt->at(idx) + std::max( 0.0, ntP->mu_pfIso04_sumNeutralHadronEt->at(idx)+ntP->mu_pfIso04_sumPhotonEt->at(idx) - 0.5*ntP->mu_pfIso04_sumPUPt->at(idx) ))/_pt : -9999;
   
   bool isLooseTTH = ( pass_pt      &&
		       pass_eta     &&
		       pass_dxy     &&
		       pass_dz      &&
		       pass_miniIso &&
		       pass_SIP     &&
		       pass_isLoose );
   
   _isLooseTTH = isLooseTTH;
   
   bool pass_fakeable_lepMVA = 0;
   bool pass_lepMVA = ( _lepMVA >= 0.90 );
   
   if( pass_lepMVA )
     if( _lepMVA_jetBTagCSV < 0.4941 ) pass_fakeable_lepMVA = 1;
   else
     if( _lepMVA_jetPtRatio > 0.6 && _lepMVA_jetBTagCSV < 0.07 && _lepMVA_mvaId > 0.3 ) pass_fakeable_lepMVA = 1;
   
   bool pass_conept = ( _conept > 10. );
   
   _isFakeableTTH = ( _isLooseTTH &&
		      pass_fakeable_lepMVA &&
		      pass_conept );
   
   _isTightTTH = ( _isFakeableTTH && 
		   pass_lepMVA && 
		   _isMedium );
}
