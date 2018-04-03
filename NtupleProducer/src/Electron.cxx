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

     _lepMVA	                        = ntP->el_lepMVA->at(idx);
     _lepMVA_TTH                        = ntP->el_lepMVA_Moriond16->at(idx);
     _mvaNonTrigV0                      = ntP->el_mvaNonTrigV0->at(idx);

     _lepMVA_miniRelIsoCharged          = ntP->el_lepMVA_miniRelIsoCharged->at(idx);
     _lepMVA_miniRelIsoNeutral          = ntP->el_lepMVA_miniRelIsoNeutral->at(idx);
     _lepMVA_jetPtRelv2                 = ntP->el_lepMVA_jetPtRelv2->at(idx);
     //_lepMVA_jetDR                    = ntP->el_lepMVA_jetDR->at(idx);
     _lepMVA_jetPtRatio                 = ntP->el_lepMVA_jetPtRatio->at(idx);
     _lepMVA_jetBTagCSV                 = ntP->el_lepMVA_jetBTagCSV->at(idx);
     _lepMVA_sip3d                      = ntP->el_lepMVA_sip3d->at(idx);
     _lepMVA_dxy                        = ntP->el_lepMVA_dxy->at(idx);
     _lepMVA_dz                         = ntP->el_lepMVA_dz->at(idx);
     _lepMVA_mvaId                      = ntP->el_lepMVA_mvaId->at(idx);
     _lepMVA_eta                        = ntP->el_lepMVA_eta->at(idx);
     _lepMVA_jetNDauChargedMVASel       = ntP->el_lepMVA_jetNDauChargedMVASel->at(idx);

     _miniIso			                = ntP->el_miniIsoTTH->at(idx);

     _sigmaIetaIeta		                = ntP->el_sigmaIetaIeta->at(idx);
     _superCluster_eta  	            = ntP->el_superCluster_eta->at(idx);
     _hadronicOverEm		            = ntP->el_hadronicOverEm->at(idx);
     _deltaEtaSuperClusterTrackAtVtx    = ntP->el_deltaEtaSuperClusterTrackAtVtx->at(idx);
     _deltaPhiSuperClusterTrackAtVtx    = ntP->el_deltaPhiSuperClusterTrackAtVtx->at(idx);
     _eSuperClusterOverP	            = ntP->el_eSuperClusterOverP->at(idx);
     _correctedEcalEnergy	            = ntP->el_correctedEcalEnergy->at(idx);
     _ecalEnergy	                    = ntP->el_ecalEnergy->at(idx);

     _trackMomentumError	            = ntP->el_trackMomentumError->at(idx);
     _tightCharge                       = ntP->el_isGsfCtfScPixChargeConsistent->at(idx) + ntP->el_isGsfScPixChargeConsistent->at(idx);
}

void Electron::init()
{

    // general informations
    _E                              = -100.;
    _pt                             = -100.;
    _ptUnc                          = -100.;
    _eta                            = -100.;
    _phi                            = -100.;
    _m                              = -100.;
    _charge                         = 0;
    _id                             = 0;

    // Id
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


    // variables for Id
    _dxy                            = -100;
    _dz                             = -100;
    _miniIso                        = -100;
    _nlosthits                      = -100;

    _passCV                         = 0;
    _isPCC                          = 0;
    _passPtEta                      = 0;
    _ip3d                           = -100.;
    _ip3dErr                        = -100.;

    // more variables
    _lepMVA                         = -100;
    _lepMVA_TTH                     = -100;


    _lepMVA_miniRelIsoCharged       = -100.;
    _lepMVA_miniRelIsoNeutral       = -100.;
    _lepMVA_jetPtRelv2              = -100.;
    //_lepMVA_jetDR                   = -100.;
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
    _mvaNonTrigV0                   = -100;
}








// ######## ######## ##     ##
//    ##       ##    ##     ##
//    ##       ##    ##     ##
//    ##       ##    #########
//    ##       ##    ##     ##
//    ##       ##    ##     ##
//    ##       ##    ##     ##


bool Electron::sel()
{

    float SIP     = fabs(_ip3d/_ip3dErr);

    //MVA ID
    bool  isLoose = false;

    if      (fabs(_eta) < 0.8  ) { isLoose = ( _mvaNonTrigV0 > -0.70 ); }
    else if (fabs(_eta) < 1.479) { isLoose = ( _mvaNonTrigV0 > -0.83 ); }
    else                         { isLoose = ( _mvaNonTrigV0 > -0.92 ); }


    // Loose
    bool pass_pt       = (_pt        > 7   );
    bool pass_eta      = (fabs(_eta) < 2.5 );
    bool pass_dxy      = (fabs(_dxy) < 0.05);
    bool pass_dz       = (fabs(_dz)  < 0.1 );
    bool pass_miniIso  = (_miniIso    < 0.4 );
    bool pass_SIP      = (SIP        < 8   );
    bool pass_isLoose  = (isLoose          );
    bool pass_losthits = (_nlosthits < 2   );

    //
    bool _passMuOverlap = 1;
    int nMuon = nt->NtMuon->size();
    for(int im=0;im<nMuon;im++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtMuon->at(im).eta(),nt->NtMuon->at(im).phi());
        if( dr < 0.05 ) _passMuOverlap = 0; //&& nt->NtMuon->at(im).isLoose() ) _passMuOverlap = 0;
    }

    bool isLooseTTH     = ( pass_pt          &&
                            pass_eta         &&
                            pass_dxy         &&
                            pass_dz          &&
                            pass_miniIso     &&
                            pass_SIP         &&
                            pass_isLoose     &&
                            pass_losthits    &&
                            _passMuOverlap   );

    _isLooseTTH = isLooseTTH; // OK

    // preselection electron for pt > 30 (aiming at closure for the fake rate method whatever)
    // from https://github.com/peruzzim/cmg-cmssw/blob/works_260116/CMGTools/TTHAnalysis/python/tools/emulateElectronTriggerCuts.py#L1-L8
    // def _susy2lss_idEmu_cuts(lep):
    //if (abs(lep.pdgId())!=11): return True OK
    //if (lep.full5x5_sigmaIetaIeta()>=(0.011 if abs(lep.superCluster().eta())<1.479 else 0.031)): return False
    //if (lep.hadronicOverEm()>=0.08): return False
    //if (abs(lep.deltaEtaSuperClusterTrackAtVtx())>=0.01): return False
    //if (abs(lep.deltaPhiSuperClusterTrackAtVtx())>=(0.04 if abs(lep.superCluster().eta())<1.479 else 0.08)): return False
    //if (abs((1.0/lep.ecalEnergy() - lep.eSuperClusterOverP()/lep.ecalEnergy()) if lep.ecalEnergy()>0. else 9e9)>=0.01): return False
    //return True

    bool _passConditions = true;
    if ( _pt > 30 )
    {

        _eInvMinusPinv  = ( _ecalEnergy > 0 ) ?  (1. / _ecalEnergy - _eSuperClusterOverP / _ecalEnergy) : 99;
        //std::cout << _eInvMinusPinv << " ";
        if ( ntP->el_hadronicOverEm->at(idx)        >= ( 0.10 - 0.03 * ( abs(_superCluster_eta) > 1.479 ) ) )   {_passConditions = false;}// std::cout << "H/E nope ";}
        if ( fabs(_deltaEtaSuperClusterTrackAtVtx)  >= ( 0.01 - 0.002 * ( abs(_superCluster_eta) > 1.479 ) ) )  {_passConditions = false;}// std::cout << "Eta nope ";}
        if ( fabs(_deltaPhiSuperClusterTrackAtVtx)  >= ( 0.04 + 0.03 * ( abs(_superCluster_eta) > 1.479  ) ) )  {_passConditions = false;}// std::cout << "Phi nope ";}
        if ( _eInvMinusPinv                          <= ( -0.05) )                                               {_passConditions = false;}// std::cout << "1/. inf nope ";}
        if ( _eInvMinusPinv                          >= ( 0.01 - 0.005 * ( abs(_superCluster_eta) > 1.479  ) ) ) {_passConditions = false;}// std::cout << "1/. sup nope ";}
        if ( _sigmaIetaIeta                         >= ( 0.011 + 0.019 * ( abs(_superCluster_eta) > 1.479 ) ) ) {_passConditions = false;}// std::cout << "See nope ";}
    }

    // Fakeable

    pass_losthits = (_nlosthits == 0 );
    pass_pt       = (_pt        >  10); // should be 0.85 * pt(jet) for fakeable object cf v4 of note

    bool pass_lepMVA_TTH  = _lepMVA_TTH > 0.90 ;

    bool pass_lepMVA_jetBTagCSV084 = _lepMVA_jetBTagCSV < 0.8484;


    bool pass_lepMVA_jetBtagCSVPtRatio = false;

    if (!pass_lepMVA_TTH && _lepMVA_jetPtRatio > 0.5 && _lepMVA_jetBTagCSV < 0.3) pass_lepMVA_jetBtagCSVPtRatio = true;
    if ( pass_lepMVA_TTH && pass_lepMVA_jetBTagCSV084) pass_lepMVA_jetBtagCSVPtRatio = true;

    bool isFakeableTTH     = ( pass_pt          &&
                               pass_eta         &&
                               pass_dxy         &&
                               pass_dz          &&
                               pass_miniIso     &&
                               pass_SIP         &&
                               pass_isLoose     &&
                               //pass_losthits    &&
                               _passConditions &&
                               _passMuOverlap   &&
            			       pass_lepMVA_jetBtagCSVPtRatio
			                 );

    _isFakeableTTH = isFakeableTTH;

    // Tight

    bool pass_CV            = (_passCV         );
    _passTightCharge        = (_tightCharge  >1);
    _noLostHits             = pass_losthits;
    _cutEventSel            = _passCV;

    bool isTightTTH     = ( pass_pt               &&
                            pass_eta              &&
                            pass_dxy              &&
                            pass_dz               &&
                            pass_miniIso          &&
                            pass_SIP              &&
                            pass_isLoose          &&
                            _passConditions      &&
                            _passMuOverlap        &&
			                pass_lepMVA_TTH       &&
			                pass_lepMVA_jetBTagCSV084
                            );

    _isTightTTH = isTightTTH;

    if(_isFakeableTTH && !_isTightTTH)
    {
        float dr_min = 0.5, new_pt = -100;
        int n_jets = ntP->jet_pt->size();
        for(int ij=0;ij<n_jets;ij++)
        {
            float dr = GetDeltaR(_eta,_phi,ntP->jet_eta->at(ij),ntP->jet_phi->at(ij));
            if( dr < dr_min ) new_pt = ntP->jet_pt->at(ij) * 0.85;
            //std::cout << "jet[" << ij << "]  dr: " << dr << "  pt: " << ntP->jet_pt->at(ij) << std::endl;
        }
        _pt = new_pt;
    }

    cout<<std::setiosflags(ios::fixed)<<setprecision(5);

    // synchronization printout
    if( false )
    //if( isLooseTTH )
    {
        std::cout   //<< "Electrons: "
                    << nt->NtEvent->at(0).id()                          << " "
                    << _pt                                              << " "
                    << _ptUnc                                           << " "
                    << _eta                                             << " "
                    << _phi                                             << " "
                    << _E                                               << " "
                    << ntP->el_id->at(idx)                              << " "
                    << ntP->el_charge->at(idx)                          << " "
                    << ntP->el_miniIsoTTH->at(idx)                      << " "
                    << ntP->el_lepMVA_miniRelIsoCharged->at(idx)        << " "
                    << ntP->el_lepMVA_miniRelIsoNeutral->at(idx)        << " "
                    << ntP->el_lepMVA_jetPtRelv2->at(idx)               << " "
                    //<< 0.0                                              << " "
                    << ntP->el_lepMVA_jetBTagCSV->at(idx)               << " "
                    << ntP->el_lepMVA_jetPtRatio->at(idx)               << " "
                    << SIP                                              << " "
                    << fabs(_dxy)                                       << " "
                    << fabs(_dz)                                        << " "
                    << fabs( ntP->el_mvaNonTrigV0->at(idx) )            << " "
                    << _lepMVA_TTH                                      << " "
                    << std::endl;


    }

    return isLooseTTH;
}














// ######## ##     ##  #######      #######    #####      ##    #######
//    ##    ##     ## ##     ##    ##     ##  ##   ##   ####   ##     ##
//    ##    ##     ## ##     ##           ## ##     ##    ##   ##
//    ##    ######### ##     ##     #######  ##     ##    ##   ########
//    ##    ##     ## ##  ## ##    ##        ##     ##    ##   ##     ##
//    ##    ##     ## ##    ##     ##         ##   ##     ##   ##     ##
//    ##    ##     ##  ##### ##    #########   #####    ######  #######

/*
bool Electron::sel()
{
    float SIP     = fabs(_ip3d/_ip3dErr);

    //MVA ID
    bool  pass_MVAid = false;
    bool  _passConditions = false;
    _eInvMinusPinv  = ( _ecalEnergy > 0 ) ?  (1. / _ecalEnergy - _eSuperClusterOverP / _ecalEnergy) : 99;

    if(fabs(_eta) < 0.8)
    {
        if(_lepMVA_mvaId > 0) {pass_MVAid = true;}

        if(_sigmaIetaIeta < 0.011 && ntP->el_hadronicOverEm->at(idx) < 0.10 && fabs(_deltaEtaSuperClusterTrackAtVtx) < 0.01 && fabs(_deltaPhiSuperClusterTrackAtVtx) < 0.04 && _eInvMinusPinv > -0.05 && _eInvMinusPinv < 0.010) {_passConditions = true;}
    }
    else if(fabs(_eta) < 1.479)
    {
        if(_lepMVA_mvaId > 0) {pass_MVAid = true;}

        if(_sigmaIetaIeta < 0.011 && ntP->el_hadronicOverEm->at(idx) < 0.10 && fabs(_deltaEtaSuperClusterTrackAtVtx) < 0.01 && fabs(_deltaPhiSuperClusterTrackAtVtx) < 0.04 && _eInvMinusPinv > -0.05 && _eInvMinusPinv < 0.010) {_passConditions = true;}
    }
    else
    {
        if(_lepMVA_mvaId > 0.7) {pass_MVAid = true;}

        if(_sigmaIetaIeta < 0.030 && ntP->el_hadronicOverEm->at(idx) < 0.07 && fabs(_deltaEtaSuperClusterTrackAtVtx) < 0.008 && fabs(_deltaPhiSuperClusterTrackAtVtx) < 0.07 && _eInvMinusPinv > -0.05 && _eInvMinusPinv < 0.005) {_passConditions = true;}
    }



    // Loose
    bool pass_eta      = (fabs(_eta) < 2.5 );
    bool pass_pt       = (_pt        > 7   );
    bool pass_dxy      = (fabs(_dxy) < 0.05);
    bool pass_dz       = (fabs(_dz)  < 0.1 );
    bool pass_SIP      = (SIP        < 8   );
    bool pass_miniIso  = (_miniIso    < 0.4 );
    bool pass_losthits = (_nlosthits < 2   );

    bool _passMuOverlap = 1;
    int nMuon = nt->NtMuon->size();
    for(int im=0;im<nMuon;im++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtMuon->at(im).eta(),nt->NtMuon->at(im).phi());
        if( dr < 0.05 ) _passMuOverlap = 0; //&& nt->NtMuon->at(im).isLoose() ) _passMuOverlap = 0;
    }

    bool isLooseTTH     = ( pass_pt          &&
                            pass_eta         &&
                            pass_dxy         &&
                            pass_dz          &&
                            pass_SIP         &&
							pass_miniIso     &&
                            pass_losthits    &&
                            pass_MVAid   );

    _isLooseTTH = isLooseTTH;


    // Fakeable
    // pass_pt = (_pt > 15); // should be 0.85 * pt(jet) for fakeable object cf v4 of note
	pass_losthits = (_nlosthits == 0 );
    bool pass_lepMVA_TTH  = _lepMVA_TTH > 0.90 ;
    bool pass_lepMVA_jetBTagCSV084 = _lepMVA_jetBTagCSV < 0.8484;

	bool pass_lepMVA_jetBtagCSVPtRatio = false;
    if (!pass_lepMVA_TTH && _lepMVA_jetPtRatio > 0.5 && _lepMVA_jetBTagCSV < 0.3) pass_lepMVA_jetBtagCSVPtRatio = true;
    if ( pass_lepMVA_TTH && pass_lepMVA_jetBTagCSV084) pass_lepMVA_jetBtagCSVPtRatio = true;

    bool isFakeableTTH     = ( pass_pt          &&
                               pass_eta         &&
                               pass_dxy         &&
                               pass_dz          &&
                               pass_SIP         &&
							   pass_miniIso     &&
                               _passConditions    &&
                               pass_lepMVA_jetBtagCSVPtRatio &&
                               pass_losthits
			                 );

    _isFakeableTTH = isFakeableTTH;

    // Tight
    _passTightCharge        = (_tightCharge  >1);
    _cutEventSel            = _passCV;

    bool isTightTTH     = ( pass_pt               &&
                            pass_eta              &&
                            pass_dxy              &&
                            pass_dz               &&
                            pass_miniIso          &&
                            pass_SIP              &&
                            _passConditions      &&
							pass_losthits        &&
                            pass_lepMVA_jetBTagCSV084 &&
                            _passTightCharge    &&
                            _passCV         &&
			                pass_lepMVA_TTH
                            );

    _isTightTTH = isTightTTH;


    if(_isFakeableTTH && !_isTightTTH)
    {
        float dr_min = 0.5, new_pt = -100;
        int n_jets = ntP->jet_pt->size();
        for(int ij=0;ij<n_jets;ij++)
        {
            float dr = GetDeltaR(_eta,_phi,ntP->jet_eta->at(ij),ntP->jet_phi->at(ij));
            if( dr < dr_min ) new_pt = ntP->jet_pt->at(ij) * 0.85;
            //std::cout << "jet[" << ij << "]  dr: " << dr << "  pt: " << ntP->jet_pt->at(ij) << std::endl;
        }
        _pt = new_pt;
    }


    return isLooseTTH;
}
*/
