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

    // general informations
    _E                              = ntP->mu_E->at(idx);
    _pt                             = ntP->mu_pt->at(idx);
    _ptUnc                          = ntP->mu_pt->at(idx);
    _eta                            = ntP->mu_eta->at(idx);
    _phi                            = ntP->mu_phi->at(idx);
    _m                              = ntP->mu_m->at(idx);
    _charge                         = ntP->mu_charge->at(idx);
    _id                             = ntP->mu_id->at(idx);

    //
    _isLoose                        = ntP->mu_isLooseMuon->at(idx);
    _isMedium	                    = ntP->mu_isMediumMuon->at(idx);
    _isTight	                    = ntP->mu_isTightMuon->at(idx);
    _isPFMuon                       = ntP->mu_isPFMuon->at(idx);

    // preselection variables
    _dxy                            = ntP->mu_innerTrack_PV_dxy->at(idx);
    _dz                             = ntP->mu_innerTrack_PV_dz->at(idx);
    _iso                            = ntP->mu_miniIsoTTH->at(idx);
    _sip3d                          = ntP->mu_ip3d->at(idx) / ntP->mu_ip3dErr->at(idx);
    _bestTrack_pt                   = ntP->mu_bestTrack_pt->at(idx);
    _bestTrack_ptError              = ntP->mu_bestTrack_ptError->at(idx);
    _tightCharge                    = ntP->mu_innerTrack_ptError->at(idx) / ntP->mu_innerTrack_pt->at(idx);

    // mva-based selection additionnal variables
    _lepMVA	                        = ntP->mu_lepMVA->at(idx);
    _lepMVA_TTH                     = ntP->mu_lepMVA_Moriond16->at(idx);

    _lepMVA_miniRelIsoCharged       = ntP->mu_lepMVA_miniRelIsoCharged->at(idx);
    _lepMVA_miniRelIsoNeutral       = ntP->mu_lepMVA_miniRelIsoNeutral->at(idx);
    _lepMVA_jetPtRelv2              = ntP->mu_lepMVA_jetPtRelv2->at(idx);
    //_lepMVA_jetDR                 = ntP->mu_lepMVA_jetDR->at(idx); // branch not available in FlatTree
    _lepMVA_jetPtRatio              = ntP->mu_lepMVA_jetPtRatio->at(idx);
    _lepMVA_jetBTagCSV              = ntP->mu_lepMVA_jetBTagCSV->at(idx);
    _lepMVA_sip3d                   = ntP->mu_lepMVA_sip3d->at(idx);
    _lepMVA_dxy                     = ntP->mu_lepMVA_dxy->at(idx);
    _lepMVA_dz                      = ntP->mu_lepMVA_dz->at(idx);
    _lepMVA_mvaId                   = ntP->mu_lepMVA_mvaId->at(idx);
    _lepMVA_eta                     = ntP->mu_lepMVA_eta->at(idx);
    _lepMVA_jetNDauChargedMVASel    = ntP->mu_lepMVA_jetNDauChargedMVASel->at(idx);

    // more variables

    //_dB3D  = ntP->mu_ip3dBS->at(idx);
    //_edB3D = ntP->mu_ip3dBSErr->at(idx);
}

void Muon::init()
{
    // general informations

    _fakeType          = -100.;

    _E                 = -100.;
    _pt                = -100.;
    _ptUnc             = -100.;
    _eta               = -100.;
    _phi               = -100.;
    _m                 = -100.;
    _charge            = 0;
    _id                = 0;
    _passPtEta         = 0;

    // Id
    _isLoose           = 0;
    _isMedium          = 0;
    _isTight           = 0;
    _isPFMuon          = 0;

    _isLooseTTH        = 0;
    _isFakeableTTH     = 0;
    _isTightTTH        = 0;

    _isLoose_tZq       = 0;
    _isTight_tZq       = 0;

    // variables for Id

    _dxy                = -100.;
    _dz                 = -100.;
    _iso                = -100.;
    _sip3d              = -100.;
    _bestTrack_pt       = -100.;
    _bestTrack_ptError  = -100.;
    //_dB3D               = -100.;
    //_edB3D              = -100.;
    _tightCharge        = 999;
    _cutEventSel        = true;

    // more variables

    _lepMVA            = -100.;
    _lepMVA_TTH  = -100.;


    _lepMVA_miniRelIsoCharged    = -100.;
    _lepMVA_miniRelIsoNeutral    = -100.;
    _lepMVA_jetPtRelv2           = -100.;
    //_lepMVA_jetDR                = -100.;
    _lepMVA_jetPtRatio           = -100.;
    _lepMVA_jetBTagCSV           = -100.;
    _lepMVA_sip3d                = -100.;
    _lepMVA_dxy                  = -100.;
    _lepMVA_dz                   = -100.;
    _lepMVA_mvaId                = -100.;
    _lepMVA_eta                  = -100.;
    _lepMVA_jetNDauChargedMVASel = -100.;

}

bool Muon::sel()
{
    // Loose
    bool pass_pt      = (_pt          >  5    );
    bool pass_eta     = (fabs(_eta)   <  2.4  );
    bool pass_dxy     = (fabs(_dxy)   <  0.05 );
    bool pass_dz      = (fabs(_dz)    <  0.1  );
    bool pass_miniIso = (_iso         <  0.4  );
    bool pass_SIP     = (fabs(_sip3d) <  8    );
    bool pass_isLoose = (_isLoose             );

    bool isLooseTTH = ( pass_pt      &&
                        pass_eta     &&
                        pass_dxy     &&
                        pass_dz      &&
                        pass_miniIso &&
                        pass_SIP     &&
                        pass_isLoose );

    _isLooseTTH = isLooseTTH;


    // Fakeable

    bool pass_lepMVA_TTH  = _lepMVA_TTH > 0.75 ;
    bool pass_lepMVA_jetBTagCSV089 = _lepMVA_jetBTagCSV < 0.89;

    bool pass_lepMVA_jetBtagCSVPtRatio = false;

    if (!pass_lepMVA_TTH && _lepMVA_jetPtRatio > 0.3 && _lepMVA_jetBTagCSV < 0.605) pass_lepMVA_jetBtagCSVPtRatio = true;
    if (pass_lepMVA_TTH && pass_lepMVA_jetBTagCSV089) pass_lepMVA_jetBtagCSVPtRatio = true;

    _isFakeableTTH = isLooseTTH && pass_lepMVA_jetBtagCSVPtRatio;

    // Tight

    _isTightTTH = isLooseTTH && pass_lepMVA_TTH && _isMedium && pass_lepMVA_jetBTagCSV089;

    _passTightCharge = (_tightCharge < 0.2);
    //_isTightTTH = isLooseTTH && pass_lepMVA_TTH && _isMedium && pass_lepMVA_jetBTagCSV089 && _pass_tightCharge;

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

    /*cout<<std::setiosflags(ios::fixed)<<setprecision(5);

    // synchronization printout
    if( true )//isLooseTTH )
    {
        std::cout                      << nt->NtEvent->at(0).id()                       << " "
                                       << _pt                                           << " "
                                       << _ptUnc                                        << " "
                                       << _eta                                          << " "
                                       << _phi                                          << " "
                                       << _E                                            << " "
                                       << ntP->mu_id->at(idx)                           << " "
                                       << ntP->mu_charge->at(idx)                       << " "
                                       << ntP->mu_lepMVA_jetNDauChargedMVASel->at(idx)  << " "
                                       << "miniIso"                                     << " "
                                       << ntP->mu_lepMVA_miniRelIsoCharged->at(idx)     << " "
                                       << ntP->mu_lepMVA_miniRelIsoNeutral->at(idx)     << " "
                                       << ntP->mu_lepMVA_jetPtRelv2->at(idx)            << " "
                                       << ntP->mu_lepMVA_jetBTagCSV->at(idx)            << " "
                                       << ntP->mu_lepMVA_jetPtRatio->at(idx)            << " "
                                       << fabs(_sip3d)                                  << " "
                                       << fabs(_dxy)                                    << " "
                                       << _dz                                           << " "
                                       << ntP->mu_segmentCompatibility->at(idx)         << " "
                                       << _lepMVA_TTH                                   << std::endl;
    }*/

    return isLooseTTH;
}


void Muon::sel_tZq()
{
    bool pass_pt = (_pt > 20.);
    bool pass_eta = (fabs(_eta) < 2.4);
    _passPtEta = (pass_pt && pass_eta);

    float sumChargedHadronPt = ntP->mu_pfIso04_sumChargedHadronPt->at(idx);
    float sumNeutralHadronEt = ntP->mu_pfIso04_sumNeutralHadronEt->at(idx);
    float sumPhotonEt = ntP->mu_pfIso04_sumPhotonEt->at(idx);
    float sumPUPt = ntP->mu_pfIso04_sumPUPt->at(idx);
    float effArea = effectiveArea(30,_eta);
    _iso = (sumChargedHadronPt + std::max(float(sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt),float(0.)))/_pt;
    //_iso = (sumChargedHadronPt + std::max(sumNeutralHadronEt+sumPhotonEt-0.5*effArea,float(0.)))/_pt;

    _isLoose_tZq =  (_isPFMuon &&
                    (ntP->mu_isGlobalMuon->at(idx) ||
                    ntP->mu_isTrackerMuon->at(idx))
                    );

    _isTight_tZq =  (
                    ntP->mu_isGlobalMuon->at(idx)                                    &&
                    _isPFMuon                                       		     &&
    	            ntP->mu_globalTrack_normalizedChi2->at(idx)              < 10    &&
    	            ntP->mu_globalTrack_numberOfValidMuonHits->at(idx)	     > 0     &&
                    ntP->mu_numberOfMatchedStations->at(idx)                 > 1     &&
                    fabs(ntP->mu_bestTrack_PV_dxy->at(idx) )		     < 0.2   &&
                    fabs(ntP->mu_bestTrack_PV_dz->at(idx) )                  < 0.5   &&
    	            ntP->mu_innerTrack_numberOfValidPixelHits->at(idx)	     > 0     &&
                    ntP->mu_innerTrack_trackerLayersWithMeasurement->at(idx) > 5
                    );


    for(int id=0;id<evdebug->size();id++)
    {
        if( nt->NtEvent->at(0).id() == evdebug->at(id) )
        {
            std::cout << "Muon #"                   << _ID                  << std::endl;
            std::cout << "   pt="                   << _pt  << " eta=" << _eta << " phi=" << _phi << std::endl;
            std::cout << "   passPtEta="            << _passPtEta           << std::endl;
            std::cout << "   isLoose="              << _isLoose             << std::endl;
            std::cout << "   isTight="              << _isTight             << std::endl;
            std::cout << "   iso="                  << _iso                 << std::endl;
            std::cout << "   effArea="              << effArea              << std::endl;
        }
    }
    //cout<<nt->NtEvent->at(0).id()<<" sumChargedHadronPt="<<sumChargedHadronPt<<", sumNeutralHadronEt="<<sumNeutralHadronEt<<", sumPhotonEt="<<sumPhotonEt<<", sumPUPt="<<sumPUPt<<endl;

}

float Muon::effectiveArea(int dr,float eta)
{
    float effArea = 0.;

    if( fabs(eta) >= 0 && fabs(eta) < 0.8 )
    {
        if( dr == 30 ) effArea = 0.0913;
        else if( dr == 40 ) effArea = 0.1564;
    }
    else if( fabs(eta) >= 0.8 && fabs(eta) < 1.3 )
    {
        if( dr == 30 ) effArea = 0.0765;
        else if( dr == 40 ) effArea = 0.1325;
    }
    else if( fabs(eta) >= 1.3 && fabs(eta) < 2.0 )
    {
        if( dr == 30 ) effArea = 0.0546;
        else if( dr == 40 ) effArea = 0.0913;
    }
    else if( fabs(eta) >= 2.0 && fabs(eta) < 2.2 )
    {
        if( dr == 30 ) effArea = 0.0728;
        else if( dr == 40 ) effArea = 0.1212;
    }
    else if( fabs(eta) >= 2.2 && fabs(eta) < 2.5 )
    {
        if( dr == 30 ) effArea = 0.1177;
        else if( dr == 40 ) effArea = 0.2085;
    }

    return effArea;
}
