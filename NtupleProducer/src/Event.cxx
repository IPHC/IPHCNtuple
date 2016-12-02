#include "include/NtupleProducer.h"
#include <string.h>

ClassImp(Event)

Event::Event()
{
}

Event::~Event()
{
}

void Event::read(bool isdata)
{
    // Event
    _id               = ntP->ev_id;
    _run              = ntP->ev_run;
    _lumi             = ntP->ev_lumi;
    _rho              = ntP->ev_rho;

    //std::cout << "id: " << _id << std::endl;

    // pv
    _pv_n             = ntP->nvertex;
    _pv_z             = ntP->pv_z;
    _pv_zError        = ntP->pv_zError;

    // MET
    _metpt            = ntP->met_pt;
    _metphi           = ntP->met_phi;
    _metsumet         = ntP->met_sumet;
    _metcov00	      = ntP->met_cov00;
    _metcov01         = ntP->met_cov01;
    _metcov10         = ntP->met_cov10;
    _metcov11         = ntP->met_cov11;

    _metNoHF_pt        = ntP->metNoHF_pt;
    _metNoHF_phi       = ntP->metNoHF_phi;
    _metNoHF_sumet     = ntP->metNoHF_sumet;


    if (!isdata)
    {

        _weight_scale_muF0p5 = ntP->weight_scale_muF0p5;
        _weight_scale_muF2   = ntP->weight_scale_muF2;
        _weight_scale_muR0p5 = ntP->weight_scale_muR0p5;
        _weight_scale_muR2   = ntP->weight_scale_muR2;

        _mc_weight           = ntP->mc_weight;
        _mc_ptHat            = ntP->mc_ptHat;
        _mc_pu_trueNumInt    = ntP->mc_pu_trueNumInt;

        _pdf_weights = *ntP->mc_pdfweights;   

        _pdf_ids = *ntP->mc_pdfweightIds;    

    }


    for( int i = 0; i < ntP->trigger->size(); i++)
    {
        std::string currentpath ("Nopathsofar");
        if( ntP->trigger_pass->at(i) == 1) { currentpath = ntP->trigger_name->at(i); }

        std::size_t eee   = currentpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"           );
        std::size_t me    = currentpath.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v" );
        std::size_t em    = currentpath.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"  );
        std::size_t ee    = currentpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"       );
        std::size_t mm    = currentpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"             );
        std::size_t mmTk  = currentpath.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"           );
        std::size_t m     = currentpath.find("HLT_IsoMu20_v"                                     );
        std::size_t mTk   = currentpath.find("HLT_IsoTkMu20_v"                                   );
        std::size_t e1    = currentpath.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v"		     );
        std::size_t e2    = currentpath.find("HLT_Ele35_WPLoose_Gsf_v"			     );
        std::size_t mme   = currentpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"                 );
        std::size_t eem   = currentpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"                );
        std::size_t mmm   = currentpath.find("HLT_TripleMu_5_3_3_v"                              );

        std::size_t ee_noDz   = currentpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"	  );
        std::size_t mm_noDz   = currentpath.find("HLT_Mu27_TrkIsoVVL_Mu8_TrkIsoVVL_v"	  );
        std::size_t mmTk_noDz = currentpath.find("HLT_Mu27_TrkIsoVVL_TkMu8_TrkIsoVVL_v" 	  ); 


        if(eee  != std::string::npos)  { _TRIGeee  = true ;}
        if(me	 != std::string::npos)  { _TRIGme   = true ;}
        if(em	 != std::string::npos)  { _TRIGem   = true ;}
        if(ee	 != std::string::npos)  { _TRIGee   = true ;}
        if(mm	 != std::string::npos)  { _TRIGmm   = true ;}
        if(mmTk != std::string::npos)  { _TRIGmmTk = true ;}
        if(m	 != std::string::npos)  { _TRIGm    = true ;}
        if(mTk  != std::string::npos)  { _TRIGmTk  = true ;}
        if(e1	 != std::string::npos || e2 != std::string::npos)  { _TRIGe   = true ;}
        if(mme  !=std::string::npos)   { _TRIGmme = true ;}
        if(eem  !=std::string::npos)   { _TRIGeem = true ;}
        if(mmm  !=std::string::npos)   { _TRIGmmm = true ;}
        if(ee_noDz   != std::string::npos) { _TRIGee_noDz   = true ;}
        if(mm_noDz   != std::string::npos) { _TRIGmm_noDz   = true ;}
        if(mmTk_noDz != std::string::npos) { _TRIGmmTk_noDz = true ;}
    }

    // discriminant vs tt
    _disc_TT = 0;
}

void Event::init()
{
    _id                    = -888;
    _run                   = -888;
    _lumi                  = -888;
    _rho                   = -888;

    _pv_n                  = -1;
    _pv_z                  = -888;
    _pv_zError             = -888;

    _metpt                 = -888;
    _metphi                = -888;
    _metsumet              = -888;
    _metcov00              = -888;
    _metcov01              = -888;
    _metcov10              = -888;
    _metcov11              = -888;

    _metNoHF_pt            = -888;
    _metNoHF_phi           = -888;
    _metNoHF_sumet         = -888; 

    _weight_scale_muF0p5   = -888;
    _weight_scale_muF2     = -888;
    _weight_scale_muR0p5   = -888;
    _weight_scale_muR2     = -888;

    _mc_weight             = -888;
    _mc_ptHat              = -888;
    _mc_pu_trueNumInt      = -888;

    _tth_channel           = -888;

    _disc_TT               = -888;

    _pdf_weights.clear();
    _pdf_ids.clear();

    _TRIGm    = false;
    _TRIGe    = false;
    _TRIGmTk  = false;  
    _TRIGee   = false;
    _TRIGmm   = false, 
              _TRIGme   = false, 
              _TRIGem   = false;
    _TRIGmmTk = false;
    _TRIGeee  = false;
    _TRIGmme  = false;
    _TRIGeem  = false;
    _TRIGmmm  = false;
}
