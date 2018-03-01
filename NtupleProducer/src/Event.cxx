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

    //std::cout << "id: "   << _id << std::endl;
    //std::cout << "run: "  << _run << std::endl;
    //std::cout << "lumi: " << _lumi << std::endl;

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

    // trigger
    
    /*
    for( int i = 0; i < ntP->trigger->size(); i++)
    {
        std::string currentpath ("Nopathsofar");
        if( ntP->trigger_pass->at(i) == 1) { currentpath = ntP->trigger_name->at(i); }

        // Single lepton
        std::size_t e1    = currentpath.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v"		           );
        std::size_t e2    = currentpath.find("HLT_Ele35_WPLoose_Gsf_v"			           );
        std::size_t e3    = currentpath.find("HLT_Ele27_WPTight_Gsf_v"                             );
        std::size_t e4    = currentpath.find("HLT_Ele25_eta2p1_WPTight_Gsf_v"                      );
        //std::size_t m1    = currentpath.find("HLT_IsoMu22_v"                                       ); //-- THQ2016
	std::size_t m1    = currentpath.find("HLT_IsoMu20_v"                                       );
        //std::size_t mTk1  = currentpath.find("HLT_IsoTkMu22_v"                                     ); //-- THQ2016
	std::size_t mTk1  = currentpath.find("HLT_IsoTkMu20_v"                                     );
        std::size_t m2    = currentpath.find("HLT_IsoMu22_eta2p1_v"                                );
        std::size_t mTk2  = currentpath.find("HLT_IsoTkMu22_eta2p1_v"                              );
        std::size_t m3    = currentpath.find("HLT_IsoMu24_v"                                       );
        std::size_t mTk3  = currentpath.find("HLT_IsoTkMu24_v"                                     );

        // Double lepton
        std::size_t ee    = currentpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"         );
        //std::size_t me    = currentpath.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"   ); //-- THQ2016
	std::size_t me    = currentpath.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"   );
        //std::size_t em    = currentpath.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"    ); //-- THQ2016
	std::size_t em    = currentpath.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"    );
        std::size_t medZ  = currentpath.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v" );
        std::size_t emdZ  = currentpath.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v" );
        std::size_t mm    = currentpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"               );
        std::size_t mmTk  = currentpath.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"             );

        // Triple lepton
        std::size_t eee   = currentpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"             );
        std::size_t mme   = currentpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"                   );
        std::size_t eem   = currentpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"                  );
        //std::size_t mmm   = currentpath.find("HLT_TripleMu_12_10_5_v"                                ); //--  THQ2016
	std::size_t mmm   = currentpath.find("HLT_TripleMu_5_3_3_v"                                );

        std::size_t ee_noDz   = currentpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"	   );
        std::size_t mm_noDz   = currentpath.find("HLT_Mu27_TrkIsoVVL_Mu8_TrkIsoVVL_v"	           );
        std::size_t mmTk_noDz = currentpath.find("HLT_Mu27_TrkIsoVVL_TkMu8_TrkIsoVVL_v" 	   ); 
        
        if(e1        != std::string::npos || e2   != std::string::npos || e3   != std::string::npos || e4 != std::string::npos )    { _TRIGe    = true ;}
        if(m1        != std::string::npos || m2   != std::string::npos || m2   != std::string::npos)                                { _TRIGm    = true ;}
        if(mTk1      != std::string::npos || mTk2 != std::string::npos || mTk3 != std::string::npos)                                { _TRIGmTk  = true ;}

        if(ee	     != std::string::npos)  { _TRIGee   = true ;}
        if(me        != std::string::npos || medZ != std::string::npos)  { _TRIGme   = true ;}
        if(em        != std::string::npos || emdZ != std::string::npos)  { _TRIGem   = true ;}
        if(mm	     != std::string::npos)  { _TRIGmm   = true ;}
        if(mmTk      != std::string::npos)  { _TRIGmmTk = true ;}
 
        if(eee       != std::string::npos)  { _TRIGeee  = true ;}
        if(mme       !=std::string::npos)   { _TRIGmme = true ;}
        if(eem       !=std::string::npos)   { _TRIGeem = true ;}
        if(mmm       !=std::string::npos)   { _TRIGmmm = true ;}

        if(ee_noDz   != std::string::npos)  { _TRIGee_noDz   = true ;}
        if(mm_noDz   != std::string::npos)  { _TRIGmm_noDz   = true ;}
        if(mmTk_noDz != std::string::npos)  { _TRIGmmTk_noDz = true ;}
    }
    */
    
    
    //  ### CHANGED -- TRIGGERS FROM 2016 THQ AN ###
    
    for( int i = 0; i < ntP->trigger->size(); i++)
    {
	    std::string currentpath ("Nopathsofar");
	    if( ntP->trigger_pass->at(i) == 1) { currentpath = ntP->trigger_name->at(i); }

	    // Single lepton
	    std::size_t e1    = currentpath.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v");
	    std::size_t e2    = currentpath.find("HLT_Ele27_WPTight_Gsf_v");
	    std::size_t e3    = currentpath.find("HLT_Ele25_eta2p1_WPTight_Gsf_v");
	    std::size_t m1    = currentpath.find("HLT_IsoMu22_v");
	    std::size_t m2    = currentpath.find("HLT_IsoMu22_eta2p1_v");
	    std::size_t m3    = currentpath.find("HLT_IsoMu24_v");
	    std::size_t mTk1    = currentpath.find("HLT_IsoTkMu22_v");
	    std::size_t mTk2    = currentpath.find("HLT_IsoTkMu22_eta2p1_v");
	    std::size_t mTk3    = currentpath.find("HLT_IsoTkMu24_v");

	    // Double lepton
	    std::size_t ee    = currentpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
	    std::size_t me1    = currentpath.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
	    std::size_t me2    = currentpath.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
	    std::size_t em1    = currentpath.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
	    std::size_t em2    = currentpath.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
	    std::size_t mm    = currentpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
	    std::size_t mmTk    = currentpath.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

	    // Triple lepton
	    std::size_t eee    = currentpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v");
	    std::size_t mme    = currentpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");
	    std::size_t eem    = currentpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v");
	    std::size_t mmm    = currentpath.find("HLT_TripleMu_12_10_5_v");


	    if(e1        != std::string::npos || e2   != std::string::npos || e3   != std::string::npos)    { _TRIGe    = true ;}
	    if(m1        != std::string::npos || m2   != std::string::npos || m2   != std::string::npos)    { _TRIGm    = true ;}
	    if(mTk1      != std::string::npos || mTk2 != std::string::npos || mTk3 != std::string::npos)    { _TRIGmTk  = true ;}

	    if(ee	     != std::string::npos)  { _TRIGee   = true ;}
	    if(me1        != std::string::npos || me2 != std::string::npos)  { _TRIGme   = true ;}
	    if(em1        != std::string::npos || em2 != std::string::npos)  { _TRIGem   = true ;}
	    if(mm	     != std::string::npos)  { _TRIGmm   = true ;}
	    if(mmTk      != std::string::npos)  { _TRIGmmTk = true ;}

	    if(eee       != std::string::npos)  { _TRIGeee  = true ;}
	    if(mme       !=std::string::npos)   { _TRIGmme = true ;}
	    if(eem       !=std::string::npos)   { _TRIGeem = true ;}
	    if(mmm       !=std::string::npos)   { _TRIGmmm = true ;}
    }

    // discriminant vs tt
    _disc_TT = 0;
}

void Event::init()
{
    _id                    = -1000;
    _run                   = -1000;
    _lumi                  = -1000;
    _rho                   = -1000;

    _pv_n                  = -1;
    _pv_z                  = -1000;
    _pv_zError             = -1000;

    _metpt                 = -1000;
    _metphi                = -1000;
    _metsumet              = -1000;
    _metcov00              = -1000;
    _metcov01              = -1000;
    _metcov10              = -1000;
    _metcov11              = -1000;

    _metNoHF_pt            = -1000;
    _metNoHF_phi           = -1000;
    _metNoHF_sumet         = -1000; 

    _weight_scale_muF0p5   = -1000;
    _weight_scale_muF2     = -1000;
    _weight_scale_muR0p5   = -1000;
    _weight_scale_muR2     = -1000;

    _mc_weight             = -1000;
    _mc_ptHat              = -1000;
    _mc_pu_trueNumInt      = -1000;

    _tth_channel           = -1000;

    _disc_TT               = -1000;

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
