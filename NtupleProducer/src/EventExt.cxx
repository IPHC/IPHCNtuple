#include "include/NtupleProducer.h"
#include <string.h>

ClassImp(EventExt)

EventExt::EventExt()
{
}

EventExt::~EventExt()
{
}

void EventExt::read(bool isdata)
{
   id               = ntP->ev_id;
   run              = ntP->ev_run;
   lumi             = ntP->ev_lumi;
   rho              = ntP->ev_rho;

   pv_n             = ntP->nvertex;
   pv_z             = ntP->pv_z;
   pv_zError        = ntP->pv_zError;

   metpt            = ntP->met_pt;
   metphi           = ntP->met_phi;
   metsumet         = ntP->met_sumet;
   metUncorrectedPt            = ntP->met_uncorrectedPt;
   metUncorrectedPhi           = ntP->met_uncorrectedPhi;
   metUncorrectedSumEt         = ntP->met_uncorrectedSumEt;
   metcov00	      = ntP->met_cov00;
   metcov01         = ntP->met_cov01;
   metcov10         = ntP->met_cov10;
   metcov11         = ntP->met_cov11;

   metNoHF_pt        = ntP->metNoHF_pt;
   metNoHF_phi       = ntP->metNoHF_phi;
   metNoHF_sumet     = ntP->metNoHF_sumet;
   
   if( !isdata )
     {
        weight_scale_muF0p5 = ntP->weight_scale_muF0p5;
        weight_scale_muF2   = ntP->weight_scale_muF2;
        weight_scale_muR0p5 = ntP->weight_scale_muR0p5;
        weight_scale_muR2   = ntP->weight_scale_muR2;

        mc_weight           = ntP->mc_weight;
        mc_ptHat            = ntP->mc_ptHat;
        mc_pu_trueNumInt    = ntP->mc_pu_trueNumInt;

        pdf_weights = *ntP->mc_pdfweights;   

        pdf_ids = *ntP->mc_pdfweightIds;    
     }

   for(int i=0;i<ntP->trigger->size();i++)
     {
        std::string tpath("");
	
        if( ntP->trigger_pass->at(i) == 1) 
	  {
	     tpath = ntP->trigger_name->at(i);

	     std::size_t e1    = tpath.find("HLT_Ele32_WPTight_Gsf_v");
	     std::size_t e2    = tpath.find("HLT_Ele35_WPTight_Gsf_v");
	     std::size_t e3    = tpath.find("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v");
	     std::size_t m1    = tpath.find("HLT_IsoMu24_v");
	     std::size_t m2    = tpath.find("HLT_IsoMu27_v");
	     std::size_t m3    = tpath.find("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v");

	     std::size_t ee    = tpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
	     std::size_t ee2   = tpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
	     std::size_t em    = tpath.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
	     std::size_t em2   = tpath.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
	     std::size_t em3   = tpath.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
	     std::size_t mm    = tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
	     std::size_t mm2   = tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");

	     std::size_t eee   = tpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v");
	     std::size_t eem   = tpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v");
	     std::size_t emm   = tpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v");
	     std::size_t mmm   = tpath.find("HLT_TripleMu_12_10_5_v");

	     trig_e = (e1!=std::string::npos || e2!=std::string::npos || e3!=std::string::npos);
	     trig_m = (m1!=std::string::npos || m2!=std::string::npos || m3!=std::string::npos);
	     
	     trig_ee = (ee!=std::string::npos || ee2!=std::string::npos);
	     trig_em = (em!=std::string::npos || em2!=std::string::npos || em3!=std::string::npos);
	     trig_mm = (mm!=std::string::npos || mm2!=std::string::npos);
	     
	     trig_eee = (eee!=std::string::npos);
	     trig_eem = (eee!=std::string::npos);
	     trig_emm = (eee!=std::string::npos);
	     trig_mmm = (eee!=std::string::npos);
	  }
     }   
}

void EventExt::init()
{
   id                    = -1000;
   run                   = -1000;
   lumi                  = -1000;
   rho                   = -1000;
   
   pv_n                  = -1;
   pv_z                  = -1000;
   pv_zError             = -1000;
   
   metpt                 = -1000;
   metphi                = -1000;
   metsumet              = -1000;
   metUncorrectedPt      = -1000;
   metUncorrectedPhi     = -1000;
   metUncorrectedSumEt   = -1000;
   metcov00              = -1000;
   metcov01              = -1000;
   metcov10              = -1000;
   metcov11              = -1000;
   
   metNoHF_pt            = -1000;
   metNoHF_phi           = -1000;
   metNoHF_sumet         = -1000; 
   
   weight_scale_muF0p5   = -1000;
   weight_scale_muF2     = -1000;
   weight_scale_muR0p5   = -1000;
   weight_scale_muR2     = -1000;
   
   mc_weight             = -1000;
   mc_ptHat              = -1000;
   mc_pu_trueNumInt      = -1000;
   
   tth_channel           = -1000;
   
   disc_TT               = -1000;
   
   pdf_weights.clear();
   pdf_ids.clear();
   
   trig_e = 0;
   trig_m = 0;
   
   trig_ee = 0;
   trig_em = 0;
   trig_mm = 0;
	     
   trig_eee = 0;
   trig_eem = 0;
   trig_emm = 0;
   trig_mmm = 0;
}
