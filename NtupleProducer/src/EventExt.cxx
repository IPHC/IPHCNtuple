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
//   metUncorrectedPt            = ntP->met_uncorrectedPt;
//   metUncorrectedPhi           = ntP->met_uncorrectedPhi;
//   metUncorrectedSumEt         = ntP->met_uncorrectedSumEt;
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

   bool e1_trig = 0;
   bool e2_trig = 0;
   bool e3_trig = 0;
   bool m1_trig = 0;
   bool m2_trig = 0;
   bool m3_trig = 0;

   bool ee_trig  = 0;
   bool ee2_trig = 0;
   bool em_trig  = 0;
   bool em2_trig = 0;
   bool em3_trig = 0;
   bool mm_trig  = 0;
   bool mm2_trig = 0;
   
   bool eee_trig = 0;
   bool eem_trig = 0;
   bool emm_trig = 0;
   bool mmm_trig = 0;
   
   for(int i=0;i<ntP->trigger->size();i++)
     {
        std::string tpath("");
	
        if( ntP->trigger_pass->at(i) == 1 )
	  {
	     tpath = ntP->trigger_name->at(i);

	     if( !e1_trig && (tpath.find("HLT_Ele32_WPTight_Gsf_v") != std::string::npos) ) {e1_trig=1;}
	     if( !e2_trig && (tpath.find("HLT_Ele35_WPTight_Gsf_v") != std::string::npos) ) {e2_trig=1;}
	     if( !e3_trig && (tpath.find("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v") != std::string::npos) ) {e3_trig=1;}
	     if( !m1_trig && (tpath.find("HLT_IsoMu24_v") != std::string::npos) ) {m1_trig=1;}
	     if( !m2_trig && (tpath.find("HLT_IsoMu27_v") != std::string::npos) ) {m2_trig=1;}
	     if( !m3_trig && (tpath.find("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v") != std::string::npos) ) {m3_trig=1;}

	     if( !ee_trig  && (tpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {ee_trig=1;};
	     if( !ee2_trig && (tpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {ee2_trig=1;};
	     if( !em_trig  && (tpath.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em_trig=1;};
	     if( !em2_trig && (tpath.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {em2_trig=1;};
	     if( !em3_trig && (tpath.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {em3_trig=1;};
	     if( !mm_trig  && (tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") != std::string::npos) ) {mm_trig=1;};
	     if( !mm2_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") != std::string::npos) ) {mm2_trig=1;};

	     if( !eee_trig && (tpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") != std::string::npos) ) {eee_trig=1;};
	     if( !eem_trig && (tpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") != std::string::npos) ) {eem_trig=1;};
	     if( !emm_trig && (tpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v") != std::string::npos) ) {emm_trig=1;};
	     if( !mmm_trig && (tpath.find("HLT_TripleMu_12_10_5_v") != std::string::npos) ) {mmm_trig=1;};
	  }
     }   

   trig_e = (e1_trig || e2_trig);
   trig_etau = (e3_trig);
   trig_m = (m1_trig || m2_trig);
   trig_mtau = (m3_trig);
   
   trig_ee = (ee_trig || ee2_trig);
   trig_em = (em_trig || em2_trig || em3_trig);
   trig_mm = (mm_trig || mm2_trig);
   
   trig_eee = (eee_trig);
   trig_eem = (eem_trig);
   trig_emm = (emm_trig);
   trig_mmm = (mmm_trig);
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
   metLD              	 = -1000;
   metsumet              = -1000;
//   metUncorrectedPt      = -1000;
//   metUncorrectedPhi     = -1000;
//   metUncorrectedSumEt   = -1000;
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
   trig_etau = 0;
   trig_m = 0;
   trig_mtau = 0;
   
   trig_ee = 0;
   trig_em = 0;
   trig_mm = 0;
	     
   trig_eee = 0;
   trig_eem = 0;
   trig_emm = 0;
   trig_mmm = 0;

   is_1l2tau             = 0;
   is_1l2tau_SR_Data     = 0;
   is_1l2tau_SR          = 0;
   is_1l2tau_Fake        = 0;
   is_2lSS               = 0;
   is_2lSS_SR_Data       = 0;
   is_2lSS_SR            = 0;
   is_2lSS_Fake          = 0;
   is_2lSS_Flip_Data     = 0;
   is_2lSS_Flip          = 0;
   is_2lSS1tau           = 0;
   is_2lSS1tau_SR_Data   = 0;
   is_2lSS1tau_SR        = 0;
   is_2lSS1tau_Fake      = 0;
   is_2lSS1tau_Flip_Data = 0;
   is_2lSS1tau_Flip      = 0;
   is_2l2tau             = 0;
   is_2l2tau_SR_Data     = 0;
   is_2l2tau_SR          = 0;
   is_2l2tau_Fake        = 0;
   is_3l                 = 0;
   is_3l_SR_Data         = 0;
   is_3l_SR              = 0;
   is_3l_Fake            = 0;
   is_3l1tau             = 0;
   is_3l1tau_SR_Data     = 0;
   is_3l1tau_SR          = 0;
   is_3l1tau_Fake        = 0;
   is_4l                 = 0;
   is_4l_SR_Data         = 0;
   is_4l_SR              = 0;
   is_4l_Fake            = 0;
   is_ttWctrl            = 0;
   is_ttWctrl_SR_Data    = 0;
   is_ttWctrl_SR         = 0;
   is_ttWctrl_Fake       = 0;
   is_ttWctrl_Flip_Data  = 0;
   is_ttWctrl_Flip       = 0;
   is_ttZctrl            = 0;
   is_ttZctrl_SR_Data    = 0;
   is_ttZctrl_SR         = 0;
   is_ttZctrl_Fake       = 0;
}
