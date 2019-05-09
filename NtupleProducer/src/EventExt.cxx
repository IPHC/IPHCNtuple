#include "include/NtupleProducer.h"
#include <string.h>

ClassImp(EventExt)

EventExt::EventExt()
{
}

EventExt::~EventExt()
{
}

void EventExt::read(bool isdata,int year)
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
     	/*
        weight_originalXWGTUP = ntP->weight_originalXWGTUP;
        weight_scale_muF0p5 = ntP->weight_scale_muF0p5;
        weight_scale_muF2   = ntP->weight_scale_muF2;
        weight_scale_muR0p5 = ntP->weight_scale_muR0p5;
        weight_scale_muR2   = ntP->weight_scale_muR2;
        weight_scale_muR2muF2   = ntP->weight_scale_muR2muF2;
        weight_scale_muR0p5muF0p5   = ntP->weight_scale_muR0p5muF0p5;*/
	
	weight_scale_index2 = ntP->weight_scale_index2;
	weight_scale_index3 = ntP->weight_scale_index3;
	weight_scale_index4 = ntP->weight_scale_index4;
	weight_scale_index5 = ntP->weight_scale_index5;
	weight_scale_index6 = ntP->weight_scale_index6;
	weight_scale_index7 = ntP->weight_scale_index7;
	weight_scale_index8 = ntP->weight_scale_index8;
	weight_scale_index9 = ntP->weight_scale_index9;

        mc_weight           = ntP->mc_weight;
	mc_weight_originalValue = ntP->mc_weight_originalValue;
        mc_ptHat            = ntP->mc_ptHat;
        mc_pu_trueNumInt    = ntP->mc_pu_trueNumInt;
	//if(mc_pu_trueNumInt <=0) {std::cout<<"mc_pu_trueNumInt = "<<mc_pu_trueNumInt<<std::endl;} //FIXME

        pdf_weights = *ntP->mc_pdfweights;   
        pdf_ids = *ntP->mc_pdfweightIds;    
     }

   bool e1_trig = 0;
   bool e2_trig = 0;
   bool e3_trig = 0;
   bool m1_trig = 0;
   bool m2_trig = 0;
   bool m3_trig = 0;
   bool m4_trig = 0;
   bool m5_trig = 0;
   bool m6_trig = 0;

   bool ee1_trig = 0;
   bool ee2_trig = 0;
   bool em1_trig = 0;
   bool em2_trig = 0;
   bool em3_trig = 0;
   bool em4_trig = 0;
   bool em5_trig = 0;
   bool em6_trig = 0;
   bool em7_trig = 0;
   bool mm1_trig = 0;
   bool mm2_trig = 0;
   bool mm3_trig = 0;
   bool mm4_trig = 0;
   bool mm5_trig = 0;
   bool mm6_trig = 0;
   
   bool et1_trig = 0;
   bool et2_trig = 0;
   bool et3_trig = 0;
   bool mt1_trig = 0;
   bool mt2_trig = 0;
   
   bool tt1_trig = 0;
   bool tt2_trig = 0;
   bool tt3_trig = 0;
   bool tt4_trig = 0;
   
   bool eee1_trig = 0;
   bool eem1_trig = 0;
   bool eem2_trig = 0;
   bool emm1_trig = 0;
   bool mmm1_trig = 0;
   
   for(int i=0;i<ntP->trigger->size();i++)
     {
        std::string tpath("");
	
        if( ntP->trigger_pass->at(i) == 1 )
	  {
	     tpath = ntP->trigger_name->at(i);

	     if( year == 2016 )
	       {
		  if( !e1_trig && (tpath.find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos) ) {e1_trig=1;}
		  if( !e2_trig && (tpath.find("HLT_Ele25_eta2p1_WPTight_Gsf_v") != std::string::npos) ) {e2_trig=1;}
		  if( !e3_trig && (tpath.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v") != std::string::npos) ) {e3_trig=1;}
		  if( !m1_trig && (tpath.find("HLT_IsoMu22_v") != std::string::npos) ) {m1_trig=1;}
		  if( !m2_trig && (tpath.find("HLT_IsoTkMu22_v") != std::string::npos) ) {m2_trig=1;}
		  if( !m3_trig && (tpath.find("HLT_IsoMu22_eta2p1_v") != std::string::npos) ) {m3_trig=1;}
		  if( !m4_trig && (tpath.find("HLT_IsoTkMu22_eta2p1_v") != std::string::npos) ) {m4_trig=1;}
		  if( !m5_trig && (tpath.find("HLT_IsoMu24_v") != std::string::npos) ) {m5_trig=1;}
		  if( !m6_trig && (tpath.find("HLT_IsoTkMu24_v") != std::string::npos) ) {m6_trig=1;}
		  
		  if( !ee1_trig && (tpath.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {ee1_trig=1;};
		  if( !ee2_trig && (tpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {ee2_trig=1;};
		  if( !em1_trig && (tpath.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em1_trig=1;};
		  if( !em2_trig && (tpath.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {em2_trig=1;};
		  if( !em3_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em3_trig=1;};
		  if( !em4_trig && (tpath.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {em4_trig=1;};
		  if( !em5_trig && (tpath.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em5_trig=1;};
		  if( !em6_trig && (tpath.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em6_trig=1;};
		  if( !em7_trig && (tpath.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em7_trig=1;};
		  if( !mm1_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") != std::string::npos) ) {mm1_trig=1;};
		  if( !mm2_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos) ) {mm2_trig=1;};
		  if( !mm3_trig && (tpath.find("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos) ) {mm3_trig=1;};
		  if( !mm4_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") != std::string::npos) ) {mm4_trig=1;};
		  if( !mm5_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") != std::string::npos) ) {mm5_trig=1;};
		  if( !mm6_trig && (tpath.find("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") != std::string::npos) ) {mm6_trig=1;};
		  
		  if( !et1_trig  && (tpath.find("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v") != std::string::npos) ) {et1_trig=1;};
		  if( !et2_trig  && (tpath.find("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v") != std::string::npos) ) {et2_trig=1;};
		  if( !et3_trig  && (tpath.find("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v") != std::string::npos) ) {et3_trig=1;};
		  if( !mt1_trig  && (tpath.find("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v") != std::string::npos) ) {mt1_trig=1;};
		  
		  if( !tt1_trig  && (tpath.find("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v") != std::string::npos) ) {tt1_trig=1;};
		  if( !tt2_trig  && (tpath.find("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v") != std::string::npos) ) {tt2_trig=1;};
		  
		  if( !eee1_trig && (tpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") != std::string::npos) ) {eee1_trig=1;};
		  if( !eem1_trig && (tpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") != std::string::npos) ) {eem1_trig=1;};
		  if( !emm1_trig && (tpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v") != std::string::npos) ) {emm1_trig=1;};
		  if( !mmm1_trig && (tpath.find("HLT_TripleMu_12_10_5_v") != std::string::npos) ) {mmm1_trig=1;};
	       }	     
	     else if( year == 2017 )
	       { 		  
		  if( !e1_trig && (tpath.find("HLT_Ele35_WPTight_Gsf_v") != std::string::npos) ) {e1_trig=1;}
		  if( !m1_trig && (tpath.find("HLT_IsoMu24_v") != std::string::npos) ) {m1_trig=1;}
		  if( !m2_trig && (tpath.find("HLT_IsoMu27_v") != std::string::npos) ) {m2_trig=1;}
		  
		  if( !ee1_trig && (tpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {ee1_trig=1;};
		  if( !em1_trig && (tpath.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {em1_trig=1;};
		  if( !em2_trig && (tpath.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {em2_trig=1;};
		  if( !em3_trig && (tpath.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em3_trig=1;};
		  if( !em4_trig && (tpath.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em4_trig=1;};
		  if( !mm1_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") != std::string::npos) ) {mm1_trig=1;};
		  if( !mm2_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v") != std::string::npos) ) {mm2_trig=1;};
		  
		  if( !et1_trig  && (tpath.find("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v") != std::string::npos) ) {et1_trig=1;};
		  if( !mt1_trig  && (tpath.find("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v") != std::string::npos) ) {mt1_trig=1;};
		  
		  if( !tt1_trig  && (tpath.find("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v") != std::string::npos) ) {tt1_trig=1;};
		  if( !tt2_trig  && (tpath.find("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v") != std::string::npos) ) {tt2_trig=1;};
		  if( !tt3_trig  && (tpath.find("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v") != std::string::npos) ) {tt3_trig=1;};
		  if( !tt4_trig  && (tpath.find("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v") != std::string::npos) ) {tt4_trig=1;};
		  
		  if( !eee1_trig && (tpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") != std::string::npos) ) {eee1_trig=1;};
		  if( !eem1_trig && (tpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v") != std::string::npos) ) {eem1_trig=1;};
		  if( !eem2_trig && (tpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") != std::string::npos) ) {eem2_trig=1;};
		  if( !emm1_trig && (tpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v") != std::string::npos) ) {emm1_trig=1;};
		  if( !mmm1_trig && (tpath.find("HLT_TripleMu_12_10_5_v") != std::string::npos) ) {mmm1_trig=1;};
	       }	     
	     else
	       { 		  
		  if( !e1_trig && (tpath.find("HLT_Ele32_WPTight_Gsf_v") != std::string::npos) ) {e1_trig=1;}
		  if( !m1_trig && (tpath.find("HLT_IsoMu24_v") != std::string::npos) ) {m1_trig=1;}
		  if( !m2_trig && (tpath.find("HLT_IsoMu27_v") != std::string::npos) ) {m2_trig=1;}
		  
		  if( !ee1_trig && (tpath.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {ee1_trig=1;};
		  if( !em1_trig && (tpath.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {em1_trig=1;};
		  if( !em2_trig && (tpath.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) ) {em2_trig=1;};
		  if( !em3_trig && (tpath.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) ) {em3_trig=1;};
		  if( !mm1_trig && (tpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") != std::string::npos) ) {mm1_trig=1;};
		  
		  if( !et1_trig  && (tpath.find("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v") != std::string::npos) ) {et1_trig=1;};
		  if( !et2_trig  && (tpath.find("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v") != std::string::npos) ) {et2_trig=1;};
		  if( !mt1_trig  && (tpath.find("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v") != std::string::npos) ) {mt1_trig=1;};
		  if( !mt2_trig  && (tpath.find("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v") != std::string::npos) ) {mt2_trig=1;};

		  if( !tt1_trig  && (tpath.find("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v") != std::string::npos) ) {tt1_trig=1;};
		  if( !tt2_trig  && (tpath.find("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v") != std::string::npos) ) {tt2_trig=1;};
		  if( !tt3_trig  && (tpath.find("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v") != std::string::npos) ) {tt3_trig=1;};
		  if( !tt4_trig  && (tpath.find("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v") != std::string::npos) ) {tt4_trig=1;};

		  if( !eee1_trig && (tpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") != std::string::npos) ) {eee1_trig=1;};
		  if( !eem1_trig && (tpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") != std::string::npos) ) {eem1_trig=1;};
		  if( !emm1_trig && (tpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v") != std::string::npos) ) {emm1_trig=1;};
		  if( !mmm1_trig && (tpath.find("HLT_TripleMu_12_10_5_v") != std::string::npos) ) {mmm1_trig=1;};
	       }	     
	  }
     }   

   if( year == 2016 )
     {	
	if( isdata )
	  {	     
	     trig_e = (e1_trig || e2_trig || e3_trig);
	     trig_et = (et1_trig || et2_trig || et3_trig);
	     trig_m = (m1_trig || m2_trig || m3_trig || m4_trig || m5_trig || m6_trig);
	     trig_mt = (mt1_trig);
	     
	     trig_ee = (ee1_trig || ee2_trig);
	     trig_em = (em1_trig || em2_trig || em3_trig || em4_trig || em5_trig);
	     trig_mm = (mm1_trig || mm2_trig || mm3_trig);
	     
	     trig_tt = (tt1_trig || tt2_trig);
	     
	     trig_eee = (eee1_trig);
	     trig_eem = (eem1_trig);
	     trig_emm = (emm1_trig);
	     trig_mmm = (mmm1_trig);
	  }
	else
	  {
	     trig_e = (e1_trig || e2_trig || e3_trig);
	     trig_et = (et1_trig || et2_trig || et3_trig);
	     trig_m = (m1_trig || m2_trig || m3_trig || m4_trig || m5_trig || m6_trig);
	     trig_mt = (mt1_trig);
	     
	     trig_ee = (ee1_trig || ee2_trig);
	     trig_em = (em1_trig || em3_trig || em5_trig || em6_trig || em7_trig);
	     trig_mm = (mm4_trig || mm5_trig || mm6_trig);
	     
	     trig_tt = (tt1_trig || tt2_trig);
	     
	     trig_eee = (eee1_trig);
	     trig_eem = (eem1_trig);
	     trig_emm = (emm1_trig);
	     trig_mmm = (mmm1_trig);
	  }	
     }   
   else if( year == 2017 )
     {	
	if( isdata )
	  {	     
	     trig_e = (e1_trig);
	     trig_et = (et1_trig);
	     trig_m = (m1_trig || m2_trig);
	     trig_mt = (mt1_trig);
	     
	     trig_ee = (ee1_trig);
	     trig_em = (em1_trig || em2_trig || em3_trig || em4_trig);
	     trig_mm = (mm1_trig || mm2_trig);
	     
	     trig_tt = (tt1_trig || tt2_trig || tt3_trig || tt4_trig);
	     
	     trig_eee = (eee1_trig);
	     trig_eem = (eem1_trig);
	     trig_emm = (emm1_trig);
	     trig_mmm = (mmm1_trig);
	  }
	else
	  {
	     trig_e = (e1_trig);
	     trig_et = (et1_trig);
	     trig_m = (m1_trig || m2_trig);
	     trig_mt = (mt1_trig);
	     
	     trig_ee = (ee1_trig);
	     trig_em = (em1_trig || em2_trig || em3_trig || em4_trig);
	     trig_mm = (mm1_trig || mm2_trig);
	     
	     trig_tt = (tt1_trig || tt2_trig || tt3_trig || tt4_trig);
	     
	     trig_eee = (eee1_trig);
	     trig_eem = (eem2_trig);
	     trig_emm = (emm1_trig);
	     trig_mmm = (mmm1_trig);
	  }	
     }   
   else
     {	
	if( isdata )
	  {	     
	     trig_e = (e1_trig);
	     trig_et = (et1_trig || et2_trig);
	     trig_m = (m1_trig || m2_trig);
	     trig_mt = (mt1_trig || mt2_trig);
	     
	     trig_ee = (ee1_trig);
	     trig_em = (em1_trig || em2_trig || em3_trig);
	     trig_mm = (mm1_trig);
	     
	     trig_tt = (tt1_trig || tt2_trig || tt3_trig || tt4_trig);
	     
	     trig_eee = (eee1_trig);
	     trig_eem = (eem1_trig);
	     trig_emm = (emm1_trig);
	     trig_mmm = (mmm1_trig);
	  }
	else
	  {
	     trig_e = (e1_trig);
	     trig_et = (et1_trig || et2_trig);
	     trig_m = (m1_trig || m2_trig);
	     trig_mt = (mt1_trig || mt2_trig);
	     
	     trig_ee = (ee1_trig);
	     trig_em = (em1_trig || em2_trig || em3_trig);
	     trig_mm = (mm1_trig);
	     
	     trig_tt = (tt1_trig || tt2_trig || tt3_trig || tt4_trig);
	     
	     trig_eee = (eee1_trig);
	     trig_eem = (eem1_trig);
	     trig_emm = (emm1_trig);
	     trig_mmm = (mmm1_trig);
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
   metLD              	 = -1000;
   metLD_JESup = -1000;
   metLD_JESdown = -1000;
   metLD_JERup = -1000;
   metLD_JERdown = -1000;
   
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
   
   /*
   weight_originalXWGTUP   = -1000;
   weight_scale_muF0p5   = -1000;
   weight_scale_muF2     = -1000;
   weight_scale_muR0p5   = -1000;
   weight_scale_muR2     = -1000;
   weight_scale_muR2muF2     = -1000;
   weight_scale_muR0p5muF0p5     = -1000;*/

   weight_scale_index2 = -1000;
   weight_scale_index3 = -1000;
   weight_scale_index4 = -1000;
   weight_scale_index5 = -1000;
   weight_scale_index6 = -1000;
   weight_scale_index7 = -1000;
   weight_scale_index8 = -1000;
   weight_scale_index9 = -1000;   
   
   mc_weight             = -1000;
   mc_weight_originalValue = -1000;
   mc_ptHat              = -1000;
   mc_pu_trueNumInt      = -1000;
   
   tth_channel           = -1000;
   
   disc_TT               = -1000;
   
   pdf_weights.clear();
   pdf_ids.clear();
   
   trig_e = 0;
   trig_et = 0;
   trig_m = 0;
   trig_mt = 0;
   
   trig_tt = 0;
   
   trig_ee = 0;
   trig_em = 0;
   trig_mm = 0;
	     
   trig_eee = 0;
   trig_eem = 0;
   trig_emm = 0;
   trig_mmm = 0;

   is_0l2tau             = 0;
   is_0l2tau_SR_Data     = 0;
   is_0l2tau_SR          = 0;
   is_0l2tau_Fake        = 0;
   is_1l1tau             = 0;
   is_1l1tau_SR_Data     = 0;
   is_1l1tau_SR          = 0;
   is_1l1tau_Fake        = 0;   
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
   is_2lOS1tau           = 0;
   is_2lOS1tau_SR_Data   = 0;
   is_2lOS1tau_SR        = 0;
   is_2lOS1tau_Fake      = 0;
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
   is_WZctrl             = 0;
   is_WZctrl_SR_Data     = 0;
   is_WZctrl_SR          = 0;
   is_WZctrl_Fake        = 0;
   is_ZZctrl             = 0;
   is_ZZctrl_SR_Data     = 0;
   is_ZZctrl_SR          = 0;
   is_ZZctrl_Fake        = 0;
}
