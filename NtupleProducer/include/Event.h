#ifndef EVENT_H
#define EVENT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Event : public Base
{
 public:
   Event();
   virtual ~Event();

   float mc_weight; //set to +/-1
   float mc_weight_originalValue; //original weight, can differ from +/-1
   Int_t id;
   int   run;
   int   lumi;
      
   float weight_originalXWGTUP;
   
   /*
   float weight_scale_muF0p5;
   float weight_scale_muF2;
   float weight_scale_muR0p5;
   float weight_scale_muR2;
   float weight_scale_muR2muF2;
   float weight_scale_muR0p5muF0p5;*/

   float weight_scale_index2;
   float weight_scale_index3;
   float weight_scale_index4;
   float weight_scale_index5;
   float weight_scale_index6;
   float weight_scale_index7;
   float weight_scale_index8;
   float weight_scale_index9;
   
   std::vector<float> pdf_weights;
   std::vector<std::string> pdf_ids;
   
   bool is_0l2tau;
   bool is_0l2tau_SR_Data;
   bool is_0l2tau_SR;
   bool is_0l2tau_Fake;
   bool is_1l1tau;
   bool is_1l1tau_SR_Data;
   bool is_1l1tau_SR;
   bool is_1l1tau_Fake;
   bool is_1l1tau_Flip;
   bool is_1l2tau;
   bool is_1l2tau_SR_Data;
   bool is_1l2tau_SR;
   bool is_1l2tau_Fake;
   bool is_2lSS;
   bool is_2lSS_ttH;
   bool is_2lSS_tHq;
   bool is_2lSS_SR_Data;
   bool is_2lSS_SR_Data_ttH;
   bool is_2lSS_SR_Data_tHq;
   bool is_2lSS_SR;
   bool is_2lSS_SR_ttH;
   bool is_2lSS_SR_tHq;
   bool is_2lSS_Fake;
   bool is_2lSS_Fake_ttH;
   bool is_2lSS_Fake_tHq;
   bool is_2lSS_Flip_Data;
   bool is_2lSS_Flip_Data_ttH;
   bool is_2lSS_Flip_Data_tHq;
   bool is_2lSS_Flip;
   bool is_2lSS_Flip_ttH;
   bool is_2lSS_Flip_tHq;
   bool is_2lSS1tau;
   bool is_2lSS1tau_ttH;
   bool is_2lSS1tau_tHq;
   bool is_2lSS1tau_SR_Data;
   bool is_2lSS1tau_SR_Data_ttH;
   bool is_2lSS1tau_SR_Data_tHq;
   bool is_2lSS1tau_SR;
   bool is_2lSS1tau_SR_ttH;
   bool is_2lSS1tau_SR_tHq;
   bool is_2lSS1tau_Fake;
   bool is_2lSS1tau_Fake_ttH;
   bool is_2lSS1tau_Fake_tHq;
   bool is_2lSS1tau_Flip_Data;
   bool is_2lSS1tau_Flip_Data_ttH;
   bool is_2lSS1tau_Flip_Data_tHq;
   bool is_2lSS1tau_Flip;
   bool is_2lSS1tau_Flip_ttH;
   bool is_2lSS1tau_Flip_tHq;
   bool is_2lOS1tau;
   bool is_2lOS1tau_SR_Data;
   bool is_2lOS1tau_SR;
   bool is_2lOS1tau_Fake;
   bool is_2l2tau;
   bool is_2l2tau_SR_Data;
   bool is_2l2tau_SR;
   bool is_2l2tau_Fake;
   bool is_3l;
   bool is_3l_SR_Data;
   bool is_3l_SR;
   bool is_3l_Fake;
   bool is_3l1tau;
   bool is_3l1tau_SR_Data;
   bool is_3l1tau_SR;
   bool is_3l1tau_Fake;
   bool is_4l;
   bool is_4l_SR_Data;
   bool is_4l_SR;
   bool is_4l_Fake;
   bool is_ttWctrl;
   bool is_ttWctrl_SR_Data;
   bool is_ttWctrl_SR;
   bool is_ttWctrl_Fake;
   bool is_ttWctrl_Flip_Data;
   bool is_ttWctrl_Flip;
   bool is_ttZctrl;
   bool is_ttZctrl_SR_Data;
   bool is_ttZctrl_SR;
   bool is_ttZctrl_Fake;
   bool is_WZctrl;
   bool is_WZctrl_SR_Data;
   bool is_WZctrl_SR;
   bool is_WZctrl_Fake;
   bool is_ZZctrl;
   bool is_ZZctrl_SR_Data;
   bool is_ZZctrl_SR;
   bool is_ZZctrl_Fake;
   bool is_2lSS_Training;
   bool is_3l_Training;
   bool is_2lSS_LooseSel;
   bool is_3l_LooseSel;

   float metphi;
   float metpt;
   float metsumet;
   float metLD;

//New -- also stored metLD recomputed for JES/JER variations
   float metLD_JESup;
   float metLD_JESdown;
   float metLD_JERup;
   float metLD_JERdown;

   float metcov00;
   float metcov01;
   float metcov10;
   float metcov11;
   
   bool trig_e;
   bool trig_et;
   bool trig_m;
   bool trig_mt;
   
   bool trig_tt;
   
   bool trig_ee;
   bool trig_em;
   bool trig_mm;
	     
   bool trig_eee;
   bool trig_eem;
   bool trig_emm;
   bool trig_mmm;
   
   int mc_pu_trueNumInt;
   int pv_n;

   ClassDef(Event,1)
};

#endif
