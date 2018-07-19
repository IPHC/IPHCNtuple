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

   float mc_weight;
   float   id;
   int   run;
   int   lumi;
   
   bool is_1l2tau;
   bool is_1l2tau_SR_Data;
   bool is_1l2tau_SR;
   bool is_1l2tau_Fake;
   bool is_2lSS;
   bool is_2lSS_SR_Data;
   bool is_2lSS_SR;
   bool is_2lSS_Fake;
   bool is_2lSS_Flip_Data;
   bool is_2lSS_Flip;
   bool is_2lSS1tau;
   bool is_2lSS1tau_SR_Data;
   bool is_2lSS1tau_SR;
   bool is_2lSS1tau_Fake;
   bool is_2lSS1tau_Flip_Data;
   bool is_2lSS1tau_Flip;
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

   float metphi;
   float metpt;
   float metsumet;

   float metcov00;
   float metcov01;
   float metcov10;
   float metcov11;
   
   bool trig_e;
   bool trig_etau;
   bool trig_m;
   bool trig_mtau;
   
   bool trig_ee;
   bool trig_em;
   bool trig_mm;
	     
   bool trig_eee;
   bool trig_eem;
   bool trig_emm;
   bool trig_mmm;
   
   ClassDef(Event,1)
};

#endif
