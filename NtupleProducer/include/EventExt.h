#ifndef EVENTEXT_H
#define EVENTEXT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Event.h"

class EventExt : public Event
{
 public:
   EventExt();
   virtual ~EventExt();
   
   void read(bool isdata);
   void init();
   
   float rho;
   
   float metphi;
   float metsumet;
   float metUncorrectedPt;
   float metUncorrectedPhi;
   float metUncorrectedSumEt;
   
   float metNoHF_pt;
   float metNoHF_phi;
   float metNoHF_sumet;
   
   int   pv_n;
   float pv_z;
   float pv_zError;
   
   float weight_scale_muF0p5;
   float weight_scale_muF2;
   float weight_scale_muR0p5;
   float weight_scale_muR2;
   
   std::vector<float> pdf_weights;
   std::vector<std::string> pdf_ids;
   
   float mc_weight;
   float mc_ptHat;
   int   mc_pu_trueNumInt;
   
   float disc_TT;
  
 public:
   int tth_channel;
   
   ClassDef(EventExt,1)
};

#endif
