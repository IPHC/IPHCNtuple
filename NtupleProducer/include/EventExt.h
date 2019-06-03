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
   
   void read(bool isdata,int year,bool=false);
   void init();
   
   float rho;
   
   float metUncorrectedPt;
   float metUncorrectedPhi;
   float metUncorrectedSumEt;
   
   float metNoHF_pt;
   float metNoHF_phi;
   float metNoHF_sumet;
   
   float pv_z;
   float pv_zError;
   
   float mc_ptHat;
   
   float disc_TT;
  
 public:
   int tth_channel;
   
   ClassDef(EventExt,1)
};

#endif
