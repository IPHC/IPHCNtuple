#ifndef HELPER_H
#define HELPER_H

#include "TObject.h"
#include "TLorentzVector.h"

extern unsigned int idx;

class Helper : public TObject
{
 public:
   Helper();
   virtual ~Helper();
   
   float GetDPhi(float phi1,float phi2);
   float GetDeltaR(float eta1,float phi1,float eta2,float phi2);
   
   ClassDef(Helper,1)
};

#endif
