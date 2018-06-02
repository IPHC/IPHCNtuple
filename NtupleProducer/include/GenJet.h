#ifndef GENJET_H
#define GENJET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class GenJet : public Base
{
 public:
   
   GenJet();
   virtual ~GenJet();
   
   int    ID;
   
   float  genJet_pt;
   float  genJet_eta;
   float  genJet_phi;
   float  genJet_E;

   ClassDef(GenJet,1)
};

#endif
