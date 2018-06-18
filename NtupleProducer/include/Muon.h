#ifndef MUON_H
#define MUON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Muon : public Base
{
 public:
  
   int id;
   
   Muon();
   virtual ~Muon();
   
   ClassDef(Muon,1)
};

#endif
