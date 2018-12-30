#ifndef ELECTRON_H
#define ELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Electron : public Base
{
 public:
   
   int id;
   int hasPhotonMCMatch;
   
   Electron();
   virtual ~Electron();
   
   ClassDef(Electron,1)
};

#endif
