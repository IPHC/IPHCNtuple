#ifndef JET_H
#define JET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Jet : public Base
{
 public:
   
   Jet();
   virtual ~Jet();

   bool isLooseBTag;
   bool isMediumBTag;
   bool isTightBTag;
   
   ClassDef(Jet,1)
};

#endif
