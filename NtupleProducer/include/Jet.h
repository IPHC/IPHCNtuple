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
   
   float CSVv2; //Also add other btaggers ?

   ClassDef(Jet,1)
};

#endif
