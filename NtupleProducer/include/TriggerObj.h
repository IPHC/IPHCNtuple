#ifndef TriggerObj_H
#define TriggerObj_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class TriggerObj : public Base
{
 public:
   
   TriggerObj();
   virtual ~TriggerObj();
   
   ClassDef(TriggerObj,1)
};

#endif
