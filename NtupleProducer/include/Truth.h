#ifndef TRUTH_H
#define TRUTH_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Truth : public Base
{
 public:
   
   Truth();
   virtual ~Truth();
   
   ClassDef(Truth,1)
};

#endif
