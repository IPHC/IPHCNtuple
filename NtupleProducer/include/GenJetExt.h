#ifndef GENJETEXT_H
#define GENJETEXT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "GenJet.h"

class GenJetExt : public GenJet
{
 public:
   
   GenJetExt();
   virtual ~GenJetExt();
   
   void read(); 
   void init();
   bool sel();
   
 protected:
   
   ClassDef(GenJetExt,1)
};

#endif
