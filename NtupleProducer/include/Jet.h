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
   float DeepCSVbtag; // <-> (jet_DeepCSVProbb + jet_DeepCSVProbbb)
   
   float jet_partonFlavour  ;
   float jet_hadronFlavour  ;
   float qg;

   ClassDef(Jet,1)
};

#endif
