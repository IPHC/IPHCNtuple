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
   float qgtag;
   float pileupJetId;

   //Moved from JetExt.h
   float JES_uncert;
   float pt_JER_down;
   float pt_JER_up;

   ClassDef(Jet,1)
};

#endif
