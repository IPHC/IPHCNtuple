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
   
   bool isLooseFwdTTH;

   float CSVv2;
   float DeepCSVbtag;
   float DeepFlavourbtag;
   float QGdiscr;

   float jet_partonFlavour  ;
   float jet_hadronFlavour  ;
   float pileupJetId;

   float JER_corr;
   float JER_corr_up;
   float JER_corr_down;

   float pt_JES_up;
   float pt_JES_down;
   float pt_JER_up;
   float pt_JER_down;
   
   float E_JES_up;
   float E_JES_down;
   float E_JER_up;
   float E_JER_down;

   ClassDef(Jet,1)
};

#endif
