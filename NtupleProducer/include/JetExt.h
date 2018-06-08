#ifndef JETEXT_H
#define JETEXT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Jet.h"

class JetExt : public Jet
{
 public:
   
   JetExt();
   virtual ~JetExt();

   void  sel(int sync);
   void  read(bool isdata);
   void  init();
   void  setJESUncertainty(float unc);
   void  JECUncertainty();

   int ID;
   
   float qg;
   
   bool isTight;
   bool tightJetID;
   
   int ntrk;
   
   //float CSVv2; //Declared in Jet.h -- do the same with other btaggers ?
   float cMVAv2;
   float deepCSV;
   float deepCSVudsg;
   float deepCSVb;
   float deepCSVbb;
   float deepCSVc;
   float deepCSVcc;
   
   float jet_partonFlavour  ;
   float jet_hadronFlavour  ;
   
   float jet_genJet_pt      ;
   float jet_genJet_E       ;
   
   float jet_genParton_pt     ;
   float jet_genParton_E      ;
   float jet_genParton_id     ;
   
   float JES_uncert;
   
   float pt_JER;
   float pt_JER_down;
   float pt_JER_up;
   
   ClassDef(JetExt,1)
};

#endif
