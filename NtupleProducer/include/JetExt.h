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

   void  sel(int sync, bool=false, int year=2016);
   void  read(bool isdata);
   void  init();
   void  setJESUncertainty(bool isdata, float unc);
   //void  JECUncertainty(); //obsolete JER
   void apply_JER_smearing(bool, float, float, float, float);

   int ID;
      
   bool isTight;
   bool tightJetID;
   bool looseJetID;
   
   int ntrk;
   
   float cMVAv2;
   float deepCSV;
   float deepCSVudsg;
   float deepCSVb;
   float deepCSVbb;
   float deepCSVc;
   float deepCSVcc;

   float deepFlavourb;
   float deepFlavourbb;
   float deepFlavourlepb;
   float deepFlavourc;
   float deepFlavouruds;
   float deepFlavourg;
   
   float jet_genJet_pt      ;
   float jet_genJet_E       ;
   
   float jet_genParton_pt     ;
   float jet_genParton_E      ;
   float jet_genParton_id     ;
   
   bool isSoftLooseTTH;
   
   ClassDef(JetExt,1)
};

#endif
