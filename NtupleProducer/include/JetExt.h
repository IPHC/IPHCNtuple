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
   void  setJESUncertainty(bool isdata, float unc);
   //void  JECUncertainty(); //obsolete JER
   void apply_JER_smearing(bool, float, float, float, float);

   int ID;
      
   bool isTight;
   bool tightJetID;
   
   int ntrk;
   
   //float CSVv2; //Declared in Jet.h (-> available in output files)
   float cMVAv2;
   float deepCSV;
   float deepCSVudsg;
   float deepCSVb;
   float deepCSVbb;
   float deepCSVc;
   float deepCSVcc;
   
   float jet_genJet_pt      ;
   float jet_genJet_E       ;
   
   float jet_genParton_pt     ;
   float jet_genParton_E      ;
   float jet_genParton_id     ;
   
   ClassDef(JetExt,1)
};

#endif
