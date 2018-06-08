#ifndef TRUTHEXT_H
#define TRUTHEXT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Truth.h"

class TruthExt : public Truth
{
 public:
   
   TruthExt();
   virtual ~TruthExt();

   void read();
   void readMultiLepton();
   void init();

 protected:
    
   float metGen_px;
   float metGen_py;
   float metGen_pt;
   float metGen_phi;
   float metGen_sumet;
   float metGen_MuonEt;
   
   std::vector<int> Bjets_id;
   std::vector<int> Leptons_id;
   std::vector<int> Jets_id;
   std::vector<int> AllJets_id;
   std::vector<int> JetsHighestPt_id;
   std::vector<int> JetsClosestMw_id;
   std::vector<int> JetsLowestMjj_id;
   std::vector<int> QuarksFromWs_id;
   std::vector<int> JetsFromWs_id;
   
   std::vector<float> Bjets_pt;
   std::vector<float> Leptons_pt;
   std::vector<float> Jets_pt;
   std::vector<float> AllJets_pt;
   std::vector<float> JetsHighestPt_pt;
   std::vector<float> JetsClosestMw_pt;
   std::vector<float> JetsLowestMjj_pt;
   std::vector<float> QuarksFromWs_pt;
   std::vector<float> JetsFromWs_pt;

   std::vector<float> Bjets_eta;
   std::vector<float> Leptons_eta;
   std::vector<float> Jets_eta;
   std::vector<float> AllJets_eta;
   std::vector<float> JetsHighestPt_eta;
   std::vector<float> JetsClosestMw_eta;
   std::vector<float> JetsLowestMjj_eta;
   std::vector<float> QuarksFromWs_eta;
   std::vector<float> JetsFromWs_eta;
   
   std::vector<float> Bjets_phi;
   std::vector<float> Leptons_phi;
   std::vector<float> Jets_phi;
   std::vector<float> AllJets_phi;
   std::vector<float> JetsHighestPt_phi;
   std::vector<float> JetsClosestMw_phi;
   std::vector<float> JetsLowestMjj_phi;
   std::vector<float> QuarksFromWs_phi;
   std::vector<float> JetsFromWs_phi;
   
   std::vector<float> Bjets_E;
   std::vector<float> Leptons_E;
   std::vector<float> Jets_E;
   std::vector<float> AllJets_E;
   std::vector<float> JetsHighestPt_E;
   std::vector<float> JetsClosestMw_E;
   std::vector<float> JetsLowestMjj_E;
   std::vector<float> QuarksFromWs_E;
   std::vector<float> JetsFromWs_E;
   
   int boson_decay;
   int ttbar_decay;

   ClassDef(TruthExt,1)
};

#endif
