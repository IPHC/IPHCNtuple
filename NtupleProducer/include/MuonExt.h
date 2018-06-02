#ifndef MUONEXT_H
#define MUONEXT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Muon.h"

class MuonExt : public Muon
{
 public:
   
   MuonExt();
   virtual ~MuonExt();
   
   void sel();
   void read();
   void init();
   
   int ID;
   
   int fakeType;
   
   float ptCor;
   float ptUnc;
   int id;
   
   bool isLoose;
   bool isMedium;
   bool isTight;
   bool isPFMuon;
   
   float dxy;
   float dz;
   float iso;
   float isoR04;
   float sip3d;
   float bestTrack_pt;
   float bestTrack_ptError;
   bool cutEventSel;
   bool noLostHits;

   float lepMVA_miniRelIsoCharged;
   float lepMVA_miniRelIsoNeutral;
   float lepMVA_jetPtRelv2;
   float lepMVA_jetPtRatio;
   float lepMVA_jetBTagCSV;
   float lepMVA_jetBTagDeepCSV;
   float lepMVA_sip3d;
   float lepMVA_dxy;
   float lepMVA_dz;
   float lepMVA_mvaId;
   float lepMVA_eta;
   float lepMVA_jetNDauChargedMVASel;

   float gen_pt;
   float gen_eta;
   float gen_phi;
   float gen_m;
   float gen_E;
   int gen_status;
   int gen_id;
   int gen_charge;
   float gen_dr;
   
   ClassDef(MuonExt,1)
};

#endif
