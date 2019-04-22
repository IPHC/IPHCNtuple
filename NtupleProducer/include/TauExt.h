#ifndef TauExt_H
#define TauExt_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Tau.h"

class TauExt : public Tau
{
 public:
   
   TauExt();
   virtual ~TauExt();
   
   void sel(bool=false,int year=2016);
   void read(bool isdata);
   void init();
   
   int   ID;
   
   int   fakeType;
   
   float ptCor;
   float ptUnc;
   
   float dxy;
   float dz;
   
   float lepMVA_TTH;
   bool passTightCharge;
   bool cutEventSel;
   bool noLostHits;

   int   decayMode;
   bool  hasLeadChargedHadrCand;
   float leadingTrackPt;
   float leadingTrackCharge;

   float chargedIsoPtSum;
   float neutralIsoPtSum;
   float puCorrPtSum;
   
   float pfEssential_jet_pt;
   float pfEssential_jet_eta;
   float pfEssential_jet_phi;
   float pfEssential_jet_m;
   float pfEssential_jetCorr_pt;
   float pfEssential_jetCorr_eta;
   float pfEssential_jetCorr_phi;
   float pfEssential_jetCorr_m;
   float pfEssential_hasSV;
   float pfEssential_sv_x;
   float pfEssential_sv_y;
   float pfEssential_sv_z;
   float pfEssential_flightLengthSig;
   float pfEssential_dxy;
   float pfEssential_dxy_error;
   float pfEssential_dxy_Sig;

   float decayModeFinding;
   float decayModeFindingNewDMs;
   bool byLooseCombinedIsolationDeltaBetaCorr3Hits;
   bool byMediumCombinedIsolationDeltaBetaCorr3Hits;
   bool byTightCombinedIsolationDeltaBetaCorr3Hits;
   bool againstElectronVLooseMVA6;
   bool againstElectronLooseMVA6;
   bool againstElectronMediumMVA6;
   bool againstElectronTightMVA6;

   bool byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   bool byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   bool byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   bool byTightIsolationMVArun2v1DBdR03oldDMwLT;
   bool byVTightIsolationMVArun2v1DBdR03oldDMwLT;

   bool byCombinedIsolationDeltaBetaCorrRaw3Hits;
   bool againstMuonLoose3;
   bool againstMuonTight3;

   float gen_pt;
   float gen_eta;
   float gen_phi;
   float gen_m;
   float gen_E;
   int gen_status;
   int gen_id;
   int gen_charge;
   float gen_dr;
   
   ClassDef(TauExt,1)
};

#endif
