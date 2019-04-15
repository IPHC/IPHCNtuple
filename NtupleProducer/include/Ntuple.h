#ifndef NTUPLE_H
#define NTUPLE_H

#include <iostream>

#include "EventExt.h"
#include "ElectronExt.h"
#include "MuonExt.h"
#include "TauExt.h"
#include "JetExt.h"
#include "TruthExt.h"
#include "GenJetExt.h"
#include "TriggerObjExt.h"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"

class Ntuple
{
 public:
   
   Ntuple(std::string fname_out);
   virtual ~Ntuple();
   
   std::vector<Event>*         NtEvent;
   std::vector<Electron>*      NtElectronLoose;
   std::vector<Electron>*      NtElectronFakeable;
   std::vector<Electron>*      NtElectronTight;
   std::vector<Muon>*          NtMuonLoose;
   std::vector<Muon>*          NtMuonFakeable;
   std::vector<Muon>*          NtMuonTight;
   std::vector<Tau>*           NtTauFakeable;
   std::vector<Tau>*           NtTauLoose;
   std::vector<Tau>*           NtTauMedium;
   std::vector<Jet>*           NtJetLoose;
   std::vector<Jet>*           NtJetLooseSoft;
   std::vector<Truth>*         NtTruth;
   std::vector<GenJet>*        NtGenJet;
   std::vector<TriggerObj>*    NtTriggerObj;

   std::vector<EventExt>*         NtEventExt;
   std::vector<ElectronExt>*      NtElectronLooseExt;
   std::vector<ElectronExt>*      NtElectronFakeableExt;
   std::vector<ElectronExt>*      NtElectronTightExt;
   std::vector<MuonExt>*          NtMuonLooseExt;
   std::vector<MuonExt>*          NtMuonFakeableExt;
   std::vector<MuonExt>*          NtMuonTightExt;
   std::vector<TauExt>*           NtTauFakeableExt;
   std::vector<TauExt>*           NtTauLooseExt;
   std::vector<TauExt>*           NtTauMediumExt;
   std::vector<JetExt>*           NtJetLooseExt;
   std::vector<JetExt>*           NtJetLooseSoftExt;
   std::vector<TruthExt>*         NtTruthExt;
   std::vector<GenJetExt>*        NtGenJetExt;
   std::vector<TriggerObjExt>*    NtTriggerObjExt;
   
   //TH1Fs storing sum of weights before preselection for scales
   /*
   TH1F* h_sumWeights_nominal;
   TH1F* h_sumWeightsScale_originalXWGTUP;
   TH1F* h_sumWeightsScale_muF0p5;
   TH1F* h_sumWeightsScale_muF2;
   TH1F* h_sumWeightsScale_muR0p5;
   TH1F* h_sumWeightsScale_muR2;
   TH1F* h_sumWeightsScale_muR2muF2;
   TH1F* h_sumWeightsScale_muR0p5muF0p5;
   */

   void Init();
   
   void setBranchAddress();
   void createVar();
   void clearVar();
   void fill();
   //void fill_histograms(EventExt);
   void write();
   void convert();
   
   TFile*  m_file;
   
 private:
   
   TTree*  m_tree;
   TChain* m_chain;
   std::string fname_out;
};

#endif
