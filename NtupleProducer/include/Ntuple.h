#ifndef NTUPLE_H
#define NTUPLE_H

#include "EventExt.h"
#include "ElectronExt.h"
#include "MuonExt.h"
#include "TauExt.h"
#include "JetExt.h"
#include "TruthExt.h"
#include "GenJetExt.h"
#include "TriggerObjExt.h"

#include "TFile.h"
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
   std::vector<TruthExt>*         NtTruthExt;
   std::vector<GenJetExt>*        NtGenJetExt;
   std::vector<TriggerObjExt>*    NtTriggerObjExt;
   
   void Init();
   
   void setBranchAddress();
   void createVar();
   void clearVar();
   void fill();
   void write();
   void convert();
   
   TFile*  m_file;
   
 private:
   
   TTree*  m_tree;
   TChain* m_chain;
   std::string fname_out;
};

#endif
