#ifndef SYNC_H
#define SYNC_H

#include "Event.h"
#include "Electron.h"
#include "Muon.h"
#include "Tau.h"
#include "Jet.h"
#include "Truth.h"
#include "GenJet.h"
#include "TriggerObj.h"
#include "Ntuple.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

class Sync
{
 public:
   
   Sync(std::string fname_out);
   virtual ~Sync();
   
   void Init();
   
   void setBranchAddress();
   void initVar();
   void get(Ntuple *nt,int n_presel_el,int n_presel_mu,int n_presel_tau,int n_presel_jet);
   void fill();
   
   TFile*  m_file;

   // https://gitlab.cern.ch/ttH_leptons/doc/blob/master/2017/taus/sync/vars.md
   
   int nEvent;
   int ls;
   int run;
   int n_presel_mu;
   int n_fakeablesel_mu;
   int n_mvasel_mu;
   int n_presel_ele;
   int n_fakeablesel_ele;
   int n_mvasel_ele;
   int n_presel_tau;
   int n_presel_jet;
   
   float mu1_pt;
   float mu1_conept;
   float mu1_eta;
   float mu1_phi;
   float mu1_E;
   float mu1_charge;
   float mu1_miniRelIso;
   float mu1_miniIsoCharged;
   float mu1_miniIsoNeutral;
   float mu1_PFRelIso04;
   float mu1_jetNDauChargedMVASel;
   float mu1_jetPtRel;
   float mu1_jetPtRatio;
   float mu1_jetCSV;
   float mu1_sip3D;
   float mu1_dxy;
   float mu1_dxyAbs;
   float mu1_dz;
   float mu1_segmentCompatibility;
   float mu1_leptonMVA;
   float mu1_mediumID;
   float mu1_dpt_div_pt;
   bool mu1_isfakeablesel;
   bool mu1_ismvasel;

   float mu2_pt;
   float mu2_conept;
   float mu2_eta;
   float mu2_phi;
   float mu2_E;
   float mu2_charge;
   float mu2_miniRelIso;
   float mu2_miniIsoCharged;
   float mu2_miniIsoNeutral;
   float mu2_PFRelIso04;
   float mu2_jetNDauChargedMVASel;
   float mu2_jetPtRel;
   float mu2_jetPtRatio;
   float mu2_jetCSV;
   float mu2_sip3D;
   float mu2_dxy;
   float mu2_dxyAbs;
   float mu2_dz;
   float mu2_segmentCompatibility;
   float mu2_leptonMVA;
   float mu2_mediumID;
   float mu2_dpt_div_pt;
   int mu2_isfakeablesel;
   int mu2_ismvasel;
   
   float ele1_pt;
   float ele1_conept;
   float ele1_eta;
   float ele1_phi;
   float ele1_E;
   float ele1_charge;
   float ele1_miniRelIso;
   float ele1_miniIsoCharged;
   float ele1_miniIsoNeutral;
   float ele1_PFRelIso04;
   float ele1_jetNDauChargedMVASel;
   float ele1_jetPtRel;
   float ele1_jetPtRatio;
   float ele1_jetCSV;
   float ele1_sip3D;
   float ele1_dxy;
   float ele1_dxyAbs;
   float ele1_dz;
   float ele1_ntMVAeleID;
   float ele1_leptonMVA;
   int ele1_isChargeConsistent;
   float ele1_passesConversionVeto;
   float ele1_nMissingHits;
   float ele1_sigmaEtaEta;
   float ele1_HoE;
   float ele1_deltaEta;
   float ele1_deltaPhi;
   float ele1_OoEminusOoP;
   int ele1_isfakeablesel;
   int ele1_ismvasel;

   float ele2_pt;
   float ele2_conept;
   float ele2_eta;
   float ele2_phi;
   float ele2_E;
   float ele2_charge;
   float ele2_miniRelIso;
   float ele2_miniIsoCharged;
   float ele2_miniIsoNeutral;
   float ele2_PFRelIso04;
   float ele2_jetNDauChargedMVASel;
   float ele2_jetPtRel;
   float ele2_jetPtRatio;
   float ele2_jetCSV;
   float ele2_sip3D;
   float ele2_dxy;
   float ele2_dxyAbs;
   float ele2_dz;
   float ele2_ntMVAeleID;
   float ele2_leptonMVA;
   int ele2_isChargeConsistent;
   float ele2_passesConversionVeto;
   float ele2_nMissingHits;
   float ele2_sigmaEtaEta;
   float ele2_HoE;
   float ele2_deltaEta;
   float ele2_deltaPhi;
   float ele2_OoEminusOoP;
   int ele2_isfakeablesel;
   int ele2_ismvasel;

   float tau1_pt;
   float tau1_eta;
   float tau1_phi;
   float tau1_E;
   float tau1_charge;
   float tau1_dxy;
   float tau1_dz;
   float tau1_decayModeFindingOldDMs;
   float tau1_decayModeFindingNewDMs;
   float tau1_byCombinedIsolationDeltaBetaCorr3Hits;
   float tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   float tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   float tau1_byTightCombinedIsolationDeltaBetaCorr3Hits;
   float tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03;
   float tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03;
   float tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03;
   float tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   float tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   float tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   float tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   float tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
   float tau1_rawMVArun2v1DBdR03oldDMwLT;
   float tau1_againstMuonLoose3;
   float tau1_againstMuonTight3;
   float tau1_againstElectronVLooseMVA6;
   float tau1_againstElectronLooseMVA6;
   float tau1_againstElectronMediumMVA6;
   float tau1_againstElectronTightMVA6;
   
   float tau2_pt;
   float tau2_eta;
   float tau2_phi;
   float tau2_E;
   float tau2_charge;
   float tau2_dxy;
   float tau2_dz;
   float tau2_decayModeFindingOldDMs;
   float tau2_decayModeFindingNewDMs;
   float tau2_byCombinedIsolationDeltaBetaCorr3Hits;
   float tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   float tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   float tau2_byTightCombinedIsolationDeltaBetaCorr3Hits;
   float tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03;
   float tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03;
   float tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03;
   float tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   float tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   float tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   float tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   float tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
   float tau2_rawMVArun2v1DBdR03oldDMwLT;
   float tau2_againstMuonLoose3;
   float tau2_againstMuonTight3;
   float tau2_againstElectronVLooseMVA6;
   float tau2_againstElectronLooseMVA6;
   float tau2_againstElectronMediumMVA6;
   float tau2_againstElectronTightMVA6;   
   
   float jet1_pt;
   float jet1_eta;
   float jet1_phi;
   float jet1_E;
   float jet1_CSV;

   float jet2_pt;
   float jet2_eta;
   float jet2_phi;
   float jet2_E;
   float jet2_CSV;

   float jet3_pt;
   float jet3_eta;
   float jet3_phi;
   float jet3_E;
   float jet3_CSV;

   float jet4_pt;
   float jet4_eta;
   float jet4_phi;
   float jet4_E;
   float jet4_CSV;
   
   float PFMET;
   float PFMETphi;
   float MHT;
   float metLD;
   int isGenMatched;
   
   float lep1_conePt;
   float lep2_conePt;
   float lep3_conePt;
   
   float mindr_lep1_jet;
   float mindr_lep2_jet;
   float mindr_lep3_jet;
   
   float mindr_tau1_jet;
   float mindr_tau2_jet;
   float mindr_tau3_jet;
   
   float avg_dr_jet;
   float avr_dr_lep_tau;
   float max_dr_lep_tau;
   float mindr_tau_jet;
   float min_dr_lep_tau;
   float min_dr_lep_jet;
   float dr_leps;
   float dr_taus;
   float dR_lep_tau_ss;
   float dr_lep1_tau;
   float dr_lep2_tau;
   float max_lep_eta;
   float mT_lep1;
   float mT_lep2;
   float mTauTauVis;
   float mTauTauVis1;
   float mTauTauVis2;
   float mbb;
   float mbb_loose;
   float cosThetaS_hadTau;
   float HTT;
   float HadTop_pt;
   float Hj_tagger;
   int nBJetLoose;
   
   float mvaOutput_plainKin_ttV;
   float mvaOutput_plainKin_ttbar;
   float mvaOutput_1l_2tau_HTT_SUM_VT;
   float mvaOutput_2l_2tau_plainKin_1B_VT;
   float mvaOutput_2l_2tau_plainKin_SUM_VT;
   float mvaOutput_2lss_ttV;
   float mvaOutput_2lss_ttbar;
   float mvaOutput_2lss_1tau_plainKin_ttbar;
   float mvaOutput_2lss_1tau_plainKin_ttV;
   float mvaOutput_2lss_1tau_plainKin_1B_M;
   float mvaOutput_2lss_1tau_plainKin_SUM_M;
   float mvaOutput_2lss_1tau_HTT_SUM_M;
   float mvaOutput_2lss_1tau_HTTMEM_SUM_M;
   float mvaOutput_3l_ttV;
   float mvaOutput_3l_ttbar;
   float mvaOutput_3l_1tau_plainKin_SUM_M;
   float mvaOutput_3l_1tau_plainKin_1B_M;
   
   float FR_weight;
   float triggerSF_weight;
   float leptonSF_weight;
   float tauSF_weight;
   float bTagSF_weight;
   float PU_weight;
   float MC_weight;
   
   float Integral_ttH;
   float Integral_ttZ;
   float Integral_ttZ_Zll;
   float Integral_ttbar;
   float integration_type;
   float memOutput_LR;
   
 private:
     
   TTree*  m_tree;
   TChain* m_chain;
   std::string _fname_out;
};

#endif
