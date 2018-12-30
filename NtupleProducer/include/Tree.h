#ifndef Tree_h
#define Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class Tree {
   
 public :
   
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   // Declaration of leaf types
   // 
   Int_t           ev_run;
   Int_t        ev_id;
   Int_t           ev_lumi;
   Float_t         ev_rho;
   
   vector<int>     *trigger;
   vector<bool>    *trigger_pass;
   vector<string>  *trigger_name;
   
   Float_t         met_pt;
   Float_t         met_phi;
   Float_t         met_sumet;
   Float_t         metPuppi_pt;
   Float_t         metPuppi_phi;
   Float_t         metPuppi_sumet;
//   Float_t         met_uncorrectedPt;
//   Float_t         met_uncorrectedPhi;
//   Float_t         met_uncorrectedSumEt;
   Double_t         met_cov00;
   Double_t         met_cov01;
   Double_t         met_cov10;
   Double_t         met_cov11;
   
   Float_t         metNoHF_pt;
   Float_t         metNoHF_phi;
   Float_t         metNoHF_sumet;
   
   Int_t           nvertex;
   Float_t         pv_x;
   Float_t         pv_y;
   Float_t         pv_z;
   Float_t         pv_zError;
   
   Float_t         weight_scale_muF0p5;
   Float_t         weight_scale_muF2;
   Float_t         weight_scale_muR0p5;
   Float_t         weight_scale_muR2;
   Int_t           mc_id;
   Int_t           mc_f1;
   Int_t           mc_f2;
   Float_t         mc_x1;
   Float_t         mc_x2;
   Float_t         mc_scale;
   Float_t         mc_ptHat;
   Float_t         mc_weight;
   vector<float>       *mc_pdfweights;
   vector<std::string> *mc_pdfweightIds;

   Int_t                 mc_pu_intime_NumInt;
   Int_t                 mc_pu_trueNumInt;
   Int_t                 mc_pu_before_npu;
   Int_t                 mc_pu_after_npu;
   Int_t                 mc_pu_Npvi;
   vector<int>           *mc_pu_Nzpositions;
   vector<int>           *mc_pu_BunchCrossing;
   vector<vector<float>> *mc_pu_zpositions;
   vector<vector<float>> *mc_pu_sumpT_lowpT;
   vector<vector<float>> *mc_pu_sumpT_highpT;
   vector<vector<int>>   *mc_pu_ntrks_lowpT;
   vector<vector<int>>   *mc_pu_ntrks_highpT;
   
   Int_t           el_n;
   vector<float>   *el_pt;
   vector<float>   *el_eta;
   vector<float>   *el_phi;
   vector<float>   *el_m;
   vector<float>   *el_E;
   vector<float>   *el_looseCBId;
   vector<float>   *el_mediumCBId;
   vector<int>     *el_numberOfLostHits;
   vector<float>   *el_gsfTrack_PV_dxy;
   vector<float>   *el_gsfTrack_PV_dz;
   vector<float>   *el_ip3d;
   vector<float>   *el_ip3dErr;
   vector<float>   *el_miniIso;
   vector<float>   *el_miniIsoTTH;
   vector<float>   *el_PFRelIso04;
   vector<float>   *el_conept;
   vector<float>   *el_ooEmooP;
   vector<float>   *el_E_postCorr;
   vector<float>   *el_pt_postCorr;
   
   vector<int>     *el_id;
   vector<int>     *el_charge;
   vector<float>   *el_neutralHadronIso;
   vector<float>   *el_chargedHadronIso;
   vector<float>   *el_puChargedHadronIso;
   vector<float>   *el_ecalIso;
   vector<float>   *el_hcalIso;
   vector<float>   *el_particleIso;
   vector<float>   *el_photonIso;
   vector<float>   *el_trackIso;
   vector<int>     *el_isLoose;
   vector<int>     *el_isTight;
   vector<int>     *el_isRobustLoose;
   vector<int>     *el_isRobustTight;
   vector<int>     *el_isRobustHighEnergy;
   vector<float>   *el_vx;
   vector<float>   *el_vy;
   vector<float>   *el_vz;
   vector<bool>    *el_isGsf;
   vector<float>   *el_dxy;
   vector<float>   *el_dz;
   vector<float>   *el_dxyError;
   vector<float>   *el_dzError;
   vector<float>   *el_mvaIso;
   vector<float>   *el_mvaNoIso;
   vector<bool>    *el_mvaPassMedium;
   vector<bool>    *el_mvaPassTight;
   vector<int>     *el_numberOfHits;
   vector<float>   *el_pfIso_sumChargedHadronPt;
   vector<float>   *el_pfIso_sumNeutralHadronEt;
   vector<float>   *el_pfIso_sumPhotonEt;
   vector<float>   *el_pfIso_sumPUPt;
   vector<float>   *el_lepMVA;
   vector<float>   *el_lepMVA_miniRelIsoCharged;
   vector<float>   *el_lepMVA_miniRelIsoNeutral;
   vector<float>   *el_lepMVA_jetPtRelv2;
   vector<float>   *el_lepMVA_jetPtRatio;
   vector<float>   *el_lepMVA_jetBTagCSV;
   vector<float>   *el_lepMVA_jetBTagDeepCSV;
   vector<float>   *el_lepMVA_sip3d;
   vector<float>   *el_lepMVA_dxy;
   vector<float>   *el_lepMVA_dz;
   vector<float>   *el_lepMVA_mvaId;
   vector<float>   *el_lepMVA_eta;
   vector<float>   *el_lepMVA_jetNDauChargedMVASel;
   vector<int>     *el_isGsfCtfScPixChargeConsistent;
   vector<int>     *el_isGsfScPixChargeConsistent;
   vector<int>     *el_passConversionVeto;
   vector<float>   *el_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *el_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *el_hadronicOverEm;
   vector<float>   *el_scleta;
   vector<bool>    *el_hasMatchedConversion;
   vector<float>   *el_sigmaIetaIeta;
   vector<float>   *el_full5x5_sigmaIetaIeta;
   vector<float>   *el_full5x5_sigmaEtaEta;
   vector<float>   *el_superCluster_eta;
   vector<float>   *el_correctedEcalEnergy;
   vector<float>   *el_ecalEnergy;
   vector<float>   *el_eSuperClusterOverP;
   vector<float>   *el_trackMomentumError;

   vector<int>     *el_hasMCMatch;
   vector<int>     *el_hasChargeMCMatch;
   vector<int>     *el_hasPhotonMCMatch;
   vector<float>   *el_gen_pt;
   vector<float>   *el_gen_eta;
   vector<float>   *el_gen_phi;
   vector<float>   *el_gen_m;
   vector<float>   *el_gen_E;
   vector<int>     *el_gen_status;
   vector<int>     *el_gen_id;
   vector<int>     *el_gen_charge;
   vector<float>   *el_gen_dr;
   
   Int_t           mu_n;
   vector<float>   *mu_pt;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<float>   *mu_m;
   vector<float>   *mu_E;
   vector<int>     *mu_id;
   vector<int>     *mu_charge;
   vector<float>   *mu_ip3d;
   vector<float>   *mu_ip3dErr;
   vector<float>   *mu_miniIso;
   vector<float>   *mu_miniIsoTTH;
   vector<float>   *mu_PFRelIso04;
   vector<bool>    *mu_isLooseMuon;
   vector<float>   *mu_conept;
   
   vector<float>   *mu_neutralHadronIso;
   vector<float>   *mu_chargedHadronIso;
   vector<float>   *mu_ecalIso;
   vector<float>   *mu_hcalIso;
   vector<float>   *mu_photonIso;
   vector<float>   *mu_trackIso;
   vector<int>     *mu_isGlobalMuon;
   vector<int>     *mu_isTrackerMuon;
   vector<int>     *mu_isStandAloneMuon;
   vector<int>     *mu_isCaloMuon;
   vector<int>     *mu_isPFMuon;
   vector<bool>    *mu_isMediumMuon;
   vector<bool>    *mu_isTightMuon;
   vector<float>   *mu_vx;
   vector<float>   *mu_vy;
   vector<float>   *mu_vz;
   vector<float>   *mu_segmentCompatibility;
   vector<int>     *mu_hasGlobalTrack;
   vector<float>   *mu_globalTrack_dxy;
   vector<float>   *mu_globalTrack_dz;
   vector<float>   *mu_globalTrack_dxyError;
   vector<float>   *mu_globalTrack_dzError;
   vector<float>   *mu_globalTrack_normalizedChi2;
   vector<float>   *mu_combinedQuality_chi2LocalPosition;
   vector<float>   *mu_combinedQuality_trkKink;
   vector<int>     *mu_hasInnerTrack;
   vector<float>   *mu_innerTrack_dxy;
   vector<float>   *mu_innerTrack_dz;
   vector<float>   *mu_innerTrack_PV_dxy;
   vector<float>   *mu_innerTrack_PV_dz;
   vector<float>   *mu_innerTrack_dxyError;
   vector<float>   *mu_innerTrack_dzError;
   vector<float>   *mu_innerTrack_validFraction;
   vector<float>   *mu_bestTrack_dxy;
   vector<float>   *mu_bestTrack_PV_dz;
   vector<float>   *mu_bestTrack_dxyError;
   vector<float>   *mu_bestTrack_dzError;
   vector<float>   *mu_bestTrack_pt;
   vector<float>   *mu_bestTrack_ptError;
   vector<int>     *mu_numberOfMatches;
   vector<int>     *mu_numberOfValidMuonHits;
   vector<float>   *mu_pfIso04_sumChargedHadronPt;
   vector<float>   *mu_pfIso04_sumNeutralHadronEt;
   vector<float>   *mu_pfIso04_sumPhotonEt;
   vector<float>   *mu_pfIso04_sumPUPt;
   vector<float>   *mu_lepMVA;
   vector<float>   *mu_lepMVA_miniRelIsoCharged;
   vector<float>   *mu_lepMVA_miniRelIsoNeutral;
   vector<float>   *mu_lepMVA_jetPtRelv2;
   vector<float>   *mu_lepMVA_jetPtRatio;
   vector<float>   *mu_lepMVA_jetBTagCSV;
   vector<float>   *mu_lepMVA_jetBTagDeepCSV;
   vector<float>   *mu_lepMVA_sip3d;
   vector<float>   *mu_lepMVA_dxy;
   vector<float>   *mu_lepMVA_dz;
   vector<float>   *mu_lepMVA_mvaId;
   vector<float>   *mu_lepMVA_eta;
   vector<float>   *mu_lepMVA_jetNDauChargedMVASel;
   vector<float>   *mu_innerTrack_pt;
   vector<float>   *mu_innerTrack_ptError;

   vector<int>     *mu_hasMCMatch;
   vector<int>     *mu_hasChargeMCMatch;
   vector<float>   *mu_gen_pt;
   vector<float>   *mu_gen_eta;
   vector<float>   *mu_gen_phi;
   vector<float>   *mu_gen_m;
   vector<float>   *mu_gen_E;
   vector<int>     *mu_gen_status;
   vector<int>     *mu_gen_id;
   vector<int>     *mu_gen_charge;
   vector<float>   *mu_gen_dr;
   
   Int_t           tau_n;
   vector<float>   *tau_E;
   vector<float>   *tau_pt;
   vector<float>   *tau_eta;
   vector<float>   *tau_phi;
   vector<float>   *tau_m;
   vector<float>   *tau_leadingTrackDxy;
   vector<float>   *tau_leadingTrackDz;
   vector<int>     *tau_charge;
   vector<int>     *tau_id;
   vector<int>     *tau_decayMode;
   vector<bool>    *tau_hasLeadChargedHadrCand;
   vector<float>   *tau_leadingTrackPt;
   vector<float>   *tau_leadingTrackCharge;
   vector<float>   *tau_decayModeFinding;
//   vector<float>   *tau_decayModeFindingOldDMs;
   vector<float>   *tau_decayModeFindingNewDMs;
   vector<float>   *tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_againstElectronVLooseMVA6;
   vector<float>   *tau_againstElectronLooseMVA6;
   vector<float>   *tau_againstElectronMediumMVA6;
   vector<float>   *tau_againstElectronTightMVA6;

   vector<float>   *tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
   
   vector<float>   *tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<float>   *tau_chargedIsoPtSum;
   vector<float>   *tau_neutralIsoPtSum;
   vector<float>   *tau_puCorrPtSum;
   vector<float>   *tau_againstMuonLoose3;
   vector<float>   *tau_againstMuonTight3;
   vector<float>   *tau_pfEssential_jet_pt;
   vector<float>   *tau_pfEssential_jet_eta;
   vector<float>   *tau_pfEssential_jet_phi;
   vector<float>   *tau_pfEssential_jet_m;
   vector<float>   *tau_pfEssential_jetCorr_pt;
   vector<float>   *tau_pfEssential_jetCorr_eta;
   vector<float>   *tau_pfEssential_jetCorr_phi;
   vector<float>   *tau_pfEssential_jetCorr_m;
   vector<float>   *tau_pfEssential_hasSV;
   vector<float>   *tau_pfEssential_sv_x;
   vector<float>   *tau_pfEssential_sv_y;
   vector<float>   *tau_pfEssential_sv_z;
   vector<float>   *tau_pfEssential_flightLengthSig;
   vector<float>   *tau_pfEssential_dxy;
   vector<float>   *tau_pfEssential_dxy_error;
   vector<float>   *tau_pfEssential_dxy_Sig;

   vector<int>     *tau_hasMCMatch;
   vector<int>     *tau_hasChargeMCMatch;
   vector<float>   *tau_gen_pt;
   vector<float>   *tau_gen_eta;
   vector<float>   *tau_gen_phi;
   vector<float>   *tau_gen_m;
   vector<float>   *tau_gen_E;
   vector<int>     *tau_gen_status;
   vector<int>     *tau_gen_id;
   vector<int>     *tau_gen_charge;
   vector<float>   *tau_gen_dr;
   
   Int_t           jet_n;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_m;
   vector<float>   *jet_E;
   vector<float>   *jet_qgtag;
   vector<int>     *jet_ntrk;
   vector<float>   *jet_CSVv2;
   vector<float>   *jet_cMVAv2;
   vector<float>   *jet_DeepCSVProbudsg;
   vector<float>   *jet_DeepCSVProbb;
   vector<float>   *jet_DeepCSVProbbb;
   vector<float>   *jet_DeepCSVProbc;
   vector<float>   *jet_DeepCSVProbcc;
   vector<bool>    *jet_tightJetID;
   vector<int>     *jet_partonFlavour;
   vector<int>     *jet_hadronFlavour;
   vector<float>   *jet_neutralHadronEnergy;
   vector<float>   *jet_neutralEmEnergy;
   vector<float>   *jet_neutralEmEnergyFraction; //CHANGED -- for debug
   vector<float>   *jet_chargedHadronEnergy;
   vector<float>   *jet_chargedEmEnergy;
   vector<float>   *jet_electronEnergy;
   vector<float>   *jet_muonEnergy;
   vector<float>   *jet_photonEnergy;
   vector<float>   *jet_genJet_pt;
   vector<float>   *jet_genJet_eta;
   vector<float>   *jet_genJet_phi;
   vector<float>   *jet_genJet_m;
   vector<float>   *jet_genJet_E;
   vector<int>     *jet_genJet_status;
   vector<int>     *jet_genJet_id;
   vector<float>   *jet_genParton_pt;
   vector<float>   *jet_genParton_eta;
   vector<float>   *jet_genParton_phi;
   //vector<float>   *jet_genParton_m;
   vector<float>   *jet_genParton_E;
   //vector<int>     *jet_genParton_status;
   vector<int>     *jet_genParton_id;
   vector<float>   *jet_pileupJetId;
   
   Int_t gen_n;
   std::vector<float> *gen_pt;
   std::vector<float> *gen_eta;
   std::vector<float> *gen_phi;
   std::vector<float> *gen_m;
   std::vector<int> *gen_id;
   std::vector<int> *gen_status;
   std::vector<int> *gen_index;
   std::vector<int> *gen_mother_index;

   Float_t gen_PVz;
   
   Float_t metGen_px;
   Float_t metGen_py;
   Float_t metGen_pt;
   Float_t metGen_phi;
   Float_t metGen_sumet;
   Float_t metGen_MuonEt;
   
   Int_t           mc_truth_h0_id;
   Int_t           mc_truth_h0W1_id;
   Int_t           mc_truth_h0Wl1_id;
   Int_t           mc_truth_h0Wtau1_id;
   Int_t           mc_truth_h0Wtaul1_id;
   Int_t           mc_truth_h0Wq11_id;
   Int_t           mc_truth_h0Wq21_id;
   Int_t           mc_truth_h0W2_id;
   Int_t           mc_truth_h0Wl2_id;
   Int_t           mc_truth_h0Wtau2_id;
   Int_t           mc_truth_h0Wtaul2_id;
   Int_t           mc_truth_h0Wq12_id;
   Int_t           mc_truth_h0Wq22_id;
   Int_t           mc_truth_h0Z1_id;
   Int_t           mc_truth_h0Zl11_id;
   Int_t           mc_truth_h0Zl21_id;
   Int_t	        mc_truth_h0Ztau11_id;
   Int_t	        mc_truth_h0Ztau21_id;
   Int_t	        mc_truth_h0Ztaul11_id;
   Int_t	        mc_truth_h0Ztaul21_id;
   Int_t           mc_truth_h0Zq11_id;
   Int_t           mc_truth_h0Zq21_id;
   Int_t           mc_truth_h0Z2_id;
   Int_t           mc_truth_h0Zl12_id;
   Int_t           mc_truth_h0Zl22_id;
   Int_t	        mc_truth_h0Ztau12_id;
   Int_t	        mc_truth_h0Ztau22_id;
   Int_t	        mc_truth_h0Ztaul12_id;
   Int_t	        mc_truth_h0Ztaul22_id;
   Int_t           mc_truth_h0Zq12_id;
   Int_t           mc_truth_h0Zq22_id;
   Int_t           mc_truth_h0tau1_id;
   Int_t           mc_truth_h0tau2_id;
   Int_t           mc_truth_h0taul1_id;
   Int_t           mc_truth_h0taul2_id;
   Int_t           mc_truth_h0b1_id;
   Int_t           mc_truth_h0b2_id;
   
   Int_t           mc_truth_t1_id;
   Int_t           mc_truth_t2_id;
   Int_t           mc_truth_tb1_id;
   Int_t           mc_truth_tb2_id;
   Int_t           mc_truth_tW1_id;
   Int_t           mc_truth_tW2_id;
   Int_t           mc_truth_tWl1_id;
   Int_t           mc_truth_tWl2_id;
   Int_t           mc_truth_tWtau1_id;
   Int_t           mc_truth_tWtau2_id;
   Int_t           mc_truth_tWtaul1_id;
   Int_t           mc_truth_tWtaul2_id;
   Int_t           mc_truth_tWq11_id;
   Int_t           mc_truth_tWq21_id;
   Int_t           mc_truth_tWq12_id;
   Int_t           mc_truth_tWq22_id;
   
   Int_t           mc_truth_t_id;
   Int_t           mc_truth_tb_id;
   Int_t           mc_truth_tW_id;
   Int_t           mc_truth_tWl_id;
   Int_t           mc_truth_tWq1_id;
   Int_t           mc_truth_tWq2_id;
   
   Int_t           mc_truth_W_id;
   Int_t           mc_truth_Wl_id;
   Int_t           mc_truth_Wtau_id;
   Int_t           mc_truth_Wtaul_id;
   Int_t           mc_truth_Wq1_id;
   Int_t           mc_truth_Wq2_id;
   Int_t           mc_truth_Z_id;
   Int_t           mc_truth_Zl1_id;
   Int_t           mc_truth_Zl2_id;
   Int_t           mc_truth_Ztau1_id;
   Int_t           mc_truth_Ztau2_id;
   Int_t           mc_truth_Ztaul1_id;
   Int_t           mc_truth_Ztaul2_id;
   Int_t           mc_truth_Zq1_id;
   Int_t           mc_truth_Zq2_id;
   Int_t           mc_truth_gammal1_id;
   Int_t           mc_truth_gammal2_id;
   Int_t           mc_truth_gammatau1_id;
   Int_t           mc_truth_gammatau2_id;
   Int_t           mc_truth_gammataul1_id;
   Int_t           mc_truth_gammataul2_id;
   Int_t           mc_truth_gamma_id;
   
   Float_t           mc_truth_t1_pt;
   Float_t           mc_truth_t2_pt;
   Float_t           mc_truth_tb1_pt;
   Float_t           mc_truth_tb2_pt;
   Float_t           mc_truth_tW1_pt;
   Float_t           mc_truth_tW2_pt;
   Float_t           mc_truth_tWl1_pt;
   Float_t           mc_truth_tWl2_pt;
   Float_t           mc_truth_tWtau1_pt;
   Float_t           mc_truth_tWtau2_pt;
   Float_t           mc_truth_tWtaul1_pt;
   Float_t           mc_truth_tWtaul2_pt;
   Float_t           mc_truth_tWq11_pt;
   Float_t           mc_truth_tWq21_pt;
   Float_t           mc_truth_tWq12_pt;
   Float_t           mc_truth_tWq22_pt;
   
   Float_t           mc_truth_t_pt;
   Float_t           mc_truth_tb_pt;
   Float_t           mc_truth_tW_pt;
   Float_t           mc_truth_tWl_pt;
   Float_t           mc_truth_tWq1_pt;
   Float_t           mc_truth_tWq2_pt;
   
   Float_t           mc_truth_W_pt;
   Float_t           mc_truth_Wl_pt;
   Float_t	          mc_truth_Wtau_pt;
   Float_t	          mc_truth_Wtaul_pt;
   Float_t           mc_truth_Wq1_pt;
   Float_t           mc_truth_Wq2_pt;
   Float_t           mc_truth_Z_pt;
   Float_t           mc_truth_Zl1_pt;
   Float_t           mc_truth_Zl2_pt;
   Float_t           mc_truth_Zq1_pt;
   Float_t           mc_truth_Zq2_pt;
   Float_t 	      mc_truth_Ztau1_pt;
   Float_t 	      mc_truth_Ztau2_pt;
   Float_t 	      mc_truth_Ztaul1_pt;
   Float_t 	      mc_truth_Ztaul2_pt;
   Float_t           mc_truth_gammal1_pt;
   Float_t           mc_truth_gammal2_pt;
   Float_t	          mc_truth_gammatau1_pt;
   Float_t	          mc_truth_gammatau2_pt;
   Float_t	          mc_truth_gammataul1_pt;
   Float_t	          mc_truth_gammataul2_pt;
   Float_t	          mc_truth_gamma_pt;
   
   Float_t           mc_truth_h0_pt;
   Float_t           mc_truth_h0W1_pt;
   Float_t           mc_truth_h0Wtau1_pt;
   Float_t           mc_truth_h0Wtaul1_pt;
   Float_t           mc_truth_h0Wl1_pt;
   Float_t           mc_truth_h0Wq11_pt;
   Float_t           mc_truth_h0Wq21_pt;
   Float_t           mc_truth_h0W2_pt;
   Float_t           mc_truth_h0Wl2_pt;
   Float_t           mc_truth_h0Wtau2_pt;
   Float_t           mc_truth_h0Wtaul2_pt;
   Float_t           mc_truth_h0Wq12_pt;
   Float_t           mc_truth_h0Wq22_pt;
   Float_t           mc_truth_h0Z1_pt;
   Float_t           mc_truth_h0Zl11_pt;
   Float_t           mc_truth_h0Zl21_pt;
   Float_t           mc_truth_h0Ztau11_pt;
   Float_t           mc_truth_h0Ztau21_pt;
   Float_t           mc_truth_h0Ztaul11_pt;
   Float_t           mc_truth_h0Ztaul21_pt;
   Float_t           mc_truth_h0Zq11_pt;
   Float_t           mc_truth_h0Zq21_pt;
   Float_t           mc_truth_h0Z2_pt;
   Float_t           mc_truth_h0Zl12_pt;
   Float_t           mc_truth_h0Zl22_pt;
   Float_t           mc_truth_h0Ztau12_pt;
   Float_t           mc_truth_h0Ztau22_pt;
   Float_t           mc_truth_h0Ztaul12_pt;
   Float_t           mc_truth_h0Ztaul22_pt;
   Float_t           mc_truth_h0Zq12_pt;
   Float_t           mc_truth_h0Zq22_pt;
   Float_t           mc_truth_h0tau1_pt;
   Float_t           mc_truth_h0tau2_pt;
   Float_t           mc_truth_h0taul1_pt;
   Float_t           mc_truth_h0taul2_pt;
   Float_t           mc_truth_h0b1_pt;
   Float_t           mc_truth_h0b2_pt;
   
   Float_t           mc_truth_t1_eta;
   Float_t           mc_truth_t2_eta;
   Float_t           mc_truth_tb1_eta;
   Float_t           mc_truth_tb2_eta;
   Float_t           mc_truth_tW1_eta;
   Float_t           mc_truth_tW2_eta;
   Float_t           mc_truth_tWl1_eta;
   Float_t           mc_truth_tWl2_eta;
   Float_t	          mc_truth_tWtau1_eta;
   Float_t	          mc_truth_tWtau2_eta;
   Float_t	          mc_truth_tWtaul1_eta;
   Float_t	          mc_truth_tWtaul2_eta;
   Float_t           mc_truth_tWq11_eta;
   Float_t           mc_truth_tWq21_eta;
   Float_t           mc_truth_tWq12_eta;
   Float_t           mc_truth_tWq22_eta;
   
   Float_t           mc_truth_t_eta;
   Float_t           mc_truth_tb_eta;
   Float_t           mc_truth_tW_eta;
   Float_t           mc_truth_tWl_eta;
   Float_t           mc_truth_tWq1_eta;
   Float_t           mc_truth_tWq2_eta;
   
   Float_t           mc_truth_W_eta;
   Float_t           mc_truth_Wl_eta;
   Float_t	          mc_truth_Wtau_eta;
   Float_t	          mc_truth_Wtaul_eta;
   Float_t           mc_truth_Wq1_eta;
   Float_t           mc_truth_Wq2_eta;
   Float_t           mc_truth_Z_eta;
   Float_t           mc_truth_Zl1_eta;
   Float_t           mc_truth_Zl2_eta;
   Float_t           mc_truth_Zq1_eta;
   Float_t           mc_truth_Zq2_eta;
   Float_t 	      mc_truth_Ztau1_eta;
   Float_t 	      mc_truth_Ztau2_eta;
   Float_t 	      mc_truth_Ztaul1_eta;
   Float_t 	      mc_truth_Ztaul2_eta;
   Float_t           mc_truth_gammal1_eta;
   Float_t           mc_truth_gammal2_eta;
   Float_t	          mc_truth_gammatau1_eta;
   Float_t	          mc_truth_gammatau2_eta;
   Float_t	          mc_truth_gammataul1_eta;
   Float_t	          mc_truth_gammataul2_eta;
   Float_t	          mc_truth_gamma_eta;

   Float_t           mc_truth_h0_eta;
   Float_t           mc_truth_h0W1_eta;
   Float_t           mc_truth_h0Wl1_eta;
   Float_t           mc_truth_h0Wtau1_eta;
   Float_t           mc_truth_h0Wtaul1_eta;
   Float_t           mc_truth_h0Wq11_eta;
   Float_t           mc_truth_h0Wq21_eta;
   Float_t           mc_truth_h0W2_eta;
   Float_t           mc_truth_h0Wl2_eta;
   Float_t           mc_truth_h0Wtau2_eta;
   Float_t           mc_truth_h0Wtaul2_eta;
   Float_t           mc_truth_h0Wq12_eta;
   Float_t           mc_truth_h0Wq22_eta;
   Float_t           mc_truth_h0Z1_eta;
   Float_t           mc_truth_h0Zl11_eta;
   Float_t           mc_truth_h0Zl21_eta;
   Float_t           mc_truth_h0Ztau11_eta;
   Float_t           mc_truth_h0Ztau21_eta;
   Float_t           mc_truth_h0Ztaul11_eta;
   Float_t           mc_truth_h0Ztaul21_eta;
   Float_t           mc_truth_h0Zq11_eta;
   Float_t           mc_truth_h0Zq21_eta;
   Float_t           mc_truth_h0Z2_eta;
   Float_t           mc_truth_h0Zl12_eta;
   Float_t           mc_truth_h0Zl22_eta;
   Float_t           mc_truth_h0Ztau12_eta;
   Float_t           mc_truth_h0Ztau22_eta;
   Float_t           mc_truth_h0Ztaul12_eta;
   Float_t           mc_truth_h0Ztaul22_eta;
   Float_t           mc_truth_h0Zq12_eta;
   Float_t           mc_truth_h0Zq22_eta;
   Float_t           mc_truth_h0tau1_eta;
   Float_t           mc_truth_h0tau2_eta;
   Float_t           mc_truth_h0taul1_eta;
   Float_t           mc_truth_h0taul2_eta;
   Float_t           mc_truth_h0b1_eta;
   Float_t           mc_truth_h0b2_eta;
   
   Float_t           mc_truth_t1_phi;
   Float_t           mc_truth_t2_phi;
   Float_t           mc_truth_tb1_phi;
   Float_t           mc_truth_tb2_phi;
   Float_t           mc_truth_tW1_phi;
   Float_t           mc_truth_tW2_phi;
   Float_t           mc_truth_tWl1_phi;
   Float_t           mc_truth_tWl2_phi;
   Float_t	          mc_truth_tWtau1_phi;
   Float_t	          mc_truth_tWtau2_phi;
   Float_t	          mc_truth_tWtaul1_phi;
   Float_t	          mc_truth_tWtaul2_phi;
   Float_t           mc_truth_tWq11_phi;
   Float_t           mc_truth_tWq21_phi;
   Float_t           mc_truth_tWq12_phi;
   Float_t           mc_truth_tWq22_phi;
   
   Float_t           mc_truth_t_phi;
   Float_t           mc_truth_tb_phi;
   Float_t           mc_truth_tW_phi;
   Float_t           mc_truth_tWl_phi;
   Float_t           mc_truth_tWq1_phi;
   Float_t           mc_truth_tWq2_phi;

   Float_t           mc_truth_W_phi;
   Float_t           mc_truth_Wl_phi;
   Float_t	          mc_truth_Wtau_phi;
   Float_t	          mc_truth_Wtaul_phi;
   Float_t           mc_truth_Wq1_phi;
   Float_t           mc_truth_Wq2_phi;
   Float_t           mc_truth_Z_phi;
   Float_t           mc_truth_Zl1_phi;
   Float_t           mc_truth_Zl2_phi;
   Float_t           mc_truth_Zq1_phi;
   Float_t           mc_truth_Zq2_phi;
   Float_t 	      mc_truth_Ztau1_phi;
   Float_t 	      mc_truth_Ztau2_phi;
   Float_t 	      mc_truth_Ztaul1_phi;
   Float_t 	      mc_truth_Ztaul2_phi;
   Float_t           mc_truth_gammal1_phi;
   Float_t           mc_truth_gammal2_phi;
   Float_t	          mc_truth_gammatau1_phi;
   Float_t	          mc_truth_gammatau2_phi;
   Float_t	          mc_truth_gammataul1_phi;
   Float_t	          mc_truth_gammataul2_phi;
   Float_t	          mc_truth_gamma_phi;
   
   Float_t           mc_truth_h0_phi;
   Float_t           mc_truth_h0W1_phi;
   Float_t           mc_truth_h0Wl1_phi;
   Float_t           mc_truth_h0Wtau1_phi;
   Float_t           mc_truth_h0Wtaul1_phi;
   Float_t           mc_truth_h0Wq11_phi;
   Float_t           mc_truth_h0Wq21_phi;
   Float_t           mc_truth_h0W2_phi;
   Float_t           mc_truth_h0Wl2_phi;
   Float_t           mc_truth_h0Wtau2_phi;
   Float_t           mc_truth_h0Wtaul2_phi;
   Float_t           mc_truth_h0Wq12_phi;
   Float_t           mc_truth_h0Wq22_phi;
   Float_t           mc_truth_h0Z1_phi;
   Float_t           mc_truth_h0Zl11_phi;
   Float_t           mc_truth_h0Zl21_phi;
   Float_t           mc_truth_h0Ztau11_phi;
   Float_t           mc_truth_h0Ztau21_phi;
   Float_t           mc_truth_h0Ztaul11_phi;
   Float_t           mc_truth_h0Ztaul21_phi;
   Float_t           mc_truth_h0Zq11_phi;
   Float_t           mc_truth_h0Zq21_phi;
   Float_t           mc_truth_h0Z2_phi;
   Float_t           mc_truth_h0Zl12_phi;
   Float_t           mc_truth_h0Zl22_phi;
   Float_t           mc_truth_h0Ztau12_phi;
   Float_t           mc_truth_h0Ztau22_phi;
   Float_t           mc_truth_h0Ztaul12_phi;
   Float_t           mc_truth_h0Ztaul22_phi;
   Float_t           mc_truth_h0Zq12_phi;
   Float_t           mc_truth_h0Zq22_phi;
   Float_t           mc_truth_h0tau1_phi;
   Float_t           mc_truth_h0tau2_phi;
   Float_t           mc_truth_h0taul1_phi;
   Float_t           mc_truth_h0taul2_phi;
   Float_t           mc_truth_h0b1_phi;
   Float_t           mc_truth_h0b2_phi;

   Float_t           mc_truth_t1_E;
   Float_t           mc_truth_t2_E;
   Float_t           mc_truth_tb1_E;
   Float_t           mc_truth_tb2_E;
   Float_t           mc_truth_tW1_E;
   Float_t           mc_truth_tW2_E;
   Float_t           mc_truth_tWl1_E;
   Float_t           mc_truth_tWl2_E;
   Float_t	          mc_truth_tWtau1_E;
   Float_t	          mc_truth_tWtau2_E;
   Float_t	          mc_truth_tWtaul1_E;
   Float_t	          mc_truth_tWtaul2_E;
   Float_t           mc_truth_tWq11_E;
   Float_t           mc_truth_tWq21_E;
   Float_t           mc_truth_tWq12_E;
   Float_t           mc_truth_tWq22_E;
   
   Float_t           mc_truth_t_E;
   Float_t           mc_truth_tb_E;
   Float_t           mc_truth_tW_E;
   Float_t           mc_truth_tWl_E;
   Float_t           mc_truth_tWq1_E;
   Float_t           mc_truth_tWq2_E;
   
   Float_t           mc_truth_W_E;
   Float_t           mc_truth_Wl_E;
   Float_t	          mc_truth_Wtau_E;
   Float_t	          mc_truth_Wtaul_E;
   Float_t           mc_truth_Wq1_E;
   Float_t           mc_truth_Wq2_E;
   Float_t           mc_truth_Z_E;
   Float_t           mc_truth_Zl1_E;
   Float_t           mc_truth_Zl2_E;
   Float_t 	      mc_truth_Ztau1_E;
   Float_t 	      mc_truth_Ztau2_E;
   Float_t 	      mc_truth_Ztaul1_E;
   Float_t 	      mc_truth_Ztaul2_E;
   Float_t           mc_truth_Zq1_E;
   Float_t           mc_truth_Zq2_E;
   Float_t           mc_truth_gammal1_E;
   Float_t           mc_truth_gammal2_E;
   Float_t	          mc_truth_gammatau1_E;
   Float_t	          mc_truth_gammatau2_E;
   Float_t	          mc_truth_gammataul1_E;
   Float_t	          mc_truth_gammataul2_E;
   Float_t	          mc_truth_gamma_E;
   
   Float_t           mc_truth_h0_E;
   Float_t           mc_truth_h0W1_E;
   Float_t           mc_truth_h0Wl1_E;
   Float_t           mc_truth_h0Wtau1_E;
   Float_t           mc_truth_h0Wtaul1_E;
   Float_t           mc_truth_h0Wq11_E;
   Float_t           mc_truth_h0Wq21_E;
   Float_t           mc_truth_h0W2_E;
   Float_t           mc_truth_h0Wl2_E;
   Float_t           mc_truth_h0Wtau2_E;
   Float_t           mc_truth_h0Wtaul2_E;
   Float_t           mc_truth_h0Wq12_E;
   Float_t           mc_truth_h0Wq22_E;
   Float_t           mc_truth_h0Z1_E;
   Float_t           mc_truth_h0Zl11_E;
   Float_t           mc_truth_h0Zl21_E;
   Float_t           mc_truth_h0Ztau11_E;
   Float_t           mc_truth_h0Ztau21_E;
   Float_t           mc_truth_h0Ztaul11_E;
   Float_t           mc_truth_h0Ztaul21_E;
   Float_t           mc_truth_h0Zq11_E;
   Float_t           mc_truth_h0Zq21_E;
   Float_t           mc_truth_h0Z2_E;
   Float_t           mc_truth_h0Zl12_E;
   Float_t           mc_truth_h0Zl22_E;
   Float_t           mc_truth_h0Ztau12_E;
   Float_t           mc_truth_h0Ztau22_E;
   Float_t           mc_truth_h0Ztaul12_E;
   Float_t           mc_truth_h0Ztaul22_E;
   Float_t           mc_truth_h0Zq12_E;
   Float_t           mc_truth_h0Zq22_E;
   Float_t           mc_truth_h0tau1_E;
   Float_t           mc_truth_h0tau2_E;
   Float_t           mc_truth_h0taul1_E;
   Float_t           mc_truth_h0taul2_E;
   Float_t           mc_truth_h0b1_E;
   Float_t           mc_truth_h0b2_E;

   Int_t           genJet_n;
   vector<float>   *genJet_pt;
   vector<float>   *genJet_eta;
   vector<float>   *genJet_phi;
   vector<float>   *genJet_m;
   vector<float>   *genJet_E;

   Int_t          triggerobject_n;
   vector<float>  *triggerobject_pt;
   vector<float>  *triggerobject_eta;
   vector<float>  *triggerobject_phi;
   vector<string> *triggerobject_collection;
   vector<int>    *triggerobject_filterIds_n;
   vector<int>    *triggerobject_filterLabels_n;
   vector<int>    *triggerobject_pathNamesAll_n;
   
   vector<string> *triggerobject_pathNamesAll;
   vector<bool>   *triggerobject_pathNamesAll_isL3;
   vector<bool>   *triggerobject_pathNamesAll_isLF;
   vector<bool>   *triggerobject_pathNamesAll_isBoth;
   vector<bool>   *triggerobject_pathNamesAll_isNone;
   
   vector<int>    *triggerobject_filterIds;
   vector<string> *triggerobject_filterLabels;
   
   // List of branches

   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_id;   //!
   TBranch        *b_ev_lumi;   //!
   TBranch        *b_ev_rho;   //!
   
   TBranch        *b_trigger;
   TBranch        *b_trigger_pass;
   TBranch        *b_trigger_name;
   
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_metPuppi_pt;   //!
   TBranch        *b_metPuppi_phi;   //!
   TBranch        *b_metPuppi_sumet;   //!
//   TBranch        *b_met_uncorrectedPt;   //!
//   TBranch        *b_met_uncorrectedPhi;   //!
//   TBranch        *b_met_uncorrectedSumEt;   //!
   TBranch        *b_met_cov00;   //!
   TBranch        *b_met_cov01;   //!
   TBranch        *b_met_cov10;   //!
   TBranch        *b_met_cov11;   //!
   
   TBranch        *b_metNoHF_pt;   //!
   TBranch        *b_metNoHF_phi;   //!
   TBranch        *b_metNoHF_sumet;   //!

   TBranch        *b_nvertex;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_zError;   //!
   
   TBranch        *b_weight_scale_muF0p5;
   TBranch        *b_weight_scale_muF2;
   TBranch        *b_weight_scale_muR0p5;
   TBranch        *b_weight_scale_muR2;
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_id;   //!
   TBranch        *b_mc_f1;   //!
   TBranch        *b_mc_f2;   //!
   TBranch        *b_mc_x1;   //!
   TBranch        *b_mc_x2;   //!
   TBranch        *b_mc_scale;   //!
   TBranch        *b_mc_ptHat;   //!
   TBranch        *b_mc_pdfweights;  //!
   TBranch        *b_mc_pdfweightIds;  //!

   TBranch        *b_mc_pu_intime_NumInt;   //!
   TBranch        *b_mc_pu_trueNumInt;   //!
   TBranch        *b_mc_pu_before_npu;   //!
   TBranch        *b_mc_pu_after_npu;   //!
   TBranch        *b_mc_pu_Npvi;   //!
   TBranch        *b_mc_pu_Nzpositions;   //!
   TBranch        *b_mc_pu_BunchCrossing;   //!
   TBranch        *b_mc_pu_zpositions;   //!
   TBranch        *b_mc_pu_sumpT_lowpT;   //!
   TBranch        *b_mc_pu_sumpT_highpT;   //!
   TBranch        *b_mc_pu_ntrks_lowpT;   //!
   TBranch        *b_mc_pu_ntrks_highpT;   //!
   
   TBranch        *b_el_n;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_m;   //!
   TBranch        *b_el_E;   //!
   TBranch        *b_el_looseCBId;   //!
   TBranch        *b_el_mediumCBId;   //!
   TBranch        *b_el_numberOfLostHits;   //!
   TBranch        *b_el_gsfTrack_PV_dxy;
   TBranch        *b_el_gsfTrack_PV_dz;
   TBranch        *b_el_ip3d;
   TBranch        *b_el_ip3dErr;
   TBranch        *b_el_miniIso;
   TBranch        *b_el_miniIsoTTH;
   TBranch        *b_el_PFRelIso04;
   TBranch        *b_el_conept;   //!
   TBranch        *b_el_ooEmooP;   //!
   TBranch        *b_el_E_postCorr;   //!
   TBranch        *b_el_pt_postCorr;   //!
   
   TBranch        *b_el_id;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_neutralHadronIso;   //!
   TBranch        *b_el_chargedHadronIso;   //!
   TBranch        *b_el_puChargedHadronIso;   //!
   TBranch        *b_el_ecalIso;   //!
   TBranch        *b_el_hcalIso;   //!
   TBranch        *b_el_particleIso;   //!
   TBranch        *b_el_photonIso;   //!
   TBranch        *b_el_trackIso;   //!
   TBranch        *b_el_isLoose;   //!
   TBranch        *b_el_isTight;   //!
   TBranch        *b_el_isRobustLoose;   //!
   TBranch        *b_el_isRobustTight;   //!
   TBranch        *b_el_isRobustHighEnergy;   //!
   TBranch        *b_el_vx;   //!
   TBranch        *b_el_vy;   //!
   TBranch        *b_el_vz;   //!
   TBranch        *b_el_isGsf;   //!
   TBranch        *b_el_dxy;   //!
   TBranch        *b_el_dz;   //!
   TBranch        *b_el_dxyError;   //!
   TBranch        *b_el_dzError;   //!
   TBranch        *b_el_mvaIso;   //!
   TBranch        *b_el_mvaNoIso;   //!
   TBranch        *b_el_mvaPassMedium;   //!
   TBranch        *b_el_mvaPassTight;   //!
   TBranch        *b_el_numberOfHits;   //!
   TBranch        *b_el_pfIso_sumChargedHadronPt;   //!
   TBranch        *b_el_pfIso_sumNeutralHadronEt;   //!
   TBranch        *b_el_pfIso_sumPhotonEt;   //!
   TBranch        *b_el_pfIso_sumPUPt;   //!
   TBranch        *b_el_lepMVA;   //!
   TBranch        *b_el_lepMVA_miniRelIsoCharged;
   TBranch        *b_el_lepMVA_miniRelIsoNeutral;
   TBranch        *b_el_lepMVA_jetPtRelv2;
   TBranch        *b_el_lepMVA_jetPtRatio;
   TBranch        *b_el_lepMVA_jetBTagCSV;
   TBranch        *b_el_lepMVA_jetBTagDeepCSV;
   TBranch        *b_el_lepMVA_sip3d;
   TBranch        *b_el_lepMVA_dxy;
   TBranch        *b_el_lepMVA_dz;
   TBranch        *b_el_lepMVA_mvaId;
   TBranch        *b_el_lepMVA_eta;   //!
   TBranch        *b_el_lepMVA_jetNDauChargedMVASel;   //!
   TBranch        *b_el_isGsfCtfScPixChargeConsistent;   //!
   TBranch        *b_el_isGsfScPixChargeConsistent;
   TBranch        *b_el_passConversionVeto;   //!
   TBranch        *b_el_deltaEtaSuperClusterTrackAtVtx;
   TBranch        *b_el_deltaPhiSuperClusterTrackAtVtx;
   TBranch        *b_el_hadronicOverEm;
   TBranch        *b_el_scleta;
   TBranch        *b_el_hasMatchedConversion;
   TBranch        *b_el_sigmaIetaIeta;
   TBranch        *b_el_full5x5_sigmaIetaIeta;
   TBranch        *b_el_full5x5_sigmaEtaEta;
   TBranch        *b_el_superCluster_eta;
   TBranch        *b_el_correctedEcalEnergy;
   TBranch        *b_el_ecalEnergy;
   TBranch        *b_el_eSuperClusterOverP;
   TBranch        *b_el_trackMomentumError;

   TBranch        *b_el_hasMCMatch;
   TBranch        *b_el_hasChargeMCMatch;
   TBranch        *b_el_hasPhotonMCMatch;
   TBranch        *b_el_gen_pt;
   TBranch        *b_el_gen_eta;
   TBranch        *b_el_gen_phi;
   TBranch        *b_el_gen_m;
   TBranch        *b_el_gen_E;
   TBranch        *b_el_gen_status;
   TBranch        *b_el_gen_id;
   TBranch        *b_el_gen_charge;
   TBranch        *b_el_gen_dr;
   
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_m;   //!
   TBranch        *b_mu_E;   //!
   TBranch        *b_mu_id;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_ip3d;   //!
   TBranch        *b_mu_ip3dErr;   //!
   TBranch        *b_mu_miniIso;
   TBranch        *b_mu_miniIsoTTH;
   TBranch	  *b_mu_PFRelIso04;
   TBranch        *b_mu_isLooseMuon;
   TBranch        *b_mu_conept;   //!

   TBranch        *b_mu_neutralHadronIso;   //!
   TBranch        *b_mu_chargedHadronIso;   //!
   TBranch        *b_mu_ecalIso;   //!
   TBranch        *b_mu_hcalIso;   //!
   TBranch        *b_mu_photonIso;   //!
   TBranch        *b_mu_trackIso;   //!
   TBranch        *b_mu_isGlobalMuon;   //!
   TBranch        *b_mu_isTrackerMuon;   //!
   TBranch        *b_mu_isStandAloneMuon;   //!
   TBranch        *b_mu_isCaloMuon;   //!
   TBranch        *b_mu_isPFMuon;   //!
   TBranch        *b_mu_isMediumMuon;
   TBranch        *b_mu_isTightMuon;   //!
   TBranch        *b_mu_vx;   //!
   TBranch        *b_mu_vy;   //!
   TBranch        *b_mu_vz;   //!
   TBranch        *b_mu_segmentCompatibility;   //!
   TBranch        *b_mu_hasGlobalTrack;   //!
   TBranch        *b_mu_globalTrack_dxy;   //!
   TBranch        *b_mu_globalTrack_dz;   //!
   TBranch        *b_mu_globalTrack_dxyError;   //!
   TBranch        *b_mu_globalTrack_dzError;   //!
   TBranch        *b_mu_globalTrack_normalizedChi2;   //!
   TBranch        *b_mu_combinedQuality_chi2LocalPosition;   //!
   TBranch        *b_mu_combinedQuality_trkKink;   //!
   TBranch        *b_mu_hasInnerTrack;   //!
   TBranch        *b_mu_innerTrack_dxy;   //!
   TBranch        *b_mu_innerTrack_PV_dz;   //!
   TBranch        *b_mu_innerTrack_PV_dxy;   //!
   TBranch        *b_mu_innerTrack_dz;   //!
   TBranch        *b_mu_innerTrack_dxyError;   //!
   TBranch        *b_mu_innerTrack_dzError;   //!
   TBranch        *b_mu_innerTrack_validFraction;   //!
   TBranch        *b_mu_bestTrack_dxy;   //!
   TBranch        *b_mu_bestTrack_PV_dz;   //!
   TBranch        *b_mu_bestTrack_dxyError;   //!
   TBranch        *b_mu_bestTrack_dzError;   //!
   TBranch        *b_mu_bestTrack_pt;
   TBranch        *b_mu_bestTrack_ptError;
   TBranch        *b_mu_numberOfMatches;   //!
   TBranch        *b_mu_numberOfValidMuonHits;   //!
   TBranch        *b_mu_pfIso04_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIso04_sumNeutralHadronEt;   //!
   TBranch        *b_mu_pfIso04_sumPhotonEt;   //!
   TBranch        *b_mu_pfIso04_sumPUPt;   //!
   TBranch        *b_mu_lepMVA;   //!
   TBranch        *b_mu_lepMVA_miniRelIsoCharged;
   TBranch        *b_mu_lepMVA_miniRelIsoNeutral;
   TBranch        *b_mu_lepMVA_jetPtRelv2;
   TBranch        *b_mu_lepMVA_jetPtRatio;
   TBranch        *b_mu_lepMVA_jetBTagCSV;
   TBranch        *b_mu_lepMVA_jetBTagDeepCSV;
   TBranch        *b_mu_lepMVA_sip3d;
   TBranch        *b_mu_lepMVA_dxy;
   TBranch        *b_mu_lepMVA_dz;
   TBranch        *b_mu_lepMVA_mvaId;
   TBranch        *b_mu_lepMVA_eta;   //!
   TBranch        *b_mu_lepMVA_jetNDauChargedMVASel;   //!
   TBranch        *b_mu_innerTrack_pt;   //!
   TBranch        *b_mu_innerTrack_ptError;   //!

   TBranch        *b_mu_hasMCMatch;
   TBranch        *b_mu_hasChargeMCMatch;
   TBranch        *b_mu_gen_pt;
   TBranch        *b_mu_gen_eta;
   TBranch        *b_mu_gen_phi;
   TBranch        *b_mu_gen_m;
   TBranch        *b_mu_gen_E;
   TBranch        *b_mu_gen_status;
   TBranch        *b_mu_gen_id;
   TBranch        *b_mu_gen_charge;
   TBranch        *b_mu_gen_dr;
   
   TBranch   *b_tau_n;
   TBranch   *b_tau_E;
   TBranch   *b_tau_pt;
   TBranch   *b_tau_eta;
   TBranch   *b_tau_phi;
   TBranch   *b_tau_m;
   TBranch   *b_tau_leadingTrackDxy;
   TBranch   *b_tau_leadingTrackDz;
   TBranch   *b_tau_charge;
   TBranch   *b_tau_id;
   TBranch   *b_tau_decayMode;
   TBranch   *b_tau_hasLeadChargedHadrCand;
   TBranch   *b_tau_leadingTrackPt;
   TBranch   *b_tau_leadingTrackCharge;
   TBranch   *b_tau_decayModeFinding;
//   TBranch   *b_tau_decayModeFindingOldDMs;
   TBranch   *b_tau_decayModeFindingNewDMs;
   TBranch   *b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   TBranch   *b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   TBranch   *b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   TBranch   *b_tau_againstElectronVLooseMVA6;
   TBranch   *b_tau_againstElectronLooseMVA6;
   TBranch   *b_tau_againstElectronMediumMVA6;
   TBranch   *b_tau_againstElectronTightMVA6;

   TBranch   *b_tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   TBranch   *b_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   TBranch   *b_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   TBranch   *b_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   TBranch   *b_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
   
   TBranch   *b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   TBranch   *b_tau_chargedIsoPtSum;
   TBranch   *b_tau_neutralIsoPtSum;
   TBranch   *b_tau_puCorrPtSum;
   TBranch   *b_tau_againstMuonLoose3;
   TBranch   *b_tau_againstMuonTight3;
   //AC8X
   TBranch   *b_tau_pfEssential_jet_pt;
   TBranch   *b_tau_pfEssential_jet_eta;
   TBranch   *b_tau_pfEssential_jet_phi;
   TBranch   *b_tau_pfEssential_jet_m;
   TBranch   *b_tau_pfEssential_jetCorr_pt;
   TBranch   *b_tau_pfEssential_jetCorr_eta;
   TBranch   *b_tau_pfEssential_jetCorr_phi;
   TBranch   *b_tau_pfEssential_jetCorr_m;
   TBranch   *b_tau_pfEssential_hasSV;
   TBranch   *b_tau_pfEssential_sv_x;
   TBranch   *b_tau_pfEssential_sv_y;
   TBranch   *b_tau_pfEssential_sv_z;
   TBranch   *b_tau_pfEssential_flightLengthSig;
   TBranch   *b_tau_pfEssential_dxy;
   TBranch   *b_tau_pfEssential_dxy_error;
   TBranch   *b_tau_pfEssential_dxy_Sig;

   TBranch        *b_tau_hasMCMatch;
   TBranch        *b_tau_hasChargeMCMatch;
   TBranch        *b_tau_gen_pt;
   TBranch        *b_tau_gen_eta;
   TBranch        *b_tau_gen_phi;
   TBranch        *b_tau_gen_m;
   TBranch        *b_tau_gen_E;
   TBranch        *b_tau_gen_status;
   TBranch        *b_tau_gen_id;
   TBranch        *b_tau_gen_charge;
   TBranch        *b_tau_gen_dr;
   
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_qgtag;
   TBranch        *b_jet_ntrk;   //!
   TBranch        *b_jet_CSVv2;   //!
   TBranch        *b_jet_cMVAv2;
   TBranch        *b_jet_DeepCSVProbudsg;
   TBranch        *b_jet_DeepCSVProbb;
   TBranch        *b_jet_DeepCSVProbbb;
   TBranch        *b_jet_DeepCSVProbc;
   TBranch        *b_jet_DeepCSVProbcc;
   TBranch        *b_jet_tightJetID;
   TBranch        *b_jet_partonFlavour;   //!
   TBranch        *b_jet_hadronFlavour;   //!
   TBranch        *b_jet_neutralHadronEnergy;   //!
   TBranch        *b_jet_neutralEmEnergy;   //!
   TBranch        *b_jet_neutralEmEnergyFraction;   //!
   TBranch        *b_jet_chargedHadronEnergy;   //!
   TBranch        *b_jet_chargedEmEnergy;   //!
   TBranch        *b_jet_electronEnergy;   //!
   TBranch        *b_jet_muonEnergy;   //!
   TBranch        *b_jet_photonEnergy;   //!
   TBranch        *b_jet_genJet_pt;   //!
   TBranch        *b_jet_genJet_eta;   //!
   TBranch        *b_jet_genJet_phi;   //!
   TBranch        *b_jet_genJet_m;   //!
   TBranch        *b_jet_genJet_E;   //!
   TBranch        *b_jet_genJet_status;   //!
   TBranch        *b_jet_genJet_id;   //!
   TBranch        *b_jet_genParton_pt; //!
   TBranch        *b_jet_genParton_eta; //!
   TBranch        *b_jet_genParton_phi; //!
   //TBranch      *b_jet_genParton_m; //!
   TBranch        *b_jet_genParton_E; //!
   //TBranch      *b_jet_genParton_status; //!
   TBranch        *b_jet_genParton_id; //!
   TBranch        *b_jet_pileupJetId;   //!
   
   TBranch        *b_gen_n;
   TBranch        *b_gen_pt;
   TBranch        *b_gen_eta;
   TBranch        *b_gen_phi;
   TBranch        *b_gen_m;
   TBranch        *b_gen_id;
   TBranch        *b_gen_status;
   TBranch        *b_gen_index;
   TBranch        *b_gen_mother_index;
   
   TBranch        *b_gen_PVz;
   
   TBranch        *b_metGen_px;
   TBranch        *b_metGen_py;
   TBranch        *b_metGen_pt;
   TBranch        *b_metGen_phi;
   TBranch        *b_metGen_sumet;
   TBranch        *b_metGen_MuonEt;
   
   TBranch        *b_mc_truth_h0_id;
   TBranch        *b_mc_truth_h0W1_id;
   TBranch        *b_mc_truth_h0Wl1_id;
   TBranch        *b_mc_truth_h0Wtau1_id;
   TBranch        *b_mc_truth_h0Wtaul1_id;
   TBranch        *b_mc_truth_h0Wq11_id;
   TBranch        *b_mc_truth_h0Wq21_id;
   TBranch        *b_mc_truth_h0W2_id;
   TBranch        *b_mc_truth_h0Wl2_id;
   TBranch        *b_mc_truth_h0Wtau2_id;
   TBranch        *b_mc_truth_h0Wtaul2_id;
   TBranch        *b_mc_truth_h0Wq12_id;
   TBranch        *b_mc_truth_h0Wq22_id;
   TBranch        *b_mc_truth_h0Z1_id;
   TBranch        *b_mc_truth_h0Zl11_id;
   TBranch        *b_mc_truth_h0Zl21_id;
   TBranch        *b_mc_truth_h0Ztau11_id;
   TBranch        *b_mc_truth_h0Ztau21_id;
   TBranch        *b_mc_truth_h0Ztaul11_id;
   TBranch        *b_mc_truth_h0Ztaul21_id;
   TBranch        *b_mc_truth_h0Zq11_id;
   TBranch        *b_mc_truth_h0Zq21_id;
   TBranch        *b_mc_truth_h0Z2_id;
   TBranch        *b_mc_truth_h0Zl12_id;
   TBranch        *b_mc_truth_h0Zl22_id;
   TBranch        *b_mc_truth_h0Ztau12_id;
   TBranch        *b_mc_truth_h0Ztau22_id;
   TBranch        *b_mc_truth_h0Ztaul12_id;
   TBranch        *b_mc_truth_h0Ztaul22_id;
   
   TBranch        *b_mc_truth_h0Zq12_id;
   TBranch        *b_mc_truth_h0Zq22_id;
   TBranch        *b_mc_truth_h0tau1_id;
   TBranch        *b_mc_truth_h0tau2_id;
   TBranch        *b_mc_truth_h0taul1_id;
   TBranch        *b_mc_truth_h0taul2_id;
   TBranch        *b_mc_truth_h0b1_id;
   TBranch        *b_mc_truth_h0b2_id;
   
   TBranch        *b_mc_truth_t1_id;
   TBranch        *b_mc_truth_t2_id;
   TBranch        *b_mc_truth_tb1_id;
   TBranch        *b_mc_truth_tb2_id;
   TBranch        *b_mc_truth_tW1_id;
   TBranch        *b_mc_truth_tW2_id;
   TBranch        *b_mc_truth_tWl1_id;
   TBranch        *b_mc_truth_tWl2_id;
   TBranch        *b_mc_truth_tWtau1_id;
   TBranch        *b_mc_truth_tWtau2_id;
   TBranch        *b_mc_truth_tWtaul1_id;
   TBranch        *b_mc_truth_tWtaul2_id;
   TBranch        *b_mc_truth_tWq11_id;
   TBranch        *b_mc_truth_tWq21_id;
   TBranch        *b_mc_truth_tWq12_id;
   TBranch        *b_mc_truth_tWq22_id;

   TBranch        *b_mc_truth_t_id;
   TBranch        *b_mc_truth_tb_id;
   TBranch        *b_mc_truth_tW_id;
   TBranch        *b_mc_truth_tWl_id;
   TBranch        *b_mc_truth_tWq1_id;
   TBranch         *b_mc_truth_tWq2_id;
   
   TBranch        *b_mc_truth_W_id;
   TBranch        *b_mc_truth_Wl_id;
   TBranch        *b_mc_truth_Wtau_id;
   TBranch        *b_mc_truth_Wtaul_id;
   TBranch        *b_mc_truth_Wq1_id;
   TBranch        *b_mc_truth_Wq2_id;
   TBranch        *b_mc_truth_Z_id;
   TBranch        *b_mc_truth_Zl1_id;
   TBranch        *b_mc_truth_Zl2_id;
   TBranch        *b_mc_truth_Zq1_id;
   TBranch        *b_mc_truth_Zq2_id;
   TBranch        *b_mc_truth_Ztau1_id;
   TBranch        *b_mc_truth_Ztau2_id;
   TBranch        *b_mc_truth_Ztaul1_id;
   TBranch        *b_mc_truth_Ztaul2_id;
   TBranch        *b_mc_truth_gammal1_id;
   TBranch        *b_mc_truth_gammal2_id;
   TBranch        *b_mc_truth_gammatau1_id;
   TBranch        *b_mc_truth_gammatau2_id;
   TBranch        *b_mc_truth_gammataul1_id;
   TBranch        *b_mc_truth_gammataul2_id;
   TBranch        *b_mc_truth_gamma_id;
   
   TBranch        *b_mc_truth_h0_pt;
   TBranch        *b_mc_truth_h0W1_pt;
   TBranch        *b_mc_truth_h0Wl1_pt;
   TBranch        *b_mc_truth_h0Wtau1_pt;
   TBranch        *b_mc_truth_h0Wtaul1_pt;
   TBranch        *b_mc_truth_h0Wq11_pt;
   TBranch        *b_mc_truth_h0Wq21_pt;
   TBranch        *b_mc_truth_h0W2_pt;
   TBranch        *b_mc_truth_h0Wl2_pt;
   TBranch        *b_mc_truth_h0Wtau2_pt;
   TBranch        *b_mc_truth_h0Wtaul2_pt;
   TBranch        *b_mc_truth_h0Wq12_pt;
   TBranch        *b_mc_truth_h0Wq22_pt;
   TBranch        *b_mc_truth_h0Z1_pt;
   TBranch        *b_mc_truth_h0Zl11_pt;
   TBranch        *b_mc_truth_h0Zl21_pt;
   TBranch        *b_mc_truth_h0Ztau11_pt;
   TBranch        *b_mc_truth_h0Ztau21_pt;
   TBranch        *b_mc_truth_h0Ztaul11_pt;
   TBranch        *b_mc_truth_h0Ztaul21_pt;
   TBranch        *b_mc_truth_h0Zq11_pt;
   TBranch        *b_mc_truth_h0Zq21_pt;
   TBranch        *b_mc_truth_h0Z2_pt;
   TBranch        *b_mc_truth_h0Zl12_pt;
   TBranch        *b_mc_truth_h0Zl22_pt;
   TBranch        *b_mc_truth_h0Ztau12_pt;
   TBranch        *b_mc_truth_h0Ztau22_pt;
   TBranch        *b_mc_truth_h0Ztaul12_pt;
   TBranch        *b_mc_truth_h0Ztaul22_pt;
   TBranch        *b_mc_truth_h0Zq12_pt;
   TBranch        *b_mc_truth_h0Zq22_pt;
   TBranch        *b_mc_truth_h0tau1_pt;
   TBranch        *b_mc_truth_h0tau2_pt;
   TBranch        *b_mc_truth_h0taul1_pt;
   TBranch        *b_mc_truth_h0taul2_pt;
   TBranch        *b_mc_truth_h0b1_pt;
   TBranch        *b_mc_truth_h0b2_pt;
   
   TBranch        *b_mc_truth_t1_pt;
   TBranch        *b_mc_truth_t2_pt;
   TBranch        *b_mc_truth_tb1_pt;
   TBranch        *b_mc_truth_tb2_pt;
   TBranch        *b_mc_truth_tW1_pt;
   TBranch        *b_mc_truth_tW2_pt;
   TBranch        *b_mc_truth_tWl1_pt;
   TBranch        *b_mc_truth_tWl2_pt;
   TBranch        *b_mc_truth_tWtau1_pt;
   TBranch        *b_mc_truth_tWtau2_pt;
   TBranch        *b_mc_truth_tWtaul1_pt;
   TBranch        *b_mc_truth_tWtaul2_pt;
   TBranch        *b_mc_truth_tWq11_pt;
   TBranch        *b_mc_truth_tWq21_pt;
   TBranch        *b_mc_truth_tWq12_pt;
   TBranch        *b_mc_truth_tWq22_pt;

   TBranch        *b_mc_truth_t_pt;
   TBranch        *b_mc_truth_tb_pt;
   TBranch        *b_mc_truth_tW_pt;
   TBranch        *b_mc_truth_tWl_pt;
   TBranch        *b_mc_truth_tWq1_pt;
   TBranch        *b_mc_truth_tWq2_pt;
   
   TBranch        *b_mc_truth_W_pt;
   TBranch        *b_mc_truth_Wl_pt;
   TBranch        *b_mc_truth_Wtau_pt;
   TBranch        *b_mc_truth_Wtaul_pt;
   TBranch        *b_mc_truth_Wq1_pt;
   TBranch        *b_mc_truth_Wq2_pt;
   TBranch        *b_mc_truth_Z_pt;
   TBranch        *b_mc_truth_Zl1_pt;
   TBranch        *b_mc_truth_Zl2_pt;
   TBranch        *b_mc_truth_Zq1_pt;
   TBranch        *b_mc_truth_Zq2_pt;
   TBranch        *b_mc_truth_Ztau1_pt;
   TBranch        *b_mc_truth_Ztau2_pt;
   TBranch        *b_mc_truth_Ztaul1_pt;
   TBranch        *b_mc_truth_Ztaul2_pt;
   TBranch        *b_mc_truth_gammal1_pt;
   TBranch        *b_mc_truth_gammal2_pt;
   TBranch        *b_mc_truth_gammatau1_pt;
   TBranch        *b_mc_truth_gammatau2_pt;
   TBranch        *b_mc_truth_gammataul1_pt;
   TBranch        *b_mc_truth_gammataul2_pt;
   TBranch        *b_mc_truth_gamma_pt;
   
   TBranch        *b_mc_truth_h0_eta;
   TBranch        *b_mc_truth_h0W1_eta;
   TBranch        *b_mc_truth_h0Wl1_eta;
   TBranch        *b_mc_truth_h0Wtau1_eta;
   TBranch        *b_mc_truth_h0Wtaul1_eta;
   TBranch        *b_mc_truth_h0Wq11_eta;
   TBranch        *b_mc_truth_h0Wq21_eta;
   TBranch        *b_mc_truth_h0W2_eta;
   TBranch        *b_mc_truth_h0Wl2_eta;
   TBranch        *b_mc_truth_h0Wtau2_eta;
   TBranch        *b_mc_truth_h0Wtaul2_eta;
   TBranch        *b_mc_truth_h0Wq12_eta;
   TBranch        *b_mc_truth_h0Wq22_eta;
   TBranch        *b_mc_truth_h0Z1_eta;
   TBranch        *b_mc_truth_h0Zl11_eta;
   TBranch        *b_mc_truth_h0Zl21_eta;
   TBranch        *b_mc_truth_h0Ztau11_eta;
   TBranch        *b_mc_truth_h0Ztau21_eta;
   TBranch        *b_mc_truth_h0Ztaul11_eta;
   TBranch        *b_mc_truth_h0Ztaul21_eta;
   TBranch        *b_mc_truth_h0Zq11_eta;
   TBranch        *b_mc_truth_h0Zq21_eta;
   TBranch        *b_mc_truth_h0Z2_eta;
   TBranch        *b_mc_truth_h0Zl12_eta;
   TBranch        *b_mc_truth_h0Zl22_eta;
   TBranch        *b_mc_truth_h0Ztau12_eta;
   TBranch        *b_mc_truth_h0Ztau22_eta;
   TBranch        *b_mc_truth_h0Ztaul12_eta;
   TBranch        *b_mc_truth_h0Ztaul22_eta;
   TBranch        *b_mc_truth_h0Zq12_eta;
   TBranch        *b_mc_truth_h0Zq22_eta;
   TBranch        *b_mc_truth_h0tau1_eta;
   TBranch        *b_mc_truth_h0tau2_eta;
   TBranch        *b_mc_truth_h0taul1_eta;
   TBranch        *b_mc_truth_h0taul2_eta;
   TBranch        *b_mc_truth_h0b1_eta;
   TBranch        *b_mc_truth_h0b2_eta;
   
   TBranch        *b_mc_truth_t1_eta;
   TBranch        *b_mc_truth_t2_eta;
   TBranch        *b_mc_truth_tb1_eta;
   TBranch        *b_mc_truth_tb2_eta;
   TBranch        *b_mc_truth_tW1_eta;
   TBranch        *b_mc_truth_tW2_eta;
   TBranch        *b_mc_truth_tWl1_eta;
   TBranch        *b_mc_truth_tWl2_eta;
   TBranch        *b_mc_truth_tWtau1_eta;
   TBranch        *b_mc_truth_tWtau2_eta;
   TBranch        *b_mc_truth_tWtaul1_eta;
   TBranch        *b_mc_truth_tWtaul2_eta;
   TBranch        *b_mc_truth_tWq11_eta;
   TBranch        *b_mc_truth_tWq21_eta;
   TBranch        *b_mc_truth_tWq12_eta;
   TBranch        *b_mc_truth_tWq22_eta;
   
   TBranch        *b_mc_truth_t_eta;
   TBranch        *b_mc_truth_tb_eta;
   TBranch        *b_mc_truth_tW_eta;
   TBranch        *b_mc_truth_tWl_eta;
   TBranch        *b_mc_truth_tWq1_eta;
   TBranch        *b_mc_truth_tWq2_eta;
   
   TBranch        *b_mc_truth_W_eta;
   TBranch        *b_mc_truth_Wl_eta;
   TBranch        *b_mc_truth_Wtau_eta;
   TBranch        *b_mc_truth_Wtaul_eta;
   TBranch        *b_mc_truth_Wq1_eta;
   TBranch        *b_mc_truth_Wq2_eta;
   TBranch        *b_mc_truth_Z_eta;
   TBranch        *b_mc_truth_Zl1_eta;
   TBranch        *b_mc_truth_Zl2_eta;
   TBranch        *b_mc_truth_Zq1_eta;
   TBranch        *b_mc_truth_Zq2_eta;
   TBranch        *b_mc_truth_Ztau1_eta;
   TBranch        *b_mc_truth_Ztau2_eta;
   TBranch        *b_mc_truth_Ztaul1_eta;
   TBranch        *b_mc_truth_Ztaul2_eta;
   TBranch        *b_mc_truth_gammal1_eta;
   TBranch        *b_mc_truth_gammal2_eta;
   TBranch        *b_mc_truth_gammatau1_eta;
   TBranch        *b_mc_truth_gammatau2_eta;
   TBranch        *b_mc_truth_gammataul1_eta;
   TBranch        *b_mc_truth_gammataul2_eta;
   TBranch        *b_mc_truth_gamma_eta;
   
   TBranch        *b_mc_truth_h0_phi;
   TBranch        *b_mc_truth_h0W1_phi;
   TBranch        *b_mc_truth_h0Wl1_phi;
   TBranch        *b_mc_truth_h0Wtau1_phi;
   TBranch        *b_mc_truth_h0Wtaul1_phi;
   TBranch        *b_mc_truth_h0Wq11_phi;
   TBranch        *b_mc_truth_h0Wq21_phi;
   TBranch        *b_mc_truth_h0W2_phi;
   TBranch        *b_mc_truth_h0Wl2_phi;
   TBranch        *b_mc_truth_h0Wtau2_phi;
   TBranch        *b_mc_truth_h0Wtaul2_phi;
   TBranch        *b_mc_truth_h0Wq12_phi;
   TBranch        *b_mc_truth_h0Wq22_phi;
   TBranch        *b_mc_truth_h0Z1_phi;
   TBranch        *b_mc_truth_h0Zl11_phi;
   TBranch        *b_mc_truth_h0Zl21_phi;
   TBranch        *b_mc_truth_h0Ztau11_phi;
   TBranch        *b_mc_truth_h0Ztau21_phi;
   TBranch        *b_mc_truth_h0Ztaul11_phi;
   TBranch        *b_mc_truth_h0Ztaul21_phi;
   TBranch        *b_mc_truth_h0Zq11_phi;
   TBranch        *b_mc_truth_h0Zq21_phi;
   TBranch        *b_mc_truth_h0Z2_phi;
   TBranch        *b_mc_truth_h0Zl12_phi;
   TBranch        *b_mc_truth_h0Zl22_phi;
   TBranch        *b_mc_truth_h0Ztau12_phi;
   TBranch        *b_mc_truth_h0Ztau22_phi;
   TBranch        *b_mc_truth_h0Ztaul12_phi;
   TBranch        *b_mc_truth_h0Ztaul22_phi;
   TBranch        *b_mc_truth_h0Zq12_phi;
   TBranch        *b_mc_truth_h0Zq22_phi;
   TBranch        *b_mc_truth_h0tau1_phi;
   TBranch        *b_mc_truth_h0tau2_phi;
   TBranch        *b_mc_truth_h0taul1_phi;
   TBranch        *b_mc_truth_h0taul2_phi;
   TBranch        *b_mc_truth_h0b1_phi;
   TBranch        *b_mc_truth_h0b2_phi;
   
   TBranch        *b_mc_truth_t1_phi;
   TBranch        *b_mc_truth_t2_phi;
   TBranch        *b_mc_truth_tb1_phi;
   TBranch        *b_mc_truth_tb2_phi;
   TBranch        *b_mc_truth_tW1_phi;
   TBranch        *b_mc_truth_tW2_phi;
   TBranch        *b_mc_truth_tWl1_phi;
   TBranch        *b_mc_truth_tWl2_phi;
   TBranch        *b_mc_truth_tWtau1_phi;
   TBranch        *b_mc_truth_tWtau2_phi;
   TBranch        *b_mc_truth_tWtaul1_phi;
   TBranch        *b_mc_truth_tWtaul2_phi;
   TBranch        *b_mc_truth_tWq11_phi;
   TBranch        *b_mc_truth_tWq21_phi;
   TBranch        *b_mc_truth_tWq12_phi;
   TBranch        *b_mc_truth_tWq22_phi;
   
   TBranch        *b_mc_truth_t_phi;
   TBranch        *b_mc_truth_tb_phi;
   TBranch        *b_mc_truth_tW_phi;
   TBranch        *b_mc_truth_tWl_phi;
   TBranch        *b_mc_truth_tWq1_phi;
   TBranch        *b_mc_truth_tWq2_phi;
   
   TBranch        *b_mc_truth_W_phi;
   TBranch        *b_mc_truth_Wl_phi;
   TBranch        *b_mc_truth_Wtau_phi;
   TBranch        *b_mc_truth_Wtaul_phi;
   TBranch        *b_mc_truth_Wq1_phi;
   TBranch        *b_mc_truth_Wq2_phi;
   TBranch        *b_mc_truth_Z_phi;
   TBranch        *b_mc_truth_Zl1_phi;
   TBranch        *b_mc_truth_Zl2_phi;
   TBranch        *b_mc_truth_Zq1_phi;
   TBranch        *b_mc_truth_Zq2_phi;
   TBranch        *b_mc_truth_Ztau1_phi;
   TBranch        *b_mc_truth_Ztau2_phi;
   TBranch        *b_mc_truth_Ztaul1_phi;
   TBranch        *b_mc_truth_Ztaul2_phi;
   TBranch        *b_mc_truth_gammal1_phi;
   TBranch        *b_mc_truth_gammal2_phi;
   TBranch        *b_mc_truth_gammatau1_phi;
   TBranch        *b_mc_truth_gammatau2_phi;
   TBranch        *b_mc_truth_gammataul1_phi;
   TBranch        *b_mc_truth_gammataul2_phi;
   TBranch        *b_mc_truth_gamma_phi;
   
   TBranch        *b_mc_truth_h0_E;
   TBranch        *b_mc_truth_h0W1_E;
   TBranch        *b_mc_truth_h0Wl1_E;
   TBranch        *b_mc_truth_h0Wtau1_E;
   TBranch        *b_mc_truth_h0Wtaul1_E;
   TBranch        *b_mc_truth_h0Wq11_E;
   TBranch        *b_mc_truth_h0Wq21_E;
   TBranch        *b_mc_truth_h0W2_E;
   TBranch        *b_mc_truth_h0Wl2_E;
   TBranch        *b_mc_truth_h0Wtau2_E;
   TBranch        *b_mc_truth_h0Wtaul2_E;
   TBranch        *b_mc_truth_h0Wq12_E;
   TBranch        *b_mc_truth_h0Wq22_E;
   TBranch        *b_mc_truth_h0Z1_E;
   TBranch        *b_mc_truth_h0Zl11_E;
   TBranch        *b_mc_truth_h0Zl21_E;
   TBranch        *b_mc_truth_h0Ztau11_E;
   TBranch        *b_mc_truth_h0Ztau21_E;
   TBranch        *b_mc_truth_h0Ztaul11_E;
   TBranch        *b_mc_truth_h0Ztaul21_E;
   TBranch        *b_mc_truth_h0Zq11_E;
   TBranch        *b_mc_truth_h0Zq21_E;
   TBranch        *b_mc_truth_h0Z2_E;
   TBranch        *b_mc_truth_h0Zl12_E;
   TBranch        *b_mc_truth_h0Zl22_E;
   TBranch        *b_mc_truth_h0Ztau12_E;
   TBranch        *b_mc_truth_h0Ztau22_E;
   TBranch        *b_mc_truth_h0Ztaul12_E;
   TBranch        *b_mc_truth_h0Ztaul22_E;
   TBranch        *b_mc_truth_h0Zq12_E;
   TBranch        *b_mc_truth_h0Zq22_E;
   TBranch        *b_mc_truth_h0tau1_E;
   TBranch        *b_mc_truth_h0tau2_E;
   TBranch        *b_mc_truth_h0taul1_E;
   TBranch        *b_mc_truth_h0taul2_E;
   TBranch        *b_mc_truth_h0b1_E;
   TBranch        *b_mc_truth_h0b2_E;

   TBranch        *b_mc_truth_t1_E;
   TBranch        *b_mc_truth_t2_E;
   TBranch        *b_mc_truth_tb1_E;
   TBranch        *b_mc_truth_tb2_E;
   TBranch        *b_mc_truth_tW1_E;
   TBranch        *b_mc_truth_tW2_E;
   TBranch        *b_mc_truth_tWl1_E;
   TBranch        *b_mc_truth_tWl2_E;
   TBranch        *b_mc_truth_tWtau1_E;
   TBranch        *b_mc_truth_tWtau2_E;
   TBranch        *b_mc_truth_tWtaul1_E;
   TBranch        *b_mc_truth_tWtaul2_E;
   TBranch        *b_mc_truth_tWq11_E;
   TBranch        *b_mc_truth_tWq21_E;
   TBranch        *b_mc_truth_tWq12_E;
   TBranch        *b_mc_truth_tWq22_E;
   
   TBranch        *b_mc_truth_t_E;
   TBranch        *b_mc_truth_tb_E;
   TBranch        *b_mc_truth_tW_E;
   TBranch        *b_mc_truth_tWl_E;
   TBranch        *b_mc_truth_tWq1_E;
   TBranch        *b_mc_truth_tWq2_E;
   
   TBranch        *b_mc_truth_W_E;
   TBranch        *b_mc_truth_Wl_E;
   TBranch        *b_mc_truth_Wtau_E;
   TBranch        *b_mc_truth_Wtaul_E;
   TBranch        *b_mc_truth_Wq1_E;
   TBranch        *b_mc_truth_Wq2_E;
   TBranch        *b_mc_truth_Z_E;
   TBranch        *b_mc_truth_Zl1_E;
   TBranch        *b_mc_truth_Zl2_E;
   TBranch        *b_mc_truth_Zq1_E;
   TBranch        *b_mc_truth_Zq2_E;
   TBranch        *b_mc_truth_Ztau1_E;
   TBranch        *b_mc_truth_Ztau2_E;
   TBranch        *b_mc_truth_Ztaul1_E;
   TBranch        *b_mc_truth_Ztaul2_E;
   TBranch        *b_mc_truth_gammal1_E;
   TBranch        *b_mc_truth_gammal2_E;
   TBranch        *b_mc_truth_gammatau1_E;
   TBranch        *b_mc_truth_gammatau2_E;
   TBranch        *b_mc_truth_gammataul1_E;
   TBranch        *b_mc_truth_gammataul2_E;
   TBranch        *b_mc_truth_gamma_E;
   
   TBranch        *b_genJet_n;
   TBranch        *b_genJet_pt;
   TBranch        *b_genJet_eta;
   TBranch        *b_genJet_phi;
   TBranch        *b_genJet_m;
   TBranch        *b_genJet_E;
   
   TBranch        *b_triggerobject_n;
   TBranch        *b_triggerobject_pt;
   TBranch        *b_triggerobject_eta;
   TBranch        *b_triggerobject_phi;
   TBranch        *b_triggerobject_collection;
   
   TBranch        *b_triggerobject_filterIds_n;
   TBranch        *b_triggerobject_filterLabels_n;
   TBranch        *b_triggerobject_pathNamesAll_n;
   
   TBranch        *b_triggerobject_pathNamesAll;
   TBranch        *b_triggerobject_pathNamesAll_isL3;
   TBranch        *b_triggerobject_pathNamesAll_isLF;
   TBranch        *b_triggerobject_pathNamesAll_isBoth;
   TBranch        *b_triggerobject_pathNamesAll_isNone;
   
   TBranch        *b_triggerobject_filterIds;
   TBranch        *b_triggerobject_filterLabels;
   
   Tree(TChain *tree=0,std::string fname="output.root",std::string treename="FlatTree/tree");
   virtual ~Tree();
   virtual Int_t    GetEntry(Long64_t entry);
   //   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *ch);
   virtual void     registerInputBranches(TChain *ch);
};

#endif
