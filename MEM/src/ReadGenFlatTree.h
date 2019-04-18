#ifndef READGENFLATTREE_H
#define READGENFLATTREE_H


class ReadGenFlatTree {

  public:

  string InputFileName;
  TFile* fInput;
  TTree* tInput;

//  string OutputFileName;
  TFile* fOutput;
  TTree* tOutput;

  ReadGenFlatTree();
  ~ReadGenFlatTree();
  void InitializeDryRun(string);
  void InitializeMEMRun(string);
  void FillGenMultilepton(Long64_t, MultiLepton*);
  int ApplyGenSelection(Long64_t, MultiLepton*);
  void WriteMultilepton(MultiLepton*);
  void ReadMultilepton(Long64_t, MultiLepton*);

  Float_t mc_weight;
  Float_t weight;
  Float_t PV_weight;
  Float_t weight_scale_muF0p5, weight_scale_muF2, weight_scale_muR0p5, weight_scale_muR2;
  Float_t weight_csv_down, weight_csv_up;

  Int_t mc_truth_h0_id;
  TLorentzVector* mc_truth_h0_p4;
  Int_t mc_truth_h0Wl1_id;
  TLorentzVector* mc_truth_h0Wl1_p4;
  Int_t mc_truth_h0Wl2_id;
  TLorentzVector* mc_truth_h0Wl2_p4;
  TLorentzVector* mc_truth_h0Wnu1_p4;
  TLorentzVector* mc_truth_h0Wnu2_p4;
  Int_t mc_truth_h0Wq11_id;
  TLorentzVector* mc_truth_h0Wq11_p4;
  Int_t mc_truth_h0Wq12_id;
  TLorentzVector* mc_truth_h0Wq12_p4;
  Int_t mc_truth_h0Wq21_id;
  TLorentzVector* mc_truth_h0Wq21_p4;
  Int_t mc_truth_h0Wq22_id;
  TLorentzVector* mc_truth_h0Wq22_p4;
  Int_t mc_truth_t1_id;
  TLorentzVector* mc_truth_t1_p4;
  Int_t mc_truth_t2_id;
  TLorentzVector* mc_truth_t2_p4;
  Int_t mc_truth_tb1_id;
  TLorentzVector* mc_truth_tb1_p4;
  Int_t mc_truth_tb2_id;
  TLorentzVector* mc_truth_tb2_p4;
  Int_t mc_truth_tWl1_id;
  TLorentzVector* mc_truth_tWl1_p4;
  TLorentzVector* mc_truth_tWnu1_p4;
  Int_t mc_truth_tWl2_id;
  TLorentzVector* mc_truth_tWl2_p4;
  TLorentzVector* mc_truth_tWnu2_p4;
  Int_t mc_truth_tWq11_id;
  TLorentzVector* mc_truth_tWq11_p4;
  Int_t mc_truth_tWq21_id;
  TLorentzVector* mc_truth_tWq21_p4;
  Int_t mc_truth_tWq12_id;
  TLorentzVector* mc_truth_tWq12_p4;
  Int_t mc_truth_tWq22_id;
  TLorentzVector* mc_truth_tWq22_p4;
  Int_t mc_truth_Z_id;
  TLorentzVector* mc_truth_Z_p4;
  Int_t mc_truth_Zl1_id;
  TLorentzVector* mc_truth_Zl1_p4;
  Int_t mc_truth_Zl2_id;
  TLorentzVector* mc_truth_Zl2_p4;
  Int_t mc_truth_W_id;
  TLorentzVector*  mc_truth_W_p4;
  Int_t mc_truth_Wl_id;
  TLorentzVector*  mc_truth_Wl_p4;
  Int_t genJet_n;
  Int_t mc_truth_gammal1_id;
  TLorentzVector*  mc_truth_gammal1_p4;
  Int_t mc_truth_gammal2_id;
  TLorentzVector*  mc_truth_gammal2_p4;
  std::vector<float>* genJet_pt;
  std::vector<float>* genJet_eta;
  std::vector<float>* genJet_phi;
  std::vector<float>* genJet_E;

  Float_t nJet25_Recl;
  Float_t max_Lep_eta, MT_met_lep1, mindr_lep1_jet, mindr_lep2_jet, LepGood_conePt0, LepGood_conePt1,  met, avg_dr_jet, mhtJet25_Recl;
  Float_t signal_2lss_TT_MVA;
  Float_t signal_2lss_TTV_MVA;
  Float_t signal_3l_TT_MVA;
  Float_t signal_3l_TTV_MVA;

  TBranch* b_mc_event;
  TBranch* b_mc_weight;
  TBranch* b_weight;
  TBranch* b_PV_weight;
  TBranch* b_weight_scale_muF0p5;
  TBranch* b_weight_scale_muF2;
  TBranch* b_weight_scale_muR0p5;
  TBranch* b_weight_scale_muR2;
  TBranch* b_weight_csv_down;
  TBranch* b_weight_csv_up;

  TBranch* b_mc_truth_h0_id;
  TBranch* b_mc_truth_h0_p4;
  TBranch* b_mc_truth_h0Wl1_id;
  TBranch* b_mc_truth_h0Wl1_p4;
  TBranch* b_mc_truth_h0Wl2_id;
  TBranch* b_mc_truth_h0Wl2_p4;
  TBranch* b_mc_truth_h0Wnu1_p4;
  TBranch* b_mc_truth_h0Wnu2_p4;
  TBranch* b_mc_truth_h0Wq11_id;
  TBranch* b_mc_truth_h0Wq11_p4;
  TBranch* b_mc_truth_h0Wq12_id;
  TBranch* b_mc_truth_h0Wq12_p4;
  TBranch* b_mc_truth_h0Wq21_id;
  TBranch* b_mc_truth_h0Wq21_p4;
  TBranch* b_mc_truth_h0Wq22_id;
  TBranch* b_mc_truth_h0Wq22_p4;
  TBranch* b_mc_truth_t1_id;
  TBranch* b_mc_truth_t1_p4;
  TBranch* b_mc_truth_t2_id;
  TBranch* b_mc_truth_t2_p4;
  TBranch* b_mc_truth_tb1_id;
  TBranch* b_mc_truth_tb1_p4;
  TBranch* b_mc_truth_tb2_id;
  TBranch* b_mc_truth_tb2_p4;
  TBranch* b_mc_truth_tWl1_id;
  TBranch* b_mc_truth_tWl1_p4;
  TBranch* b_mc_truth_tWl2_id;
  TBranch* b_mc_truth_tWl2_p4;
  TBranch* b_mc_truth_tWnu1_p4;
  TBranch* b_mc_truth_tWnu2_p4;
  TBranch* b_mc_truth_tWq11_id;
  TBranch* b_mc_truth_tWq11_p4;
  TBranch* b_mc_truth_tWq21_id;
  TBranch* b_mc_truth_tWq21_p4;
  TBranch* b_mc_truth_tWq12_id;
  TBranch* b_mc_truth_tWq12_p4;
  TBranch* b_mc_truth_tWq22_id;
  TBranch* b_mc_truth_tWq22_p4;
  TBranch* b_mc_truth_Z_id;
  TBranch* b_mc_truth_Z_p4;
  TBranch* b_mc_truth_Zl1_id;
  TBranch* b_mc_truth_Zl1_p4;
  TBranch* b_mc_truth_Zl2_id;
  TBranch* b_mc_truth_Zl2_p4;
  TBranch* b_mc_truth_W_id;
  TBranch* b_mc_truth_W_p4;
  TBranch* b_mc_truth_Wl_id;
  TBranch* b_mc_truth_Wl_p4;
  TBranch* b_mc_truth_gammal1_id;
  TBranch* b_mc_truth_gammal1_p4;
  TBranch* b_mc_truth_gammal2_id;
  TBranch* b_mc_truth_gammal2_p4;
  TBranch* b_genJet_n;
  TBranch* b_genJet_pt;
  TBranch* b_genJet_eta;
  TBranch* b_genJet_phi;
  TBranch* b_genJet_E;

  TBranch* b_max_Lep_eta;
  TBranch* b_MT_met_lep1;
  TBranch* b_mindr_lep1_jet;
  TBranch* b_mindr_lep2_jet;
  TBranch* b_LepGood_conePt0;
  TBranch* b_LepGood_conePt1;
  TBranch* b_nJet25_Recl;
  TBranch* b_met;
  TBranch* b_avg_dr_jet;
  TBranch* b_mhtJet25_Recl;

  TBranch* b_signal_2lss_TT_MVA;
  TBranch* b_signal_2lss_TTV_MVA;
  TBranch* b_signal_3l_TT_MVA;
  TBranch* b_signal_3l_TTV_MVA;

  Int_t mc_event;
  Float_t mc_totp4_px;
  Float_t mc_totp4_py;
  Float_t mc_totp4_pt;
  Float_t mc_thad_pt;
  Float_t mc_thad_b_pt;
  Float_t mc_thad_b_eta;
  Float_t mc_thad_j1_pt;
  Float_t mc_thad_j1_eta;
  Float_t mc_thad_j2_pt;
  Float_t mc_thad_j2_eta;
  Float_t mc_tlep_pt;
  Float_t mc_tlep_b_pt;
  Float_t mc_tlep_b_eta;
  Float_t mc_tlep_l_pt;
  Float_t mc_tlep_l_eta;
  Float_t mc_tlep2_pt;
  Float_t mc_tlep2_b_pt;
  Float_t mc_tlep2_b_eta;
  Float_t mc_tlep2_l_pt;
  Float_t mc_tlep2_l_eta;
  Int_t   mc_ttbar_decay;
  Int_t   mc_boson_decay;
  Float_t mc_boson_pt;
  Float_t mc_boson_l1_pt;
  Float_t mc_boson_l1_eta;
  Float_t mc_boson_l2_pt;
  Float_t mc_boson_l2_eta;
  Float_t mc_boson_ll_mass;
  Float_t mc_boson_ll_pt;
  Float_t mc_boson_ll_dphi;
  Float_t mc_boson_j1_pt;
  Float_t mc_boson_j1_eta;
  Float_t mc_boson_j2_pt;
  Float_t mc_boson_j2_eta;
  Float_t mc_boson_jj_mass;
  Float_t mc_boson_jj_pt;
  Float_t mc_boson_jj_dphi;
  Float_t mc_met;
  Int_t mc_njets25;

  Int_t mc_3l_category;

  Int_t mc_ttZhypAllowed;
  Int_t mc_hasLLcombZpeak;
  Int_t mc_passMllGt12;
  Int_t mc_passLepPresel;
  Int_t mc_passJetPresel25;
  Int_t mc_passBjetPresel25;

  Int_t catJets;
  Char_t is_2lss_TTH_SR;
  Char_t is_3l_TTH_SR;
  Char_t is_emu_TT_CR;
  //Char_t is_Zl_CR;
  Char_t is_3l_TTZ_CR;
  Char_t is_3l_WZrel_CR;
  Char_t is_3l_TZQ_SR;

  Int_t is_2bTight;
  Float_t is_2bTight_float;

  Char_t cat_HtoWW;
  Char_t cat_HtoZZ;
  Char_t cat_Htott;

  Int_t 		multilepton_Lepton1_Id, 		multilepton_Lepton2_Id, 		multilepton_Lepton3_Id, 		multilepton_Lepton4_Id;
  TLorentzVector 	multilepton_Lepton1_P4, 		multilepton_Lepton2_P4, 		multilepton_Lepton3_P4, 		multilepton_Lepton4_P4;

  Float_t           	multilepton_Lepton1_DeltaR_Matched,  	multilepton_Lepton2_DeltaR_Matched,  	multilepton_Lepton3_DeltaR_Matched,  	multilepton_Lepton4_DeltaR_Matched;
  Int_t           	multilepton_Lepton1_Label_Matched,   	multilepton_Lepton2_Label_Matched,   	multilepton_Lepton3_Label_Matched,   	multilepton_Lepton4_Label_Matched;
  Int_t           	multilepton_Lepton1_Id_Matched,      	multilepton_Lepton2_Id_Matched,      	multilepton_Lepton3_Id_Matched,      	multilepton_Lepton4_Id_Matched;
  TLorentzVector  	multilepton_Lepton1_P4_Matched,      	multilepton_Lepton2_P4_Matched,      	multilepton_Lepton3_P4_Matched,      	multilepton_Lepton4_P4_Matched;

  Int_t 		multilepton_Bjet1_Id, 								multilepton_Bjet2_Id;
  TLorentzVector 	multilepton_Bjet1_P4, 								multilepton_Bjet2_P4;
  Float_t 		multilepton_Bjet1_CSV, 								multilepton_Bjet2_CSV;
  Float_t		multilepton_Bjet1_JEC_Up, 		multilepton_Bjet1_JEC_Down, 		multilepton_Bjet2_JEC_Up, 		multilepton_Bjet2_JEC_Down;
  Float_t 		multilepton_Bjet1_JER_Up, 		multilepton_Bjet1_JER_Down, 		multilepton_Bjet2_JER_Up, 		multilepton_Bjet2_JER_Down;

  Float_t           	multilepton_Bjet1_DeltaR_Matched,   						multilepton_Bjet2_DeltaR_Matched;
  Int_t           	multilepton_Bjet1_Label_Matched,    						multilepton_Bjet2_Label_Matched;
  Int_t           	multilepton_Bjet1_Id_Matched,       						multilepton_Bjet2_Id_Matched;
  TLorentzVector  	multilepton_Bjet1_P4_Matched,       						multilepton_Bjet2_P4_Matched;

  Int_t 		multilepton_h0_Id, 			multilepton_t1_Id,			multilepton_t2_Id;
  Int_t                 multilepton_h0_Label,                   multilepton_t1_Label,                   multilepton_t2_Label;
  TLorentzVector 	multilepton_h0_P4,			multilepton_t1_P4,			multilepton_t2_P4;

  Int_t multilepton_JetHighestPt1_Id, multilepton_JetHighestPt2_Id, multilepton_JetClosestMw1_Id, multilepton_JetClosestMw2_Id, multilepton_JetLowestMjj1_Id, multilepton_JetLowestMjj2_Id;
  TLorentzVector multilepton_JetHighestPt1_P4, multilepton_JetHighestPt2_P4, multilepton_JetClosestMw1_P4, multilepton_JetClosestMw2_P4, multilepton_JetLowestMjj1_P4, multilepton_JetLowestMjj2_P4;
  Float_t multilepton_JetHighestPt1_CSV, multilepton_JetHighestPt2_CSV, multilepton_JetClosestMw1_CSV, multilepton_JetClosestMw2_CSV, multilepton_JetLowestMjj1_CSV, multilepton_JetLowestMjj2_CSV;
  Float_t multilepton_JetHighestPt1_JEC_Up, multilepton_JetHighestPt2_JEC_Up, multilepton_JetClosestMw1_JEC_Up, multilepton_JetClosestMw2_JEC_Up, multilepton_JetLowestMjj1_JEC_Up, multilepton_JetLowestMjj2_JEC_Up;
  Float_t multilepton_JetHighestPt1_JEC_Down, multilepton_JetHighestPt2_JEC_Down, multilepton_JetClosestMw1_JEC_Down, multilepton_JetClosestMw2_JEC_Down, multilepton_JetLowestMjj1_JEC_Down, multilepton_JetLowestMjj2_JEC_Down;
  Float_t multilepton_JetHighestPt1_JER_Up, multilepton_JetHighestPt2_JER_Up, multilepton_JetClosestMw1_JER_Up, multilepton_JetClosestMw2_JER_Up, multilepton_JetLowestMjj1_JER_Up, multilepton_JetLowestMjj2_JER_Up;
  Float_t multilepton_JetHighestPt1_JER_Down, multilepton_JetHighestPt2_JER_Down, multilepton_JetClosestMw1_JER_Down, multilepton_JetClosestMw2_JER_Down, multilepton_JetLowestMjj1_JER_Down, multilepton_JetLowestMjj2_JER_Down;

  Int_t multilepton_JetHighestPt1_2ndPair_Id, multilepton_JetHighestPt2_2ndPair_Id, multilepton_JetClosestMw1_2ndPair_Id, multilepton_JetClosestMw2_2ndPair_Id, multilepton_JetLowestMjj1_2ndPair_Id, multilepton_JetLowestMjj2_2ndPair_Id;
  TLorentzVector  multilepton_JetHighestPt1_2ndPair_P4, multilepton_JetHighestPt2_2ndPair_P4, multilepton_JetClosestMw1_2ndPair_P4, multilepton_JetClosestMw2_2ndPair_P4, multilepton_JetLowestMjj1_2ndPair_P4, multilepton_JetLowestMjj2_2ndPair_P4;
  Float_t multilepton_JetHighestPt1_2ndPair_CSV, multilepton_JetHighestPt2_2ndPair_CSV, multilepton_JetClosestMw1_2ndPair_CSV, multilepton_JetClosestMw2_2ndPair_CSV, multilepton_JetLowestMjj1_2ndPair_CSV, multilepton_JetLowestMjj2_2ndPair_CSV;
  Float_t multilepton_JetHighestPt1_2ndPair_JEC_Up, multilepton_JetHighestPt2_2ndPair_JEC_Up, multilepton_JetClosestMw1_2ndPair_JEC_Up, multilepton_JetClosestMw2_2ndPair_JEC_Up, multilepton_JetLowestMjj1_2ndPair_JEC_Up, multilepton_JetLowestMjj2_2ndPair_JEC_Up;
  Float_t multilepton_JetHighestPt1_2ndPair_JEC_Down, multilepton_JetHighestPt2_2ndPair_JEC_Down, multilepton_JetClosestMw1_2ndPair_JEC_Down, multilepton_JetClosestMw2_2ndPair_JEC_Down, multilepton_JetLowestMjj1_2ndPair_JEC_Down, multilepton_JetLowestMjj2_2ndPair_JEC_Down;
  Float_t multilepton_JetHighestPt1_2ndPair_JER_Up, multilepton_JetHighestPt2_2ndPair_JER_Up, multilepton_JetClosestMw1_2ndPair_JER_Up, multilepton_JetClosestMw2_2ndPair_JER_Up, multilepton_JetLowestMjj1_2ndPair_JER_Up, multilepton_JetLowestMjj2_2ndPair_JER_Up;
  Float_t multilepton_JetHighestPt1_2ndPair_JER_Down, multilepton_JetHighestPt2_2ndPair_JER_Down, multilepton_JetClosestMw1_2ndPair_JER_Down, multilepton_JetClosestMw2_2ndPair_JER_Down, multilepton_JetLowestMjj1_2ndPair_JER_Down, multilepton_JetLowestMjj2_2ndPair_JER_Down;

/*
  Int_t multilepton_Bjet1_Id;
  TLorentzVector multilepton_Bjet1_P4;
  Int_t multilepton_Bjet2_Id;
  TLorentzVector multilepton_Bjet2_P4;
  Int_t multilepton_Lepton1_Id;
  TLorentzVector multilepton_Lepton1_P4;
  Int_t multilepton_Lepton2_Id;
  TLorentzVector multilepton_Lepton2_P4;
  Int_t multilepton_Lepton3_Id;
  TLorentzVector multilepton_Lepton3_P4;
  Int_t multilepton_Lepton4_Id;
  TLorentzVector multilepton_Lepton4_P4;
  Int_t multilepton_JetHighestPt1_Id;
  TLorentzVector multilepton_JetHighestPt1_P4;
  Int_t multilepton_JetHighestPt2_Id;
  TLorentzVector multilepton_JetHighestPt2_P4;
  Float_t multilepton_JetHighestPt_Mjj;
  Int_t multilepton_JetClosestMw1_Id;
  TLorentzVector multilepton_JetClosestMw1_P4;
  Int_t multilepton_JetClosestMw2_Id;
  TLorentzVector multilepton_JetClosestMw2_P4;
  Float_t multilepton_JetClosestMw_Mjj;
  Int_t multilepton_JetLowestMjj1_Id;
  TLorentzVector multilepton_JetLowestMjj1_P4;
  Int_t multilepton_JetLowestMjj2_Id;
  TLorentzVector multilepton_JetLowestMjj2_P4;
  Float_t multilepton_JetLowestMjj_Mjj;
  Int_t multilepton_JetHighestPt1_2ndPair_Id, multilepton_JetHighestPt2_2ndPair_Id, multilepton_JetClosestMw1_2ndPair_Id, multilepton_JetClosestMw2_2ndPair_Id, multilepton_JetLowestMjj1_2ndPair_Id, multilepton_JetLowestMjj2_2ndPair_Id;
  TLorentzVector  multilepton_JetHighestPt1_2ndPair_P4, multilepton_JetHighestPt2_2ndPair_P4, multilepton_JetClosestMw1_2ndPair_P4, multilepton_JetClosestMw2_2ndPair_P4, multilepton_JetLowestMjj1_2ndPair_P4, multilepton_JetLowestMjj2_2ndPair_P4;
*/
  TLorentzVector multilepton_mET;
  Double_t multilepton_mETcov00, multilepton_mETcov01, multilepton_mETcov10, multilepton_mETcov11;
  Float_t multilepton_mHT;
  TLorentzVector multilepton_Ptot;

  TLorentzVector* multilepton_Bjet1_P4_ptr;
  TLorentzVector* multilepton_Bjet1_P4_Matched_ptr;
  TLorentzVector* multilepton_Bjet2_P4_ptr;
  TLorentzVector* multilepton_Bjet2_P4_Matched_ptr;
  TLorentzVector* multilepton_Lepton1_P4_ptr;
  TLorentzVector* multilepton_Lepton1_P4_Matched_ptr;
  TLorentzVector* multilepton_Lepton2_P4_ptr;
  TLorentzVector* multilepton_Lepton2_P4_Matched_ptr;
  TLorentzVector* multilepton_Lepton3_P4_ptr;
  TLorentzVector* multilepton_Lepton3_P4_Matched_ptr;
  TLorentzVector* multilepton_Lepton4_P4_ptr;
  TLorentzVector* multilepton_Lepton4_P4_Matched_ptr;
  TLorentzVector* multilepton_JetHighestPt1_P4_ptr;
  TLorentzVector* multilepton_JetHighestPt2_P4_ptr;
  TLorentzVector* multilepton_JetClosestMw1_P4_ptr;
  TLorentzVector* multilepton_JetClosestMw2_P4_ptr;
  TLorentzVector* multilepton_JetLowestMjj1_P4_ptr;
  TLorentzVector* multilepton_JetLowestMjj2_P4_ptr;
  TLorentzVector* multilepton_JetHighestPt1_2ndPair_P4_ptr;
  TLorentzVector* multilepton_JetHighestPt2_2ndPair_P4_ptr;
  TLorentzVector* multilepton_JetClosestMw1_2ndPair_P4_ptr;
  TLorentzVector* multilepton_JetClosestMw2_2ndPair_P4_ptr;
  TLorentzVector* multilepton_JetLowestMjj1_2ndPair_P4_ptr;
  TLorentzVector* multilepton_JetLowestMjj2_2ndPair_P4_ptr;
  TLorentzVector* multilepton_h0_P4_ptr;
  TLorentzVector* multilepton_t1_P4_ptr;
  TLorentzVector* multilepton_t2_P4_ptr;

  TLorentzVector* multilepton_mET_ptr;
  TLorentzVector* multilepton_Ptot_ptr;

  Double_t mc_mem_tthfl_weight;
  Double_t mc_mem_tthfl_weight_JEC_up, mc_mem_tthfl_weight_JEC_down, mc_mem_tthfl_weight_JER_up, mc_mem_tthfl_weight_JER_down;
  Double_t mc_mem_tthfl_weight_log;
  Double_t mc_mem_tthfl_weight_err;
  Float_t mc_mem_tthfl_weight_chi2;
  Float_t mc_mem_tthfl_weight_time;
  Double_t mc_mem_tthfl_weight_max;
  Double_t mc_mem_tthfl_weight_avg;
  Double_t mc_mem_tthfl_weight_logmean;
  Double_t mc_kin_tthfl_weight_logmax;
  Double_t mc_kin_tthfl_weight_logmaxint;
  Double_t mc_mem_tthfl_weight_kinmax;
  Double_t mc_mem_tthfl_weight_kinmaxint;

  TLorentzVector mc_kin_tthfl_tophad_P4;
  TLorentzVector mc_kin_tthfl_tophad_Bjet_P4;
  TLorentzVector mc_kin_tthfl_tophad_W_P4;
  TLorentzVector mc_kin_tthfl_tophad_Jet1_P4;
  TLorentzVector mc_kin_tthfl_tophad_Jet2_P4;
  Float_t mc_kin_tthfl_tophad_Pt;
  Float_t mc_kin_tthfl_tophad_Wmass;
  Float_t mc_kin_tthfl_tophad_Benergy;
  Float_t mc_kin_tthfl_tophad_Jet1energy;
  Float_t mc_kin_tthfl_tophad_Jet2energy;
  TLorentzVector mc_kin_tthfl_toplep_P4;
  TLorentzVector mc_kin_tthfl_toplep_Bjet_P4;
  TLorentzVector mc_kin_tthfl_toplep_W_P4;
  TLorentzVector mc_kin_tthfl_toplep_Lep_P4;
  TLorentzVector mc_kin_tthfl_toplep_Neut_P4;
  Float_t mc_kin_tthfl_toplep_Pt;
  Float_t mc_kin_tthfl_toplep_Wmass;
  Float_t mc_kin_tthfl_toplep_Benergy;
  Float_t mc_kin_tthfl_toplep_Neutenergy;
  TLorentzVector mc_kin_tthfl_toplep2_P4;
  TLorentzVector mc_kin_tthfl_toplep2_Bjet_P4;
  TLorentzVector mc_kin_tthfl_toplep2_W_P4;
  TLorentzVector mc_kin_tthfl_toplep2_Lep_P4;
  TLorentzVector mc_kin_tthfl_toplep2_Neut_P4;
  Float_t mc_kin_tthfl_toplep2_Pt;
  Float_t mc_kin_tthfl_toplep2_Wmass;
  Float_t mc_kin_tthfl_toplep2_Benergy;
  Float_t mc_kin_tthfl_toplep2_Neutenergy;
  TLorentzVector mc_kin_tthfl_h2l2nu_P4;
  TLorentzVector mc_kin_tthfl_h2l2nu_W1_P4;
  TLorentzVector mc_kin_tthfl_h2l2nu_W2_P4;
  TLorentzVector mc_kin_tthfl_h2l2nu_Lep1_P4;
  TLorentzVector mc_kin_tthfl_h2l2nu_Neut1_P4;
  TLorentzVector mc_kin_tthfl_h2l2nu_Lep2_P4;
  TLorentzVector mc_kin_tthfl_h2l2nu_Neut2_P4;
  Float_t mc_kin_tthfl_h2l2nu_Pt;
  Float_t mc_kin_tthfl_h2l2nu_W1mass;
  Float_t mc_kin_tthfl_h2l2nu_Neut1energy;
  Float_t mc_kin_tthfl_h2l2nu_W2mass;
  Float_t mc_kin_tthfl_h2l2nu_Neut2energy;

  Double_t mc_mem_tthsl_weight;
  Double_t mc_mem_tthsl_weight_JEC_up, mc_mem_tthsl_weight_JEC_down, mc_mem_tthsl_weight_JER_up, mc_mem_tthsl_weight_JER_down;
  Double_t mc_mem_tthsl_weight_log;
  Double_t mc_mem_tthsl_weight_err;
  Float_t mc_mem_tthsl_weight_chi2;
  Float_t mc_mem_tthsl_weight_time;
  Double_t mc_mem_tthsl_weight_max;
  Double_t mc_mem_tthsl_weight_avg;
  Double_t mc_mem_tthsl_weight_logmean;
  Double_t mc_kin_tthsl_weight_logmax;
  Double_t mc_kin_tthsl_weight_logmaxint;
  Double_t mc_mem_tthsl_weight_kinmax;
  Double_t mc_mem_tthsl_weight_kinmaxint;

  TLorentzVector mc_kin_tthsl_tophad_P4;
  TLorentzVector mc_kin_tthsl_tophad_Bjet_P4;
  TLorentzVector mc_kin_tthsl_tophad_W_P4;
  TLorentzVector mc_kin_tthsl_tophad_Jet1_P4;
  TLorentzVector mc_kin_tthsl_tophad_Jet2_P4;
  Float_t mc_kin_tthsl_tophad_Pt;
  Float_t mc_kin_tthsl_tophad_Wmass;
  Float_t mc_kin_tthsl_tophad_Benergy;
  Float_t mc_kin_tthsl_tophad_Jet1energy;
  Float_t mc_kin_tthsl_tophad_Jet2energy;
  TLorentzVector mc_kin_tthsl_toplep_P4;
  TLorentzVector mc_kin_tthsl_toplep_Bjet_P4;
  TLorentzVector mc_kin_tthsl_toplep_W_P4;
  TLorentzVector mc_kin_tthsl_toplep_Lep_P4;
  TLorentzVector mc_kin_tthsl_toplep_Neut_P4;
  Float_t mc_kin_tthsl_toplep_Pt;
  Float_t mc_kin_tthsl_toplep_Wmass;
  Float_t mc_kin_tthsl_toplep_Benergy;
  Float_t mc_kin_tthsl_toplep_Neutenergy;
  TLorentzVector mc_kin_tthsl_toplep2_P4;
  TLorentzVector mc_kin_tthsl_toplep2_Bjet_P4;
  TLorentzVector mc_kin_tthsl_toplep2_W_P4;
  TLorentzVector mc_kin_tthsl_toplep2_Lep_P4;
  TLorentzVector mc_kin_tthsl_toplep2_Neut_P4;
  Float_t mc_kin_tthsl_toplep2_Pt;
  Float_t mc_kin_tthsl_toplep2_Wmass;
  Float_t mc_kin_tthsl_toplep2_Benergy;
  Float_t mc_kin_tthsl_toplep2_Neutenergy;
  TLorentzVector mc_kin_tthsl_hlnujj_P4;
  TLorentzVector mc_kin_tthsl_hlnujj_W1_P4;
  TLorentzVector mc_kin_tthsl_hlnujj_W2_P4;
  TLorentzVector mc_kin_tthsl_hlnujj_Lep_P4;
  TLorentzVector mc_kin_tthsl_hlnujj_Neut_P4;
  TLorentzVector mc_kin_tthsl_hlnujj_Jet1_P4;
  TLorentzVector mc_kin_tthsl_hlnujj_Jet2_P4;
  Float_t mc_kin_tthsl_hlnujj_Pt;
  Float_t mc_kin_tthsl_hlnujj_W1mass;
  Float_t mc_kin_tthsl_hlnujj_Neut1energy;
  Float_t mc_kin_tthsl_hlnujj_W2mass;
  Float_t mc_kin_tthsl_hlnujj_Jet1energy;
  Float_t mc_kin_tthsl_hlnujj_Jet2energy;

  Double_t mc_mem_tth_weight;
  Double_t mc_mem_tth_weight_JEC_up, mc_mem_tth_weight_JEC_down, mc_mem_tth_weight_JER_up, mc_mem_tth_weight_JER_down;
  Double_t mc_mem_tth_weight_log;
  Double_t mc_mem_tth_weight_err;
  Float_t mc_mem_tth_weight_chi2;
  Float_t mc_mem_tth_weight_time;
  Double_t mc_mem_tth_weight_max;
  Double_t mc_mem_tth_weight_avg;
  Double_t mc_mem_tth_weight_logmean;
  Double_t mc_kin_tth_weight_logmax;
  Double_t mc_kin_tth_weight_logmaxint;
  Double_t mc_mem_tth_weight_kinmax;
  Double_t mc_mem_tth_weight_kinmaxint;

  Double_t mc_mem_ttz_weight;
  Double_t mc_mem_ttz_weight_JEC_up, mc_mem_ttz_weight_JEC_down, mc_mem_ttz_weight_JER_up, mc_mem_ttz_weight_JER_down;
  Double_t mc_mem_ttz_weight_log;
  Double_t mc_mem_ttz_weight_err;
  Float_t mc_mem_ttz_weight_chi2;
  Float_t mc_mem_ttz_weight_time;
  Double_t mc_mem_ttz_weight_max;
  Double_t mc_mem_ttz_weight_avg;
  Double_t mc_mem_ttz_weight_logmean;
  Double_t mc_kin_ttz_weight_logmax;
  Double_t mc_kin_ttz_weight_logmaxint;
  Double_t mc_mem_ttz_weight_kinmax;
  Double_t mc_mem_ttz_weight_kinmaxint;

  TLorentzVector mc_kin_ttz_tophad_P4;
  TLorentzVector mc_kin_ttz_tophad_Bjet_P4;
  TLorentzVector mc_kin_ttz_tophad_W_P4;
  TLorentzVector mc_kin_ttz_tophad_Jet1_P4;
  TLorentzVector mc_kin_ttz_tophad_Jet2_P4;
  Float_t mc_kin_ttz_tophad_Pt;
  Float_t mc_kin_ttz_tophad_Wmass;
  Float_t mc_kin_ttz_tophad_Benergy;
  Float_t mc_kin_ttz_tophad_Jet1energy;
  Float_t mc_kin_ttz_tophad_Jet2energy;
  TLorentzVector mc_kin_ttz_toplep_P4;
  TLorentzVector mc_kin_ttz_toplep_Bjet_P4;
  TLorentzVector mc_kin_ttz_toplep_W_P4;
  TLorentzVector mc_kin_ttz_toplep_Lep_P4;
  TLorentzVector mc_kin_ttz_toplep_Neut_P4;
  Float_t mc_kin_ttz_toplep_Pt;
  Float_t mc_kin_ttz_toplep_Wmass;
  Float_t mc_kin_ttz_toplep_Benergy;
  Float_t mc_kin_ttz_toplep_Neutenergy;
  TLorentzVector mc_kin_ttz_toplep2_P4;
  TLorentzVector mc_kin_ttz_toplep2_Bjet_P4;
  TLorentzVector mc_kin_ttz_toplep2_W_P4;
  TLorentzVector mc_kin_ttz_toplep2_Lep_P4;
  TLorentzVector mc_kin_ttz_toplep2_Neut_P4;
  Float_t mc_kin_ttz_toplep2_Pt;
  Float_t mc_kin_ttz_toplep2_Wmass;
  Float_t mc_kin_ttz_toplep2_Benergy;
  Float_t mc_kin_ttz_toplep2_Neutenergy;
  TLorentzVector mc_kin_ttz_zll_P4;
  TLorentzVector mc_kin_ttz_zll_Lep1_P4;
  TLorentzVector mc_kin_ttz_zll_Lep2_P4;
  Float_t mc_kin_ttz_zll_Pt;
  Float_t mc_kin_ttz_zll_Zmass;

  Double_t mc_mem_ttw_weight;
  Double_t mc_mem_ttw_weight_JEC_up, mc_mem_ttw_weight_JEC_down, mc_mem_ttw_weight_JER_up, mc_mem_ttw_weight_JER_down;
  Double_t mc_mem_ttw_weight_log;
  Double_t mc_mem_ttw_weight_err;
  Float_t mc_mem_ttw_weight_chi2;
  Float_t mc_mem_ttw_weight_time;
  Double_t mc_mem_ttw_weight_max;
  Double_t mc_mem_ttw_weight_avg;
  Double_t mc_mem_ttw_weight_logmean;
  Double_t mc_kin_ttw_weight_logmax;
  Double_t mc_kin_ttw_weight_logmaxint;
  Double_t mc_mem_ttw_weight_kinmax;
  Double_t mc_mem_ttw_weight_kinmaxint;

  TLorentzVector mc_kin_ttw_tophad_P4;
  TLorentzVector mc_kin_ttw_tophad_Bjet_P4;
  TLorentzVector mc_kin_ttw_tophad_W_P4;
  TLorentzVector mc_kin_ttw_tophad_Jet1_P4;
  TLorentzVector mc_kin_ttw_tophad_Jet2_P4;
  Float_t mc_kin_ttw_tophad_Pt;
  Float_t mc_kin_ttw_tophad_Wmass;
  Float_t mc_kin_ttw_tophad_Benergy;
  Float_t mc_kin_ttw_tophad_Jet1energy;
  Float_t mc_kin_ttw_tophad_Jet2energy;
  TLorentzVector mc_kin_ttw_toplep_P4;
  TLorentzVector mc_kin_ttw_toplep_Bjet_P4;
  TLorentzVector mc_kin_ttw_toplep_W_P4;
  TLorentzVector mc_kin_ttw_toplep_Lep_P4;
  TLorentzVector mc_kin_ttw_toplep_Neut_P4;
  Float_t mc_kin_ttw_toplep_Pt;
  Float_t mc_kin_ttw_toplep_Wmass;
  Float_t mc_kin_ttw_toplep_Benergy;
  Float_t mc_kin_ttw_toplep_Neutenergy;
  TLorentzVector mc_kin_ttw_toplep2_P4;
  TLorentzVector mc_kin_ttw_toplep2_Bjet_P4;
  TLorentzVector mc_kin_ttw_toplep2_W_P4;
  TLorentzVector mc_kin_ttw_toplep2_Lep_P4;
  TLorentzVector mc_kin_ttw_toplep2_Neut_P4;
  Float_t mc_kin_ttw_toplep2_Pt;
  Float_t mc_kin_ttw_toplep2_Wmass;
  Float_t mc_kin_ttw_toplep2_Benergy;
  Float_t mc_kin_ttw_toplep2_Neutenergy;
  TLorentzVector mc_kin_ttw_wlnu_W_P4;
  TLorentzVector mc_kin_ttw_wlnu_Lep_P4;
  TLorentzVector mc_kin_ttw_wlnu_Neut_P4;
  Float_t mc_kin_ttw_wlnu_Pt;
  Float_t mc_kin_ttw_wlnu_Wmass;
  Float_t mc_kin_ttw_wlnu_Neutenergy;

  Double_t mc_mem_ttwjj_weight;
  Double_t mc_mem_ttwjj_weight_JEC_up, mc_mem_ttwjj_weight_JEC_down, mc_mem_ttwjj_weight_JER_up, mc_mem_ttwjj_weight_JER_down;
  Double_t mc_mem_ttwjj_weight_log;
  Double_t mc_mem_ttwjj_weight_err;
  Float_t mc_mem_ttwjj_weight_chi2;
  Float_t mc_mem_ttwjj_weight_time;
  Double_t mc_mem_ttwjj_weight_max;
  Double_t mc_mem_ttwjj_weight_avg;
  Double_t mc_mem_ttwjj_weight_logmean;
  Double_t mc_kin_ttwjj_weight_logmax;
  Double_t mc_kin_ttwjj_weight_logmaxint;
  Double_t mc_mem_ttwjj_weight_kinmax;
  Double_t mc_mem_ttwjj_weight_kinmaxint;

  Double_t mc_mem_ttbarfl_weight;
  Double_t mc_mem_ttbarfl_weight_JEC_up, mc_mem_ttbarfl_weight_JEC_down, mc_mem_ttbarfl_weight_JER_up, mc_mem_ttbarfl_weight_JER_down;
  Double_t mc_mem_ttbarfl_weight_log;
  Double_t mc_mem_ttbarfl_weight_err;
  Float_t mc_mem_ttbarfl_weight_chi2;
  Float_t mc_mem_ttbarfl_weight_time;
  Double_t mc_mem_ttbarfl_weight_max;
  Double_t mc_mem_ttbarfl_weight_avg;
  Double_t mc_mem_ttbarfl_weight_logmean;
  Double_t mc_kin_ttbarfl_weight_logmax;
  Double_t mc_kin_ttbarfl_weight_logmaxint;
  Double_t mc_mem_ttbarfl_weight_kinmax;
  Double_t mc_mem_ttbarfl_weight_kinmaxint;

  Double_t mc_mem_ttbarsl_weight;
  Double_t mc_mem_ttbarsl_weight_JEC_up, mc_mem_ttbarsl_weight_JEC_down, mc_mem_ttbarsl_weight_JER_up, mc_mem_ttbarsl_weight_JER_down;
  Double_t mc_mem_ttbarsl_weight_log;
  Double_t mc_mem_ttbarsl_weight_err;
  Float_t mc_mem_ttbarsl_weight_chi2;
  Float_t mc_mem_ttbarsl_weight_time;
  Double_t mc_mem_ttbarsl_weight_max;
  Double_t mc_mem_ttbarsl_weight_avg;
  Double_t mc_mem_ttbarsl_weight_logmean;
  Double_t mc_kin_ttbarsl_weight_logmax;
  Double_t mc_kin_ttbarsl_weight_logmaxint;
  Double_t mc_mem_ttbarsl_weight_kinmax;
  Double_t mc_mem_ttbarsl_weight_kinmaxint;

  Double_t mc_mem_ttbar_weight;
  Double_t mc_mem_ttbar_weight_JEC_up, mc_mem_ttbar_weight_JEC_down, mc_mem_ttbar_weight_JER_up, mc_mem_ttbar_weight_JER_down;
  Double_t mc_mem_ttbar_weight_log;
  Double_t mc_mem_ttbar_weight_err;
  Float_t mc_mem_ttbar_weight_chi2;
  Float_t mc_mem_ttbar_weight_time;
  Double_t mc_mem_ttbar_weight_max;
  Double_t mc_mem_ttbar_weight_avg;
  Double_t mc_mem_ttbar_weight_logmean;
  Double_t mc_kin_ttbar_weight_logmax;
  Double_t mc_kin_ttbar_weight_logmaxint;
  Double_t mc_mem_ttbar_weight_kinmax;
  Double_t mc_mem_ttbar_weight_kinmaxint;

  Double_t mc_mem_tllj_weight;
  Double_t mc_mem_tllj_weight_JEC_up, mc_mem_tllj_weight_JEC_down, mc_mem_tllj_weight_JER_up, mc_mem_tllj_weight_JER_down;
  Double_t mc_mem_tllj_weight_log;
  Double_t mc_mem_tllj_weight_err;
  Float_t mc_mem_tllj_weight_chi2;
  Float_t mc_mem_tllj_weight_time;
  Double_t mc_mem_tllj_weight_max;
  Double_t mc_mem_tllj_weight_avg;
  Double_t mc_mem_tllj_weight_logmean;
  Double_t mc_kin_tllj_weight_logmax;
  Double_t mc_kin_tllj_weight_logmaxint;
  Double_t mc_mem_tllj_weight_kinmax;
  Double_t mc_mem_tllj_weight_kinmaxint;

  //Double_t mc_mem_ttz_tthfl_likelihood;
  //Double_t mc_mem_ttz_tthsl_likelihood;
  //Double_t mc_mem_ttw_tthfl_likelihood;
  //Double_t mc_mem_ttw_tthsl_likelihood;
  Double_t mc_mem_ttz_tth_likelihood;
  Double_t mc_mem_ttz_tth_likelihood_nlog;
  Double_t mc_mem_ttz_tth_likelihood_max;
  Double_t mc_mem_ttz_tth_likelihood_avg;

  Double_t mc_mem_ttw_tth_likelihood;
  Double_t mc_mem_ttw_tth_likelihood_nlog;
  Double_t mc_mem_ttw_tth_likelihood_max;
  Double_t mc_mem_ttw_tth_likelihood_avg;

  Double_t mc_mem_ttwjj_tth_likelihood;
  Double_t mc_mem_ttwjj_tth_likelihood_nlog;
  Double_t mc_mem_ttwjj_tth_likelihood_max;
  Double_t mc_mem_ttwjj_tth_likelihood_avg;

  Double_t mc_mem_ttbar_tth_likelihood;
  Double_t mc_mem_ttbar_tth_likelihood_nlog;
  Double_t mc_mem_ttbar_tth_likelihood_max;
  Double_t mc_mem_ttbar_tth_likelihood_avg;

  Double_t mc_mem_ttz_tllj_likelihood;
  Double_t mc_mem_ttz_tllj_likelihood_nlog;
  Double_t mc_mem_ttz_tllj_likelihood_max;
  Double_t mc_mem_ttz_tllj_likelihood_avg;

  Double_t mc_mem_ttv_tth_likelihood;
  Double_t mc_mem_ttv_tth_likelihood_nlog;
  Double_t mc_mem_ttv_tth_likelihood_max;
  Double_t mc_mem_ttv_tth_likelihood_avg;

  Double_t mc_mem_ttvjj_tth_likelihood;
  Double_t mc_mem_ttvjj_tth_likelihood_nlog;
  Double_t mc_mem_ttvjj_tth_likelihood_max;
  Double_t mc_mem_ttvjj_tth_likelihood_avg;
  /*
  std::vector<double>* MEAllWeights_TTLL;
  std::vector<double>* MEAllWeights_TTHfl;
  std::vector<double>* MEAllWeights_TTHsl;
  std::vector<double>* MEAllWeights_TTH;
  std::vector<double>* MEAllWeights_TTW;
  std::vector<double>* MEAllWeights_TTWJJ;
  std::vector<double>* MEAllWeights_TTbarfl;
  std::vector<double>* MEAllWeights_TTbarsl;
  std::vector<double>* MEAllWeights_TTbar;
  std::vector<double>* MEAllWeights_TLLJ;

  std::vector<float>* MEAllWeights_TTLL_log;
  std::vector<float>* MEAllWeights_TTHfl_log;
  std::vector<float>* MEAllWeights_TTHsl_log;
  std::vector<float>* MEAllWeights_TTH_log;
  std::vector<float>* MEAllWeights_TTW_log;
  std::vector<float>* MEAllWeights_TTWJJ_log;
  std::vector<float>* MEAllWeights_TTbarfl_log;
  std::vector<float>* MEAllWeights_TTbarsl_log;
  std::vector<float>* MEAllWeights_TTbar_log;
  std::vector<float>* MEAllWeights_TLLJ_log;
  */
  TBranch* b_catJets;
  TBranch* b_is_2lss_TTH_SR;
  TBranch* b_is_3l_TTH_SR;
  TBranch* b_is_emu_TT_CR;
  //TBranch* b_is_Zl_CR;
  TBranch* b_is_3l_TTZ_CR;
  //TBranch* b_is_3l_WZrel_CR;
  TBranch* b_is_3l_TZQ_SR;

  TBranch* b_is_2bTight;
  TBranch* b_is_2bTight_float;

  TBranch* b_mc_3l_category;
  TBranch* b_mc_ttbar_decay;
  TBranch* b_mc_boson_decay;
  TBranch* b_mc_ttZhypAllowed;

  TBranch* b_cat_HtoWW;
  TBranch* b_cat_HtoZZ;
  TBranch* b_cat_Htott;

  TBranch* b_multilepton_Lepton1_Id;
  TBranch* b_multilepton_Lepton1_P4;
  TBranch* b_multilepton_Lepton1_DeltaR_Matched;
  TBranch* b_multilepton_Lepton1_Label_Matched;
  TBranch* b_multilepton_Lepton1_Id_Matched;
  TBranch* b_multilepton_Lepton1_P4_Matched;

  TBranch* b_multilepton_Lepton2_Id;
  TBranch* b_multilepton_Lepton2_P4;
  TBranch* b_multilepton_Lepton2_DeltaR_Matched;
  TBranch* b_multilepton_Lepton2_Label_Matched;
  TBranch* b_multilepton_Lepton2_Id_Matched;
  TBranch* b_multilepton_Lepton2_P4_Matched;

  TBranch* b_multilepton_Lepton3_Id;
  TBranch* b_multilepton_Lepton3_P4;
  TBranch* b_multilepton_Lepton3_DeltaR_Matched;
  TBranch* b_multilepton_Lepton3_Label_Matched;
  TBranch* b_multilepton_Lepton3_Id_Matched;
  TBranch* b_multilepton_Lepton3_P4_Matched;

  TBranch* b_multilepton_Lepton4_Id;
  TBranch* b_multilepton_Lepton4_P4;
  TBranch* b_multilepton_Lepton4_DeltaR_Matched;
  TBranch* b_multilepton_Lepton4_Label_Matched;
  TBranch* b_multilepton_Lepton4_Id_Matched;
  TBranch* b_multilepton_Lepton4_P4_Matched;

  TBranch* b_multilepton_Bjet1_Id;
  TBranch* b_multilepton_Bjet1_P4;
  TBranch* b_multilepton_Bjet1_CSV;
  TBranch* b_multilepton_Bjet1_JEC_Up;
  TBranch* b_multilepton_Bjet1_JEC_Down;
  TBranch* b_multilepton_Bjet1_JER_Up;
  TBranch* b_multilepton_Bjet1_JER_Down;
  TBranch* b_multilepton_Bjet1_DeltaR_Matched;
  TBranch* b_multilepton_Bjet1_Label_Matched;
  TBranch* b_multilepton_Bjet1_Id_Matched;
  TBranch* b_multilepton_Bjet1_P4_Matched;

  TBranch* b_multilepton_Bjet2_Id;
  TBranch* b_multilepton_Bjet2_P4;
  TBranch* b_multilepton_Bjet2_CSV;
  TBranch* b_multilepton_Bjet2_JEC_Up;
  TBranch* b_multilepton_Bjet2_JEC_Down;
  TBranch* b_multilepton_Bjet2_JER_Up;
  TBranch* b_multilepton_Bjet2_JER_Down;
  TBranch* b_multilepton_Bjet2_DeltaR_Matched;
  TBranch* b_multilepton_Bjet2_Label_Matched;
  TBranch* b_multilepton_Bjet2_Id_Matched;
  TBranch* b_multilepton_Bjet2_P4_Matched;

  TBranch* b_multilepton_JetHighestPt1_Id;
  TBranch* b_multilepton_JetHighestPt1_P4;
  TBranch* b_multilepton_JetHighestPt1_CSV;
  TBranch* b_multilepton_JetHighestPt1_JEC_Up;
  TBranch* b_multilepton_JetHighestPt1_JEC_Down;
  TBranch* b_multilepton_JetHighestPt1_JER_Up;
  TBranch* b_multilepton_JetHighestPt1_JER_Down;

  TBranch* b_multilepton_JetHighestPt2_Id;
  TBranch* b_multilepton_JetHighestPt2_P4;
  TBranch* b_multilepton_JetHighestPt2_CSV;
  TBranch* b_multilepton_JetHighestPt2_JEC_Up;
  TBranch* b_multilepton_JetHighestPt2_JEC_Down;
  TBranch* b_multilepton_JetHighestPt2_JER_Up;
  TBranch* b_multilepton_JetHighestPt2_JER_Down;

  TBranch* b_multilepton_JetClosestMw1_Id;
  TBranch* b_multilepton_JetClosestMw1_P4;
  TBranch* b_multilepton_JetClosestMw1_CSV;
  TBranch* b_multilepton_JetClosestMw1_JEC_Up;
  TBranch* b_multilepton_JetClosestMw1_JEC_Down;
  TBranch* b_multilepton_JetClosestMw1_JER_Up;
  TBranch* b_multilepton_JetClosestMw1_JER_Down;

  TBranch* b_multilepton_JetClosestMw2_Id;
  TBranch* b_multilepton_JetClosestMw2_P4;
  TBranch* b_multilepton_JetClosestMw2_CSV;
  TBranch* b_multilepton_JetClosestMw2_JEC_Up;
  TBranch* b_multilepton_JetClosestMw2_JEC_Down;
  TBranch* b_multilepton_JetClosestMw2_JER_Up;
  TBranch* b_multilepton_JetClosestMw2_JER_Down;

  TBranch* b_multilepton_JetLowestMjj1_Id;
  TBranch* b_multilepton_JetLowestMjj1_P4;
  TBranch* b_multilepton_JetLowestMjj1_CSV;
  TBranch* b_multilepton_JetLowestMjj1_JEC_Up;
  TBranch* b_multilepton_JetLowestMjj1_JEC_Down;
  TBranch* b_multilepton_JetLowestMjj1_JER_Up;
  TBranch* b_multilepton_JetLowestMjj1_JER_Down;

  TBranch* b_multilepton_JetLowestMjj2_Id;
  TBranch* b_multilepton_JetLowestMjj2_P4;
  TBranch* b_multilepton_JetLowestMjj2_CSV;
  TBranch* b_multilepton_JetLowestMjj2_JEC_Up;
  TBranch* b_multilepton_JetLowestMjj2_JEC_Down;
  TBranch* b_multilepton_JetLowestMjj2_JER_Up;
  TBranch* b_multilepton_JetLowestMjj2_JER_Down;

  TBranch* b_multilepton_JetHighestPt_Mjj;
  TBranch* b_multilepton_JetClosestMw_Mjj;
  TBranch* b_multilepton_JetLowestMjj_Mjj;

  TBranch* b_multilepton_JetHighestPt1_2ndPair_P4;
  TBranch* b_multilepton_JetHighestPt1_2ndPair_Id;
  TBranch* b_multilepton_JetHighestPt1_2ndPair_CSV;
  TBranch* b_multilepton_JetHighestPt1_2ndPair_JEC_Up;
  TBranch* b_multilepton_JetHighestPt1_2ndPair_JEC_Down;
  TBranch* b_multilepton_JetHighestPt1_2ndPair_JER_Up;
  TBranch* b_multilepton_JetHighestPt1_2ndPair_JER_Down;

  TBranch* b_multilepton_JetHighestPt2_2ndPair_P4;
  TBranch* b_multilepton_JetHighestPt2_2ndPair_Id;
  TBranch* b_multilepton_JetHighestPt2_2ndPair_CSV;
  TBranch* b_multilepton_JetHighestPt2_2ndPair_JEC_Up;
  TBranch* b_multilepton_JetHighestPt2_2ndPair_JEC_Down;
  TBranch* b_multilepton_JetHighestPt2_2ndPair_JER_Up;
  TBranch* b_multilepton_JetHighestPt2_2ndPair_JER_Down;

  TBranch* b_multilepton_JetClosestMw1_2ndPair_Id;
  TBranch* b_multilepton_JetClosestMw1_2ndPair_P4;
  TBranch* b_multilepton_JetClosestMw1_2ndPair_CSV;
  TBranch* b_multilepton_JetClosestMw1_2ndPair_JEC_Up;
  TBranch* b_multilepton_JetClosestMw1_2ndPair_JEC_Down;
  TBranch* b_multilepton_JetClosestMw1_2ndPair_JER_Up;
  TBranch* b_multilepton_JetClosestMw1_2ndPair_JER_Down;

  TBranch* b_multilepton_JetClosestMw2_2ndPair_Id;
  TBranch* b_multilepton_JetClosestMw2_2ndPair_P4;
  TBranch* b_multilepton_JetClosestMw2_2ndPair_CSV;
  TBranch* b_multilepton_JetClosestMw2_2ndPair_JEC_Up;
  TBranch* b_multilepton_JetClosestMw2_2ndPair_JEC_Down;
  TBranch* b_multilepton_JetClosestMw2_2ndPair_JER_Up;
  TBranch* b_multilepton_JetClosestMw2_2ndPair_JER_Down;

  TBranch* b_multilepton_JetLowestMjj1_2ndPair_Id;
  TBranch* b_multilepton_JetLowestMjj1_2ndPair_P4;
  TBranch* b_multilepton_JetLowestMjj1_2ndPair_CSV;
  TBranch* b_multilepton_JetLowestMjj1_2ndPair_JEC_Up;
  TBranch* b_multilepton_JetLowestMjj1_2ndPair_JEC_Down;
  TBranch* b_multilepton_JetLowestMjj1_2ndPair_JER_Up;
  TBranch* b_multilepton_JetLowestMjj1_2ndPair_JER_Down;

  TBranch* b_multilepton_JetLowestMjj2_2ndPair_Id;
  TBranch* b_multilepton_JetLowestMjj2_2ndPair_P4;
  TBranch* b_multilepton_JetLowestMjj2_2ndPair_CSV;
  TBranch* b_multilepton_JetLowestMjj2_2ndPair_JEC_Up;
  TBranch* b_multilepton_JetLowestMjj2_2ndPair_JEC_Down;
  TBranch* b_multilepton_JetLowestMjj2_2ndPair_JER_Up;
  TBranch* b_multilepton_JetLowestMjj2_2ndPair_JER_Down;

  TBranch* b_multilepton_h0_P4;
  TBranch* b_multilepton_h0_Id;
  TBranch* b_multilepton_h0_Label;
  TBranch* b_multilepton_t1_P4;
  TBranch* b_multilepton_t1_Id;
  TBranch* b_multilepton_t1_Label;
  TBranch* b_multilepton_t2_P4;
  TBranch* b_multilepton_t2_Id;
  TBranch* b_multilepton_t2_Label;

  TBranch* b_multilepton_mET;
  TBranch* b_multilepton_mETcov00;
  TBranch* b_multilepton_mETcov01;
  TBranch* b_multilepton_mETcov10;
  TBranch* b_multilepton_mETcov11;
  TBranch* b_multilepton_mHT;
  TBranch* b_multilepton_Ptot;

  TH1F** hMEPhaseSpace_Error;
  TH1F** hMEPhaseSpace_ErrorTot;
  TH1F** hMEPhaseSpace_ErrorTot_Fail;
  TH1F** hMEPhaseSpace_ErrorTot_Pass;

  TH1D* hMEAllWeights[8][12];
  TH1D* hMEAllWeights_nlog[8][12];

  private:
};

ReadGenFlatTree::ReadGenFlatTree(){

}

ReadGenFlatTree::~ReadGenFlatTree(){

}

void ReadGenFlatTree::InitializeDryRun(string InputFileName){

  cout << "Opening input tree"<<endl;

  fInput = TFile::Open(InputFileName.c_str(),"READ");
  //tInput = (TTree*)fInput->Get("FlatTree/tree");
  tInput = (TTree*)fInput->Get("tree");

  mc_truth_h0_p4 = 0;
  mc_truth_h0Wl1_p4 = 0;
  mc_truth_h0Wl2_p4 = 0;
  mc_truth_h0Wnu1_p4 = 0;
  mc_truth_h0Wnu2_p4 = 0;
  mc_truth_h0Wq11_p4 = 0;
  mc_truth_h0Wq12_p4 = 0;
  mc_truth_h0Wq21_p4 = 0;
  mc_truth_h0Wq22_p4 = 0;
  mc_truth_t1_p4 = 0;
  mc_truth_t2_p4 = 0;
  mc_truth_tb1_p4 = 0;
  mc_truth_tb2_p4 = 0;
  mc_truth_tWl1_p4 = 0;
  mc_truth_tWnu1_p4 = 0;
  mc_truth_tWl2_p4 = 0;
  mc_truth_tWnu2_p4 = 0;
  mc_truth_tWq11_p4 = 0;
  mc_truth_tWq21_p4 = 0;
  mc_truth_tWq12_p4 = 0;
  mc_truth_tWq22_p4 = 0;
  mc_truth_Z_p4 = 0;
  mc_truth_Zl1_p4 = 0;
  mc_truth_Zl2_p4 = 0;
  mc_truth_W_p4 = 0;
  mc_truth_Wl_p4 = 0;
  mc_truth_gammal1_p4 = 0;
  mc_truth_gammal2_p4 = 0;
  genJet_pt = 0;
  genJet_eta = 0;
  genJet_phi = 0;
  genJet_E = 0;

  tInput->SetBranchAddress("mc_weight",&mc_weight,&b_mc_weight);
  tInput->SetBranchAddress("mc_truth_h0_id",&mc_truth_h0_id,&b_mc_truth_h0_id);
  tInput->SetBranchAddress("mc_truth_h0_p4",&mc_truth_h0_p4,&b_mc_truth_h0_p4);
  tInput->SetBranchAddress("mc_truth_h0Wl1_id",&mc_truth_h0Wl1_id,&b_mc_truth_h0Wl1_id);
  tInput->SetBranchAddress("mc_truth_h0Wl1_p4",&mc_truth_h0Wl1_p4,&b_mc_truth_h0Wl1_p4);
  tInput->SetBranchAddress("mc_truth_h0Wl2_id",&mc_truth_h0Wl2_id,&b_mc_truth_h0Wl2_id);
  tInput->SetBranchAddress("mc_truth_h0Wl2_p4",&mc_truth_h0Wl2_p4,&b_mc_truth_h0Wl2_p4);
  tInput->SetBranchAddress("mc_truth_h0Wnu1_p4",&mc_truth_h0Wnu1_p4,&b_mc_truth_h0Wnu1_p4);
  tInput->SetBranchAddress("mc_truth_h0Wnu2_p4",&mc_truth_h0Wnu2_p4,&b_mc_truth_h0Wnu2_p4);
  tInput->SetBranchAddress("mc_truth_h0Wq11_id",&mc_truth_h0Wq11_id,&b_mc_truth_h0Wq11_id);
  tInput->SetBranchAddress("mc_truth_h0Wq11_p4",&mc_truth_h0Wq11_p4,&b_mc_truth_h0Wq11_p4);
  tInput->SetBranchAddress("mc_truth_h0Wq12_id",&mc_truth_h0Wq12_id,&b_mc_truth_h0Wq12_id);
  tInput->SetBranchAddress("mc_truth_h0Wq12_p4",&mc_truth_h0Wq12_p4,&b_mc_truth_h0Wq12_p4);
  tInput->SetBranchAddress("mc_truth_h0Wq21_id",&mc_truth_h0Wq21_id,&b_mc_truth_h0Wq21_id);
  tInput->SetBranchAddress("mc_truth_h0Wq21_p4",&mc_truth_h0Wq21_p4,&b_mc_truth_h0Wq21_p4);
  tInput->SetBranchAddress("mc_truth_h0Wq22_id",&mc_truth_h0Wq22_id,&b_mc_truth_h0Wq22_id);
  tInput->SetBranchAddress("mc_truth_h0Wq22_p4",&mc_truth_h0Wq22_p4,&b_mc_truth_h0Wq22_p4);
  tInput->SetBranchAddress("mc_truth_t1_id",&mc_truth_t1_id,&b_mc_truth_t1_id);
  tInput->SetBranchAddress("mc_truth_t1_p4",&mc_truth_t1_p4,&b_mc_truth_t1_p4);
  tInput->SetBranchAddress("mc_truth_t2_id",&mc_truth_t2_id,&b_mc_truth_t2_id);
  tInput->SetBranchAddress("mc_truth_t2_p4",&mc_truth_t2_p4,&b_mc_truth_t2_p4);
  tInput->SetBranchAddress("mc_truth_tb1_id",&mc_truth_tb1_id,&b_mc_truth_tb1_id);
  tInput->SetBranchAddress("mc_truth_tb1_p4",&mc_truth_tb1_p4,&b_mc_truth_tb1_p4);
  tInput->SetBranchAddress("mc_truth_tb2_id",&mc_truth_tb2_id,&b_mc_truth_tb2_id);
  tInput->SetBranchAddress("mc_truth_tb2_p4",&mc_truth_tb2_p4,&b_mc_truth_tb2_p4);
  tInput->SetBranchAddress("mc_truth_tWl1_id",&mc_truth_tWl1_id,&b_mc_truth_tWl1_id);
  tInput->SetBranchAddress("mc_truth_tWl1_p4",&mc_truth_tWl1_p4,&b_mc_truth_tWl1_p4);
  tInput->SetBranchAddress("mc_truth_tWl2_id",&mc_truth_tWl2_id,&b_mc_truth_tWl2_id);
  tInput->SetBranchAddress("mc_truth_tWl2_p4",&mc_truth_tWl2_p4,&b_mc_truth_tWl2_p4);
  tInput->SetBranchAddress("mc_truth_tWnu1_p4",&mc_truth_tWnu1_p4,&b_mc_truth_tWnu1_p4);
  tInput->SetBranchAddress("mc_truth_tWnu2_p4",&mc_truth_tWnu2_p4,&b_mc_truth_tWnu2_p4);
  tInput->SetBranchAddress("mc_truth_tWq11_id",&mc_truth_tWq11_id,&b_mc_truth_tWq11_id);
  tInput->SetBranchAddress("mc_truth_tWq11_p4",&mc_truth_tWq11_p4,&b_mc_truth_tWq11_p4);
  tInput->SetBranchAddress("mc_truth_tWq21_id",&mc_truth_tWq21_id,&b_mc_truth_tWq21_id);
  tInput->SetBranchAddress("mc_truth_tWq21_p4",&mc_truth_tWq21_p4,&b_mc_truth_tWq21_p4);
  tInput->SetBranchAddress("mc_truth_tWq12_id",&mc_truth_tWq12_id,&b_mc_truth_tWq12_id);
  tInput->SetBranchAddress("mc_truth_tWq12_p4",&mc_truth_tWq12_p4,&b_mc_truth_tWq12_p4);
  tInput->SetBranchAddress("mc_truth_tWq22_id",&mc_truth_tWq22_id,&b_mc_truth_tWq22_id);
  tInput->SetBranchAddress("mc_truth_tWq22_p4",&mc_truth_tWq22_p4,&b_mc_truth_tWq22_p4);
  tInput->SetBranchAddress("mc_truth_Z_id",&mc_truth_Z_id,&b_mc_truth_Z_id);
  tInput->SetBranchAddress("mc_truth_Z_p4",&mc_truth_Z_p4,&b_mc_truth_Z_p4);
  tInput->SetBranchAddress("mc_truth_Zl1_id",&mc_truth_Zl1_id,&b_mc_truth_Zl1_id);
  tInput->SetBranchAddress("mc_truth_Zl1_p4",&mc_truth_Zl1_p4,&b_mc_truth_Zl1_p4);
  tInput->SetBranchAddress("mc_truth_Zl2_id",&mc_truth_Zl2_id,&b_mc_truth_Zl2_id);
  tInput->SetBranchAddress("mc_truth_Zl2_p4",&mc_truth_Zl2_p4,&b_mc_truth_Zl2_p4);
  tInput->SetBranchAddress("mc_truth_W_id",&mc_truth_W_id,&b_mc_truth_W_id);
  tInput->SetBranchAddress("mc_truth_W_p4",&mc_truth_W_p4,&b_mc_truth_W_p4);
  tInput->SetBranchAddress("mc_truth_Wl_id",&mc_truth_Wl_id,&b_mc_truth_Wl_id);
  tInput->SetBranchAddress("mc_truth_Wl_p4",&mc_truth_Wl_p4,&b_mc_truth_Wl_p4);
  tInput->SetBranchAddress("mc_truth_gammal1_id",&mc_truth_gammal1_id,&b_mc_truth_gammal1_id);
  tInput->SetBranchAddress("mc_truth_gammal1_p4",&mc_truth_gammal1_p4,&b_mc_truth_gammal1_p4);
  tInput->SetBranchAddress("mc_truth_gammal2_id",&mc_truth_gammal2_id,&b_mc_truth_gammal2_id);
  tInput->SetBranchAddress("mc_truth_gammal2_p4",&mc_truth_gammal2_p4,&b_mc_truth_gammal2_p4);
  tInput->SetBranchAddress("genJet_n",&genJet_n,&b_genJet_n);
  tInput->SetBranchAddress("genJet_pt",&genJet_pt,&b_genJet_pt);
  tInput->SetBranchAddress("genJet_eta",&genJet_eta,&b_genJet_eta);
  tInput->SetBranchAddress("genJet_phi",&genJet_phi,&b_genJet_phi);
  tInput->SetBranchAddress("genJet_E",&genJet_E,&b_genJet_E);

  cout << "Creating output tree"<<endl;

  fOutput = new TFile("output.root", "RECREATE");
  tOutput = new TTree("Tree", "Tree");

  tOutput->Branch("mc_event",&mc_event,"mc_event/I");
  tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");
  tOutput->Branch("mc_totp4_px",&mc_totp4_px,"mc_totp4_px/F");
  tOutput->Branch("mc_totp4_py",&mc_totp4_py,"mc_totp4_py/F");
  tOutput->Branch("mc_totp4_pt",&mc_totp4_pt,"mc_totp4_pt/F");
  tOutput->Branch("mc_thad_pt",&mc_thad_pt,"mc_thad_pt/F");
  tOutput->Branch("mc_thad_b_pt",&mc_thad_b_pt,"mc_thad_b_pt/F");
  tOutput->Branch("mc_thad_b_eta",&mc_thad_b_eta,"mc_thad_b_eta/F");
  tOutput->Branch("mc_thad_j1_pt",&mc_thad_j1_pt,"mc_thad_j1_pt/F");
  tOutput->Branch("mc_thad_j1_eta",&mc_thad_j1_eta,"mc_thad_j1_eta/F");
  tOutput->Branch("mc_thad_j2_pt",&mc_thad_j2_pt,"mc_thad_j2_pt/F");
  tOutput->Branch("mc_thad_j2_eta",&mc_thad_j2_eta,"mc_thad_j2_eta/F");
  tOutput->Branch("mc_tlep_pt",&mc_tlep_pt,"mc_tlep_pt/F");
  tOutput->Branch("mc_tlep_b_pt",&mc_tlep_b_pt,"mc_tlep_b_pt/F");
  tOutput->Branch("mc_tlep_b_eta",&mc_tlep_b_eta,"mc_tlep_b_eta/F");
  tOutput->Branch("mc_tlep_l_pt",&mc_tlep_l_pt,"mc_tlep_l_pt/F");
  tOutput->Branch("mc_tlep_l_eta",&mc_tlep_l_eta,"mc_tlep_l_eta/F");
  tOutput->Branch("mc_tlep2_pt",&mc_tlep2_pt,"mc_tlep2_pt/F");
  tOutput->Branch("mc_tlep2_b_pt",&mc_tlep2_b_pt,"mc_tlep2_b_pt/F");
  tOutput->Branch("mc_tlep2_b_eta",&mc_tlep2_b_eta,"mc_tlep2_b_eta/F");
  tOutput->Branch("mc_tlep2_l_pt",&mc_tlep2_l_pt,"mc_tlep2_l_pt/F");
  tOutput->Branch("mc_tlep2_l_eta",&mc_tlep2_l_eta,"mc_tlep2_l_eta/F");
  tOutput->Branch("mc_ttbar_decay",&mc_ttbar_decay,"mc_ttbar_decay/I");
  tOutput->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
  tOutput->Branch("mc_boson_pt",&mc_boson_pt,"mc_boson_pt/F");
  tOutput->Branch("mc_boson_l1_pt",&mc_boson_l1_pt,"mc_boson_l1_pt/F");
  tOutput->Branch("mc_boson_l1_eta",&mc_boson_l1_eta,"mc_boson_l1_eta/F");
  tOutput->Branch("mc_boson_l2_pt",&mc_boson_l2_pt,"mc_boson_l2_pt/F");
  tOutput->Branch("mc_boson_l2_eta",&mc_boson_l2_eta,"mc_boson_l2_eta/F");
  tOutput->Branch("mc_boson_ll_mass",&mc_boson_ll_mass,"mc_boson_ll_mass/F");
  tOutput->Branch("mc_boson_ll_pt",&mc_boson_ll_pt,"mc_boson_ll_pt/F");
  tOutput->Branch("mc_boson_ll_dphi",&mc_boson_ll_dphi,"mc_boson_ll_dphi/F");
  tOutput->Branch("mc_boson_j1_pt",&mc_boson_j1_pt,"mc_boson_j1_pt/F");
  tOutput->Branch("mc_boson_j1_eta",&mc_boson_j1_eta,"mc_boson_j1_eta/F");
  tOutput->Branch("mc_boson_j2_pt",&mc_boson_j2_pt,"mc_boson_j2_pt/F");
  tOutput->Branch("mc_boson_j2_eta",&mc_boson_j2_eta,"mc_boson_j2_eta/F");
  tOutput->Branch("mc_boson_jj_mass",&mc_boson_jj_mass,"mc_boson_jj_mass/F");
  tOutput->Branch("mc_boson_jj_pt",&mc_boson_jj_pt,"mc_boson_jj_pt/F");
  tOutput->Branch("mc_boson_jj_dphi",&mc_boson_jj_dphi,"mc_boson_jj_dphi/F");
  tOutput->Branch("mc_met",&mc_met,"mc_met/F");
  tOutput->Branch("mc_njets25",&mc_njets25,"mc_njets25/I");
  tOutput->Branch("mc_3l_category",&mc_3l_category,"mc_3l_category/I");
  tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");
  tOutput->Branch("mc_hasLLcombZpeak",&mc_hasLLcombZpeak,"mc_hasLLcombZpeak/I");
  tOutput->Branch("mc_passMllGt12",&mc_passMllGt12,"mc_passMllGt12/I");
  tOutput->Branch("mc_passLepPresel",&mc_passLepPresel,"mc_passLepPresel/I");
  tOutput->Branch("mc_passJetPresel25",&mc_passJetPresel25,"mc_passJetPresel25/I");
  tOutput->Branch("mc_passBjetPresel25",&mc_passBjetPresel25,"mc_passBjetPresel25/I");

  cout << "multilepton variables"<<endl;
  tOutput->Branch("multilepton_Bjet1_Id",			        &multilepton_Bjet1_Id,			        "multilepton_Bjet1_Id/I");
  tOutput->Branch("multilepton_Bjet1_P4",			        "TLorentzVector",			            &multilepton_Bjet1_P4);
  tOutput->Branch("multilepton_Bjet1_CSV",                  &multilepton_Bjet1_CSV,                 "multilepton_Bjet1_CSV/F");
  tOutput->Branch("multilepton_Bjet1_DeltaR_Matched",     	&multilepton_Bjet1_DeltaR_Matched,  	"multilepton_Bjet1_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Bjet1_Label_Matched",      	&multilepton_Bjet1_Label_Matched,   	"multilepton_Bjet1_Label_Matched/I");
  tOutput->Branch("multilepton_Bjet1_Id_Matched",         	&multilepton_Bjet1_Id_Matched,      	"multilepton_Bjet1_Id_Matched/I");
  tOutput->Branch("multilepton_Bjet1_P4_Matched",         	"TLorentzVector",                   	&multilepton_Bjet1_P4_Matched);
  tOutput->Branch("multilepton_Bjet2_Id",			        &multilepton_Bjet2_Id,			        "multilepton_Bjet2_Id/I");
  tOutput->Branch("multilepton_Bjet2_P4",			        "TLorentzVector",			            &multilepton_Bjet2_P4);
  tOutput->Branch("multilepton_Bjet2_CSV",                  &multilepton_Bjet2_CSV,                 "multilepton_Bjet2_CSV/F");
  tOutput->Branch("multilepton_Bjet2_DeltaR_Matched",     	&multilepton_Bjet2_DeltaR_Matched,  	"multilepton_Bjet2_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Bjet2_Label_Matched",      	&multilepton_Bjet2_Label_Matched,   	"multilepton_Bjet2_Label_Matched/I");
  tOutput->Branch("multilepton_Bjet2_Id_Matched",         	&multilepton_Bjet2_Id_Matched,      	"multilepton_Bjet2_Id_Matched/I");
  tOutput->Branch("multilepton_Bjet2_P4_Matched",         	"TLorentzVector",                   	&multilepton_Bjet2_P4_Matched);
  tOutput->Branch("multilepton_Lepton1_Id",			        &multilepton_Lepton1_Id,		        "multilepton_Lepton1_Id/I");
  tOutput->Branch("multilepton_Lepton1_P4",			        "TLorentzVector",			            &multilepton_Lepton1_P4);
  tOutput->Branch("multilepton_Lepton1_DeltaR_Matched",    	&multilepton_Lepton1_DeltaR_Matched, 	"multilepton_Lepton1_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Lepton1_Label_Matched",    	&multilepton_Lepton1_Label_Matched, 	"multilepton_Lepton1_Label_Matched/I");
  tOutput->Branch("multilepton_Lepton1_Id_Matched",       	&multilepton_Lepton1_Id_Matched,    	"multilepton_Lepton1_Id_Matched/I");
  tOutput->Branch("multilepton_Lepton1_P4_Matched",       	"TLorentzVector",                   	&multilepton_Lepton1_P4_Matched);
  tOutput->Branch("multilepton_Lepton2_Id",			        &multilepton_Lepton2_Id,		        "multilepton_Lepton2_Id/I");
  tOutput->Branch("multilepton_Lepton2_P4",			        "TLorentzVector",			            &multilepton_Lepton2_P4);
  tOutput->Branch("multilepton_Lepton2_DeltaR_Matched",    	&multilepton_Lepton2_DeltaR_Matched, 	"multilepton_Lepton2_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Lepton2_Label_Matched",    	&multilepton_Lepton2_Label_Matched, 	"multilepton_Lepton2_Label_Matched/I");
  tOutput->Branch("multilepton_Lepton2_Id_Matched",       	&multilepton_Lepton2_Id_Matched,    	"multilepton_Lepton2_Id_Matched/I");
  tOutput->Branch("multilepton_Lepton2_P4_Matched",       	"TLorentzVector",                   	&multilepton_Lepton2_P4_Matched);
  tOutput->Branch("multilepton_Lepton3_Id",			        &multilepton_Lepton3_Id,		        "multilepton_Lepton3_Id/I");
  tOutput->Branch("multilepton_Lepton3_P4",			        "TLorentzVector",			            &multilepton_Lepton3_P4);
  tOutput->Branch("multilepton_Lepton3_DeltaR_Matched",    	&multilepton_Lepton3_DeltaR_Matched, 	"multilepton_Lepton3_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Lepton3_Label_Matched",   	&multilepton_Lepton3_Label_Matched, 	"multilepton_Lepton3_Label_Matched/I");
  tOutput->Branch("multilepton_Lepton3_Id_Matched",       	&multilepton_Lepton3_Id_Matched,    	"multilepton_Lepton3_Id_Matched/I");
  tOutput->Branch("multilepton_Lepton3_P4_Matched",       	"TLorentzVector",                   	&multilepton_Lepton3_P4_Matched);
  tOutput->Branch("multilepton_Lepton4_Id",			        &multilepton_Lepton4_Id,		        "multilepton_Lepton4_Id/I");
  tOutput->Branch("multilepton_Lepton4_P4",			        "TLorentzVector",			            &multilepton_Lepton4_P4);
  tOutput->Branch("multilepton_Lepton4_DeltaR_Matched",    	&multilepton_Lepton4_DeltaR_Matched, 	"multilepton_Lepton4_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Lepton4_Label_Matched",    	&multilepton_Lepton4_Label_Matched, 	"multilepton_Lepton4_Label_Matched/I");
  tOutput->Branch("multilepton_Lepton4_Id_Matched",       	&multilepton_Lepton4_Id_Matched,    	"multilepton_Lepton4_Id_Matched/I");
  tOutput->Branch("multilepton_Lepton4_P4_Matched",       	"TLorentzVector",                   	&multilepton_Lepton4_P4_Matched);

  tOutput->Branch("multilepton_h0_Label",	            	&multilepton_h0_Label,       		    "multilepton_h0_Label/I");
  tOutput->Branch("multilepton_h0_Id",               		&multilepton_h0_Id,          		    "multilepton_h0_Id/I");
  tOutput->Branch("multilepton_h0_P4",               		"TLorentzVector",                       &multilepton_h0_P4);
  tOutput->Branch("multilepton_t1_Label",                       &multilepton_t1_Label,                  "multilepton_t1_Label/I");
  tOutput->Branch("multilepton_t1_Id",                          &multilepton_t1_Id,                     "multilepton_t1_Id/I");
  tOutput->Branch("multilepton_t1_P4",                          "TLorentzVector",                       &multilepton_t1_P4);
  tOutput->Branch("multilepton_t2_Label",                       &multilepton_t2_Label,                  "multilepton_t2_Label/I");
  tOutput->Branch("multilepton_t2_Id",                          &multilepton_t2_Id,                     "multilepton_t2_Id/I");
  tOutput->Branch("multilepton_t2_P4",                          "TLorentzVector",                       &multilepton_t2_P4);

  tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
  tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
  tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
  tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
  tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
  tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
  tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
  tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
  tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
  tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
  tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
  tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

  tOutput->Branch("multilepton_h0_Id",			      	&multilepton_h0_Id,        		"multilepton_h0_Id/I");
  tOutput->Branch("multilepton_h0_P4",             		"TLorentzVector",                       &multilepton_h0_P4);
  tOutput->Branch("multilepton_t1_Id",                          &multilepton_t1_Id,                     "multilepton_t1_Id/I");
  tOutput->Branch("multilepton_t1_P4",                          "TLorentzVector",                       &multilepton_t1_P4);
  tOutput->Branch("multilepton_t2_Id",                          &multilepton_t2_Id,                     "multilepton_t2_Id/I");
  tOutput->Branch("multilepton_t2_P4",                          "TLorentzVector",                       &multilepton_t2_P4);

  cout << "Tree created"<<endl;

}

void ReadGenFlatTree::InitializeMEMRun(string InputFileName){

  cout << "Opening input tree"<<endl;

  fInput = TFile::Open(InputFileName.c_str(),"READ");
  tInput = (TTree*)fInput->Get("Tree");

   multilepton_Bjet1_P4_ptr 		= 0;
   multilepton_Bjet1_P4_Matched_ptr 	= 0;
   multilepton_Bjet2_P4_ptr 		= 0;
   multilepton_Bjet2_P4_Matched_ptr     = 0;
   multilepton_Lepton1_P4_ptr 		= 0;
   multilepton_Lepton1_P4_Matched_ptr 	= 0;
   multilepton_Lepton2_P4_ptr 		= 0;
   multilepton_Lepton2_P4_Matched_ptr   = 0;
   multilepton_Lepton3_P4_ptr 		= 0;
   multilepton_Lepton3_P4_Matched_ptr   = 0;
   multilepton_Lepton4_P4_ptr 		= 0;
   multilepton_Lepton4_P4_Matched_ptr   = 0;
   multilepton_h0_P4_ptr		= 0;
   multilepton_t1_P4_ptr                = 0;
   multilepton_t2_P4_ptr                = 0;
   multilepton_JetHighestPt1_P4_ptr = 0;
   multilepton_JetHighestPt2_P4_ptr = 0;
   multilepton_JetClosestMw1_P4_ptr = 0;
   multilepton_JetClosestMw2_P4_ptr = 0;
   multilepton_JetLowestMjj1_P4_ptr = 0;
   multilepton_JetLowestMjj2_P4_ptr = 0;
   multilepton_JetHighestPt1_2ndPair_P4_ptr = 0;
   multilepton_JetHighestPt2_2ndPair_P4_ptr = 0;
   multilepton_JetClosestMw1_2ndPair_P4_ptr = 0;
   multilepton_JetClosestMw2_2ndPair_P4_ptr = 0;
   multilepton_JetLowestMjj1_2ndPair_P4_ptr = 0;
   multilepton_JetLowestMjj2_2ndPair_P4_ptr = 0;
   multilepton_mET_ptr = 0;
   multilepton_Ptot_ptr = 0;
   /*
   MEAllWeights_TTLL = new std::vector<double>;
   MEAllWeights_TTHfl = new std::vector<double>;
   MEAllWeights_TTHsl = new std::vector<double>;
   MEAllWeights_TTH = new std::vector<double>;
   MEAllWeights_TTW = new std::vector<double>;
   MEAllWeights_TTWJJ = new std::vector<double>;
   MEAllWeights_TTbarfl = new std::vector<double>;
   MEAllWeights_TTbarsl = new std::vector<double>;
   MEAllWeights_TTbar = new std::vector<double>;

   MEAllWeights_TTLL_log = new std::vector<float>;
   MEAllWeights_TTHfl_log = new std::vector<float>;
   MEAllWeights_TTHsl_log = new std::vector<float>;
   MEAllWeights_TTH_log = new std::vector<float>;
   MEAllWeights_TTW_log = new std::vector<float>;
   MEAllWeights_TTWJJ_log = new std::vector<float>;
   MEAllWeights_TTbarfl_log = new std::vector<float>;
   MEAllWeights_TTbarsl_log = new std::vector<float>;
   MEAllWeights_TTbar_log = new std::vector<float>;
   */
  tInput->SetBranchAddress("mc_event",&mc_event,&b_mc_event);

  tInput->SetBranchAddress("mc_weight",&mc_weight,&b_mc_weight);
  tInput->SetBranchAddress("weight_scale_muF0p5",&weight_scale_muF0p5,&b_weight_scale_muF0p5);
  tInput->SetBranchAddress("weight_scale_muF2",&weight_scale_muF2,&b_weight_scale_muF2);
  tInput->SetBranchAddress("weight_scale_muR0p5",&weight_scale_muR0p5,&b_weight_scale_muR0p5);
  tInput->SetBranchAddress("weight_scale_muR2",&weight_scale_muR2,&b_weight_scale_muR2);
  tInput->SetBranchAddress("weight_csv_down",&weight_csv_down,&b_weight_csv_down);
  tInput->SetBranchAddress("weight_csv_up",&weight_csv_up,&b_weight_csv_up);

  tInput->SetBranchAddress("weight",&weight,&b_weight);
  tInput->SetBranchAddress("PV_weight",&PV_weight,&b_PV_weight);
  tInput->SetBranchAddress("catJets",&catJets,&b_catJets);

  tInput->SetBranchAddress("is_2lss_TTH_SR",&is_2lss_TTH_SR,&b_is_2lss_TTH_SR);
  tInput->SetBranchAddress("is_3l_TTH_SR",&is_3l_TTH_SR,&b_is_3l_TTH_SR);
  tInput->SetBranchAddress("is_emu_TT_CR",&is_emu_TT_CR,&b_is_emu_TT_CR);
  //tInput->SetBranchAddress("is_Zl_CR",&is_Zl_CR,&b_is_Zl_CR);
  tInput->SetBranchAddress("is_3l_TTZ_CR",&is_3l_TTZ_CR,&b_is_3l_TTZ_CR);
  //tInput->SetBranchAddress("is_3l_WZrel_CR ",&is_3l_WZrel_CR,&b_is_3l_WZrel_CR);
  tInput->SetBranchAddress("is_3l_TZQ_SR",&is_3l_TZQ_SR,&b_is_3l_TZQ_SR);

  tInput->SetBranchAddress("is_2bTight",        &is_2bTight,        &b_is_2bTight);
  tInput->SetBranchAddress("is_2bTight_float",    &is_2bTight_float,  &b_is_2bTight_float);

  tInput->SetBranchAddress("cat_HtoWW",&cat_HtoWW,&b_cat_HtoWW);
  tInput->SetBranchAddress("cat_HtoZZ",&cat_HtoZZ,&b_cat_HtoZZ);
  tInput->SetBranchAddress("cat_Htott",&cat_Htott,&b_cat_Htott);

  tInput->SetBranchAddress("mc_3l_category",&mc_3l_category,&b_mc_3l_category);
  tInput->SetBranchAddress("mc_ttbar_decay",&mc_ttbar_decay,&b_mc_ttbar_decay);
  tInput->SetBranchAddress("mc_boson_decay",&mc_boson_decay,&b_mc_boson_decay);
  tInput->SetBranchAddress("mc_ttZhypAllowed",&mc_ttZhypAllowed,&b_mc_ttZhypAllowed);

  tInput->SetBranchAddress("multilepton_Lepton1_Id",			&multilepton_Lepton1_Id,			&b_multilepton_Lepton1_Id);
  tInput->SetBranchAddress("multilepton_Lepton1_P4",			&multilepton_Lepton1_P4_ptr,			&b_multilepton_Lepton1_P4);
  tInput->SetBranchAddress("multilepton_Lepton1_DeltaR_Matched",        &multilepton_Lepton1_DeltaR_Matched,            &b_multilepton_Lepton1_DeltaR_Matched);
  tInput->SetBranchAddress("multilepton_Lepton1_Label_Matched",         &multilepton_Lepton1_Label_Matched,             &b_multilepton_Lepton1_Label_Matched);
  tInput->SetBranchAddress("multilepton_Lepton1_Id_Matched",            &multilepton_Lepton1_Id_Matched,                &b_multilepton_Lepton1_Id_Matched);
  tInput->SetBranchAddress("multilepton_Lepton1_P4_Matched",            &multilepton_Lepton1_P4_Matched_ptr,            &b_multilepton_Lepton1_P4_Matched);
  tInput->SetBranchAddress("multilepton_Lepton2_Id",			&multilepton_Lepton2_Id,			&b_multilepton_Lepton2_Id);
  tInput->SetBranchAddress("multilepton_Lepton2_DeltaR_Matched",        &multilepton_Lepton2_DeltaR_Matched,            &b_multilepton_Lepton2_DeltaR_Matched);
  tInput->SetBranchAddress("multilepton_Lepton2_Label_Matched",         &multilepton_Lepton2_Label_Matched,             &b_multilepton_Lepton2_Label_Matched);
  tInput->SetBranchAddress("multilepton_Lepton2_Id_Matched",            &multilepton_Lepton2_Id_Matched,                &b_multilepton_Lepton2_Id_Matched);
  tInput->SetBranchAddress("multilepton_Lepton2_P4_Matched",            &multilepton_Lepton2_P4_Matched_ptr,            &b_multilepton_Lepton2_P4_Matched);
  tInput->SetBranchAddress("multilepton_Lepton2_P4",			&multilepton_Lepton2_P4_ptr,			&b_multilepton_Lepton2_P4);
  tInput->SetBranchAddress("multilepton_Lepton3_Id",			&multilepton_Lepton3_Id,			&b_multilepton_Lepton3_Id);
  tInput->SetBranchAddress("multilepton_Lepton3_P4",			&multilepton_Lepton3_P4_ptr,			&b_multilepton_Lepton3_P4);
  tInput->SetBranchAddress("multilepton_Lepton3_DeltaR_Matched",        &multilepton_Lepton3_DeltaR_Matched,            &b_multilepton_Lepton3_DeltaR_Matched);
  tInput->SetBranchAddress("multilepton_Lepton3_Label_Matched",         &multilepton_Lepton3_Label_Matched,             &b_multilepton_Lepton3_Label_Matched);
  tInput->SetBranchAddress("multilepton_Lepton3_Id_Matched",            &multilepton_Lepton3_Id_Matched,                &b_multilepton_Lepton3_Id_Matched);
  tInput->SetBranchAddress("multilepton_Lepton3_P4_Matched",            &multilepton_Lepton3_P4_Matched_ptr,            &b_multilepton_Lepton3_P4_Matched);
  tInput->SetBranchAddress("multilepton_Lepton4_Id",			&multilepton_Lepton4_Id,			&b_multilepton_Lepton4_Id);
  tInput->SetBranchAddress("multilepton_Lepton4_P4",			&multilepton_Lepton4_P4_ptr,			&b_multilepton_Lepton4_P4);
  tInput->SetBranchAddress("multilepton_Lepton4_DeltaR_Matched",        &multilepton_Lepton4_DeltaR_Matched,            &b_multilepton_Lepton4_DeltaR_Matched);
  tInput->SetBranchAddress("multilepton_Lepton4_Label_Matched",         &multilepton_Lepton4_Label_Matched,             &b_multilepton_Lepton4_Label_Matched);
  tInput->SetBranchAddress("multilepton_Lepton4_Id_Matched",            &multilepton_Lepton4_Id_Matched,                &b_multilepton_Lepton4_Id_Matched);
  tInput->SetBranchAddress("multilepton_Lepton4_P4_Matched",            &multilepton_Lepton4_P4_Matched_ptr,            &b_multilepton_Lepton4_P4_Matched);

  tInput->SetBranchAddress("multilepton_h0_Label",         		&multilepton_h0_Label,             		&b_multilepton_h0_Label);
  tInput->SetBranchAddress("multilepton_h0_Id",            		&multilepton_h0_Id,                		&b_multilepton_h0_Id);
  tInput->SetBranchAddress("multilepton_h0_P4",            		&multilepton_h0_P4_ptr,            		&b_multilepton_h0_P4);
  tInput->SetBranchAddress("multilepton_t1_Label",                      &multilepton_t1_Label,                          &b_multilepton_t1_Label);
  tInput->SetBranchAddress("multilepton_t1_Id",                         &multilepton_t1_Id,                             &b_multilepton_t1_Id);
  tInput->SetBranchAddress("multilepton_t1_P4",                         &multilepton_t1_P4_ptr,                         &b_multilepton_t1_P4);
  tInput->SetBranchAddress("multilepton_t2_Label",                      &multilepton_t2_Label,                          &b_multilepton_t2_Label);
  tInput->SetBranchAddress("multilepton_t2_Id",                         &multilepton_t2_Id,                             &b_multilepton_t2_Id);
  tInput->SetBranchAddress("multilepton_t2_P4",                         &multilepton_t2_P4_ptr,                         &b_multilepton_t2_P4);

  tInput->SetBranchAddress("multilepton_Bjet1_Id",			&multilepton_Bjet1_Id,				&b_multilepton_Bjet1_Id);
  tInput->SetBranchAddress("multilepton_Bjet1_P4",			&multilepton_Bjet1_P4_ptr,			&b_multilepton_Bjet1_P4);
  tInput->SetBranchAddress("multilepton_Bjet1_CSV",			&multilepton_Bjet1_CSV,				&b_multilepton_Bjet1_CSV);
  tInput->SetBranchAddress("multilepton_Bjet1_JEC_Up",			&multilepton_Bjet1_JEC_Up,			&b_multilepton_Bjet1_JEC_Up);
  tInput->SetBranchAddress("multilepton_Bjet1_JEC_Down",		&multilepton_Bjet1_JEC_Down,			&b_multilepton_Bjet1_JEC_Down);
  tInput->SetBranchAddress("multilepton_Bjet1_JER_Up",			&multilepton_Bjet1_JER_Up,			&b_multilepton_Bjet1_JER_Up);
  tInput->SetBranchAddress("multilepton_Bjet1_JER_Down",		&multilepton_Bjet1_JER_Down,			&b_multilepton_Bjet1_JER_Down);
  tInput->SetBranchAddress("multilepton_Bjet1_DeltaR_Matched",          &multilepton_Bjet1_DeltaR_Matched,              &b_multilepton_Bjet1_DeltaR_Matched);
  tInput->SetBranchAddress("multilepton_Bjet1_Label_Matched",           &multilepton_Bjet1_Label_Matched,               &b_multilepton_Bjet1_Label_Matched);
  tInput->SetBranchAddress("multilepton_Bjet1_Id_Matched",         	&multilepton_Bjet1_Id_Matched,          	&b_multilepton_Bjet1_Id_Matched);
  tInput->SetBranchAddress("multilepton_Bjet1_P4_Matched",         	&multilepton_Bjet1_P4_Matched_ptr,              &b_multilepton_Bjet1_P4_Matched);

  tInput->SetBranchAddress("multilepton_Bjet2_Id",			&multilepton_Bjet2_Id,				&b_multilepton_Bjet2_Id);
  tInput->SetBranchAddress("multilepton_Bjet2_P4",			&multilepton_Bjet2_P4_ptr,			&b_multilepton_Bjet2_P4);
  tInput->SetBranchAddress("multilepton_Bjet2_CSV",			&multilepton_Bjet2_CSV,				&b_multilepton_Bjet2_CSV);
  tInput->SetBranchAddress("multilepton_Bjet2_JEC_Up",			&multilepton_Bjet2_JEC_Up,			&b_multilepton_Bjet2_JEC_Up);
  tInput->SetBranchAddress("multilepton_Bjet2_JEC_Down",		&multilepton_Bjet2_JEC_Down,			&b_multilepton_Bjet2_JEC_Down);
  tInput->SetBranchAddress("multilepton_Bjet2_JER_Up",			&multilepton_Bjet2_JER_Up,			&b_multilepton_Bjet2_JER_Up);
  tInput->SetBranchAddress("multilepton_Bjet2_JER_Down",		&multilepton_Bjet2_JER_Down,			&b_multilepton_Bjet2_JER_Down);
  tInput->SetBranchAddress("multilepton_Bjet2_DeltaR_Matched",          &multilepton_Bjet2_DeltaR_Matched,              &b_multilepton_Bjet2_DeltaR_Matched);
  tInput->SetBranchAddress("multilepton_Bjet2_Label_Matched",           &multilepton_Bjet2_Label_Matched,               &b_multilepton_Bjet2_Label_Matched);
  tInput->SetBranchAddress("multilepton_Bjet2_Id_Matched",              &multilepton_Bjet2_Id_Matched,                  &b_multilepton_Bjet2_Id_Matched);
  tInput->SetBranchAddress("multilepton_Bjet2_P4_Matched",              &multilepton_Bjet2_P4_Matched_ptr,              &b_multilepton_Bjet2_P4_Matched);

  tInput->SetBranchAddress("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,&b_multilepton_JetHighestPt1_Id);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_P4",&multilepton_JetHighestPt1_P4_ptr,&b_multilepton_JetHighestPt1_P4);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_CSV",&multilepton_JetHighestPt1_CSV,&b_multilepton_JetHighestPt1_CSV);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_JEC_Up",&multilepton_JetHighestPt1_JEC_Up,&b_multilepton_JetHighestPt1_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_JEC_Down",&multilepton_JetHighestPt1_JEC_Down,&b_multilepton_JetHighestPt1_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_JER_Up",&multilepton_JetHighestPt1_JER_Up,&b_multilepton_JetHighestPt1_JER_Up);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_JER_Down",&multilepton_JetHighestPt1_JER_Down,&b_multilepton_JetHighestPt1_JER_Down);

  tInput->SetBranchAddress("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,&b_multilepton_JetHighestPt2_Id);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_P4",&multilepton_JetHighestPt2_P4_ptr,&b_multilepton_JetHighestPt2_P4);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_CSV",&multilepton_JetHighestPt2_CSV,&b_multilepton_JetHighestPt2_CSV);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_JEC_Up",&multilepton_JetHighestPt2_JEC_Up,&b_multilepton_JetHighestPt2_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_JEC_Down",&multilepton_JetHighestPt2_JEC_Down,&b_multilepton_JetHighestPt2_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_JER_Up",&multilepton_JetHighestPt2_JER_Up,&b_multilepton_JetHighestPt2_JER_Up);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_JER_Down",&multilepton_JetHighestPt2_JER_Down,&b_multilepton_JetHighestPt2_JER_Down);

  tInput->SetBranchAddress("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,&b_multilepton_JetClosestMw1_Id);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_P4",&multilepton_JetClosestMw1_P4_ptr,&b_multilepton_JetClosestMw1_P4);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_CSV",&multilepton_JetClosestMw1_CSV,&b_multilepton_JetClosestMw1_CSV);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_JEC_Up",&multilepton_JetClosestMw1_JEC_Up,&b_multilepton_JetClosestMw1_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_JEC_Down",&multilepton_JetClosestMw1_JEC_Down,&b_multilepton_JetClosestMw1_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_JER_Up",&multilepton_JetClosestMw1_JER_Up,&b_multilepton_JetClosestMw1_JER_Up);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_JER_Down",&multilepton_JetClosestMw1_JER_Down,&b_multilepton_JetClosestMw1_JER_Down);

  tInput->SetBranchAddress("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,&b_multilepton_JetClosestMw2_Id);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_P4",&multilepton_JetClosestMw2_P4_ptr,&b_multilepton_JetClosestMw2_P4);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_CSV",&multilepton_JetClosestMw2_CSV,&b_multilepton_JetClosestMw2_CSV);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_JEC_Up",&multilepton_JetClosestMw2_JEC_Up,&b_multilepton_JetClosestMw2_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_JEC_Down",&multilepton_JetClosestMw2_JEC_Down,&b_multilepton_JetClosestMw2_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_JER_Up",&multilepton_JetClosestMw2_JER_Up,&b_multilepton_JetClosestMw2_JER_Up);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_JER_Down",&multilepton_JetClosestMw2_JER_Down,&b_multilepton_JetClosestMw2_JER_Down);

  tInput->SetBranchAddress("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,&b_multilepton_JetLowestMjj1_Id);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_P4",&multilepton_JetLowestMjj1_P4_ptr,&b_multilepton_JetLowestMjj1_P4);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_CSV",&multilepton_JetLowestMjj1_CSV,&b_multilepton_JetLowestMjj1_CSV);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_JEC_Up",&multilepton_JetLowestMjj1_JEC_Up,&b_multilepton_JetLowestMjj1_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_JEC_Down",&multilepton_JetLowestMjj1_JEC_Down,&b_multilepton_JetLowestMjj1_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_JER_Up",&multilepton_JetLowestMjj1_JER_Up,&b_multilepton_JetLowestMjj1_JER_Up);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_JER_Down",&multilepton_JetLowestMjj1_JER_Down,&b_multilepton_JetLowestMjj1_JER_Down);

  tInput->SetBranchAddress("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,&b_multilepton_JetLowestMjj2_Id);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_P4",&multilepton_JetLowestMjj2_P4_ptr,&b_multilepton_JetLowestMjj2_P4);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_CSV",&multilepton_JetLowestMjj2_CSV,&b_multilepton_JetLowestMjj2_CSV);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_JEC_Up",&multilepton_JetLowestMjj2_JEC_Up,&b_multilepton_JetLowestMjj2_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_JEC_Down",&multilepton_JetLowestMjj2_JEC_Down,&b_multilepton_JetLowestMjj2_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_JER_Up",&multilepton_JetLowestMjj2_JER_Up,&b_multilepton_JetLowestMjj2_JER_Up);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_JER_Down",&multilepton_JetLowestMjj2_JER_Down,&b_multilepton_JetLowestMjj2_JER_Down);

  tInput->SetBranchAddress("multilepton_JetHighestPt1_2ndPair_Id",&multilepton_JetHighestPt1_2ndPair_Id,&b_multilepton_JetHighestPt1_2ndPair_Id);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_2ndPair_P4",&multilepton_JetHighestPt1_2ndPair_P4_ptr,&b_multilepton_JetHighestPt1_2ndPair_P4);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_2ndPair_CSV",&multilepton_JetHighestPt1_2ndPair_CSV,&b_multilepton_JetHighestPt1_2ndPair_CSV);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_2ndPair_JEC_Up",&multilepton_JetHighestPt1_2ndPair_JEC_Up,&b_multilepton_JetHighestPt1_2ndPair_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_2ndPair_JEC_Down",&multilepton_JetHighestPt1_2ndPair_JEC_Down,&b_multilepton_JetHighestPt1_2ndPair_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_2ndPair_JER_Up",&multilepton_JetHighestPt1_2ndPair_JER_Up,&b_multilepton_JetHighestPt1_2ndPair_JER_Up);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_2ndPair_JER_Down",&multilepton_JetHighestPt1_2ndPair_JER_Down,&b_multilepton_JetHighestPt1_2ndPair_JER_Down);

  tInput->SetBranchAddress("multilepton_JetHighestPt2_2ndPair_Id",&multilepton_JetHighestPt2_2ndPair_Id,&b_multilepton_JetHighestPt2_2ndPair_Id);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_2ndPair_P4",&multilepton_JetHighestPt2_2ndPair_P4_ptr,&b_multilepton_JetHighestPt2_2ndPair_P4);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_2ndPair_CSV",&multilepton_JetHighestPt2_2ndPair_CSV,&b_multilepton_JetHighestPt2_2ndPair_CSV);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_2ndPair_JEC_Up",&multilepton_JetHighestPt2_2ndPair_JEC_Up,&b_multilepton_JetHighestPt2_2ndPair_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_2ndPair_JEC_Down",&multilepton_JetHighestPt2_2ndPair_JEC_Down,&b_multilepton_JetHighestPt2_2ndPair_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_2ndPair_JER_Up",&multilepton_JetHighestPt2_2ndPair_JER_Up,&b_multilepton_JetHighestPt2_2ndPair_JER_Up);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_2ndPair_JER_Down",&multilepton_JetHighestPt2_2ndPair_JER_Down,&b_multilepton_JetHighestPt2_2ndPair_JER_Down);

  tInput->SetBranchAddress("multilepton_JetClosestMw1_2ndPair_Id",&multilepton_JetClosestMw1_2ndPair_Id,&b_multilepton_JetClosestMw1_2ndPair_Id);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_2ndPair_P4",&multilepton_JetClosestMw1_2ndPair_P4_ptr,&b_multilepton_JetClosestMw1_2ndPair_P4);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_2ndPair_CSV",&multilepton_JetClosestMw1_2ndPair_CSV,&b_multilepton_JetClosestMw1_2ndPair_CSV);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_2ndPair_JEC_Up",&multilepton_JetClosestMw1_2ndPair_JEC_Up,&b_multilepton_JetClosestMw1_2ndPair_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_2ndPair_JEC_Down",&multilepton_JetClosestMw1_2ndPair_JEC_Down,&b_multilepton_JetClosestMw1_2ndPair_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_2ndPair_JER_Up",&multilepton_JetClosestMw1_2ndPair_JER_Up,&b_multilepton_JetClosestMw1_2ndPair_JER_Up);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_2ndPair_JER_Down",&multilepton_JetClosestMw1_2ndPair_JER_Down,&b_multilepton_JetClosestMw1_2ndPair_JER_Down);

  tInput->SetBranchAddress("multilepton_JetClosestMw2_2ndPair_Id",&multilepton_JetClosestMw2_2ndPair_Id,&b_multilepton_JetClosestMw2_2ndPair_Id);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_2ndPair_P4",&multilepton_JetClosestMw2_2ndPair_P4_ptr,&b_multilepton_JetClosestMw2_2ndPair_P4);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_2ndPair_CSV",&multilepton_JetClosestMw2_2ndPair_CSV,&b_multilepton_JetClosestMw2_2ndPair_CSV);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_2ndPair_JEC_Up",&multilepton_JetClosestMw2_2ndPair_JEC_Up,&b_multilepton_JetClosestMw2_2ndPair_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_2ndPair_JEC_Down",&multilepton_JetClosestMw2_2ndPair_JEC_Down,&b_multilepton_JetClosestMw2_2ndPair_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_2ndPair_JER_Up",&multilepton_JetClosestMw2_2ndPair_JER_Up,&b_multilepton_JetClosestMw2_2ndPair_JER_Up);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_2ndPair_JER_Down",&multilepton_JetClosestMw2_2ndPair_JER_Down,&b_multilepton_JetClosestMw2_2ndPair_JER_Down);

  tInput->SetBranchAddress("multilepton_JetLowestMjj1_2ndPair_Id",&multilepton_JetLowestMjj1_2ndPair_Id,&b_multilepton_JetLowestMjj1_2ndPair_Id);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_2ndPair_P4",&multilepton_JetLowestMjj1_2ndPair_P4_ptr,&b_multilepton_JetLowestMjj1_2ndPair_P4);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_2ndPair_CSV",&multilepton_JetLowestMjj1_2ndPair_CSV,&b_multilepton_JetLowestMjj1_2ndPair_CSV);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_2ndPair_JEC_Up",&multilepton_JetLowestMjj1_2ndPair_JEC_Up,&b_multilepton_JetLowestMjj1_2ndPair_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_2ndPair_JEC_Down",&multilepton_JetLowestMjj1_2ndPair_JEC_Down,&b_multilepton_JetLowestMjj1_2ndPair_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_2ndPair_JER_Up",&multilepton_JetLowestMjj1_2ndPair_JER_Up,&b_multilepton_JetLowestMjj1_2ndPair_JER_Up);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_2ndPair_JER_Down",&multilepton_JetLowestMjj1_2ndPair_JER_Down,&b_multilepton_JetLowestMjj1_2ndPair_JER_Down);

  tInput->SetBranchAddress("multilepton_JetLowestMjj2_2ndPair_Id",&multilepton_JetLowestMjj2_2ndPair_Id,&b_multilepton_JetLowestMjj2_2ndPair_Id);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_2ndPair_P4",&multilepton_JetLowestMjj2_2ndPair_P4_ptr,&b_multilepton_JetLowestMjj2_2ndPair_P4);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_2ndPair_CSV",&multilepton_JetLowestMjj2_2ndPair_CSV,&b_multilepton_JetLowestMjj2_2ndPair_CSV);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_2ndPair_JEC_Up",&multilepton_JetLowestMjj2_2ndPair_JEC_Up,&b_multilepton_JetLowestMjj2_2ndPair_JEC_Up);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_2ndPair_JEC_Down",&multilepton_JetLowestMjj2_2ndPair_JEC_Down,&b_multilepton_JetLowestMjj2_2ndPair_JEC_Down);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_2ndPair_JER_Up",&multilepton_JetLowestMjj2_2ndPair_JER_Up,&b_multilepton_JetLowestMjj2_2ndPair_JER_Up);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_2ndPair_JER_Down",&multilepton_JetLowestMjj2_2ndPair_JER_Down,&b_multilepton_JetLowestMjj2_2ndPair_JER_Down);


  tInput->SetBranchAddress("multilepton_mET",&multilepton_mET_ptr,&b_multilepton_mET);
  tInput->SetBranchAddress("multilepton_mHT",&multilepton_mHT,&b_multilepton_mHT);
  tInput->SetBranchAddress("multilepton_mETcov00",&multilepton_mETcov00,&b_multilepton_mETcov00);
  tInput->SetBranchAddress("multilepton_mETcov01",&multilepton_mETcov01,&b_multilepton_mETcov01);
  tInput->SetBranchAddress("multilepton_mETcov10",&multilepton_mETcov10,&b_multilepton_mETcov10);
  tInput->SetBranchAddress("multilepton_mETcov11",&multilepton_mETcov11,&b_multilepton_mETcov11);
  tInput->SetBranchAddress("multilepton_Ptot",&multilepton_Ptot_ptr,&b_multilepton_Ptot);

  tInput->SetBranchAddress("nJet25_Recl",       &nJet25_Recl,       &b_nJet25_Recl);
  tInput->SetBranchAddress("max_Lep_eta",       &max_Lep_eta,       &b_max_Lep_eta);
  tInput->SetBranchAddress("MT_met_lep1",       &MT_met_lep1,       &b_MT_met_lep1);
  tInput->SetBranchAddress("mindr_lep1_jet",    &mindr_lep1_jet,    &b_mindr_lep1_jet);
  tInput->SetBranchAddress("mindr_lep2_jet",    &mindr_lep2_jet,    &b_mindr_lep2_jet);
  tInput->SetBranchAddress("LepGood_conePt0",   &LepGood_conePt0,   &b_LepGood_conePt0);
  tInput->SetBranchAddress("LepGood_conePt1",   &LepGood_conePt1,   &b_LepGood_conePt1);
  tInput->SetBranchAddress("met",               &met,               &b_met              );
  tInput->SetBranchAddress("avg_dr_jet",        &avg_dr_jet,        &b_avg_dr_jet       );
  tInput->SetBranchAddress("mhtJet25_Recl",     &mhtJet25_Recl,     &b_mhtJet25_Recl    );

  tInput->SetBranchAddress("signal_2lss_TT_MVA",&signal_2lss_TT_MVA,&b_signal_2lss_TT_MVA);
  tInput->SetBranchAddress("signal_2lss_TTV_MVA",&signal_2lss_TTV_MVA,&b_signal_2lss_TTV_MVA);
  tInput->SetBranchAddress("signal_3l_TT_MVA",&signal_3l_TT_MVA,&b_signal_3l_TT_MVA);
  tInput->SetBranchAddress("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,&b_signal_3l_TTV_MVA);

  cout << "Creating output tree"<<endl;

  fOutput = new TFile("output.root", "RECREATE");
  tOutput = new TTree("Tree", "Tree");

  tOutput->Branch("mc_event",&mc_event,"mc_event/I");
  tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");
  tOutput->Branch("weight_scale_muF0p5",&weight_scale_muF0p5,"weight_scale_muF0p5/F");
  tOutput->Branch("weight_scale_muF2",&weight_scale_muF2,"weight_scale_muF2/F");
  tOutput->Branch("weight_scale_muR0p5",&weight_scale_muR0p5,"weight_scale_muR0p5/F");
  tOutput->Branch("weight_scale_muR2",&weight_scale_muR2,"weight_scale_muR2/F");
  tOutput->Branch("weight_csv_down",&weight_csv_down,"weight_csv_down/F");
  tOutput->Branch("weight_csv_up",&weight_csv_up,"weight_csv_up/F");

  tOutput->Branch("weight",&weight,"weight/F");
  tOutput->Branch("PV_weight",&PV_weight,"PV_weight/F");

  tOutput->Branch("catJets",&catJets,"catJets/I");
  tOutput->Branch("is_2lss_TTH_SR",&is_2lss_TTH_SR,"is_2lss_TTH_SR/B");
  tOutput->Branch("is_3l_TTH_SR",&is_3l_TTH_SR,"is_3l_TTH_SR/B");
  tOutput->Branch("is_emu_TT_CR",&is_emu_TT_CR,"is_emu_TT_CR/B");
  //tOutput->Branch("is_Zl_CR",&is_Zl_CR,"is_Zl_CR/B");
  tOutput->Branch("is_3l_TTZ_CR",&is_3l_TTZ_CR,"is_3l_TTZ_CR/B");
  //tOutput->Branch("is_3l_WZrel_CR",&is_3l_WZrel_CR,"is_3l_WZrel_CR/B");
  tOutput->Branch("is_3l_TZQ_SR",&is_3l_TZQ_SR,"is_3l_TZQ_SR/B");

  tOutput->Branch("is_2bTight",       &is_2bTight,      "is_2bTight/I");
  tOutput->Branch("is_2bTight_float", &is_2bTight_float,  "is_2bTight_float/F");

  tOutput->Branch("cat_HtoWW",&cat_HtoWW,"cat_HtoWW/B");
  tOutput->Branch("cat_HtoZZ",&cat_HtoZZ,"cat_HtoZZ/B");
  tOutput->Branch("cat_Htott",&cat_Htott,"cat_Htott/B");

  tOutput->Branch("mc_3l_category",&mc_3l_category,"mc_3l_category/I");
  tOutput->Branch("mc_ttbar_decay",&mc_ttbar_decay,"mc_ttbar_decay/I");
  tOutput->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
  tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");

  tOutput->Branch("multilepton_Bjet1_Id",                       &multilepton_Bjet1_Id,                  "multilepton_Bjet1_Id/I");
  tOutput->Branch("multilepton_Bjet1_P4",                       "TLorentzVector",                       &multilepton_Bjet1_P4);
  tOutput->Branch("multilepton_Bjet1_CSV",            &multilepton_Bjet1_CSV,         "multilepton_Bjet1_CSV/F");
  tOutput->Branch("multilepton_Bjet1_DeltaR_Matched",           &multilepton_Bjet1_DeltaR_Matched,      "multilepton_Bjet1_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Bjet1_Label_Matched",            &multilepton_Bjet1_Label_Matched,       "multilepton_Bjet1_Label_Matched/I");
  tOutput->Branch("multilepton_Bjet1_Id_Matched",               &multilepton_Bjet1_Id_Matched,          "multilepton_Bjet1_Id_Matched/I");
  tOutput->Branch("multilepton_Bjet1_P4_Matched",               "TLorentzVector",                       &multilepton_Bjet1_P4_Matched);
  tOutput->Branch("multilepton_Bjet2_Id",                       &multilepton_Bjet2_Id,                  "multilepton_Bjet2_Id/I");
  tOutput->Branch("multilepton_Bjet2_P4",                       "TLorentzVector",                       &multilepton_Bjet2_P4);
  tOutput->Branch("multilepton_Bjet2_CSV",            &multilepton_Bjet2_CSV,         "multilepton_Bjet2_CSV/F");
  tOutput->Branch("multilepton_Bjet2_DeltaR_Matched",           &multilepton_Bjet2_DeltaR_Matched,      "multilepton_Bjet2_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Bjet2_Label_Matched",            &multilepton_Bjet2_Label_Matched,       "multilepton_Bjet2_Label_Matched/I");
  tOutput->Branch("multilepton_Bjet2_Id_Matched",               &multilepton_Bjet2_Id_Matched,          "multilepton_Bjet2_Id_Matched/I");
  tOutput->Branch("multilepton_Bjet2_P4_Matched",               "TLorentzVector",                       &multilepton_Bjet2_P4_Matched);
  tOutput->Branch("multilepton_Lepton1_Id",                     &multilepton_Lepton1_Id,                "multilepton_Lepton1_Id/I");
  tOutput->Branch("multilepton_Lepton1_P4",                     "TLorentzVector",                       &multilepton_Lepton1_P4);
  tOutput->Branch("multilepton_Lepton1_DeltaR_Matched",         &multilepton_Lepton1_DeltaR_Matched,    "multilepton_Lepton1_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Lepton1_Label_Matched",          &multilepton_Lepton1_Label_Matched,     "multilepton_Lepton1_Label_Matched/I");
  tOutput->Branch("multilepton_Lepton1_Id_Matched",             &multilepton_Lepton1_Id_Matched,        "multilepton_Lepton1_Id_Matched/I");
  tOutput->Branch("multilepton_Lepton1_P4_Matched",             "TLorentzVector",                       &multilepton_Lepton1_P4_Matched);
  tOutput->Branch("multilepton_Lepton2_Id",                     &multilepton_Lepton2_Id,                "multilepton_Lepton2_Id/I");
  tOutput->Branch("multilepton_Lepton2_P4",                     "TLorentzVector",                       &multilepton_Lepton2_P4);
  tOutput->Branch("multilepton_Lepton2_DeltaR_Matched",         &multilepton_Lepton2_DeltaR_Matched,    "multilepton_Lepton2_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Lepton2_Label_Matched",          &multilepton_Lepton2_Label_Matched,     "multilepton_Lepton2_Label_Matched/I");
  tOutput->Branch("multilepton_Lepton2_Id_Matched",             &multilepton_Lepton2_Id_Matched,        "multilepton_Lepton2_Id_Matched/I");
  tOutput->Branch("multilepton_Lepton2_P4_Matched",             "TLorentzVector",                       &multilepton_Lepton2_P4_Matched);
  tOutput->Branch("multilepton_Lepton3_Id",                     &multilepton_Lepton3_Id,                "multilepton_Lepton3_Id/I");
  tOutput->Branch("multilepton_Lepton3_P4",                     "TLorentzVector",                       &multilepton_Lepton3_P4);
  tOutput->Branch("multilepton_Lepton3_DeltaR_Matched",         &multilepton_Lepton3_DeltaR_Matched,    "multilepton_Lepton3_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Lepton3_Label_Matched",          &multilepton_Lepton3_Label_Matched,     "multilepton_Lepton3_Label_Matched/I");
  tOutput->Branch("multilepton_Lepton3_Id_Matched",             &multilepton_Lepton3_Id_Matched,        "multilepton_Lepton3_Id_Matched/I");
  tOutput->Branch("multilepton_Lepton3_P4_Matched",             "TLorentzVector",                       &multilepton_Lepton3_P4_Matched);
  tOutput->Branch("multilepton_Lepton4_Id",                     &multilepton_Lepton4_Id,                "multilepton_Lepton4_Id/I");
  tOutput->Branch("multilepton_Lepton4_P4",                     "TLorentzVector",                       &multilepton_Lepton4_P4);
  tOutput->Branch("multilepton_Lepton4_DeltaR_Matched",         &multilepton_Lepton4_DeltaR_Matched,    "multilepton_Lepton4_DeltaR_Matched/F");
  tOutput->Branch("multilepton_Lepton4_Label_Matched",          &multilepton_Lepton4_Label_Matched,     "multilepton_Lepton4_Label_Matched/I");
  tOutput->Branch("multilepton_Lepton4_Id_Matched",             &multilepton_Lepton4_Id_Matched,        "multilepton_Lepton4_Id_Matched/I");
  tOutput->Branch("multilepton_Lepton4_P4_Matched",             "TLorentzVector",                       &multilepton_Lepton4_P4_Matched);

  tOutput->Branch("multilepton_h0_Label",                       &multilepton_h0_Label,                  "multilepton_h0_Label/I");
  tOutput->Branch("multilepton_h0_Id",                          &multilepton_h0_Id,                     "multilepton_h0_Id/I");
  tOutput->Branch("multilepton_h0_P4",                          "TLorentzVector",                       &multilepton_h0_P4);
  tOutput->Branch("multilepton_t1_Label",                       &multilepton_t1_Label,                  "multilepton_t1_Label/I");
  tOutput->Branch("multilepton_t1_Id",                          &multilepton_t1_Id,                     "multilepton_t1_Id/I");
  tOutput->Branch("multilepton_t1_P4",                          "TLorentzVector",                       &multilepton_t1_P4);
  tOutput->Branch("multilepton_t2_Label",                       &multilepton_t2_Label,                  "multilepton_t2_Label/I");
  tOutput->Branch("multilepton_t2_Id",                          &multilepton_t2_Id,                     "multilepton_t2_Id/I");
  tOutput->Branch("multilepton_t2_P4",                          "TLorentzVector",                       &multilepton_t2_P4);

  tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
  tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
  tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
  tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
  tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
  tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
  tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
  tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
  tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
  tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);

  tOutput->Branch("multilepton_JetHighestPt1_2ndPair_Id",&multilepton_JetHighestPt1_2ndPair_Id,"multilepton_JetHighestPt1_2ndPair_Id/I");
  tOutput->Branch("multilepton_JetHighestPt1_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt1_2ndPair_P4);
  tOutput->Branch("multilepton_JetHighestPt2_2ndPair_Id",&multilepton_JetHighestPt2_2ndPair_Id,"multilepton_JetHighestPt2_2ndPair_Id/I");
  tOutput->Branch("multilepton_JetHighestPt2_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt2_2ndPair_P4);
  tOutput->Branch("multilepton_JetClosestMw1_2ndPair_Id",&multilepton_JetClosestMw1_2ndPair_Id,"multilepton_JetClosestMw1_2ndPair_Id/I");
  tOutput->Branch("multilepton_JetClosestMw1_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw1_2ndPair_P4);
  tOutput->Branch("multilepton_JetClosestMw2_2ndPair_Id",&multilepton_JetClosestMw2_2ndPair_Id,"multilepton_JetClosestMw2_2ndPair_Id/I");
  tOutput->Branch("multilepton_JetClosestMw2_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw2_2ndPair_P4);
  tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_Id",&multilepton_JetLowestMjj1_2ndPair_Id,"multilepton_JetLowestMjj1_2ndPair_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj1_2ndPair_P4);
  tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_Id",&multilepton_JetLowestMjj2_2ndPair_Id,"multilepton_JetLowestMjj2_2ndPair_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj2_2ndPair_P4);

  tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
  tOutput->Branch("multilepton_mHT",&multilepton_mHT,"multilepton_mHT/F");
  tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

  tOutput->Branch("mc_mem_tthfl_weight",&mc_mem_tthfl_weight,"mc_mem_tthfl_weight/D");
  tOutput->Branch("mc_mem_tthfl_weight_JEC_up",&mc_mem_tthfl_weight_JEC_up,"mc_mem_tthfl_weight_JEC_up/D");
  tOutput->Branch("mc_mem_tthfl_weight_JEC_down",&mc_mem_tthfl_weight_JEC_down,"mc_mem_tthfl_weight_JEC_down/D");
  tOutput->Branch("mc_mem_tthfl_weight_JER_up",&mc_mem_tthfl_weight_JER_up,"mc_mem_tthfl_weight_JER_up/D");
  tOutput->Branch("mc_mem_tthfl_weight_JER_down",&mc_mem_tthfl_weight_JER_down,"mc_mem_tthfl_weight_JER_down/D");
  tOutput->Branch("mc_mem_tthfl_weight_log",&mc_mem_tthfl_weight_log,"mc_mem_tthfl_weight_log/D");
  tOutput->Branch("mc_mem_tthfl_weight_err",&mc_mem_tthfl_weight_err,"mc_mem_tthfl_weight_err/D");
  tOutput->Branch("mc_mem_tthfl_weight_chi2",&mc_mem_tthfl_weight_chi2,"mc_mem_tthfl_weight_chi2/F");
  tOutput->Branch("mc_mem_tthfl_weight_time",&mc_mem_tthfl_weight_time,"mc_mem_tthfl_weight_time/F");
  tOutput->Branch("mc_mem_tthfl_weight_max",&mc_mem_tthfl_weight_max,"mc_mem_tthfl_weight_max/D");
  tOutput->Branch("mc_mem_tthfl_weight_avg",&mc_mem_tthfl_weight_avg,"mc_mem_tthfl_weight_avg/D");
  tOutput->Branch("mc_mem_tthfl_weight_logmean",&mc_mem_tthfl_weight_logmean,"mc_mem_tthfl_weight_logmean/D");
  tOutput->Branch("mc_mem_tthfl_weight_kinmax",&mc_mem_tthfl_weight_kinmax,"mc_mem_tthfl_weight_kinmax/D");
  tOutput->Branch("mc_mem_tthfl_weight_kinmaxint",&mc_mem_tthfl_weight_kinmaxint,"mc_mem_tthfl_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_tthfl_weight_logmax",&mc_kin_tthfl_weight_logmax,"mc_kin_tthfl_weight_logmax/D");
  tOutput->Branch("mc_kin_tthfl_weight_logmaxint",&mc_kin_tthfl_weight_logmaxint,"mc_kin_tthfl_weight_logmaxint/D");

  tOutput->Branch("mc_kin_tthfl_tophad_P4","TLorentzVector",&mc_kin_tthfl_tophad_P4);
  tOutput->Branch("mc_kin_tthfl_tophad_Bjet_P4","TLorentzVector",&mc_kin_tthfl_tophad_Bjet_P4);
  tOutput->Branch("mc_kin_tthfl_tophad_W_P4","TLorentzVector",&mc_kin_tthfl_tophad_W_P4);
  tOutput->Branch("mc_kin_tthfl_tophad_Jet1_P4","TLorentzVector",&mc_kin_tthfl_tophad_Jet1_P4);
  tOutput->Branch("mc_kin_tthfl_tophad_Jet2_P4","TLorentzVector",&mc_kin_tthfl_tophad_Jet2_P4);
  tOutput->Branch("mc_kin_tthfl_tophad_Pt",&mc_kin_tthfl_tophad_Pt,"mc_kin_tthfl_tophad_Pt/F");
  tOutput->Branch("mc_kin_tthfl_tophad_Wmass",&mc_kin_tthfl_tophad_Wmass,"mc_kin_tthfl_tophad_Wmass/F");
  tOutput->Branch("mc_kin_tthfl_tophad_Benergy",&mc_kin_tthfl_tophad_Benergy,"mc_kin_tthfl_tophad_Benergy/F");
  tOutput->Branch("mc_kin_tthfl_tophad_Jet1energy",&mc_kin_tthfl_tophad_Jet1energy,"mc_kin_tthfl_tophad_Jet1energy/F");
  tOutput->Branch("mc_kin_tthfl_tophad_Jet2energy",&mc_kin_tthfl_tophad_Jet2energy,"mc_kin_tthfl_tophad_Jet2energy/F");
  tOutput->Branch("mc_kin_tthfl_toplep_P4","TLorentzVector",&mc_kin_tthfl_toplep_P4);
  tOutput->Branch("mc_kin_tthfl_toplep_Bjet_P4","TLorentzVector",&mc_kin_tthfl_toplep_Bjet_P4);
  tOutput->Branch("mc_kin_tthfl_toplep_W_P4","TLorentzVector",&mc_kin_tthfl_toplep_W_P4);
  tOutput->Branch("mc_kin_tthfl_toplep_Lep_P4","TLorentzVector",&mc_kin_tthfl_toplep_Lep_P4);
  tOutput->Branch("mc_kin_tthfl_toplep_Neut_P4","TLorentzVector",&mc_kin_tthfl_toplep_Neut_P4);
  tOutput->Branch("mc_kin_tthfl_toplep_Pt",&mc_kin_tthfl_toplep_Pt,"mc_kin_tthfl_toplep_Pt/F");
  tOutput->Branch("mc_kin_tthfl_toplep_Wmass",&mc_kin_tthfl_toplep_Wmass,"mc_kin_tthfl_toplep_Wmass/F");
  tOutput->Branch("mc_kin_tthfl_toplep_Benergy",&mc_kin_tthfl_toplep_Benergy,"mc_kin_tthfl_toplep_Benergy/F");
  tOutput->Branch("mc_kin_tthfl_toplep_Neutenergy",&mc_kin_tthfl_toplep_Neutenergy,"mc_kin_tthfl_toplep_Neutenergy/F");
  tOutput->Branch("mc_kin_tthfl_toplep2_P4","TLorentzVector",&mc_kin_tthfl_toplep2_P4);
  tOutput->Branch("mc_kin_tthfl_toplep2_Bjet_P4","TLorentzVector",&mc_kin_tthfl_toplep2_Bjet_P4);
  tOutput->Branch("mc_kin_tthfl_toplep2_W_P4","TLorentzVector",&mc_kin_tthfl_toplep2_W_P4);
  tOutput->Branch("mc_kin_tthfl_toplep2_Lep_P4","TLorentzVector",&mc_kin_tthfl_toplep2_Lep_P4);
  tOutput->Branch("mc_kin_tthfl_toplep2_Neut_P4","TLorentzVector",&mc_kin_tthfl_toplep2_Neut_P4);
  tOutput->Branch("mc_kin_tthfl_toplep2_Pt",&mc_kin_tthfl_toplep2_Pt,"mc_kin_tthfl_toplep2_Pt/F");
  tOutput->Branch("mc_kin_tthfl_toplep2_Wmass",&mc_kin_tthfl_toplep2_Wmass,"mc_kin_tthfl_toplep2_Wmass/F");
  tOutput->Branch("mc_kin_tthfl_toplep2_Benergy",&mc_kin_tthfl_toplep2_Benergy,"mc_kin_tthfl_toplep2_Benergy/F");
  tOutput->Branch("mc_kin_tthfl_toplep2_Neutenergy",&mc_kin_tthfl_toplep2_Neutenergy,"mc_kin_tthfl_toplep2_Neutenergy/F");
  tOutput->Branch("mc_kin_tthfl_h2l2nu_P4","TLorentzVector",&mc_kin_tthfl_h2l2nu_P4);
  tOutput->Branch("mc_kin_tthfl_h2l2nu_W1_P4","TLorentzVector",&mc_kin_tthfl_h2l2nu_W1_P4);
  tOutput->Branch("mc_kin_tthfl_h2l2nu_W2_P4","TLorentzVector",&mc_kin_tthfl_h2l2nu_W2_P4);
  tOutput->Branch("mc_kin_tthfl_h2l2nu_Lep1_P4","TLorentzVector",&mc_kin_tthfl_h2l2nu_Lep1_P4);
  tOutput->Branch("mc_kin_tthfl_h2l2nu_Neut1_P4","TLorentzVector",&mc_kin_tthfl_h2l2nu_Neut1_P4);
  tOutput->Branch("mc_kin_tthfl_h2l2nu_Lep2_P4","TLorentzVector",&mc_kin_tthfl_h2l2nu_Lep2_P4);
  tOutput->Branch("mc_kin_tthfl_h2l2nu_Neut2_P4","TLorentzVector",&mc_kin_tthfl_h2l2nu_Neut2_P4);
  tOutput->Branch("mc_kin_tthfl_h2l2nu_Pt",&mc_kin_tthfl_h2l2nu_Pt,"mc_kin_tthfl_h2l2nu_Pt/F");
  tOutput->Branch("mc_kin_tthfl_h2l2nu_W1mass",&mc_kin_tthfl_h2l2nu_W1mass,"mc_kin_tthfl_h2l2nu_W1mass/F");
  tOutput->Branch("mc_kin_tthfl_h2l2nu_Neut1energy",&mc_kin_tthfl_h2l2nu_Neut1energy,"mc_kin_tthfl_h2l2nu_Neut1energy/F");
  tOutput->Branch("mc_kin_tthfl_h2l2nu_W2mass",&mc_kin_tthfl_h2l2nu_W2mass,"mc_kin_tthfl_h2l2nu_W2mass/F");
  tOutput->Branch("mc_kin_tthfl_h2l2nu_Neut2energy",&mc_kin_tthfl_h2l2nu_Neut2energy,"mc_kin_tthfl_h2l2nu_Neut2energy/F");

  tOutput->Branch("mc_mem_tthsl_weight",&mc_mem_tthsl_weight,"mc_mem_tthsl_weight/D");
  tOutput->Branch("mc_mem_tthsl_weight_JEC_up",&mc_mem_tthsl_weight_JEC_up,"mc_mem_tthsl_weight_JEC_up/D");
  tOutput->Branch("mc_mem_tthsl_weight_JEC_down",&mc_mem_tthsl_weight_JEC_down,"mc_mem_tthsl_weight_JEC_down/D");
  tOutput->Branch("mc_mem_tthsl_weight_JER_up",&mc_mem_tthsl_weight_JER_up,"mc_mem_tthsl_weight_JER_up/D");
  tOutput->Branch("mc_mem_tthsl_weight_JER_down",&mc_mem_tthsl_weight_JER_down,"mc_mem_tthsl_weight_JER_down/D");
  tOutput->Branch("mc_mem_tthsl_weight_log",&mc_mem_tthsl_weight_log,"mc_mem_tthsl_weight_log/D");
  tOutput->Branch("mc_mem_tthsl_weight_err",&mc_mem_tthsl_weight_err,"mc_mem_tthsl_weight_err/D");
  tOutput->Branch("mc_mem_tthsl_weight_chi2",&mc_mem_tthsl_weight_chi2,"mc_mem_tthsl_weight_chi2/F");
  tOutput->Branch("mc_mem_tthsl_weight_time",&mc_mem_tthsl_weight_time,"mc_mem_tthsl_weight_time/F");
  tOutput->Branch("mc_mem_tthsl_weight_max",&mc_mem_tthsl_weight_max,"mc_mem_tthsl_weight_max/D");
  tOutput->Branch("mc_mem_tthsl_weight_avg",&mc_mem_tthsl_weight_avg,"mc_mem_tthsl_weight_avg/D");
  tOutput->Branch("mc_mem_tthsl_weight_logmean",&mc_mem_tthsl_weight_logmean,"mc_mem_tthsl_weight_logmean/D");
  tOutput->Branch("mc_mem_tthsl_weight_kinmax",&mc_mem_tthsl_weight_kinmax,"mc_mem_tthsl_weight_kinmax/D");
  tOutput->Branch("mc_mem_tthsl_weight_kinmaxint",&mc_mem_tthsl_weight_kinmaxint,"mc_mem_tthsl_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_tthsl_weight_logmax",&mc_kin_tthsl_weight_logmax,"mc_kin_tthsl_weight_logmax/D");
  tOutput->Branch("mc_kin_tthsl_weight_logmaxint",&mc_kin_tthsl_weight_logmaxint,"mc_kin_tthsl_weight_logmaxint/D");

  tOutput->Branch("mc_kin_tthsl_tophad_P4","TLorentzVector",&mc_kin_tthsl_tophad_P4);
  tOutput->Branch("mc_kin_tthsl_tophad_Bjet_P4","TLorentzVector",&mc_kin_tthsl_tophad_Bjet_P4);
  tOutput->Branch("mc_kin_tthsl_tophad_W_P4","TLorentzVector",&mc_kin_tthsl_tophad_W_P4);
  tOutput->Branch("mc_kin_tthsl_tophad_Jet1_P4","TLorentzVector",&mc_kin_tthsl_tophad_Jet1_P4);
  tOutput->Branch("mc_kin_tthsl_tophad_Jet2_P4","TLorentzVector",&mc_kin_tthsl_tophad_Jet2_P4);
  tOutput->Branch("mc_kin_tthsl_tophad_Pt",&mc_kin_tthsl_tophad_Pt,"mc_kin_tthsl_tophad_Pt/F");
  tOutput->Branch("mc_kin_tthsl_tophad_Wmass",&mc_kin_tthsl_tophad_Wmass,"mc_kin_tthsl_tophad_Wmass/F");
  tOutput->Branch("mc_kin_tthsl_tophad_Benergy",&mc_kin_tthsl_tophad_Benergy,"mc_kin_tthsl_tophad_Benergy/F");
  tOutput->Branch("mc_kin_tthsl_tophad_Jet1energy",&mc_kin_tthsl_tophad_Jet1energy,"mc_kin_tthsl_tophad_Jet1energy/F");
  tOutput->Branch("mc_kin_tthsl_tophad_Jet2energy",&mc_kin_tthsl_tophad_Jet2energy,"mc_kin_tthsl_tophad_Jet2energy/F");
  tOutput->Branch("mc_kin_tthsl_toplep_P4","TLorentzVector",&mc_kin_tthsl_toplep_P4);
  tOutput->Branch("mc_kin_tthsl_toplep_Bjet_P4","TLorentzVector",&mc_kin_tthsl_toplep_Bjet_P4);
  tOutput->Branch("mc_kin_tthsl_toplep_W_P4","TLorentzVector",&mc_kin_tthsl_toplep_W_P4);
  tOutput->Branch("mc_kin_tthsl_toplep_Lep_P4","TLorentzVector",&mc_kin_tthsl_toplep_Lep_P4);
  tOutput->Branch("mc_kin_tthsl_toplep_Neut_P4","TLorentzVector",&mc_kin_tthsl_toplep_Neut_P4);
  tOutput->Branch("mc_kin_tthsl_toplep_Pt",&mc_kin_tthsl_toplep_Pt,"mc_kin_tthsl_toplep_Pt/F");
  tOutput->Branch("mc_kin_tthsl_toplep_Wmass",&mc_kin_tthsl_toplep_Wmass,"mc_kin_tthsl_toplep_Wmass/F");
  tOutput->Branch("mc_kin_tthsl_toplep_Benergy",&mc_kin_tthsl_toplep_Benergy,"mc_kin_tthsl_toplep_Benergy/F");
  tOutput->Branch("mc_kin_tthsl_toplep_Neutenergy",&mc_kin_tthsl_toplep_Neutenergy,"mc_kin_tthsl_toplep_Neutenergy/F");
  tOutput->Branch("mc_kin_tthsl_toplep2_P4","TLorentzVector",&mc_kin_tthsl_toplep2_P4);
  tOutput->Branch("mc_kin_tthsl_toplep2_Bjet_P4","TLorentzVector",&mc_kin_tthsl_toplep2_Bjet_P4);
  tOutput->Branch("mc_kin_tthsl_toplep2_W_P4","TLorentzVector",&mc_kin_tthsl_toplep2_W_P4);
  tOutput->Branch("mc_kin_tthsl_toplep2_Lep_P4","TLorentzVector",&mc_kin_tthsl_toplep2_Lep_P4);
  tOutput->Branch("mc_kin_tthsl_toplep2_Neut_P4","TLorentzVector",&mc_kin_tthsl_toplep2_Neut_P4);
  tOutput->Branch("mc_kin_tthsl_toplep2_Pt",&mc_kin_tthsl_toplep2_Pt,"mc_kin_tthsl_toplep2_Pt/F");
  tOutput->Branch("mc_kin_tthsl_toplep2_Wmass",&mc_kin_tthsl_toplep2_Wmass,"mc_kin_tthsl_toplep2_Wmass/F");
  tOutput->Branch("mc_kin_tthsl_toplep2_Benergy",&mc_kin_tthsl_toplep2_Benergy,"mc_kin_tthsl_toplep2_Benergy/F");
  tOutput->Branch("mc_kin_tthsl_toplep2_Neutenergy",&mc_kin_tthsl_toplep2_Neutenergy,"mc_kin_tthsl_toplep2_Neutenergy/F");
  tOutput->Branch("mc_kin_tthsl_hlnujj_P4","TLorentzVector",&mc_kin_tthsl_hlnujj_P4);
  tOutput->Branch("mc_kin_tthsl_hlnujj_W1_P4","TLorentzVector",&mc_kin_tthsl_hlnujj_W1_P4);
  tOutput->Branch("mc_kin_tthsl_hlnujj_W2_P4","TLorentzVector",&mc_kin_tthsl_hlnujj_W2_P4);
  tOutput->Branch("mc_kin_tthsl_hlnujj_Lep_P4","TLorentzVector",&mc_kin_tthsl_hlnujj_Lep_P4);
  tOutput->Branch("mc_kin_tthsl_hlnujj_Neut_P4","TLorentzVector",&mc_kin_tthsl_hlnujj_Neut_P4);
  tOutput->Branch("mc_kin_tthsl_hlnujj_Jet1_P4","TLorentzVector",&mc_kin_tthsl_hlnujj_Jet1_P4);
  tOutput->Branch("mc_kin_tthsl_hlnujj_Jet2_P4","TLorentzVector",&mc_kin_tthsl_hlnujj_Jet2_P4);
  tOutput->Branch("mc_kin_tthsl_hlnujj_Pt",&mc_kin_tthsl_hlnujj_Pt,"mc_kin_tthsl_hlnujj_Pt/F");
  tOutput->Branch("mc_kin_tthsl_hlnujj_W1mass",&mc_kin_tthsl_hlnujj_W1mass,"mc_kin_tthsl_hlnujj_W1mass/F");
  tOutput->Branch("mc_kin_tthsl_hlnujj_Neut1energy",&mc_kin_tthsl_hlnujj_Neut1energy,"mc_kin_tthsl_hlnujj_Neut1energy/F");
  tOutput->Branch("mc_kin_tthsl_hlnujj_W2mass",&mc_kin_tthsl_hlnujj_W2mass,"mc_kin_tthsl_hlnujj_W2mass/F");
  tOutput->Branch("mc_kin_tthsl_hlnujj_Jet1energy",&mc_kin_tthsl_hlnujj_Jet1energy,"mc_kin_tthsl_hlnujj_Jet1energy/F");
  tOutput->Branch("mc_kin_tthsl_hlnujj_Jet2energy",&mc_kin_tthsl_hlnujj_Jet2energy,"mc_kin_tthsl_hlnujj_Jet2energy/F");

  tOutput->Branch("mc_mem_tth_weight",&mc_mem_tth_weight,"mc_mem_tth_weight/D");
  tOutput->Branch("mc_mem_tth_weight_JEC_up",&mc_mem_tth_weight_JEC_up,"mc_mem_tth_weight_JEC_up/D");
  tOutput->Branch("mc_mem_tth_weight_JEC_down",&mc_mem_tth_weight_JEC_down,"mc_mem_tth_weight_JEC_down/D");
  tOutput->Branch("mc_mem_tth_weight_JER_up",&mc_mem_tth_weight_JER_up,"mc_mem_tth_weight_JER_up/D");
  tOutput->Branch("mc_mem_tth_weight_JER_down",&mc_mem_tth_weight_JER_down,"mc_mem_tth_weight_JER_down/D");
  tOutput->Branch("mc_mem_tth_weight_log",&mc_mem_tth_weight_log,"mc_mem_tth_weight_log/D");
  tOutput->Branch("mc_mem_tth_weight_err",&mc_mem_tth_weight_err,"mc_mem_tth_weight_err/D");
  tOutput->Branch("mc_mem_tth_weight_chi2",&mc_mem_tth_weight_chi2,"mc_mem_tth_weight_chi2/F");
  tOutput->Branch("mc_mem_tth_weight_time",&mc_mem_tth_weight_time,"mc_mem_tth_weight_time/F");
  tOutput->Branch("mc_mem_tth_weight_max",&mc_mem_tth_weight_max,"mc_mem_tth_weight_max/D");
  tOutput->Branch("mc_mem_tth_weight_avg",&mc_mem_tth_weight_avg,"mc_mem_tth_weight_avg/D");
  tOutput->Branch("mc_mem_tth_weight_logmean",&mc_mem_tth_weight_logmean,"mc_mem_tth_weight_logmean/D");
  tOutput->Branch("mc_mem_tth_weight_kinmax",&mc_mem_tth_weight_kinmax,"mc_mem_tth_weight_kinmax/D");
  tOutput->Branch("mc_mem_tth_weight_kinmaxint",&mc_mem_tth_weight_kinmaxint,"mc_mem_tth_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_tth_weight_logmax",&mc_kin_tth_weight_logmax,"mc_kin_tth_weight_logmax/D");
  tOutput->Branch("mc_kin_tth_weight_logmaxint",&mc_kin_tth_weight_logmaxint,"mc_kin_tth_weight_logmaxint/D");

  tOutput->Branch("mc_mem_ttz_weight",&mc_mem_ttz_weight,"mc_mem_ttz_weight/D");
  tOutput->Branch("mc_mem_ttz_weight_JEC_up",&mc_mem_ttz_weight_JEC_up,"mc_mem_ttz_weight_JEC_up/D");
  tOutput->Branch("mc_mem_ttz_weight_JEC_down",&mc_mem_ttz_weight_JEC_down,"mc_mem_ttz_weight_JEC_down/D");
  tOutput->Branch("mc_mem_ttz_weight_JER_up",&mc_mem_ttz_weight_JER_up,"mc_mem_ttz_weight_JER_up/D");
  tOutput->Branch("mc_mem_ttz_weight_JER_down",&mc_mem_ttz_weight_JER_down,"mc_mem_ttz_weight_JER_down/D");
  tOutput->Branch("mc_mem_ttz_weight_log",&mc_mem_ttz_weight_log,"mc_mem_ttz_weight_log/D");
  tOutput->Branch("mc_mem_ttz_weight_err",&mc_mem_ttz_weight_err,"mc_mem_ttz_weight_err/D");
  tOutput->Branch("mc_mem_ttz_weight_chi2",&mc_mem_ttz_weight_chi2,"mc_mem_ttz_weight_chi2/F");
  tOutput->Branch("mc_mem_ttz_weight_time",&mc_mem_ttz_weight_time,"mc_mem_ttz_weight_time/F");
  tOutput->Branch("mc_mem_ttz_weight_max",&mc_mem_ttz_weight_max,"mc_mem_ttz_weight_max/D");
  tOutput->Branch("mc_mem_ttz_weight_avg",&mc_mem_ttz_weight_avg,"mc_mem_ttz_weight_avg/D");
  tOutput->Branch("mc_mem_ttz_weight_logmean",&mc_mem_ttz_weight_logmean,"mc_mem_ttz_weight_logmean/D");
  tOutput->Branch("mc_mem_ttz_weight_kinmax",&mc_mem_ttz_weight_kinmax,"mc_mem_ttz_weight_kinmax/D");
  tOutput->Branch("mc_mem_ttz_weight_kinmaxint",&mc_mem_ttz_weight_kinmaxint,"mc_mem_ttz_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_ttz_weight_logmax",&mc_kin_ttz_weight_logmax,"mc_kin_ttz_weight_logmax/D");
  tOutput->Branch("mc_kin_ttz_weight_logmaxint",&mc_kin_ttz_weight_logmaxint,"mc_kin_ttz_weight_logmaxint/D");

  tOutput->Branch("mc_kin_ttz_tophad_P4","TLorentzVector",&mc_kin_ttz_tophad_P4);
  tOutput->Branch("mc_kin_ttz_tophad_Bjet_P4","TLorentzVector",&mc_kin_ttz_tophad_Bjet_P4);
  tOutput->Branch("mc_kin_ttz_tophad_W_P4","TLorentzVector",&mc_kin_ttz_tophad_W_P4);
  tOutput->Branch("mc_kin_ttz_tophad_Jet1_P4","TLorentzVector",&mc_kin_ttz_tophad_Jet1_P4);
  tOutput->Branch("mc_kin_ttz_tophad_Jet2_P4","TLorentzVector",&mc_kin_ttz_tophad_Jet2_P4);
  tOutput->Branch("mc_kin_ttz_tophad_Pt",&mc_kin_ttz_tophad_Pt,"mc_kin_ttz_tophad_Pt/F");
  tOutput->Branch("mc_kin_ttz_tophad_Wmass",&mc_kin_ttz_tophad_Wmass,"mc_kin_ttz_tophad_Wmass/F");
  tOutput->Branch("mc_kin_ttz_tophad_Benergy",&mc_kin_ttz_tophad_Benergy,"mc_kin_ttz_tophad_Benergy/F");
  tOutput->Branch("mc_kin_ttz_tophad_Jet1energy",&mc_kin_ttz_tophad_Jet1energy,"mc_kin_ttz_tophad_Jet1energy/F");
  tOutput->Branch("mc_kin_ttz_tophad_Jet2energy",&mc_kin_ttz_tophad_Jet2energy,"mc_kin_ttz_tophad_Jet2energy/F");
  tOutput->Branch("mc_kin_ttz_toplep_P4","TLorentzVector",&mc_kin_ttz_toplep_P4);
  tOutput->Branch("mc_kin_ttz_toplep_Bjet_P4","TLorentzVector",&mc_kin_ttz_toplep_Bjet_P4);
  tOutput->Branch("mc_kin_ttz_toplep_W_P4","TLorentzVector",&mc_kin_ttz_toplep_W_P4);
  tOutput->Branch("mc_kin_ttz_toplep_Lep_P4","TLorentzVector",&mc_kin_ttz_toplep_Lep_P4);
  tOutput->Branch("mc_kin_ttz_toplep_Neut_P4","TLorentzVector",&mc_kin_ttz_toplep_Neut_P4);
  tOutput->Branch("mc_kin_ttz_toplep_Pt",&mc_kin_ttz_toplep_Pt,"mc_kin_ttz_toplep_Pt/F");
  tOutput->Branch("mc_kin_ttz_toplep_Wmass",&mc_kin_ttz_toplep_Wmass,"mc_kin_ttz_toplep_Wmass/F");
  tOutput->Branch("mc_kin_ttz_toplep_Benergy",&mc_kin_ttz_toplep_Benergy,"mc_kin_ttz_toplep_Benergy/F");
  tOutput->Branch("mc_kin_ttz_toplep_Neutenergy",&mc_kin_ttz_toplep_Neutenergy,"mc_kin_ttz_toplep_Neutenergy/F");
  tOutput->Branch("mc_kin_ttz_toplep2_P4","TLorentzVector",&mc_kin_ttz_toplep2_P4);
  tOutput->Branch("mc_kin_ttz_toplep2_Bjet_P4","TLorentzVector",&mc_kin_ttz_toplep2_Bjet_P4);
  tOutput->Branch("mc_kin_ttz_toplep2_W_P4","TLorentzVector",&mc_kin_ttz_toplep2_W_P4);
  tOutput->Branch("mc_kin_ttz_toplep2_Lep_P4","TLorentzVector",&mc_kin_ttz_toplep2_Lep_P4);
  tOutput->Branch("mc_kin_ttz_toplep2_Neut_P4","TLorentzVector",&mc_kin_ttz_toplep2_Neut_P4);
  tOutput->Branch("mc_kin_ttz_toplep2_Pt",&mc_kin_ttz_toplep2_Pt,"mc_kin_ttz_toplep2_Pt/F");
  tOutput->Branch("mc_kin_ttz_toplep2_Wmass",&mc_kin_ttz_toplep2_Wmass,"mc_kin_ttz_toplep2_Wmass/F");
  tOutput->Branch("mc_kin_ttz_toplep2_Benergy",&mc_kin_ttz_toplep2_Benergy,"mc_kin_ttz_toplep2_Benergy/F");
  tOutput->Branch("mc_kin_ttz_toplep2_Neutenergy",&mc_kin_ttz_toplep2_Neutenergy,"mc_kin_ttz_toplep2_Neutenergy/F");
  tOutput->Branch("mc_kin_ttz_zll_P4","TLorentzVector",&mc_kin_ttz_zll_P4);
  tOutput->Branch("mc_kin_ttz_zll_Lep1_P4","TLorentzVector",&mc_kin_ttz_zll_Lep1_P4);
  tOutput->Branch("mc_kin_ttz_zll_Lep2_P4","TLorentzVector",&mc_kin_ttz_zll_Lep2_P4);
  tOutput->Branch("mc_kin_ttz_zll_Pt",&mc_kin_ttz_zll_Pt,"mc_kin_ttz_zll_Pt/F");
  tOutput->Branch("mc_kin_ttz_zll_Zmass",&mc_kin_ttz_zll_Zmass,"mc_kin_ttz_zll_Zmass/F");

  tOutput->Branch("mc_mem_ttw_weight",&mc_mem_ttw_weight,"mc_mem_ttw_weight/D");
  tOutput->Branch("mc_mem_ttw_weight_JEC_up",&mc_mem_ttw_weight_JEC_up,"mc_mem_ttw_weight_JEC_up/D");
  tOutput->Branch("mc_mem_ttw_weight_JEC_down",&mc_mem_ttw_weight_JEC_down,"mc_mem_ttw_weight_JEC_down/D");
  tOutput->Branch("mc_mem_ttw_weight_JER_up",&mc_mem_ttw_weight_JER_up,"mc_mem_ttw_weight_JER_up/D");
  tOutput->Branch("mc_mem_ttw_weight_JER_down",&mc_mem_ttw_weight_JER_down,"mc_mem_ttw_weight_JER_down/D");
  tOutput->Branch("mc_mem_ttw_weight_log",&mc_mem_ttw_weight_log,"mc_mem_ttw_weight_log/D");
  tOutput->Branch("mc_mem_ttw_weight_err",&mc_mem_ttw_weight_err,"mc_mem_ttw_weight_err/D");
  tOutput->Branch("mc_mem_ttw_weight_chi2",&mc_mem_ttw_weight_chi2,"mc_mem_ttw_weight_chi2/F");
  tOutput->Branch("mc_mem_ttw_weight_time",&mc_mem_ttw_weight_time,"mc_mem_ttw_weight_time/F");
  tOutput->Branch("mc_mem_ttw_weight_max",&mc_mem_ttw_weight_max,"mc_mem_ttw_weight_max/D");
  tOutput->Branch("mc_mem_ttw_weight_avg",&mc_mem_ttw_weight_avg,"mc_mem_ttw_weight_avg/D");
  tOutput->Branch("mc_mem_ttw_weight_logmean",&mc_mem_ttw_weight_logmean,"mc_mem_ttw_weight_logmean/D");
  tOutput->Branch("mc_mem_ttw_weight_kinmax",&mc_mem_ttw_weight_kinmax,"mc_mem_ttw_weight_kinmax/D");
  tOutput->Branch("mc_mem_ttw_weight_kinmaxint",&mc_mem_ttw_weight_kinmaxint,"mc_mem_ttw_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_ttw_weight_logmax",&mc_kin_ttw_weight_logmax,"mc_kin_ttw_weight_logmax/D");
  tOutput->Branch("mc_kin_ttw_weight_logmaxint",&mc_kin_ttw_weight_logmaxint,"mc_kin_ttw_weight_logmaxint/D");

  tOutput->Branch("mc_kin_ttw_tophad_P4","TLorentzVector",&mc_kin_ttw_tophad_P4);
  tOutput->Branch("mc_kin_ttw_tophad_Bjet_P4","TLorentzVector",&mc_kin_ttw_tophad_Bjet_P4);
  tOutput->Branch("mc_kin_ttw_tophad_W_P4","TLorentzVector",&mc_kin_ttw_tophad_W_P4);
  tOutput->Branch("mc_kin_ttw_tophad_Jet1_P4","TLorentzVector",&mc_kin_ttw_tophad_Jet1_P4);
  tOutput->Branch("mc_kin_ttw_tophad_Jet2_P4","TLorentzVector",&mc_kin_ttw_tophad_Jet2_P4);
  tOutput->Branch("mc_kin_ttw_tophad_Pt",&mc_kin_ttw_tophad_Pt,"mc_kin_ttw_tophad_Pt/F");
  tOutput->Branch("mc_kin_ttw_tophad_Wmass",&mc_kin_ttw_tophad_Wmass,"mc_kin_ttw_tophad_Wmass/F");
  tOutput->Branch("mc_kin_ttw_tophad_Benergy",&mc_kin_ttw_tophad_Benergy,"mc_kin_ttw_tophad_Benergy/F");
  tOutput->Branch("mc_kin_ttw_tophad_Jet1energy",&mc_kin_ttw_tophad_Jet1energy,"mc_kin_ttw_tophad_Jet1energy/F");
  tOutput->Branch("mc_kin_ttw_tophad_Jet2energy",&mc_kin_ttw_tophad_Jet2energy,"mc_kin_ttw_tophad_Jet2energy/F");
  tOutput->Branch("mc_kin_ttw_toplep_P4","TLorentzVector",&mc_kin_ttw_toplep_P4);
  tOutput->Branch("mc_kin_ttw_toplep_Bjet_P4","TLorentzVector",&mc_kin_ttw_toplep_Bjet_P4);
  tOutput->Branch("mc_kin_ttw_toplep_W_P4","TLorentzVector",&mc_kin_ttw_toplep_W_P4);
  tOutput->Branch("mc_kin_ttw_toplep_Lep_P4","TLorentzVector",&mc_kin_ttw_toplep_Lep_P4);
  tOutput->Branch("mc_kin_ttw_toplep_Neut_P4","TLorentzVector",&mc_kin_ttw_toplep_Neut_P4);
  tOutput->Branch("mc_kin_ttw_toplep_Pt",&mc_kin_ttw_toplep_Pt,"mc_kin_ttw_toplep_Pt/F");
  tOutput->Branch("mc_kin_ttw_toplep_Wmass",&mc_kin_ttw_toplep_Wmass,"mc_kin_ttw_toplep_Wmass/F");
  tOutput->Branch("mc_kin_ttw_toplep_Benergy",&mc_kin_ttw_toplep_Benergy,"mc_kin_ttw_toplep_Benergy/F");
  tOutput->Branch("mc_kin_ttw_toplep_Neutenergy",&mc_kin_ttw_toplep_Neutenergy,"mc_kin_ttw_toplep_Neutenergy/F");
  tOutput->Branch("mc_kin_ttw_toplep2_P4","TLorentzVector",&mc_kin_ttw_toplep2_P4);
  tOutput->Branch("mc_kin_ttw_toplep2_Bjet_P4","TLorentzVector",&mc_kin_ttw_toplep2_Bjet_P4);
  tOutput->Branch("mc_kin_ttw_toplep2_W_P4","TLorentzVector",&mc_kin_ttw_toplep2_W_P4);
  tOutput->Branch("mc_kin_ttw_toplep2_Lep_P4","TLorentzVector",&mc_kin_ttw_toplep2_Lep_P4);
  tOutput->Branch("mc_kin_ttw_toplep2_Neut_P4","TLorentzVector",&mc_kin_ttw_toplep2_Neut_P4);
  tOutput->Branch("mc_kin_ttw_toplep2_Pt",&mc_kin_ttw_toplep2_Pt,"mc_kin_ttw_toplep2_Pt/F");
  tOutput->Branch("mc_kin_ttw_toplep2_Wmass",&mc_kin_ttw_toplep2_Wmass,"mc_kin_ttw_toplep2_Wmass/F");
  tOutput->Branch("mc_kin_ttw_toplep2_Benergy",&mc_kin_ttw_toplep2_Benergy,"mc_kin_ttw_toplep2_Benergy/F");
  tOutput->Branch("mc_kin_ttw_toplep2_Neutenergy",&mc_kin_ttw_toplep2_Neutenergy,"mc_kin_ttw_toplep2_Neutenergy/F");
  tOutput->Branch("mc_kin_ttw_wlnu_W_P4","TLorentzVector",&mc_kin_ttw_wlnu_W_P4);
  tOutput->Branch("mc_kin_ttw_wlnu_Lep_P4","TLorentzVector",&mc_kin_ttw_wlnu_Lep_P4);
  tOutput->Branch("mc_kin_ttw_wlnu_Neut_P4","TLorentzVector",&mc_kin_ttw_wlnu_Neut_P4);
  tOutput->Branch("mc_kin_ttw_wlnu_Pt",&mc_kin_ttw_wlnu_Pt,"mc_kin_ttw_wlnu_Pt/F");
  tOutput->Branch("mc_kin_ttw_wlnu_Wmass",&mc_kin_ttw_wlnu_Wmass,"mc_kin_ttw_wlnu_Wmass/F");
  tOutput->Branch("mc_kin_ttw_wlnu_Neutenergy",&mc_kin_ttw_wlnu_Neutenergy,"mc_kin_ttw_wlnu_Neutenergy/F");

  tOutput->Branch("mc_mem_ttwjj_weight",&mc_mem_ttwjj_weight,"mc_mem_ttwjj_weight/D");
  tOutput->Branch("mc_mem_ttwjj_weight_JEC_up",&mc_mem_ttwjj_weight_JEC_up,"mc_mem_ttwjj_weight_JEC_up/D");
  tOutput->Branch("mc_mem_ttwjj_weight_JEC_down",&mc_mem_ttwjj_weight_JEC_down,"mc_mem_ttwjj_weight_JEC_down/D");
  tOutput->Branch("mc_mem_ttwjj_weight_JER_up",&mc_mem_ttwjj_weight_JER_up,"mc_mem_ttwjj_weight_JER_up/D");
  tOutput->Branch("mc_mem_ttwjj_weight_JER_down",&mc_mem_ttwjj_weight_JER_down,"mc_mem_ttwjj_weight_JER_down/D");
  tOutput->Branch("mc_mem_ttwjj_weight_log",&mc_mem_ttwjj_weight_log,"mc_mem_ttwjj_weight_log/D");
  tOutput->Branch("mc_mem_ttwjj_weight_err",&mc_mem_ttwjj_weight_err,"mc_mem_ttwjj_weight_err/D");
  tOutput->Branch("mc_mem_ttwjj_weight_chi2",&mc_mem_ttwjj_weight_chi2,"mc_mem_ttwjj_weight_chi2/F");
  tOutput->Branch("mc_mem_ttwjj_weight_time",&mc_mem_ttwjj_weight_time,"mc_mem_ttwjj_weight_time/F");
  tOutput->Branch("mc_mem_ttwjj_weight_max",&mc_mem_ttwjj_weight_max,"mc_mem_ttwjj_weight_max/D");
  tOutput->Branch("mc_mem_ttwjj_weight_avg",&mc_mem_ttwjj_weight_avg,"mc_mem_ttwjj_weight_avg/D");
  tOutput->Branch("mc_mem_ttwjj_weight_logmean",&mc_mem_ttwjj_weight_logmean,"mc_mem_ttwjj_weight_logmean/D");
  tOutput->Branch("mc_mem_ttwjj_weight_kinmax",&mc_mem_ttwjj_weight_kinmax,"mc_mem_ttwjj_weight_kinmax/D");
  tOutput->Branch("mc_mem_ttwjj_weight_kinmaxint",&mc_mem_ttwjj_weight_kinmaxint,"mc_mem_ttwjj_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_ttwjj_weight_logmax",&mc_kin_ttwjj_weight_logmax,"mc_kin_ttwjj_weight_logmax/D");
  tOutput->Branch("mc_kin_ttwjj_weight_logmaxint",&mc_kin_ttwjj_weight_logmaxint,"mc_kin_ttwjj_weight_logmaxint/D");

  tOutput->Branch("mc_mem_ttbarfl_weight",&mc_mem_ttbarfl_weight,"mc_mem_ttbarfl_weight/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_JEC_up",&mc_mem_ttbarfl_weight_JEC_up,"mc_mem_ttbarfl_weight_JEC_up/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_JEC_down",&mc_mem_ttbarfl_weight_JEC_down,"mc_mem_ttbarfl_weight_JEC_down/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_JER_up",&mc_mem_ttbarfl_weight_JER_up,"mc_mem_ttbarfl_weight_JER_up/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_JER_down",&mc_mem_ttbarfl_weight_JER_down,"mc_mem_ttbarfl_weight_JER_down/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_log",&mc_mem_ttbarfl_weight_log,"mc_mem_ttbarfl_weight_log/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_err",&mc_mem_ttbarfl_weight_err,"mc_mem_ttbarfl_weight_err/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_chi2",&mc_mem_ttbarfl_weight_chi2,"mc_mem_ttbarfl_weight_chi2/F");
  tOutput->Branch("mc_mem_ttbarfl_weight_time",&mc_mem_ttbarfl_weight_time,"mc_mem_ttbarfl_weight_time/F");
  tOutput->Branch("mc_mem_ttbarfl_weight_max",&mc_mem_ttbarfl_weight_max,"mc_mem_ttbarfl_weight_max/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_avg",&mc_mem_ttbarfl_weight_avg,"mc_mem_ttbarfl_weight_avg/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_logmean",&mc_mem_ttbarfl_weight_logmean,"mc_mem_ttbarfl_weight_logmean/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_kinmax",&mc_mem_ttbarfl_weight_kinmax,"mc_mem_ttbarfl_weight_kinmax/D");
  tOutput->Branch("mc_mem_ttbarfl_weight_kinmaxint",&mc_mem_ttbarfl_weight_kinmaxint,"mc_mem_ttbarfl_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_ttbarfl_weight_logmax",&mc_kin_ttbarfl_weight_logmax,"mc_kin_ttbarfl_weight_logmax/D");
  tOutput->Branch("mc_kin_ttbarfl_weight_logmaxint",&mc_kin_ttbarfl_weight_logmaxint,"mc_kin_ttbarfl_weight_logmaxint/D");

  tOutput->Branch("mc_mem_ttbarsl_weight",&mc_mem_ttbarsl_weight,"mc_mem_ttbarsl_weight/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_JEC_up",&mc_mem_ttbarsl_weight_JEC_up,"mc_mem_ttbarsl_weight_JEC_up/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_JEC_down",&mc_mem_ttbarsl_weight_JEC_down,"mc_mem_ttbarsl_weight_JEC_down/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_JER_up",&mc_mem_ttbarsl_weight_JER_up,"mc_mem_ttbarsl_weight_JER_up/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_JER_down",&mc_mem_ttbarsl_weight_JER_down,"mc_mem_ttbarsl_weight_JER_down/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_log",&mc_mem_ttbarsl_weight_log,"mc_mem_ttbarsl_weight_log/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_err",&mc_mem_ttbarsl_weight_err,"mc_mem_ttbarsl_weight_err/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_chi2",&mc_mem_ttbarsl_weight_chi2,"mc_mem_ttbarsl_weight_chi2/F");
  tOutput->Branch("mc_mem_ttbarsl_weight_time",&mc_mem_ttbarsl_weight_time,"mc_mem_ttbarsl_weight_time/F");
  tOutput->Branch("mc_mem_ttbarsl_weight_max",&mc_mem_ttbarsl_weight_max,"mc_mem_ttbarsl_weight_max/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_avg",&mc_mem_ttbarsl_weight_avg,"mc_mem_ttbarsl_weight_avg/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_logmean",&mc_mem_ttbarsl_weight_logmean,"mc_mem_ttbarsl_weight_logmean/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_kinmax",&mc_mem_ttbarsl_weight_kinmax,"mc_mem_ttbarsl_weight_kinmax/D");
  tOutput->Branch("mc_mem_ttbarsl_weight_kinmaxint",&mc_mem_ttbarsl_weight_kinmaxint,"mc_mem_ttbarsl_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_ttbarsl_weight_logmax",&mc_kin_ttbarsl_weight_logmax,"mc_kin_ttbarsl_weight_logmax/D");
  tOutput->Branch("mc_kin_ttbarsl_weight_logmaxint",&mc_kin_ttbarsl_weight_logmaxint,"mc_kin_ttbarsl_weight_logmaxint/D");

  tOutput->Branch("mc_mem_ttbar_weight",&mc_mem_ttbar_weight,"mc_mem_ttbar_weight/D");
  tOutput->Branch("mc_mem_ttbar_weight_JEC_up",&mc_mem_ttbar_weight_JEC_up,"mc_mem_ttbar_weight_JEC_up/D");
  tOutput->Branch("mc_mem_ttbar_weight_JEC_down",&mc_mem_ttbar_weight_JEC_down,"mc_mem_ttbar_weight_JEC_down/D");
  tOutput->Branch("mc_mem_ttbar_weight_JER_up",&mc_mem_ttbar_weight_JER_up,"mc_mem_ttbar_weight_JER_up/D");
  tOutput->Branch("mc_mem_ttbar_weight_JER_down",&mc_mem_ttbar_weight_JER_down,"mc_mem_ttbar_weight_JER_down/D");
  tOutput->Branch("mc_mem_ttbar_weight_log",&mc_mem_ttbar_weight_log,"mc_mem_ttbar_weight_log/D");
  tOutput->Branch("mc_mem_ttbar_weight_err",&mc_mem_ttbar_weight_err,"mc_mem_ttbar_weight_err/D");
  tOutput->Branch("mc_mem_ttbar_weight_chi2",&mc_mem_ttbar_weight_chi2,"mc_mem_ttbar_weight_chi2/F");
  tOutput->Branch("mc_mem_ttbar_weight_time",&mc_mem_ttbar_weight_time,"mc_mem_ttbar_weight_time/F");
  tOutput->Branch("mc_mem_ttbar_weight_max",&mc_mem_ttbar_weight_max,"mc_mem_ttbar_weight_max/D");
  tOutput->Branch("mc_mem_ttbar_weight_avg",&mc_mem_ttbar_weight_avg,"mc_mem_ttbar_weight_avg/D");
  tOutput->Branch("mc_mem_ttbar_weight_logmean",&mc_mem_ttbar_weight_logmean,"mc_mem_ttbar_weight_logmean/D");
  tOutput->Branch("mc_mem_ttbar_weight_kinmax",&mc_mem_ttbar_weight_kinmax,"mc_mem_ttbar_weight_kinmax/D");
  tOutput->Branch("mc_mem_ttbar_weight_kinmaxint",&mc_mem_ttbar_weight_kinmaxint,"mc_mem_ttbar_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_ttbar_weight_logmax",&mc_kin_ttbar_weight_logmax,"mc_kin_ttbar_weight_logmax/D");
  tOutput->Branch("mc_kin_ttbar_weight_logmaxint",&mc_kin_ttbar_weight_logmaxint,"mc_kin_ttbar_weight_logmaxint/D");

  tOutput->Branch("mc_mem_tllj_weight",&mc_mem_tllj_weight,"mc_mem_tllj_weight/D");
  tOutput->Branch("mc_mem_tllj_weight_JEC_up",&mc_mem_tllj_weight_JEC_up,"mc_mem_tllj_weight_JEC_up/D");
  tOutput->Branch("mc_mem_tllj_weight_JEC_down",&mc_mem_tllj_weight_JEC_down,"mc_mem_tllj_weight_JEC_down/D");
  tOutput->Branch("mc_mem_tllj_weight_JER_up",&mc_mem_tllj_weight_JER_up,"mc_mem_tllj_weight_JER_up/D");
  tOutput->Branch("mc_mem_tllj_weight_JER_down",&mc_mem_tllj_weight_JER_down,"mc_mem_tllj_weight_JER_down/D");
  tOutput->Branch("mc_mem_tllj_weight_log",&mc_mem_tllj_weight_log,"mc_mem_tllj_weight_log/D");
  tOutput->Branch("mc_mem_tllj_weight_err",&mc_mem_tllj_weight_err,"mc_mem_tllj_weight_err/D");
  tOutput->Branch("mc_mem_tllj_weight_chi2",&mc_mem_tllj_weight_chi2,"mc_mem_tllj_weight_chi2/F");
  tOutput->Branch("mc_mem_tllj_weight_time",&mc_mem_tllj_weight_time,"mc_mem_tllj_weight_time/F");
  tOutput->Branch("mc_mem_tllj_weight_max",&mc_mem_tllj_weight_max,"mc_mem_tllj_weight_max/D");
  tOutput->Branch("mc_mem_tllj_weight_avg",&mc_mem_tllj_weight_avg,"mc_mem_tllj_weight_avg/D");
  tOutput->Branch("mc_mem_tllj_weight_logmean",&mc_mem_tllj_weight_logmean,"mc_mem_tllj_weight_logmean/D");
  tOutput->Branch("mc_mem_tllj_weight_kinmax",&mc_mem_tllj_weight_kinmax,"mc_mem_tllj_weight_kinmax/D");
  tOutput->Branch("mc_mem_tllj_weight_kinmaxint",&mc_mem_tllj_weight_kinmaxint,"mc_mem_tllj_weight_kinmaxint/D");
  tOutput->Branch("mc_kin_tllj_weight_logmax",&mc_kin_tllj_weight_logmax,"mc_kin_tllj_weight_logmax/D");
  tOutput->Branch("mc_kin_tllj_weight_logmaxint",&mc_kin_tllj_weight_logmaxint,"mc_kin_tllj_weight_logmaxint/D");


  tOutput->Branch("mc_mem_ttz_tth_likelihood",&mc_mem_ttz_tth_likelihood,"mc_mem_ttz_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttz_tth_likelihood_nlog",&mc_mem_ttz_tth_likelihood_nlog,"mc_mem_ttz_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttz_tth_likelihood_max",&mc_mem_ttz_tth_likelihood_max,"mc_mem_ttz_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttz_tth_likelihood_avg",&mc_mem_ttz_tth_likelihood_avg,"mc_mem_ttz_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttw_tth_likelihood",&mc_mem_ttw_tth_likelihood,"mc_mem_ttw_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttw_tth_likelihood_nlog",&mc_mem_ttw_tth_likelihood_nlog,"mc_mem_ttw_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttw_tth_likelihood_max",&mc_mem_ttw_tth_likelihood_max,"mc_mem_ttw_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttw_tth_likelihood_avg",&mc_mem_ttw_tth_likelihood_avg,"mc_mem_ttw_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttwjj_tth_likelihood",&mc_mem_ttwjj_tth_likelihood,"mc_mem_ttwjj_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttwjj_tth_likelihood_nlog",&mc_mem_ttwjj_tth_likelihood_nlog,"mc_mem_ttwjj_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttwjj_tth_likelihood_max",&mc_mem_ttwjj_tth_likelihood_max,"mc_mem_ttwjj_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttwjj_tth_likelihood_avg",&mc_mem_ttwjj_tth_likelihood_avg,"mc_mem_ttwjj_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttbar_tth_likelihood",&mc_mem_ttbar_tth_likelihood,"mc_mem_ttbar_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttbar_tth_likelihood_nlog",&mc_mem_ttbar_tth_likelihood_nlog,"mc_mem_ttbar_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttbar_tth_likelihood_max",&mc_mem_ttbar_tth_likelihood_max,"mc_mem_ttbar_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttbar_tth_likelihood_avg",&mc_mem_ttbar_tth_likelihood_avg,"mc_mem_ttbar_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttz_tllj_likelihood",&mc_mem_ttz_tllj_likelihood,"mc_mem_ttz_tllj_likelihood/D");
  tOutput->Branch("mc_mem_ttz_tllj_likelihood_nlog",&mc_mem_ttz_tllj_likelihood_nlog,"mc_mem_ttz_tllj_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttz_tllj_likelihood_max",&mc_mem_ttz_tllj_likelihood_max,"mc_mem_ttz_tllj_likelihood_max/D");
  tOutput->Branch("mc_mem_ttz_tllj_likelihood_avg",&mc_mem_ttz_tllj_likelihood_avg,"mc_mem_ttz_tllj_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttv_tth_likelihood",&mc_mem_ttv_tth_likelihood,"mc_mem_ttv_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttv_tth_likelihood_nlog",&mc_mem_ttv_tth_likelihood_nlog,"mc_mem_ttv_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttv_tth_likelihood_max",&mc_mem_ttv_tth_likelihood_max,"mc_mem_ttv_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttv_tth_likelihood_avg",&mc_mem_ttv_tth_likelihood_avg,"mc_mem_ttv_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttvjj_tth_likelihood",&mc_mem_ttvjj_tth_likelihood,"mc_mem_ttvjj_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttvjj_tth_likelihood_nlog",&mc_mem_ttvjj_tth_likelihood_nlog,"mc_mem_ttvjj_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttvjj_tth_likelihood_max",&mc_mem_ttvjj_tth_likelihood_max,"mc_mem_ttvjj_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttvjj_tth_likelihood_avg",&mc_mem_ttvjj_tth_likelihood_avg,"mc_mem_ttvjj_tth_likelihood_avg/D");

  tOutput->Branch("nJet25_Recl",     &nJet25_Recl,     "nJet25_Recl/F"     );
  tOutput->Branch("max_Lep_eta",     &max_Lep_eta,     "max_Lep_eta/F"     );
  tOutput->Branch("MT_met_lep1",     &MT_met_lep1,     "MT_met_lep1/F"     );
  tOutput->Branch("mindr_lep1_jet",  &mindr_lep1_jet,  "mindr_lep1_jet/F"  );
  tOutput->Branch("mindr_lep2_jet",  &mindr_lep2_jet,  "mindr_lep2_jet/F"  );
  tOutput->Branch("LepGood_conePt0", &LepGood_conePt0, "LepGood_conePt0/F" );
  tOutput->Branch("LepGood_conePt1", &LepGood_conePt1, "LepGood_conePt1/F" );
  tOutput->Branch("met",             &met,             "met/F"             );
  tOutput->Branch("avg_dr_jet",      &avg_dr_jet,      "avg_dr_jet/F"      );
  tOutput->Branch("mhtJet25_Recl",   &mhtJet25_Recl,   "mhtJet25_Recl/F"   );

  tOutput->Branch("signal_2lss_TT_MVA",&signal_2lss_TT_MVA,"signal_2lss_TT_MVA/F");
  tOutput->Branch("signal_2lss_TTV_MVA",&signal_2lss_TTV_MVA,"signal_2lss_TTV_MVA/F");
  tOutput->Branch("signal_3l_TT_MVA",&signal_3l_TT_MVA,"signal_3l_TT_MVA/F");
  tOutput->Branch("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,"signal_3l_TTV_MVA/F");
  /*
  tOutput->Branch("MEAllWeights_TTLL","vector<double>",&MEAllWeights_TTLL);
  tOutput->Branch("MEAllWeights_TTHfl","vector<double>",&MEAllWeights_TTHfl);
  tOutput->Branch("MEAllWeights_TTHsl","vector<double>",&MEAllWeights_TTHsl);
  tOutput->Branch("MEAllWeights_TTH","vector<double>",&MEAllWeights_TTH);
  tOutput->Branch("MEAllWeights_TTW","vector<double>",&MEAllWeights_TTW);
  tOutput->Branch("MEAllWeights_TTWJJ","vector<double>",&MEAllWeights_TTWJJ);
  tOutput->Branch("MEAllWeights_TTbarfl","vector<double>",&MEAllWeights_TTbarfl);
  tOutput->Branch("MEAllWeights_TTbarsl","vector<double>",&MEAllWeights_TTbarsl);
  tOutput->Branch("MEAllWeights_TTbar","vector<double>",&MEAllWeights_TTbar);
  tOutput->Branch("MEAllWeights_TLLJ","vector<double>",&MEAllWeights_TLLJ);

  tOutput->Branch("MEAllWeights_TTLL_log","vector<float>",&MEAllWeights_TTLL_log);
  tOutput->Branch("MEAllWeights_TTHfl_log","vector<float>",&MEAllWeights_TTHfl_log);
  tOutput->Branch("MEAllWeights_TTHsl_log","vector<float>",&MEAllWeights_TTHsl_log);
  tOutput->Branch("MEAllWeights_TTH_log","vector<float>",&MEAllWeights_TTH_log);
  tOutput->Branch("MEAllWeights_TTW_log","vector<float>",&MEAllWeights_TTW_log);
  tOutput->Branch("MEAllWeights_TTWJJ_log","vector<float>",&MEAllWeights_TTWJJ_log);
  tOutput->Branch("MEAllWeights_TTbarfl_log","vector<float>",&MEAllWeights_TTbarfl_log);
  tOutput->Branch("MEAllWeights_TTbarsl_log","vector<float>",&MEAllWeights_TTbarsl_log);
  tOutput->Branch("MEAllWeights_TTbar_log","vector<float>",&MEAllWeights_TTbar_log);
  tOutput->Branch("MEAllWeights_TLLJ_log","vector<float>",&MEAllWeights_TLLJ_log);
  */
  tOutput->Branch("multilepton_h0_Id",                          &multilepton_h0_Id,                     "multilepton_h0_Id/I");
  tOutput->Branch("multilepton_h0_P4",                          "TLorentzVector",                       &multilepton_h0_P4);
  tOutput->Branch("multilepton_t1_Id",                          &multilepton_t1_Id,                     "multilepton_t1_Id/I");
  tOutput->Branch("multilepton_t1_P4",                          "TLorentzVector",                       &multilepton_t1_P4);
  tOutput->Branch("multilepton_t2_Id",                          &multilepton_t2_Id,                     "multilepton_t2_Id/I");
  tOutput->Branch("multilepton_t2_P4",                          "TLorentzVector",                       &multilepton_t2_P4);
}


void ReadGenFlatTree::FillGenMultilepton(Long64_t iEvent, MultiLepton* multiLepton){

  cout << "FillGenMultilepton"<<endl;

  TLorentzVector Ptot, Pll, Pl1, PtotNeut(0,0,0,0), Pjj;

  tInput->LoadTree(iEvent);
  tInput->GetEntry(iEvent);

  mc_ttbar_decay = -1;
  mc_boson_decay = -1;

  if (mc_truth_t1_id==-666 || mc_truth_t2_id==-666) return;

  (*multiLepton).Leptons.clear();
  (*multiLepton).Jets.clear();
  (*multiLepton).Bjets.clear();
  (*multiLepton).AllJets.clear();
  (*multiLepton).JetsHighestPt.clear();
  (*multiLepton).JetsClosestMw.clear();
  (*multiLepton).JetsLowestMjj.clear();

    if (mc_truth_h0_id>-666 && mc_truth_h0Wl1_id>-666 && mc_truth_h0Wl2_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_h0Wl1_id, *mc_truth_h0Wl1_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_h0Wl2_id, *mc_truth_h0Wl2_p4);
      Ptot = (*mc_truth_h0_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pll = (*mc_truth_h0Wl1_p4) + (*mc_truth_h0Wl2_p4);
      mc_boson_decay = 0;
      mc_boson_pt = mc_truth_h0_p4->Pt();
      mc_boson_l1_pt = mc_truth_h0Wl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_h0Wl1_p4->Eta();
      mc_boson_l2_pt = mc_truth_h0Wl2_p4->Pt();
      mc_boson_l2_eta = mc_truth_h0Wl2_p4->Eta();
      mc_boson_ll_mass = Pll.M();
      mc_boson_ll_pt = Pll.Pt();
      mc_boson_ll_dphi = TMath::Abs(mc_truth_h0Wl1_p4->DeltaPhi(*mc_truth_h0Wl2_p4));
      PtotNeut += (*mc_truth_h0Wnu1_p4);
      PtotNeut += (*mc_truth_h0Wnu2_p4);
    }
    else if (mc_truth_h0_id>-666 && mc_truth_h0Wq11_id>-666 && mc_truth_h0Wl2_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_h0Wl2_id, *mc_truth_h0Wl2_p4);
      mc_boson_decay = 1;
      Ptot = (*mc_truth_h0_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pjj = (*mc_truth_h0Wq11_p4) + (*mc_truth_h0Wq12_p4);
      mc_boson_pt = mc_truth_h0_p4->Pt();
      mc_boson_l1_pt = mc_truth_h0Wl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_h0Wl1_p4->Eta();
      mc_boson_l2_pt = -666;
      mc_boson_l2_eta = -666;
      mc_boson_ll_mass = -666;
      mc_boson_ll_pt = -666;
      mc_boson_j1_pt = mc_truth_h0Wq11_p4->Pt();
      mc_boson_j1_eta = mc_truth_h0Wq11_p4->Eta();
      mc_boson_j2_pt = mc_truth_h0Wq12_p4->Pt();
      mc_boson_j2_eta = mc_truth_h0Wq12_p4->Eta();
      mc_boson_jj_mass = Pjj.M();
      mc_boson_jj_pt = Pjj.Pt();
      mc_boson_jj_dphi = TMath::Abs(mc_truth_h0Wq11_p4->DeltaPhi(*mc_truth_h0Wq12_p4));
      PtotNeut += (*mc_truth_h0Wnu1_p4);
    }
    else if (mc_truth_h0_id>-666 && mc_truth_h0Wq21_id>-666 && mc_truth_h0Wl1_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_h0Wl1_id, *mc_truth_h0Wl1_p4);
      mc_boson_decay = 1;
      Ptot = (*mc_truth_h0_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pjj = (*mc_truth_h0Wq21_p4) + (*mc_truth_h0Wq22_p4);
      mc_boson_pt = mc_truth_h0_p4->Pt();
      mc_boson_l1_pt = mc_truth_h0Wl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_h0Wl1_p4->Eta();
      mc_boson_l2_pt = -666;
      mc_boson_l2_eta = -666;
      mc_boson_ll_mass = -666;
      mc_boson_ll_pt = -666;
      mc_boson_j1_pt = mc_truth_h0Wq21_p4->Pt();
      mc_boson_j1_eta = mc_truth_h0Wq21_p4->Eta();
      mc_boson_j2_pt = mc_truth_h0Wq22_p4->Pt();
      mc_boson_j2_eta = mc_truth_h0Wq22_p4->Eta();
      mc_boson_jj_mass = Pjj.M();
      mc_boson_jj_pt = Pjj.Pt();
      mc_boson_jj_dphi = TMath::Abs(mc_truth_h0Wq21_p4->DeltaPhi(*mc_truth_h0Wq22_p4));
      PtotNeut += (*mc_truth_h0Wnu1_p4);
    }
    else if (mc_truth_W_id>-666 && mc_truth_Wl_id>-666){
      //cout << "TTW"<<endl;
      (*multiLepton).FillParticle("lepton", mc_truth_Wl_id, *mc_truth_Wl_p4);
      mc_boson_decay = 3;
      Ptot = (*mc_truth_W_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      mc_boson_pt = mc_truth_W_p4->Pt();
      mc_boson_l1_pt = mc_truth_Wl_p4->Pt();
      mc_boson_l1_eta = mc_truth_Wl_p4->Eta();
    }
    else if (mc_truth_Z_id>-666 && mc_truth_Zl1_id>-666 && mc_truth_Zl2_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_Zl1_id, *mc_truth_Zl1_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_Zl2_id, *mc_truth_Zl2_p4);
      mc_boson_decay = 2;
      Ptot = (*mc_truth_Zl1_p4) + (*mc_truth_Zl2_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pll = (*mc_truth_Zl1_p4) + (*mc_truth_Zl2_p4);
      mc_boson_pt = Pll.Pt();
      mc_boson_l1_pt = mc_truth_Zl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_Zl1_p4->Eta();
      mc_boson_l2_pt = mc_truth_Zl2_p4->Pt();
      mc_boson_l2_eta = mc_truth_Zl2_p4->Eta();
      mc_boson_ll_mass = Pll.M();
      mc_boson_ll_pt = Pll.Pt();
      mc_boson_ll_dphi = TMath::Abs(mc_truth_Zl1_p4->DeltaPhi(*mc_truth_Zl2_p4));
    }
    else if (mc_truth_gammal1_id>-666 && mc_truth_gammal2_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_gammal1_id, *mc_truth_gammal1_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_gammal2_id, *mc_truth_gammal2_p4);
      mc_boson_decay = 2;
      Ptot = (*mc_truth_gammal1_p4) + (*mc_truth_gammal2_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pll = (*mc_truth_gammal1_p4) + (*mc_truth_gammal2_p4);
      mc_boson_pt = Pll.Pt();
      mc_boson_l1_pt = mc_truth_gammal1_p4->Pt();
      mc_boson_l1_eta = mc_truth_gammal1_p4->Eta();
      mc_boson_l2_pt = mc_truth_gammal2_p4->Pt();
      mc_boson_l2_eta = mc_truth_gammal2_p4->Eta();
      mc_boson_ll_mass = Pll.M();
      mc_boson_ll_pt = Pll.Pt();
      mc_boson_ll_dphi = TMath::Abs(mc_truth_gammal1_p4->DeltaPhi(*mc_truth_gammal2_p4));
    }
    else return;

    if (mc_truth_tWq11_id>-666 && mc_truth_tWl2_id>-666){
      //cout << "Semileptonic ttbar"<<endl;
      (*multiLepton).FillParticle("lepton", mc_truth_tWl2_id, *mc_truth_tWl2_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb1_id, *mc_truth_tb1_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb2_id, *mc_truth_tb2_p4);
      mc_ttbar_decay = 1;
      mc_thad_pt = mc_truth_t1_p4->Pt();
      mc_thad_b_pt = mc_truth_tb1_p4->Pt();
      mc_thad_b_eta = mc_truth_tb1_p4->Eta();
      mc_thad_j1_pt = mc_truth_tWq11_p4->Pt();
      mc_thad_j1_eta = mc_truth_tWq11_p4->Eta();
      mc_thad_j2_pt = mc_truth_tWq12_p4->Pt();
      mc_thad_j2_eta = mc_truth_tWq12_p4->Eta();
      mc_tlep_pt = mc_truth_t2_p4->Pt();
      mc_tlep_b_pt = mc_truth_tb2_p4->Pt();
      mc_tlep_b_eta = mc_truth_tb2_p4->Eta();
      mc_tlep_l_pt = mc_truth_tWl2_p4->Pt();
      mc_tlep_l_eta = mc_truth_tWl2_p4->Eta();
      PtotNeut += (*mc_truth_tWnu2_p4);
   }
   else if (mc_truth_tWq21_id>-666 && mc_truth_tWl1_id>-666){
      //cout << "Semileptonic ttbar"<<endl;
      (*multiLepton).FillParticle("lepton", mc_truth_tWl1_id, *mc_truth_tWl1_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb1_id, *mc_truth_tb1_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb2_id, *mc_truth_tb2_p4);
      mc_ttbar_decay = 1;
      mc_thad_pt = mc_truth_t2_p4->Pt();
      mc_thad_b_pt = mc_truth_tb2_p4->Pt();
      mc_thad_b_eta = mc_truth_tb2_p4->Eta();
      mc_thad_j1_pt = mc_truth_tWq21_p4->Pt();
      mc_thad_j1_eta = mc_truth_tWq21_p4->Eta();
      mc_thad_j2_pt = mc_truth_tWq22_p4->Pt();
      mc_thad_j2_eta = mc_truth_tWq22_p4->Eta();
      mc_tlep_pt = mc_truth_t1_p4->Pt();
      mc_tlep_b_pt = mc_truth_tb1_p4->Pt();
      mc_tlep_b_eta = mc_truth_tb1_p4->Eta();
      mc_tlep_l_pt = mc_truth_tWl1_p4->Pt();
      mc_tlep_l_eta = mc_truth_tWl1_p4->Eta();
      PtotNeut += (*mc_truth_tWnu1_p4);
   }
   else if (mc_truth_tWl1_id>-666 && mc_truth_tWl2_id>-666 && (mc_boson_decay==1 || mc_boson_decay==3)){
      //cout << "Fully leptonic ttbar"<<endl;
      (*multiLepton).FillParticle("bjet", mc_truth_tb1_id, *mc_truth_tb1_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb2_id, *mc_truth_tb2_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_tWl1_id, *mc_truth_tWl1_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_tWl2_id, *mc_truth_tWl2_p4);
      mc_ttbar_decay = 2;
      mc_thad_pt = -666;
      mc_thad_b_pt = -666;
      mc_thad_b_eta = -666;
      mc_thad_j1_pt = -666;
      mc_thad_j1_eta = -666;
      mc_thad_j2_pt = -666;
      mc_thad_j2_eta = -666;
      mc_tlep_pt = mc_truth_t1_p4->Pt();
      mc_tlep_b_pt = mc_truth_tb1_p4->Pt();
      mc_tlep_b_eta = mc_truth_tb1_p4->Eta();
      mc_tlep_l_pt = mc_truth_tWl1_p4->Pt();
      mc_tlep_l_eta = mc_truth_tWl1_p4->Eta();
      mc_tlep2_pt = mc_truth_t2_p4->Pt();
      mc_tlep2_b_eta = mc_truth_tb2_p4->Eta();
      mc_tlep2_l_pt = mc_truth_tWl2_p4->Pt();
      mc_tlep2_l_eta = mc_truth_tWl2_p4->Eta();
      PtotNeut =PtotNeut + (*mc_truth_tWnu1_p4)+ (*mc_truth_tWnu2_p4);
   }
   else return;


   TLorentzVector Pjet;
   for (int i=0; i<genJet_n; i++){
     //cout << "genJet "<<i<<" pt="<<genJet_pt->at(i)<<endl;
     if (TMath::Abs(genJet_eta->at(i))>2.5) continue;
     if (genJet_pt->at(i)<25) continue;
     Pjet.SetPtEtaPhiE(genJet_pt->at(i), genJet_eta->at(i), genJet_phi->at(i), genJet_E->at(i));
     if ((*multiLepton).Bjets[0].P4.DeltaR(Pjet)<0.4 || (*multiLepton).Bjets[1].P4.DeltaR(Pjet)<0.4 || (*multiLepton).Leptons[0].P4.DeltaR(Pjet)<0.4 || (*multiLepton).Leptons[1].P4.DeltaR(Pjet)<0.4 || (*multiLepton).Leptons[2].P4.DeltaR(Pjet)<0.4) continue;
     (*multiLepton).FillParticle("alljet", 0, Pjet);
   }

   if ((*multiLepton).AllJets.size()==0) mc_3l_category = 0;
   else if ((*multiLepton).AllJets.size()==1) {
     mc_3l_category = 1;
     //(*multiLepton).FillParticle("jet", 0, Pjet);
   }
   else if ((*multiLepton).AllJets.size()>=2) {
      mc_3l_category = 2;

      TLorentzVector Pjet2;
      float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
      float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
      for (unsigned int ij=0; ij<(*multiLepton).AllJets.size(); ij++){
        if ((*multiLepton).AllJets[ij].P4.Pt() > pt_max ) {
           pt_max2 = pt_max;
           ij2 = ij1;
           pt_max = (*multiLepton).AllJets[ij].P4.Pt();
           ij1 = ij;
         }
         if ((*multiLepton).AllJets[ij].P4.Pt() < pt_max && (*multiLepton).AllJets[ij].P4.Pt() > pt_max2){
           pt_max2 = (*multiLepton).AllJets[ij].P4.Pt();
           ij2 = ij;
         }
         for (unsigned int ik=0; ik<(*multiLepton).AllJets.size(); ik++){
           if (ik==ij) continue;
           if (TMath::Abs(((*multiLepton).AllJets[ij].P4+(*multiLepton).AllJets[ik].P4).M()-80.419)<diffmass_min){
             ik1=ij;
             ik2=ik;
             diffmass_min = TMath::Abs(((*multiLepton).AllJets[ij].P4+(*multiLepton).AllJets[ik].P4).M()-80.419);
           }
	   if (((*multiLepton).AllJets[ij].P4+(*multiLepton).AllJets[ik].P4).M()<mass_min){
	     il1=ij;
             il2=ik;
	     mass_min = ((*multiLepton).AllJets[ij].P4+(*multiLepton).AllJets[ik].P4).M();
	   }
         }
      }
      if (ij1!=-1 && ij2!=-1) {
        (*multiLepton).FillParticle("jetHighestPt", 1, (*multiLepton).AllJets[ij1].P4);
        (*multiLepton).FillParticle("jetHighestPt", 1, (*multiLepton).AllJets[ij2].P4);
      }
      if (ik1!=-1 && ik2!=-1){
        (*multiLepton).FillParticle("jetClosestMw", 2, (*multiLepton).AllJets[ik1].P4);
        (*multiLepton).FillParticle("jetClosestMw", 2, (*multiLepton).AllJets[ik2].P4);
      }
      if (il1!=-1 && il2!=-1){
        (*multiLepton).FillParticle("jetLowestMjj", 3, (*multiLepton).AllJets[il1].P4);
        (*multiLepton).FillParticle("jetLowestMjj", 3, (*multiLepton).AllJets[il2].P4);
      }
   }

   mc_totp4_px = Ptot.Px();
   mc_totp4_py = Ptot.Py();
   mc_totp4_pt = Ptot.Pt();
   mc_met = PtotNeut.Pt();

   mc_njets25 = (*multiLepton).AllJets.size();

   (*multiLepton).Ptot = Ptot;
   (*multiLepton).mET = PtotNeut;

}

int ReadGenFlatTree::ApplyGenSelection(Long64_t iEvent, MultiLepton* multiLepton){

  cout << "iEvent="<<iEvent<<" BosonDecay="<<mc_boson_decay<< " TTBarDecay="<< mc_ttbar_decay<<" mc_weigth="<<mc_weight<<" PtTot="<< (*multiLepton).Ptot.Pt()<<endl;
  //cout << "MultiLepton nLepton="<<(*multiLepton).Leptons.size()<<" nJets="<<(*multiLepton).Jets.size()<<" nBjets="<<(*multiLepton).Bjets.size()<<endl;


  if (!(((mc_boson_decay==0 || mc_boson_decay==2) && mc_ttbar_decay==1) //ttHfl, ttZ
      || (mc_boson_decay==1 && mc_ttbar_decay==2) //ttHsl
      || (mc_boson_decay==3 && mc_ttbar_decay==2) //ttW
      )) return 0;

  if (!((*multiLepton).Leptons.size()==3 && (*multiLepton).Bjets.size()==2 )) return 0;

  if (mc_3l_category!=2) return 0;

   mc_ttZhypAllowed = 0;
   mc_hasLLcombZpeak = -1;
   mc_passMllGt12 = -1;

   mc_passMllGt12 = 1;
   if (((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(1).P4).M()<12) mc_passMllGt12 = 0;
   if (((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(2).P4).M()<12) mc_passMllGt12 = 0;
   if (((*multiLepton).Leptons.at(1).P4+(*multiLepton).Leptons.at(2).P4).M()<12) mc_passMllGt12 = 0;

   if ((*multiLepton).Leptons.at(0).Id==-(*multiLepton).Leptons.at(1).Id || (*multiLepton).Leptons.at(0).Id==-(*multiLepton).Leptons.at(2).Id || (*multiLepton).Leptons.at(1).Id==-(*multiLepton).Leptons.at(2).Id){

     mc_ttZhypAllowed = 1;
     mc_hasLLcombZpeak = 0;
     mc_passMllGt12 = 1;
     if ((*multiLepton).Leptons.at(0).Id==-(*multiLepton).Leptons.at(1).Id){
        if (TMath::Abs(((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(1).P4).M()-91.188)<10) mc_hasLLcombZpeak = 1;
        //if (((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(1).P4).M()<12) mc_passMllGt12 = 0;
     }
     if ((*multiLepton).Leptons.at(0).Id==-(*multiLepton).Leptons.at(2).Id){
        if (TMath::Abs(((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(2).P4).M()-91.188)<10) mc_hasLLcombZpeak = 1;
        //if (((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(2).P4).M()<12) mc_passMllGt12 = 0;
     }
     if ((*multiLepton).Leptons.at(1).Id==-(*multiLepton).Leptons.at(2).Id){
        if (TMath::Abs(((*multiLepton).Leptons.at(1).P4+(*multiLepton).Leptons.at(2).P4).M()-91.188)<10) mc_hasLLcombZpeak = 1;
        //if (((*multiLepton).Leptons.at(1).P4+(*multiLepton).Leptons.at(2).P4).M()<12) mc_passMllGt12 = 0;
     }
   }

   //cout << "Sort particles"<<endl;
   (*multiLepton).DoSort(&(*multiLepton).Leptons);
   (*multiLepton).DoSort(&(*multiLepton).Jets);
   (*multiLepton).DoSort(&(*multiLepton).Bjets);

   cout << "After sorting, Lepton0Pt="<<(*multiLepton).Leptons.at(0).P4.Pt()<<" Lepton1Pt="<<(*multiLepton).Leptons.at(1).P4.Pt() << " Lepton2Pt="<<(*multiLepton).Leptons.at(2).P4.Pt()<<endl;

   mc_passLepPresel = true;
   if (TMath::Abs((*multiLepton).Leptons.at(0).P4.Eta())>2.5 || TMath::Abs((*multiLepton).Leptons.at(1).P4.Eta())>2.5 || TMath::Abs((*multiLepton).Leptons.at(2).P4.Eta())>2.5) mc_passLepPresel=false;
   if ((*multiLepton).Leptons.at(0).P4.Pt()<10 || (*multiLepton).Leptons.at(1).P4.Pt()<10 || (*multiLepton).Leptons.at(2).P4.Pt()<10) mc_passLepPresel=false;
   if ((*multiLepton).Leptons.at(0).P4.Pt()<20) mc_passLepPresel=false;

   mc_passJetPresel25 = true;
   if (mc_3l_category!=2) mc_passJetPresel25 = false;

   mc_passBjetPresel25 = true;
   //if (TMath::Abs((*multiLepton).Jets.at(0).P4.Eta())>2.5 || TMath::Abs((*multiLepton).Jets.at(1).P4.Eta())>2.5) mc_passJetPresel25 = false;
   //if ((*multiLepton).Jets.at(0).P4.Pt()<25 || (*multiLepton).Jets.at(1).P4.Pt()<25) mc_passJetPresel25 = false;
   if (TMath::Abs((*multiLepton).Bjets.at(0).P4.Eta())>2.5 || TMath::Abs((*multiLepton).Bjets.at(1).P4.Eta())>2.5) mc_passBjetPresel25 = false;
   if ((*multiLepton).Bjets.at(0).P4.Pt()<25 || (*multiLepton).Bjets.at(1).P4.Pt()<25) mc_passBjetPresel25 = false;

   if (mc_passLepPresel && mc_passBjetPresel25 && mc_passMllGt12 && (mc_ttZhypAllowed==0 || (mc_ttZhypAllowed==1 && mc_hasLLcombZpeak!=1)))
     return 1;
   else return 0;

}

void ReadGenFlatTree::WriteMultilepton(MultiLepton* multiLepton){

  cout << "WriteMultilepton"<<endl;

  multilepton_Bjet1_Id 			= (*multiLepton).Bjets[0].Id;
  multilepton_Bjet1_P4 			= (*multiLepton).Bjets[0].P4;
  multilepton_Bjet1_CSV           	= (*multiLepton).Bjets[0].CSV;
  multilepton_Bjet1_DeltaR_Matched      = (*multiLepton).BjetsMatched[0].DeltaR;
  multilepton_Bjet1_Label_Matched       = (*multiLepton).BjetsMatched[0].Label;
  multilepton_Bjet1_Id_Matched 		= (*multiLepton).BjetsMatched[0].Id;
  multilepton_Bjet1_P4_Matched          = (*multiLepton).BjetsMatched[0].P4;

  multilepton_Bjet2_Id 			= (*multiLepton).Bjets[1].Id;
  multilepton_Bjet2_P4 			= (*multiLepton).Bjets[1].P4;
  multilepton_Bjet2_CSV           	= (*multiLepton).Bjets[1].CSV;
  multilepton_Bjet2_DeltaR_Matched      = (*multiLepton).BjetsMatched[1].DeltaR;
  multilepton_Bjet2_Label_Matched       = (*multiLepton).BjetsMatched[1].Label;
  multilepton_Bjet2_Id_Matched          = (*multiLepton).BjetsMatched[1].Id;
  multilepton_Bjet2_P4_Matched          = (*multiLepton).BjetsMatched[1].P4;

  multilepton_Lepton1_Id 		= (*multiLepton).Leptons[0].Id;
  multilepton_Lepton1_P4 		= (*multiLepton).Leptons[0].P4;
  multilepton_Lepton1_DeltaR_Matched	= (*multiLepton).LeptonsMatched[0].DeltaR;
  multilepton_Lepton1_Label_Matched    	= (*multiLepton).LeptonsMatched[0].Label;
  multilepton_Lepton1_Id_Matched    	= (*multiLepton).LeptonsMatched[0].Id;
  multilepton_Lepton1_P4_Matched    	= (*multiLepton).LeptonsMatched[0].P4;

  multilepton_Lepton2_Id 		= (*multiLepton).Leptons[1].Id;
  multilepton_Lepton2_P4 		= (*multiLepton).Leptons[1].P4;
  multilepton_Lepton2_DeltaR_Matched    = (*multiLepton).LeptonsMatched[1].DeltaR;
  multilepton_Lepton2_Label_Matched     = (*multiLepton).LeptonsMatched[1].Label;
  multilepton_Lepton2_Id_Matched        = (*multiLepton).LeptonsMatched[1].Id;
  multilepton_Lepton2_P4_Matched        = (*multiLepton).LeptonsMatched[1].P4;

  multilepton_Lepton3_Id 		= (*multiLepton).Leptons[2].Id;
  multilepton_Lepton3_P4 		= (*multiLepton).Leptons[2].P4;
  multilepton_Lepton3_DeltaR_Matched    = (*multiLepton).LeptonsMatched[2].DeltaR;
  multilepton_Lepton3_Label_Matched     = (*multiLepton).LeptonsMatched[2].Label;
  multilepton_Lepton3_Id_Matched        = (*multiLepton).LeptonsMatched[2].Id;
  multilepton_Lepton3_P4_Matched        = (*multiLepton).LeptonsMatched[2].P4;

  multilepton_Lepton4_Id 		= (*multiLepton).Leptons[3].Id;
  multilepton_Lepton4_P4 		= (*multiLepton).Leptons[3].P4;
  multilepton_Lepton4_DeltaR_Matched    = (*multiLepton).LeptonsMatched[3].DeltaR;
  multilepton_Lepton4_Label_Matched     = (*multiLepton).LeptonsMatched[3].Label;
  multilepton_Lepton4_Id_Matched        = (*multiLepton).LeptonsMatched[3].Id;
  multilepton_Lepton4_P4_Matched        = (*multiLepton).LeptonsMatched[3].P4;

  multilepton_h0_Label                  = (*multiLepton).ParticleGen[0].Label;
  multilepton_h0_Id			= (*multiLepton).ParticleGen[0].Id;
  multilepton_h0_P4			= (*multiLepton).ParticleGen[0].P4;
  multilepton_t1_Label                  = (*multiLepton).ParticleGen[1].Label;
  multilepton_t1_Id                     = (*multiLepton).ParticleGen[1].Id;
  multilepton_t1_P4                     = (*multiLepton).ParticleGen[1].P4;
  multilepton_t2_Label                  = (*multiLepton).ParticleGen[2].Label;
  multilepton_t2_Id                     = (*multiLepton).ParticleGen[2].Id;
  multilepton_t2_P4                     = (*multiLepton).ParticleGen[2].P4;

  multilepton_JetHighestPt1_Id 		= (*multiLepton).JetsHighestPt[0].Id;
  multilepton_JetHighestPt1_P4 		= (*multiLepton).JetsHighestPt[0].P4;
  multilepton_JetHighestPt2_Id 		= (*multiLepton).JetsHighestPt[1].Id;
  multilepton_JetHighestPt2_P4 		= (*multiLepton).JetsHighestPt[1].P4;
  multilepton_JetClosestMw1_Id 		= (*multiLepton).JetsClosestMw[0].Id;
  multilepton_JetClosestMw1_P4 		= (*multiLepton).JetsClosestMw[0].P4;
  multilepton_JetClosestMw2_Id 		= (*multiLepton).JetsClosestMw[1].Id;
  multilepton_JetClosestMw2_P4 		= (*multiLepton).JetsClosestMw[1].P4;
  multilepton_JetLowestMjj1_Id 		= (*multiLepton).JetsLowestMjj[0].Id;
  multilepton_JetLowestMjj1_P4 		= (*multiLepton).JetsLowestMjj[0].P4;
  multilepton_JetLowestMjj2_Id 		= (*multiLepton).JetsLowestMjj[1].Id;
  multilepton_JetLowestMjj2_P4 		= (*multiLepton).JetsLowestMjj[1].P4;
  multilepton_mET 			= (*multiLepton).mET;
  multilepton_Ptot 			= (*multiLepton).Ptot;

  //multilepton_JetHighestPt_Mjj = (multilepton_JetHighestPt1_P4+multilepton_JetHighestPt2_P4).M();
  //multilepton_JetClosestMw_Mjj = (multilepton_JetClosestMw1_P4+multilepton_JetClosestMw2_P4).M();
  //multilepton_JetLowestMjj_Mjj = (multilepton_JetLowestMjj1_P4+multilepton_JetLowestMjj2_P4).M();

  return;
}

void ReadGenFlatTree::ReadMultilepton(Long64_t iEvent, MultiLepton* multiLepton){

  cout << "ReadMultilepton"<<endl;

  tInput->LoadTree(iEvent);
  tInput->GetEntry(iEvent);

  (*multiLepton).kCatJets = catJets;

  (*multiLepton).Leptons.clear();
  if (multilepton_Lepton1_Id!=-999) (*multiLepton).FillParticle("lepton", multilepton_Lepton1_Id, *multilepton_Lepton1_P4_ptr);
  if (multilepton_Lepton2_Id!=-999) (*multiLepton).FillParticle("lepton", multilepton_Lepton2_Id, *multilepton_Lepton2_P4_ptr);
  if (multilepton_Lepton3_Id!=-999) (*multiLepton).FillParticle("lepton", multilepton_Lepton3_Id, *multilepton_Lepton3_P4_ptr);
  if (multilepton_Lepton4_Id!=-999) (*multiLepton).FillParticle("lepton", multilepton_Lepton4_Id, *multilepton_Lepton4_P4_ptr);

  (*multiLepton).LeptonsMatched.clear();
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton1_DeltaR_Matched, multilepton_Lepton1_Label_Matched, multilepton_Lepton1_Id_Matched, *multilepton_Lepton1_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton2_DeltaR_Matched, multilepton_Lepton2_Label_Matched, multilepton_Lepton2_Id_Matched, *multilepton_Lepton2_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton3_DeltaR_Matched, multilepton_Lepton3_Label_Matched, multilepton_Lepton3_Id_Matched, *multilepton_Lepton3_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton4_DeltaR_Matched, multilepton_Lepton4_Label_Matched, multilepton_Lepton4_Id_Matched, *multilepton_Lepton4_P4_Matched_ptr);

  (*multiLepton).Bjets.clear();
  //if (multilepton_Bjet1_Id!=-999) (*multiLepton).FillParticle("bjet", multilepton_Bjet1_Id, *multilepton_Bjet1_P4_ptr);
  if (multilepton_Bjet1_Id!=-999) (*multiLepton).FillParticle("bjet", multilepton_Bjet1_Id, multilepton_Bjet1_CSV, multilepton_Bjet1_JEC_Up, multilepton_Bjet1_JEC_Down, multilepton_Bjet1_JER_Up, multilepton_Bjet1_JER_Down, *multilepton_Bjet1_P4_ptr);
  //if (multilepton_Bjet2_Id!=-999) (*multiLepton).FillParticle("bjet", multilepton_Bjet2_Id, *multilepton_Bjet2_P4_ptr);
  if (multilepton_Bjet2_Id!=-999) (*multiLepton).FillParticle("bjet", multilepton_Bjet2_Id, multilepton_Bjet2_CSV, multilepton_Bjet2_JEC_Up, multilepton_Bjet2_JEC_Down, multilepton_Bjet2_JER_Up, multilepton_Bjet2_JER_Down, *multilepton_Bjet2_P4_ptr);

  (*multiLepton).BjetsMatched.clear();
  (*multiLepton).FillParticleMatched("jet", multilepton_Bjet1_DeltaR_Matched, multilepton_Bjet1_Label_Matched, multilepton_Bjet1_Id_Matched, *multilepton_Bjet1_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("jet", multilepton_Bjet2_DeltaR_Matched, multilepton_Bjet2_Label_Matched, multilepton_Bjet2_Id_Matched, *multilepton_Bjet2_P4_Matched_ptr);

  (*multiLepton).ParticleGen.clear();
  (*multiLepton).FillParticleGen("whatever", multilepton_h0_Label, multilepton_h0_Id, *multilepton_h0_P4_ptr);
  (*multiLepton).FillParticleGen("whatever", multilepton_t1_Label, multilepton_t1_Id, *multilepton_t1_P4_ptr);
  (*multiLepton).FillParticleGen("whatever", multilepton_t2_Label, multilepton_t2_Id, *multilepton_t2_P4_ptr);

  (*multiLepton).Jets.clear();
  (*multiLepton).JetsHighestPt.clear();
  (*multiLepton).JetsClosestMw.clear();
  (*multiLepton).JetsLowestMjj.clear();

  //if (multilepton_JetHighestPt1_Id!=-999) (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt1_Id, *multilepton_JetHighestPt1_P4_ptr);
  if (multilepton_JetHighestPt1_Id!=-999) (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt1_Id, multilepton_JetHighestPt1_CSV, multilepton_JetHighestPt1_JEC_Up, multilepton_JetHighestPt1_JEC_Down, multilepton_JetHighestPt1_JER_Up, multilepton_JetHighestPt1_JER_Down, *multilepton_JetHighestPt1_P4_ptr);
  //if (multilepton_JetHighestPt2_Id!=-999) (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt2_Id, *multilepton_JetHighestPt2_P4_ptr);
  if (multilepton_JetHighestPt2_Id!=-999) (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt2_Id, multilepton_JetHighestPt2_CSV, multilepton_JetHighestPt2_JEC_Up, multilepton_JetHighestPt2_JEC_Down, multilepton_JetHighestPt2_JER_Up, multilepton_JetHighestPt2_JER_Down, *multilepton_JetHighestPt2_P4_ptr);
  //if (multilepton_JetClosestMw1_Id!=-999) (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw1_Id, *multilepton_JetClosestMw1_P4_ptr);
  if (multilepton_JetClosestMw1_Id!=-999) (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw1_Id, multilepton_JetClosestMw1_CSV, multilepton_JetClosestMw1_JEC_Up, multilepton_JetClosestMw1_JEC_Down, multilepton_JetClosestMw1_JER_Up, multilepton_JetClosestMw1_JER_Down, *multilepton_JetClosestMw1_P4_ptr);
  //if (multilepton_JetClosestMw2_Id!=-999) (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw2_Id, *multilepton_JetClosestMw2_P4_ptr);
  if (multilepton_JetClosestMw2_Id!=-999) (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw2_Id, multilepton_JetClosestMw2_CSV, multilepton_JetClosestMw2_JEC_Up, multilepton_JetClosestMw2_JEC_Down, multilepton_JetClosestMw2_JER_Up, multilepton_JetClosestMw2_JER_Down, *multilepton_JetClosestMw2_P4_ptr);
  //if (multilepton_JetLowestMjj1_Id!=-999) (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj1_Id, *multilepton_JetLowestMjj1_P4_ptr);
  if (multilepton_JetLowestMjj1_Id!=-999) (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj1_Id, multilepton_JetLowestMjj1_CSV, multilepton_JetLowestMjj1_JEC_Up, multilepton_JetLowestMjj1_JEC_Down, multilepton_JetLowestMjj1_JER_Up, multilepton_JetLowestMjj1_JER_Down, *multilepton_JetLowestMjj1_P4_ptr);
  //if (multilepton_JetLowestMjj2_Id!=-999) (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj2_Id, *multilepton_JetLowestMjj2_P4_ptr);
  if (multilepton_JetLowestMjj2_Id!=-999) (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj2_Id, multilepton_JetLowestMjj2_CSV, multilepton_JetLowestMjj2_JEC_Up, multilepton_JetLowestMjj2_JEC_Down, multilepton_JetLowestMjj2_JER_Up, multilepton_JetLowestMjj2_JER_Down, *multilepton_JetLowestMjj2_P4_ptr);
  //if (multilepton_JetHighestPt1_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt1_2ndPair_Id, *multilepton_JetHighestPt1_2ndPair_P4_ptr);
  if (multilepton_JetHighestPt1_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt1_2ndPair_Id, multilepton_JetHighestPt1_2ndPair_CSV, multilepton_JetHighestPt1_2ndPair_JEC_Up, multilepton_JetHighestPt1_2ndPair_JEC_Down, multilepton_JetHighestPt1_2ndPair_JER_Up, multilepton_JetHighestPt1_2ndPair_JER_Down, *multilepton_JetHighestPt1_2ndPair_P4_ptr);
  //if (multilepton_JetHighestPt2_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt2_2ndPair_Id, *multilepton_JetHighestPt2_2ndPair_P4_ptr);
  if (multilepton_JetHighestPt2_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt2_2ndPair_Id, multilepton_JetHighestPt2_2ndPair_CSV, multilepton_JetHighestPt2_2ndPair_JEC_Up, multilepton_JetHighestPt2_2ndPair_JEC_Down, multilepton_JetHighestPt2_2ndPair_JER_Up, multilepton_JetHighestPt2_2ndPair_JER_Down, *multilepton_JetHighestPt2_2ndPair_P4_ptr);
  //if (multilepton_JetClosestMw1_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw1_2ndPair_Id, *multilepton_JetClosestMw1_2ndPair_P4_ptr);
  if (multilepton_JetClosestMw1_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw1_2ndPair_Id, multilepton_JetClosestMw1_2ndPair_CSV, multilepton_JetClosestMw1_2ndPair_JEC_Up, multilepton_JetClosestMw1_2ndPair_JEC_Down, multilepton_JetClosestMw1_2ndPair_JER_Up, multilepton_JetClosestMw1_2ndPair_JER_Down, *multilepton_JetClosestMw1_2ndPair_P4_ptr);
  //if (multilepton_JetClosestMw2_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw2_2ndPair_Id, *multilepton_JetClosestMw2_2ndPair_P4_ptr);
  if (multilepton_JetClosestMw2_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw2_2ndPair_Id, multilepton_JetClosestMw2_2ndPair_CSV, multilepton_JetClosestMw2_2ndPair_JEC_Up, multilepton_JetClosestMw2_2ndPair_JEC_Down, multilepton_JetClosestMw2_2ndPair_JER_Up, multilepton_JetClosestMw2_2ndPair_JER_Down, *multilepton_JetClosestMw2_2ndPair_P4_ptr);
  //if (multilepton_JetLowestMjj1_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj1_2ndPair_Id, *multilepton_JetLowestMjj1_2ndPair_P4_ptr);
  if (multilepton_JetLowestMjj1_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj1_2ndPair_Id, multilepton_JetLowestMjj1_2ndPair_CSV, multilepton_JetLowestMjj1_2ndPair_JEC_Up, multilepton_JetLowestMjj1_2ndPair_JEC_Down, multilepton_JetLowestMjj1_2ndPair_JER_Up, multilepton_JetLowestMjj1_2ndPair_JER_Down, *multilepton_JetLowestMjj1_2ndPair_P4_ptr);
  //if (multilepton_JetLowestMjj2_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj2_2ndPair_Id, *multilepton_JetLowestMjj2_2ndPair_P4_ptr);
  if (multilepton_JetLowestMjj2_2ndPair_Id!=-999) (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj2_2ndPair_Id, multilepton_JetLowestMjj2_2ndPair_CSV, multilepton_JetLowestMjj2_2ndPair_JEC_Up, multilepton_JetLowestMjj2_2ndPair_JEC_Down, multilepton_JetLowestMjj2_2ndPair_JER_Up, multilepton_JetLowestMjj2_2ndPair_JER_Down, *multilepton_JetLowestMjj2_2ndPair_P4_ptr);

  //(*multiLepton).Ptot = *multilepton_Ptot_ptr;
  (*multiLepton).mET = *multilepton_mET_ptr;
  (*multiLepton).mET_cov00 = multilepton_mETcov00;
  (*multiLepton).mET_cov01 = multilepton_mETcov01;
  (*multiLepton).mET_cov10 = multilepton_mETcov10;
  (*multiLepton).mET_cov11 = multilepton_mETcov11;
  (*multiLepton).mHT = multilepton_mHT;

  //cout << "Lepton0Pt="<<(*multiLepton).Leptons.at(0).P4.Pt()<<" Lepton1Pt="<<(*multiLepton).Leptons.at(1).P4.Pt() << " Lepton2Pt="<<(*multiLepton).Leptons.at(2).P4.Pt()<<endl;
  //cout << "Bjet0Pt="<<(*multiLepton).Bjets.at(0).P4.Pt()<<" Bjet1Pt="<<(*multiLepton).Bjets.at(1).P4.Pt() << endl;
  //cout << "JetHighestPt0Pt="<<(*multiLepton).JetsHighestPt.at(0).P4.Pt() << " JetHighestPt1Pt="<<(*multiLepton).JetsHighestPt.at(1).P4.Pt() << endl;
  //cout << "JetClosestMw0Pt="<<(*multiLepton).JetsClosestMw.at(0).P4.Pt() << " JetClosestMw1Pt="<<(*multiLepton).JetsClosestMw.at(1).P4.Pt() << endl;
  //cout << "JetLowestMjj0Pt="<<(*multiLepton).JetsLowestMjj.at(0).P4.Pt() << " JetLowestMjj1Pt="<<(*multiLepton).JetsLowestMjj.at(1).P4.Pt() << endl;

  cout << "MultiLepton loaded"<<endl;

  return;
}


void ReadGenFlatTree::ReadMultileptonUserDefined(Long64_t iEvent, MultiLepton* multiLepton){

  cout << "ReadMultilepton"<<endl;

  tInput->LoadTree(iEvent);
  tInput->GetEntry(iEvent);

  // EVENT #1 (run 275376 lumi  331 event 578533240)

  if(false)
  {

  // ===========================================================================================================
  (*multiLepton).kCatJets = 4;

  (*multiLepton).Leptons.clear();
  (*multiLepton).FillParticle("lepton", 13, TLorentzVector(-63.17969382,10.20649212,-176.4611832,187.7083017));
  (*multiLepton).FillParticle("lepton", -13, TLorentzVector(10.29204053,24.36307254,-56.49449953,62.37888306));
  (*multiLepton).FillParticle("lepton", -11, TLorentzVector(12.32398555,-0.1728198492,-5.328197284,13.42759014));

  (*multiLepton).LeptonsMatched.clear();
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton1_DeltaR_Matched, multilepton_Lepton1_Label_Matched, multilepton_Lepton1_Id_Matched, *multilepton_Lepton1_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton2_DeltaR_Matched, multilepton_Lepton2_Label_Matched, multilepton_Lepton2_Id_Matched, *multilepton_Lepton2_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton3_DeltaR_Matched, multilepton_Lepton3_Label_Matched, multilepton_Lepton3_Id_Matched, *multilepton_Lepton3_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton4_DeltaR_Matched, multilepton_Lepton4_Label_Matched, multilepton_Lepton4_Id_Matched, *multilepton_Lepton4_P4_Matched_ptr);

  (*multiLepton).Bjets.clear();
  (*multiLepton).FillParticle("bjet", 0, 0.8973406553, 0,0,0,0, TLorentzVector(25.81664968,-63.80021097,-279.2412525,287.7334337));
  (*multiLepton).FillParticle("bjet", 0, 0.5502513051, 0,0,0,0, TLorentzVector(29.0699138,1.73146953,143.3043164,146.3926179));

  (*multiLepton).BjetsMatched.clear();
  (*multiLepton).FillParticleMatched("jet", multilepton_Bjet1_DeltaR_Matched, multilepton_Bjet1_Label_Matched, multilepton_Bjet1_Id_Matched, *multilepton_Bjet1_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("jet", multilepton_Bjet2_DeltaR_Matched, multilepton_Bjet2_Label_Matched, multilepton_Bjet2_Id_Matched, *multilepton_Bjet2_P4_Matched_ptr);

  (*multiLepton).ParticleGen.clear();
  (*multiLepton).FillParticleGen("whatever", multilepton_h0_Label, multilepton_h0_Id, *multilepton_h0_P4_ptr);
  (*multiLepton).FillParticleGen("whatever", multilepton_t1_Label, multilepton_t1_Id, *multilepton_t1_P4_ptr);
  (*multiLepton).FillParticleGen("whatever", multilepton_t2_Label, multilepton_t2_Id, *multilepton_t2_P4_ptr);

  (*multiLepton).Jets.clear();
  (*multiLepton).JetsHighestPt.clear();
  (*multiLepton).JetsClosestMw.clear();
  (*multiLepton).JetsLowestMjj.clear();

  (*multiLepton).mET = TLorentzVector(16.17305402,16.68213381,0,23.23491478);
  (*multiLepton).mET_cov00 = 1;
  (*multiLepton).mET_cov01 = 0;
  (*multiLepton).mET_cov10 = 0;
  (*multiLepton).mET_cov11 = 1;
  (*multiLepton).mHT = 1572.618652;

  }

  // EVENT #2 (run 276776 lumi  881 event 1524218683)

  if(true)
  {

  // ===========================================================================================================
  (*multiLepton).kCatJets = 3;

  (*multiLepton).Leptons.clear();
  (*multiLepton).FillParticle("lepton", 11, TLorentzVector(32.88572615,24.57231569,-46.06872267,61.70573037));
  (*multiLepton).FillParticle("lepton", -11, TLorentzVector(14.44073434,-31.24464316,-121.2751663,126.0651749));
  (*multiLepton).FillParticle("lepton", -11, TLorentzVector(-10.16623148,-1.23847927,-25.05356718,27.0659809));

  (*multiLepton).LeptonsMatched.clear();
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton1_DeltaR_Matched, multilepton_Lepton1_Label_Matched, multilepton_Lepton1_Id_Matched, *multilepton_Lepton1_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton2_DeltaR_Matched, multilepton_Lepton2_Label_Matched, multilepton_Lepton2_Id_Matched, *multilepton_Lepton2_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton3_DeltaR_Matched, multilepton_Lepton3_Label_Matched, multilepton_Lepton3_Id_Matched, *multilepton_Lepton3_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("lepton", multilepton_Lepton4_DeltaR_Matched, multilepton_Lepton4_Label_Matched, multilepton_Lepton4_Id_Matched, *multilepton_Lepton4_P4_Matched_ptr);

  (*multiLepton).Bjets.clear();
  (*multiLepton).FillParticle("bjet", 0, 0.930102706, 0,0,0,0, TLorentzVector(-143.3783458,8.338624118,41.79349018,150.625403));

  (*multiLepton).BjetsMatched.clear();
  (*multiLepton).FillParticleMatched("jet", multilepton_Bjet1_DeltaR_Matched, multilepton_Bjet1_Label_Matched, multilepton_Bjet1_Id_Matched, *multilepton_Bjet1_P4_Matched_ptr);
  (*multiLepton).FillParticleMatched("jet", multilepton_Bjet2_DeltaR_Matched, multilepton_Bjet2_Label_Matched, multilepton_Bjet2_Id_Matched, *multilepton_Bjet2_P4_Matched_ptr);

  (*multiLepton).ParticleGen.clear();
  (*multiLepton).FillParticleGen("whatever", multilepton_h0_Label, multilepton_h0_Id, *multilepton_h0_P4_ptr);
  (*multiLepton).FillParticleGen("whatever", multilepton_t1_Label, multilepton_t1_Id, *multilepton_t1_P4_ptr);
  (*multiLepton).FillParticleGen("whatever", multilepton_t2_Label, multilepton_t2_Id, *multilepton_t2_P4_ptr);

  (*multiLepton).Jets.clear();
  (*multiLepton).JetsHighestPt.clear();
  (*multiLepton).JetsClosestMw.clear();
  (*multiLepton).JetsLowestMjj.clear();

  (*multiLepton).FillParticle("jetHighestPt", 0, 0.1687322706, 0,0,0,0, TLorentzVector(162.643124,12.89360473,-309.6293432,350.2472345));

  (*multiLepton).mET = TLorentzVector(-17.86693354,27.48508335,0,32.78196335);
  (*multiLepton).mET_cov00 = 1;
  (*multiLepton).mET_cov01 = 0;
  (*multiLepton).mET_cov10 = 0;
  (*multiLepton).mET_cov11 = 1;
  (*multiLepton).mHT = 1302.917236;

  }

  //cout << "Lepton0Pt="<<(*multiLepton).Leptons.at(0).P4.Pt()<<" Lepton1Pt="<<(*multiLepton).Leptons.at(1).P4.Pt() << " Lepton2Pt="<<(*multiLepton).Leptons.at(2).P4.Pt()<<endl;
  //cout << "Bjet0Pt="<<(*multiLepton).Bjets.at(0).P4.Pt()<<" Bjet1Pt="<<(*multiLepton).Bjets.at(1).P4.Pt() << endl;
  //cout << "JetHighestPt0Pt="<<(*multiLepton).JetsHighestPt.at(0).P4.Pt() << " JetHighestPt1Pt="<<(*multiLepton).JetsHighestPt.at(1).P4.Pt() << endl;
  //cout << "JetClosestMw0Pt="<<(*multiLepton).JetsClosestMw.at(0).P4.Pt() << " JetClosestMw1Pt="<<(*multiLepton).JetsClosestMw.at(1).P4.Pt() << endl;
  //cout << "JetLowestMjj0Pt="<<(*multiLepton).JetsLowestMjj.at(0).P4.Pt() << " JetLowestMjj1Pt="<<(*multiLepton).JetsLowestMjj.at(1).P4.Pt() << endl;

  cout << "MultiLepton loaded"<<endl;

  return;
}


#endif
