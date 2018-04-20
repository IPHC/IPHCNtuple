{

  //TFile* f = TFile::Open("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/output_ttZDelphes_testMinimizer_Simplex_new_all.root");
  TFile* f = TFile::Open("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v3/output_ttZDelphes_testDelphesEvalGen_4000k.root");
  TTree* t = (TTree*) f->Get("Tree");
  t->SetBranchStatus("*",0);
  t->SetBranchStatus("multilepton_Bjet1_P4",1);
  t->SetBranchStatus("multilepton_Bjet2_P4",1);
  t->SetBranchStatus("multilepton_JetClosestMw1_P4",1);
  t->SetBranchStatus("multilepton_JetClosestMw2_P4",1);
  t->SetBranchStatus("multilepton_Lepton1_P4",1);
  t->SetBranchStatus("multilepton_Lepton2_P4",1);
  t->SetBranchStatus("multilepton_Lepton3_P4",1);
  t->SetBranchStatus("weight",1);
  t->SetBranchStatus("is_3l_TTZ_CR",1);
  t->SetBranchStatus("mc_ttZhypAllowed",1);
  t->SetBranchStatus("catJets",1);
  t->SetBranchStatus("multilepton_mET",1);
  t->SetBranchStatus("multilepton_mHT",1);
  t->SetBranchStatus("mc_kin_ttz_inputvars",1);
  t->SetBranchStatus("mc_kinInt_ttz_inputvars",1);
  t->SetBranchStatus("mc_kin_ttz_weight_logmax",1);
  t->SetBranchStatus("mc_kin_ttz_weight_logmaxint",1);
  t->SetBranchStatus("mc_mem_ttz_weight_evalgenmax",1);
  t->SetBranchStatus("multilepton_Bjet1_P4_Matched",1);
  t->SetBranchStatus("multilepton_Bjet2_P4_Matched",1);
  t->SetBranchStatus("multilepton_mET_Matched",1);
  t->SetBranchStatus("multilepton_W1_P4_Matched",1);
  t->SetBranchStatus("multilepton_W2_P4_Matched",1);

  TFile* fNew = new TFile("input_TTZ_DelphesEvalGen_4000k.root", "RECREATE");
  TTree* tNew = new TTree("Tree","Tree");

  TBranch* b_multilepton_Bjet1_P4;
  TLorentzVector* multilepton_Bjet1_P4;
  t->SetBranchAddress("multilepton_Bjet1_P4",&multilepton_Bjet1_P4,&b_multilepton_Bjet1_P4);
  tNew->Branch("multilepton_Bjet1_P4","TLorentzVector",&multilepton_Bjet1_P4);

  TBranch* b_multilepton_Bjet2_P4;
  TLorentzVector* multilepton_Bjet2_P4;
  t->SetBranchAddress("multilepton_Bjet2_P4",&multilepton_Bjet2_P4,&b_multilepton_Bjet2_P4);
  tNew->Branch("multilepton_Bjet2_P4","TLorentzVector",&multilepton_Bjet2_P4);

  TBranch* b_multilepton_JetClosestMw1_P4;
  TLorentzVector* multilepton_JetClosestMw1_P4;
  t->SetBranchAddress("multilepton_JetClosestMw1_P4",&multilepton_JetClosestMw1_P4,&b_multilepton_JetClosestMw1_P4);
  tNew->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);

  TBranch* b_multilepton_JetClosestMw2_P4;
  TLorentzVector* multilepton_JetClosestMw2_P4;
  t->SetBranchAddress("multilepton_JetClosestMw2_P4",&multilepton_JetClosestMw2_P4,&b_multilepton_JetClosestMw2_P4);
  tNew->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);

  TBranch* b_multilepton_Lepton1_P4;
  TLorentzVector* multilepton_Lepton1_P4;
  t->SetBranchAddress("multilepton_Lepton1_P4",&multilepton_Lepton1_P4,&b_multilepton_Lepton1_P4);
  tNew->Branch("multilepton_Lepton1_P4","TLorentzVector",&multilepton_Lepton1_P4);

  TBranch* b_multilepton_Lepton2_P4;
  TLorentzVector* multilepton_Lepton2_P4;
  t->SetBranchAddress("multilepton_Lepton2_P4",&multilepton_Lepton2_P4,&b_multilepton_Lepton2_P4);
  tNew->Branch("multilepton_Lepton2_P4","TLorentzVector",&multilepton_Lepton2_P4);

  TBranch* b_multilepton_Lepton3_P4;
  TLorentzVector* multilepton_Lepton3_P4;
  t->SetBranchAddress("multilepton_Lepton3_P4",&multilepton_Lepton3_P4,&b_multilepton_Lepton3_P4);
  tNew->Branch("multilepton_Lepton3_P4","TLorentzVector",&multilepton_Lepton3_P4);

  TBranch* b_multilepton_mET;
  TLorentzVector* multilepton_mET;
  t->SetBranchAddress("multilepton_mET",&multilepton_mET,&b_multilepton_mET);
  tNew->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);

  TBranch* b_multilepton_mHT;
  Float_t multilepton_mHT;
  t->SetBranchAddress("multilepton_mHT",&multilepton_mHT,&b_multilepton_mHT);
  tNew->Branch("multilepton_mHT",&multilepton_mHT,"multilepton_mHT/F");

  TBranch* b_is_3l_TTZ_CR;
  Char_t is_3l_TTZ_CR;
  t->SetBranchAddress("is_3l_TTZ_CR",&is_3l_TTZ_CR,&b_is_3l_TTZ_CR);
  tNew->Branch("is_3l_TTZ_CR",&is_3l_TTZ_CR,"is_3l_TTZ_CR/B");

  TBranch* b_mc_ttZhypAllowed;
  Int_t mc_ttZhypAllowed;
  t->SetBranchAddress("mc_ttZhypAllowed",&mc_ttZhypAllowed,&b_mc_ttZhypAllowed);
  tNew->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");

  TBranch* b_catJets;
  Int_t catJets;
  t->SetBranchAddress("catJets",&catJets,&b_catJets);
  tNew->Branch("catJets",&catJets,"catJets/I");

  TBranch* b_mc_kin_ttz_weight_logmax;
  Double_t mc_kin_ttz_weight_logmax;
  t->SetBranchAddress("mc_kin_ttz_weight_logmax",&mc_kin_ttz_weight_logmax, &b_mc_kin_ttz_weight_logmax);
  tNew->Branch("mc_kin_ttz_weight_logmax",&mc_kin_ttz_weight_logmax,"mc_kin_ttz_weight_logmax/D");

  TBranch* b_mc_kin_ttz_weight_logmaxint;
  Double_t mc_kin_ttz_weight_logmaxint;
  t->SetBranchAddress("mc_kin_ttz_weight_logmaxint",&mc_kin_ttz_weight_logmaxint, &b_mc_kin_ttz_weight_logmaxint);
  tNew->Branch("mc_kin_ttz_weight_logmaxint",&mc_kin_ttz_weight_logmaxint,"mc_kin_ttz_weight_logmaxint/D");


  TBranch* b_kin_ttz_inputvars;
  std::vector<double>* kin_ttz_inputvars;
  t->SetBranchAddress("mc_kin_ttz_inputvars", &kin_ttz_inputvars, &b_kin_ttz_inputvars);
  tNew->Branch("mc_kin_ttz_inputvars","std::vector<double>", &kin_ttz_inputvars);

  TBranch* b_kinInt_ttz_inputvars;
  std::vector<double>* kinInt_ttz_inputvars;
  t->SetBranchAddress("mc_kinInt_ttz_inputvars", &kinInt_ttz_inputvars, &b_kinInt_ttz_inputvars);
  tNew->Branch("mc_kinInt_ttz_inputvars","std::vector<double>", &kinInt_ttz_inputvars);

  TBranch* b_mc_mem_ttz_weight_evalgenmax;
  Double_t mc_mem_ttz_weight_evalgenmax, mc_mem_ttz_weight_evalgenmax_log;
  t->SetBranchAddress("mc_mem_ttz_weight_evalgenmax",&mc_mem_ttz_weight_evalgenmax,&b_mc_mem_ttz_weight_evalgenmax);
  tNew->Branch("mc_mem_ttz_weight_evalgenmax_log",&mc_mem_ttz_weight_evalgenmax_log,"mc_mem_ttz_weight_evalgenmax_log/D");

  TBranch* b_multilepton_Bjet1_P4_Matched;
  TLorentzVector* multilepton_Bjet1_P4_Matched;
  multilepton_Bjet1_P4_Matched = 0;
  t->SetBranchAddress("multilepton_Bjet1_P4_Matched", &multilepton_Bjet1_P4_Matched, &b_multilepton_Bjet1_P4_Matched);

  TBranch* b_multilepton_Bjet2_P4_Matched;
  TLorentzVector* multilepton_Bjet2_P4_Matched = 0;
  t->SetBranchAddress("multilepton_Bjet2_P4_Matched", &multilepton_Bjet2_P4_Matched, &b_multilepton_Bjet2_P4_Matched);

  TBranch* b_multilepton_mET_Matched;
  TLorentzVector* multilepton_mET_Matched = 0;
  t->SetBranchAddress("multilepton_mET_Matched",&multilepton_mET_Matched,&b_multilepton_mET_Matched);

  TBranch* b_multilepton_W1_P4_Matched;
  TLorentzVector* multilepton_W1_P4_Matched = 0;
  t->SetBranchAddress("multilepton_W1_P4_Matched",&multilepton_W1_P4_Matched,&b_multilepton_W1_P4_Matched);

  TBranch* b_multilepton_W2_P4_Matched;
  TLorentzVector* multilepton_W2_P4_Matched = 0;
  t->SetBranchAddress("multilepton_W2_P4_Matched",&multilepton_W2_P4_Matched,&b_multilepton_W2_P4_Matched);

  float Gen_BjetTopHad_E, Gen_BjetTopLep_E;
  float Gen_WTopHad_tW, Gen_WTopHad_mW;
  float Gen_WTopLep_tW, Gen_WTopLep_mW;
  float Gen_NeutTopLep_Phi;
  tNew->Branch("Gen_BjetTopHad_E",&Gen_BjetTopHad_E,"Gen_BjetTopHad_E/F");
  tNew->Branch("Gen_WTopHad_tW",&Gen_WTopHad_tW,"Gen_WTopHad_tW/F");
  tNew->Branch("Gen_WTopHad_mW",&Gen_WTopHad_mW,"Gen_WTopHad_mW/F");
  tNew->Branch("Gen_BjetTopLep_E",&Gen_BjetTopLep_E,"Gen_BjetTopLep_E/F");
  tNew->Branch("Gen_NeutTopLep_Phi",&Gen_NeutTopLep_Phi,"Gen_NeutTopLep_Phi/F");
  tNew->Branch("Gen_WTopLep_tW",&Gen_WTopLep_tW,"Gen_WTopLep_tW/F");
  tNew->Branch("Gen_WTopLep_mW",&Gen_WTopLep_mW,"Gen_WTopLep_mW/F");


  float Kin_BjetTopHad_E, Kin_WTopHad_tW, Kin_WTopHad_mW;
  float Kin_BjetTopLep_E, Kin_NeutTopLep_Phi, Kin_WTopLep_tW, Kin_WTopLep_mW;
  tNew->Branch("Kin_BjetTopHad_E", &Kin_BjetTopHad_E, "Kin_BjetTopHad_E/F");
  tNew->Branch("Kin_WTopHad_tW", &Kin_WTopHad_tW, "Kin_WTopHad_tW/F");
  tNew->Branch("Kin_WTopHad_mW", &Kin_WTopHad_mW, "Kin_WTopHad_mW/F");
  tNew->Branch("Kin_BjetTopLep_E",&Kin_BjetTopLep_E,"Kin_BjetTopLep_E/F");
  tNew->Branch("Kin_NeutTopLep_Phi",&Kin_NeutTopLep_Phi,"Kin_NeutTopLep_Phi/F");
  tNew->Branch("Kin_WTopLep_tW",&Kin_WTopLep_tW,"Kin_WTopLep_tW/F");
  tNew->Branch("Kin_WTopLep_mW",&Kin_WTopLep_mW,"Kin_WTopLep_mW/F");

  TBranch* b_weight;
  float weight;
  t->SetBranchAddress("weight",&weight,&b_weight);
  tNew->Branch("weight",&weight,"weight/F");


  double mW = 80.419;
  double gammaW = 2.0476;
  unsigned int nEntries = t->GetEntries();
  for (Long64_t jentry=0; jentry<nEntries;jentry++){
    t->GetEntry(jentry);
    //if (jentry==280 || jentry==478) cout << "Kin_BjetTopHad_E="<<Kin_BjetTopHad_E<<endl;
    //cout << "size="<<kin_ttz_inputvars->size()<<endl;
    //if (kin_ttz_inputvars->size()!=5) continue;
    if (kin_ttz_inputvars->size()!=5 || kinInt_ttz_inputvars->size()!=5) continue;
    //if (!(kin_ttz_inputvars->at(0)>0)) continue;
    weight = 1;

    Kin_BjetTopHad_E = kin_ttz_inputvars->at(0);
    Kin_WTopHad_tW = kin_ttz_inputvars->at(1);
    Kin_WTopHad_mW = sqrt(mW*mW + mW*gammaW*tan(Kin_WTopHad_tW));
    Kin_BjetTopLep_E = kin_ttz_inputvars->at(2);
    Kin_NeutTopLep_Phi = kin_ttz_inputvars->at(3);
    Kin_WTopLep_tW = kin_ttz_inputvars->at(4);
    Kin_WTopLep_mW = sqrt(mW*mW + mW*gammaW*tan(Kin_WTopLep_tW));

    Gen_BjetTopHad_E = multilepton_Bjet1_P4_Matched->E();
    Gen_WTopHad_mW = multilepton_W1_P4_Matched->M();
    Gen_WTopHad_tW = TMath::ATan((Gen_WTopHad_mW*Gen_WTopHad_mW-mW*mW)/(mW*gammaW));
    Gen_BjetTopLep_E = multilepton_Bjet2_P4_Matched->E();
    Gen_NeutTopLep_Phi = multilepton_mET_Matched->Phi();
    Gen_WTopLep_mW = multilepton_W2_P4_Matched->M();
    Gen_WTopLep_tW = TMath::ATan((Gen_WTopLep_mW*Gen_WTopLep_mW-mW*mW)/(mW*gammaW));

    mc_mem_ttz_weight_evalgenmax_log = log(mc_mem_ttz_weight_evalgenmax);
    //cout << "mc_kin_ttz_weight_logmax="<<mc_kin_ttz_weight_logmax<<endl;
    //if (!isfinite(Kin_BjetTopHad_E)) cout << "Kin_BjetTopHad_E nan"<<endl;
    //if (jentry==280 || jentry==478) cout << "Kin_BjetTopHad_E="<<Kin_BjetTopHad_E<<endl;
    tNew->Fill();
  }

  tNew->Write();
  fNew->Write();
  fNew->Close();

}
