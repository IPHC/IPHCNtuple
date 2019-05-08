{
   gROOT->SetBatch();

   float v_muon_1;
   float v_elec_1;
   float v_tau_1;
   float v_jet_1;

   float v_muon_2;
   float v_elec_2;
   float v_tau_2;
   float v_jet_2;
   
   int n_presel_mu_1;
   int n_presel_ele_1;
   int n_presel_tau_1;
   int n_presel_jet_1;

   float n_presel_mu_2;
   float n_presel_ele_2;
   float n_presel_tau_2;
   float n_presel_jet_2;
   
   int nEvent_1;
   Long64_t nEvent_2;

   std::string treeName = "syncTree";

   std::string varName_muon = "mu1_pt";
   std::string varName_elec = "ele1_pt";
   std::string varName_tau = "tau1_pt";
   std::string varName_jet = "jet1_pt";

   TFile *f1 = TFile::Open("tHq_syncTree_Object.root");
   TTree *tr1 = (TTree*)f1->Get(treeName.c_str());
   tr1->SetBranchAddress(varName_muon.c_str(),&v_muon_1);
   tr1->SetBranchAddress(varName_elec.c_str(),&v_elec_1);
   tr1->SetBranchAddress(varName_tau.c_str(),&v_tau_1);
   tr1->SetBranchAddress(varName_jet.c_str(),&v_jet_1);
   tr1->SetBranchAddress("n_presel_mu",&n_presel_mu_1);
   tr1->SetBranchAddress("n_presel_ele",&n_presel_ele_1);
   tr1->SetBranchAddress("n_presel_tau",&n_presel_tau_1);
   tr1->SetBranchAddress("n_presel_jet",&n_presel_jet_1);
   tr1->SetBranchAddress("nEvent",&nEvent_1);

   TFile *f2 = TFile::Open("IHEP_ttHsync_2016_V5.root");
   TTree *tr2 = (TTree*)f2->Get(treeName.c_str());
   tr2->SetBranchAddress(varName_muon.c_str(),&v_muon_2);
   tr2->SetBranchAddress(varName_elec.c_str(),&v_elec_2);
   tr2->SetBranchAddress(varName_tau.c_str(),&v_tau_2);
   tr2->SetBranchAddress(varName_jet.c_str(),&v_jet_2);
   tr2->SetBranchAddress("n_presel_mu",&n_presel_mu_2);
   tr2->SetBranchAddress("n_presel_ele",&n_presel_ele_2);
   tr2->SetBranchAddress("n_presel_tau",&n_presel_tau_2);
   tr2->SetBranchAddress("n_presel_jet",&n_presel_jet_2);
   tr2->SetBranchAddress("nEvent",&nEvent_2);

   int nEnt1 = tr1->GetEntries();
   int nEnt2 = tr2->GetEntries();

   std::ofstream f1_muon("muon_1.txt");
   std::ofstream f1_elec("elec_1.txt");
   std::ofstream f1_tau("tau_1.txt");
   std::ofstream f1_jet("jet_1.txt");
   
   std::ofstream f2_muon("muon_2.txt");
   std::ofstream f2_elec("elec_2.txt");
   std::ofstream f2_tau("tau_2.txt");
   std::ofstream f2_jet("jet_2.txt");
   
   for(int i1=0;i1<nEnt1;i1++)
     {
	tr1->GetEntry(i1);
	
	if( n_presel_mu_1 > 0 ) f1_muon << nEvent_1 << " " << v_muon_1 << "\n";
	if( n_presel_ele_1 > 0 ) f1_elec << nEvent_1 << " " << v_elec_1 << "\n";
	if( n_presel_tau_1 > 0 ) f1_tau << nEvent_1 << " " << v_tau_1 << "\n";
	if( n_presel_jet_1 > 0 ) f1_jet << nEvent_1 << " " << v_jet_1 << "\n";
     }
   
   f1_muon.close();
   f1_elec.close();
   f1_tau.close();
   f1_jet.close();

   for(int i2=0;i2<nEnt2;i2++)
     {
	tr2->GetEntry(i2);
	
	if( n_presel_mu_2 > 0 ) f2_muon << nEvent_2 << " " << v_muon_2 << "\n";
	if( n_presel_ele_2 > 0 ) f2_elec << nEvent_2 << " " << v_elec_2 << "\n";
	if( n_presel_tau_2 > 0 ) f2_tau << nEvent_2 << " " << v_tau_2 << "\n";
	if( n_presel_jet_2 > 0 ) f2_jet << nEvent_2 << " " << v_jet_2 << "\n";
     }
   
   f2_muon.close();
   f2_elec.close();
   f2_tau.close();
   f2_jet.close();
   
   gApplication->Terminate();
}
