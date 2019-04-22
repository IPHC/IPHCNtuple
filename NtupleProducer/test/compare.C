{
   gROOT->SetBatch();

//   float v1;
//   int v2;
   float v1;
   float v2;
//   int v1;
//   bool v2;
   int nEvent1;
   Long64_t nEvent2;

   std::string treeName = "syncTree";
//   std::string treeName = "syncTree_1l2tau_Fake";
//   std::string treeName = "syncTree_ttZctrl_Fake";

   std::string varName = "mu1_pt";
//   std::string varName = "ele1_isfakeablesel";
//   std::string varName = "tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT";
//   std::string varName = "isGenMatched";
//   std::string varName = "";
//   std::string varName = "MHT";
//   std::string varName = "jet3_pt";

   TFile *f1 = TFile::Open("tHq_syncTree_Object.root");
   TTree *tr1 = (TTree*)f1->Get(treeName.c_str());
   tr1->SetBranchAddress(varName.c_str(),&v1);
   tr1->SetBranchAddress("nEvent",&nEvent1);

   TFile *f2 = TFile::Open("IHEP_ttHsync_2016_V4.root");
   TTree *tr2 = (TTree*)f2->Get(treeName.c_str());
   tr2->SetBranchAddress(varName.c_str(),&v2);
   tr2->SetBranchAddress("nEvent",&nEvent2);

   int nEnt1 = tr1->GetEntries();
   int nEnt2 = tr2->GetEntries();

   // others
/*   for(int i1=0;i1<nEnt1;i1++)
     {
//	if( i2 > 10 ) break;

	tr1->GetEntry(i1);
	
//	if( nEvent2 != 15302912 ) continue;

	bool found = 0;
	
	for(int i2=0;i2<nEnt2;i2++)
	  {
	     tr2->GetEntry(i2);

	     if( nEvent1 == nEvent2 )
	       {		  
		  std::cout << nEvent1 << ": " << v1 << " " << v2 << std::endl;
		  found = 1;
	       }	     
	  }
	
	if( !found )
	  {
	     std::cout << "Others reject " << nEvent1 << ": " << v1 << std::endl;
	  }	
     }*/

   // me
   for(int i2=0;i2<nEnt2;i2++)
     {
	tr2->GetEntry(i2);

//	if( i2 > 10 ) break;

//	if( nEvent2 != 14140553 ) continue;
	
	bool found = 0;
	
	for(int i1=0;i1<nEnt1;i1++)
	  {
	     tr1->GetEntry(i1);

	     if( nEvent1 == nEvent2 )
	       {
		  std::cout << nEvent1 << ": " << v1 << " " << v2 << std::endl;
		  found = 1;
	       }	     
	  }
	
	if( !found )
	  {
	     std::cout << "I reject " << nEvent2 << ": " << v2 << std::endl;
	  }
     }

   gApplication->Terminate();
}
