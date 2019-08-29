{
   gROOT->SetBatch();

   bool isEventSync = 1;
   
   TFile *f1 = TFile::Open("tHq_syncTree_Selection.root");
   
   int isTHQ;
   int isTTH;
   
   TTree *syncTree_2lSS_SR = (TTree*)f1->Get("syncTree_2lSS_SR");
   TTree *syncTree_2lSS_Fake = (TTree*)f1->Get("syncTree_2lSS_Fake");
   TTree *syncTree_2lSS_Flip = (TTree*)f1->Get("syncTree_2lSS_Flip");
   TTree *syncTree_3l_SR = (TTree*)f1->Get("syncTree_3l_SR");
   TTree *syncTree_3l_Fake = (TTree*)f1->Get("syncTree_3l_Fake");
   TTree *syncTree_4l_SR = (TTree*)f1->Get("syncTree_4l_SR");
   TTree *syncTree_4l_Fake = (TTree*)f1->Get("syncTree_4l_Fake");
   TTree *syncTree_0l2tau_SR = (TTree*)f1->Get("syncTree_0l2tau_SR");
   TTree *syncTree_0l2tau_Fake = (TTree*)f1->Get("syncTree_0l2tau_Fake");
   TTree *syncTree_1l1tau_SR = (TTree*)f1->Get("syncTree_1l1tau_SR");
   TTree *syncTree_1l1tau_Fake = (TTree*)f1->Get("syncTree_1l1tau_Fake");
   TTree *syncTree_1l1tau_Flip = (TTree*)f1->Get("syncTree_1l1tau_Flip");
   TTree *syncTree_1l2tau_SR = (TTree*)f1->Get("syncTree_1l2tau_SR");
   TTree *syncTree_1l2tau_Fake = (TTree*)f1->Get("syncTree_1l2tau_Fake");
   TTree *syncTree_2lSS1tau_SR = (TTree*)f1->Get("syncTree_2lSS1tau_SR");
   TTree *syncTree_2lSS1tau_Fake = (TTree*)f1->Get("syncTree_2lSS1tau_Fake");
   TTree *syncTree_2lSS1tau_Flip = (TTree*)f1->Get("syncTree_2lSS1tau_Flip");
   TTree *syncTree_2lOS1tau_SR = (TTree*)f1->Get("syncTree_2lOS1tau_SR");
   TTree *syncTree_2lOS1tau_Fake = (TTree*)f1->Get("syncTree_2lOS1tau_Fake");
   TTree *syncTree_2l2tau_SR = (TTree*)f1->Get("syncTree_2l2tau_SR");
   TTree *syncTree_2l2tau_Fake = (TTree*)f1->Get("syncTree_2l2tau_Fake");
   TTree *syncTree_3l1tau_SR = (TTree*)f1->Get("syncTree_3l1tau_SR");
   TTree *syncTree_3l1tau_Fake = (TTree*)f1->Get("syncTree_3l1tau_Fake");
   TTree *syncTree_ttWctrl_SR = (TTree*)f1->Get("syncTree_ttWctrl_SR");
   TTree *syncTree_ttWctrl_Fake = (TTree*)f1->Get("syncTree_ttWctrl_Fake");
   TTree *syncTree_ttWctrl_Flip = (TTree*)f1->Get("syncTree_ttWctrl_Flip");
   TTree *syncTree_ttZctrl_SR = (TTree*)f1->Get("syncTree_ttZctrl_SR");
   TTree *syncTree_ttZctrl_Fake = (TTree*)f1->Get("syncTree_ttZctrl_Fake");
   TTree *syncTree_WZctrl_SR = (TTree*)f1->Get("syncTree_WZctrl_SR");
   TTree *syncTree_WZctrl_Fake = (TTree*)f1->Get("syncTree_WZctrl_Fake");
   TTree *syncTree_ZZctrl_SR = (TTree*)f1->Get("syncTree_ZZctrl_SR");
   TTree *syncTree_ZZctrl_Fake = (TTree*)f1->Get("syncTree_ZZctrl_Fake");

   if( isEventSync )
     {	
	syncTree_2lSS_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2lSS_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2lSS_Flip->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_3l_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_3l_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_4l_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_4l_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_0l2tau_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_0l2tau_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_1l1tau_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_1l1tau_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_1l1tau_Flip->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_1l2tau_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_1l2tau_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2lSS1tau_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2lSS1tau_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2lSS1tau_Flip->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2lOS1tau_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2lOS1tau_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2l2tau_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_2l2tau_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_3l1tau_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_3l1tau_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_ttWctrl_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_ttWctrl_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_ttWctrl_Flip->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_ttZctrl_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_ttZctrl_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_WZctrl_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_WZctrl_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_ZZctrl_SR->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);
	syncTree_ZZctrl_Fake->SetBranchAddress("is_tH_like_and_not_ttH_like",&isTHQ);

	syncTree_2lSS_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2lSS_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2lSS_Flip->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_3l_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_3l_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_4l_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_4l_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_0l2tau_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_0l2tau_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_1l1tau_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_1l1tau_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_1l1tau_Flip->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_1l2tau_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_1l2tau_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2lSS1tau_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2lSS1tau_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2lSS1tau_Flip->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2lOS1tau_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2lOS1tau_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2l2tau_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_2l2tau_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_3l1tau_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_3l1tau_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_ttWctrl_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_ttWctrl_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_ttWctrl_Flip->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_ttZctrl_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_ttZctrl_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_WZctrl_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_WZctrl_Fake->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_ZZctrl_SR->SetBranchAddress("is_ttH_like",&isTTH);
	syncTree_ZZctrl_Fake->SetBranchAddress("is_ttH_like",&isTTH);
     }   
   
   int n_syncTree_2lSS_SR = syncTree_2lSS_SR->GetEntries();
   int n_syncTree_2lSS_Fake = syncTree_2lSS_Fake->GetEntries();
   int n_syncTree_2lSS_Flip = syncTree_2lSS_Flip->GetEntries();
   int n_syncTree_3l_SR = syncTree_3l_SR->GetEntries();
   int n_syncTree_3l_Fake = syncTree_3l_Fake->GetEntries();
   int n_syncTree_4l_SR = syncTree_4l_SR->GetEntries();
   int n_syncTree_4l_Fake = syncTree_4l_Fake->GetEntries();
   int n_syncTree_0l2tau_SR = syncTree_0l2tau_SR->GetEntries();
   int n_syncTree_0l2tau_Fake = syncTree_0l2tau_Fake->GetEntries();
   int n_syncTree_1l1tau_SR = syncTree_1l1tau_SR->GetEntries();
   int n_syncTree_1l1tau_Fake = syncTree_1l1tau_Fake->GetEntries();
   int n_syncTree_1l1tau_Flip = syncTree_1l1tau_Flip->GetEntries();
   int n_syncTree_1l2tau_SR = syncTree_1l2tau_SR->GetEntries();
   int n_syncTree_1l2tau_Fake = syncTree_1l2tau_Fake->GetEntries();
   int n_syncTree_2lSS1tau_SR = syncTree_2lSS1tau_SR->GetEntries();
   int n_syncTree_2lSS1tau_Fake = syncTree_2lSS1tau_Fake->GetEntries();
   int n_syncTree_2lSS1tau_Flip = syncTree_2lSS1tau_Flip->GetEntries();
   int n_syncTree_2lOS1tau_SR = syncTree_2lOS1tau_SR->GetEntries();
   int n_syncTree_2lOS1tau_Fake = syncTree_2lOS1tau_Fake->GetEntries();
   int n_syncTree_2l2tau_SR = syncTree_2l2tau_SR->GetEntries();
   int n_syncTree_2l2tau_Fake = syncTree_2l2tau_Fake->GetEntries();
   int n_syncTree_3l1tau_SR = syncTree_3l1tau_SR->GetEntries();
   int n_syncTree_3l1tau_Fake = syncTree_3l1tau_Fake->GetEntries();
   int n_syncTree_ttWctrl_SR = syncTree_ttWctrl_SR->GetEntries();
   int n_syncTree_ttWctrl_Fake = syncTree_ttWctrl_Fake->GetEntries();
   int n_syncTree_ttWctrl_Flip = syncTree_ttWctrl_Flip->GetEntries();
   int n_syncTree_ttZctrl_SR = syncTree_ttZctrl_SR->GetEntries();
   int n_syncTree_ttZctrl_Fake = syncTree_ttZctrl_Fake->GetEntries();
   int n_syncTree_WZctrl_SR = syncTree_WZctrl_SR->GetEntries();
   int n_syncTree_WZctrl_Fake = syncTree_WZctrl_Fake->GetEntries();
   int n_syncTree_ZZctrl_SR = syncTree_ZZctrl_SR->GetEntries();
   int n_syncTree_ZZctrl_Fake = syncTree_ZZctrl_Fake->GetEntries();

   int n_syncTree_2lSS_SR_ttH = 0;
   int n_syncTree_2lSS_Fake_ttH = 0;
   int n_syncTree_2lSS_Flip_ttH = 0;
   int n_syncTree_3l_SR_ttH = 0;
   int n_syncTree_3l_Fake_ttH = 0;
   int n_syncTree_2lSS1tau_SR_ttH = 0;
   int n_syncTree_2lSS1tau_Fake_ttH = 0;
   int n_syncTree_2lSS1tau_Flip_ttH = 0;
   
   int n_syncTree_2lSS_SR_tHq = 0;
   int n_syncTree_2lSS_Fake_tHq = 0;
   int n_syncTree_2lSS_Flip_tHq = 0;
   int n_syncTree_3l_SR_tHq = 0;
   int n_syncTree_3l_Fake_tHq = 0;
   int n_syncTree_2lSS1tau_SR_tHq = 0;
   int n_syncTree_2lSS1tau_Fake_tHq = 0;
   int n_syncTree_2lSS1tau_Flip_tHq = 0;
   
   if( isEventSync )
     {	   
	for(int i=0;i<n_syncTree_2lSS_SR;i++) { 
	   syncTree_2lSS_SR->GetEntry(i);
	   if( isTHQ && !isTTH ) n_syncTree_2lSS_SR_tHq++; 
	   if( isTTH ) n_syncTree_2lSS_SR_ttH++;
	}
	for(int i=0;i<n_syncTree_2lSS_Fake;i++) { 
	   syncTree_2lSS_Fake->GetEntry(i);
	   if( isTHQ && !isTTH ) n_syncTree_2lSS_Fake_tHq++; 
	   if( isTTH ) n_syncTree_2lSS_Fake_ttH++;
	}
	for(int i=0;i<n_syncTree_2lSS_Flip;i++) { 
	   syncTree_2lSS_Flip->GetEntry(i);
	   if( isTHQ && !isTTH ) n_syncTree_2lSS_Flip_tHq++; 
	   if( isTTH ) n_syncTree_2lSS_Flip_ttH++;
	}
	for(int i=0;i<n_syncTree_3l_SR;i++) {
	   syncTree_3l_SR->GetEntry(i);
	   if( isTHQ && !isTTH ) n_syncTree_3l_SR_tHq++; 
	   if( isTTH ) n_syncTree_3l_SR_ttH++;
	}
	for(int i=0;i<n_syncTree_3l_Fake;i++) { 
	   syncTree_3l_Fake->GetEntry(i);
	   if( isTHQ && !isTTH ) n_syncTree_3l_Fake_tHq++; 
	   if( isTTH ) n_syncTree_3l_Fake_ttH++;
	}
	for(int i=0;i<n_syncTree_2lSS1tau_SR;i++) {
	   syncTree_2lSS1tau_SR->GetEntry(i);
	   if( isTHQ && !isTTH ) n_syncTree_2lSS1tau_SR_tHq++; 
	   if( isTTH ) n_syncTree_2lSS1tau_SR_ttH++;
	}
	for(int i=0;i<n_syncTree_2lSS1tau_Fake;i++) { 
	   syncTree_2lSS1tau_Fake->GetEntry(i);
	   if( isTHQ && !isTTH ) n_syncTree_2lSS1tau_Fake_tHq++; 
	   if( isTTH ) n_syncTree_2lSS1tau_Fake_ttH++;
	}
	for(int i=0;i<n_syncTree_2lSS1tau_Flip;i++) { 
	   syncTree_2lSS1tau_Flip->GetEntry(i);
	   if( isTHQ && !isTTH ) n_syncTree_2lSS1tau_Flip_tHq++; 
	   if( isTTH ) n_syncTree_2lSS1tau_Flip_ttH++;
	}
     }

   if( isEventSync )
     {
	std::cout << "2lSS_SR = " << n_syncTree_2lSS_SR_ttH << "/" << n_syncTree_2lSS_SR_tHq << std::endl;
	std::cout << "2lSS_Fake = " << n_syncTree_2lSS_Fake_ttH << "/" << n_syncTree_2lSS_Fake_tHq << std::endl;
	std::cout << "2lSS_Flip = " << n_syncTree_2lSS_Flip_ttH << "/" << n_syncTree_2lSS_Flip_tHq << std::endl;
	std::cout << "3l_SR = " << n_syncTree_3l_SR_ttH << "/" << n_syncTree_3l_SR_tHq << std::endl;
	std::cout << "3l_Fake = " << n_syncTree_3l_Fake_ttH << "/" << n_syncTree_3l_Fake_tHq << std::endl;
	std::cout << "4l_SR = " << n_syncTree_4l_SR << std::endl;
	std::cout << "4l_Fake = " << n_syncTree_4l_Fake << std::endl;
	std::cout << "0l2tau_SR = " << n_syncTree_0l2tau_SR << std::endl;
	std::cout << "0l2tau_Fake = " << n_syncTree_0l2tau_Fake << std::endl;
	std::cout << "1l1tau_SR = " << n_syncTree_1l1tau_SR << std::endl;
	std::cout << "1l1tau_Fake = " << n_syncTree_1l1tau_Fake << std::endl;
	std::cout << "1l1tau_Flip = " << n_syncTree_1l1tau_Flip << std::endl;
	std::cout << "1l2tau_SR = " << n_syncTree_1l2tau_SR << std::endl;
	std::cout << "1l2tau_Fake = " << n_syncTree_1l2tau_Fake << std::endl;
	std::cout << "2lSS1tau_SR = " << n_syncTree_2lSS1tau_SR_ttH << "/" << n_syncTree_2lSS1tau_SR_tHq << std::endl;
	std::cout << "2lSS1tau_Fake = " << n_syncTree_2lSS1tau_Fake_ttH << "/" << n_syncTree_2lSS1tau_Fake_tHq << std::endl;
	std::cout << "2lSS1tau_Flip = " << n_syncTree_2lSS1tau_Flip_ttH << "/" << n_syncTree_2lSS1tau_Flip_tHq << std::endl;
	std::cout << "2lOS1tau_SR = " << n_syncTree_2lOS1tau_SR << std::endl;
	std::cout << "2lOS1tau_Fake = " << n_syncTree_2lOS1tau_Fake << std::endl;
	std::cout << "2l2tau_SR = " << n_syncTree_2l2tau_SR << std::endl;
	std::cout << "2l2tau_Fake = " << n_syncTree_2l2tau_Fake << std::endl;
	std::cout << "3l1tau_SR = " << n_syncTree_3l1tau_SR << std::endl;
	std::cout << "3l1tau_Fake = " << n_syncTree_3l1tau_Fake << std::endl;
	std::cout << "ttWctrl_SR = " << n_syncTree_ttWctrl_SR << std::endl;
	std::cout << "ttWctrl_Fake = " << n_syncTree_ttWctrl_Fake << std::endl;
	std::cout << "ttWctrl_Flip = " << n_syncTree_ttWctrl_Flip << std::endl;
	std::cout << "ttZctrl_SR = " << n_syncTree_ttZctrl_SR << std::endl;
	std::cout << "ttZctrl_Fake = " << n_syncTree_ttZctrl_Fake << std::endl;
	std::cout << "WZctrl_SR = " << n_syncTree_WZctrl_SR << std::endl;
	std::cout << "WZctrl_Fake = " << n_syncTree_WZctrl_Fake << std::endl;
	std::cout << "ZZctrl_SR = " << n_syncTree_ZZctrl_SR << std::endl;
	std::cout << "ZZctrl_Fake = " << n_syncTree_ZZctrl_Fake << std::endl;	
     }
   else
     {	
	std::cout << "2lSS_SR = " << n_syncTree_2lSS_SR << std::endl;
	std::cout << "2lSS_Fake = " << n_syncTree_2lSS_Fake << std::endl;
	std::cout << "2lSS_Flip = " << n_syncTree_2lSS_Flip << std::endl;
	std::cout << "3l_SR = " << n_syncTree_3l_SR << std::endl;
	std::cout << "3l_Fake = " << n_syncTree_3l_Fake << std::endl;
	std::cout << "4l_SR = " << n_syncTree_4l_SR << std::endl;
	std::cout << "4l_Fake = " << n_syncTree_4l_Fake << std::endl;
	std::cout << "0l2tau_SR = " << n_syncTree_0l2tau_SR << std::endl;
	std::cout << "0l2tau_Fake = " << n_syncTree_0l2tau_Fake << std::endl;
	std::cout << "1l1tau_SR = " << n_syncTree_1l1tau_SR << std::endl;
	std::cout << "1l1tau_Fake = " << n_syncTree_1l1tau_Fake << std::endl;
	std::cout << "1l1tau_Flip = " << n_syncTree_1l1tau_Flip << std::endl;
	std::cout << "1l2tau_SR = " << n_syncTree_1l2tau_SR << std::endl;
	std::cout << "1l2tau_Fake = " << n_syncTree_1l2tau_Fake << std::endl;
	std::cout << "2lSS1tau_SR = " << n_syncTree_2lSS1tau_SR << std::endl;
	std::cout << "2lSS1tau_Fake = " << n_syncTree_2lSS1tau_Fake << std::endl;
	std::cout << "2lSS1tau_Flip = " << n_syncTree_2lSS1tau_Flip << std::endl;
	std::cout << "2lOS1tau_SR = " << n_syncTree_2lOS1tau_SR << std::endl;
	std::cout << "2lOS1tau_Fake = " << n_syncTree_2lOS1tau_Fake << std::endl;
	std::cout << "2l2tau_SR = " << n_syncTree_2l2tau_SR << std::endl;
	std::cout << "2l2tau_Fake = " << n_syncTree_2l2tau_Fake << std::endl;
	std::cout << "3l1tau_SR = " << n_syncTree_3l1tau_SR << std::endl;
	std::cout << "3l1tau_Fake = " << n_syncTree_3l1tau_Fake << std::endl;
	std::cout << "ttWctrl_SR = " << n_syncTree_ttWctrl_SR << std::endl;
	std::cout << "ttWctrl_Fake = " << n_syncTree_ttWctrl_Fake << std::endl;
	std::cout << "ttWctrl_Flip = " << n_syncTree_ttWctrl_Flip << std::endl;
	std::cout << "ttZctrl_SR = " << n_syncTree_ttZctrl_SR << std::endl;
	std::cout << "ttZctrl_Fake = " << n_syncTree_ttZctrl_Fake << std::endl;
	std::cout << "WZctrl_SR = " << n_syncTree_WZctrl_SR << std::endl;
	std::cout << "WZctrl_Fake = " << n_syncTree_WZctrl_Fake << std::endl;
	std::cout << "ZZctrl_SR = " << n_syncTree_ZZctrl_SR << std::endl;
	std::cout << "ZZctrl_Fake = " << n_syncTree_ZZctrl_Fake << std::endl;
     }   
   
   gApplication->Terminate();
}
