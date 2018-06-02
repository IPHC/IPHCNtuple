{
   gROOT->SetBatch();

   TFile *f1 = TFile::Open("output_sync.root");
   
   TTree *syncTree_2lSS_SR = (TTree*)f1->Get("syncTree_2lSS_SR");
   TTree *syncTree_2lSS_Fake = (TTree*)f1->Get("syncTree_2lSS_Fake");
   TTree *syncTree_2lSS_Flip = (TTree*)f1->Get("syncTree_2lSS_Flip");
   TTree *syncTree_3l_SR = (TTree*)f1->Get("syncTree_3l_SR");
   TTree *syncTree_3l_Fake = (TTree*)f1->Get("syncTree_3l_Fake");
   TTree *syncTree_4l_SR = (TTree*)f1->Get("syncTree_4l_SR");
   TTree *syncTree_4l_Fake = (TTree*)f1->Get("syncTree_4l_Fake");
   TTree *syncTree_1l2tau_SR = (TTree*)f1->Get("syncTree_1l2tau_SR");
   TTree *syncTree_1l2tau_Fake = (TTree*)f1->Get("syncTree_1l2tau_Fake");
   TTree *syncTree_2lSS1tau_SR = (TTree*)f1->Get("syncTree_2lSS1tau_SR");
   TTree *syncTree_2lSS1tau_Fake = (TTree*)f1->Get("syncTree_2lSS1tau_Fake");
   TTree *syncTree_2lSS1tau_Flip = (TTree*)f1->Get("syncTree_2lSS1tau_Flip");
   TTree *syncTree_2l2tau_SR = (TTree*)f1->Get("syncTree_2l2tau_SR");
   TTree *syncTree_2l2tau_Fake = (TTree*)f1->Get("syncTree_2l2tau_Fake");
   TTree *syncTree_3l1tau_SR = (TTree*)f1->Get("syncTree_3l1tau_SR");
   TTree *syncTree_3l1tau_Fake = (TTree*)f1->Get("syncTree_3l1tau_Fake");
   TTree *syncTree_ttWctrl_SR = (TTree*)f1->Get("syncTree_ttWctrl_SR");
   TTree *syncTree_ttWctrl_Fake = (TTree*)f1->Get("syncTree_ttWctrl_Fake");
   TTree *syncTree_ttWctrl_Flip = (TTree*)f1->Get("syncTree_ttWctrl_Flip");
   TTree *syncTree_ttZctrl_SR = (TTree*)f1->Get("syncTree_ttZctrl_SR");
   TTree *syncTree_ttZctrl_Fake = (TTree*)f1->Get("syncTree_ttZctrl_Fake");

   int n_syncTree_2lSS_SR = syncTree_2lSS_SR->GetEntries();
   int n_syncTree_2lSS_Fake = syncTree_2lSS_Fake->GetEntries();
   int n_syncTree_2lSS_Flip = syncTree_2lSS_Flip->GetEntries();
   int n_syncTree_3l_SR = syncTree_3l_SR->GetEntries();
   int n_syncTree_3l_Fake = syncTree_3l_Fake->GetEntries();
   int n_syncTree_4l_SR = syncTree_4l_SR->GetEntries();
   int n_syncTree_4l_Fake = syncTree_4l_Fake->GetEntries();
   int n_syncTree_1l2tau_SR = syncTree_1l2tau_SR->GetEntries();
   int n_syncTree_1l2tau_Fake = syncTree_1l2tau_Fake->GetEntries();
   int n_syncTree_2lSS1tau_SR = syncTree_2lSS1tau_SR->GetEntries();
   int n_syncTree_2lSS1tau_Fake = syncTree_2lSS1tau_Fake->GetEntries();
   int n_syncTree_2lSS1tau_Flip = syncTree_2lSS1tau_Flip->GetEntries();
   int n_syncTree_2l2tau_SR = syncTree_2l2tau_SR->GetEntries();
   int n_syncTree_2l2tau_Fake = syncTree_2l2tau_Fake->GetEntries();
   int n_syncTree_3l1tau_SR = syncTree_3l1tau_SR->GetEntries();
   int n_syncTree_3l1tau_Fake = syncTree_3l1tau_Fake->GetEntries();
   int n_syncTree_ttWctrl_SR = syncTree_ttWctrl_SR->GetEntries();
   int n_syncTree_ttWctrl_Fake = syncTree_ttWctrl_Fake->GetEntries();
   int n_syncTree_ttWctrl_Flip = syncTree_ttWctrl_Flip->GetEntries();
   int n_syncTree_ttZctrl_SR = syncTree_ttZctrl_SR->GetEntries();
   int n_syncTree_ttZctrl_Fake = syncTree_ttZctrl_Fake->GetEntries();
   
   std::cout << "2lSS_SR = " << n_syncTree_2lSS_SR << std::endl;
   std::cout << "2lSS_Fake = " << n_syncTree_2lSS_Fake << std::endl;
   std::cout << "2lSS_Flip = " << n_syncTree_2lSS_Flip << std::endl;
   std::cout << "3l_SR = " << n_syncTree_3l_SR << std::endl;
   std::cout << "3l_Fake = " << n_syncTree_3l_Fake << std::endl;
   std::cout << "4l_SR = " << n_syncTree_4l_SR << std::endl;
   std::cout << "4l_Fake = " << n_syncTree_4l_Fake << std::endl;
   std::cout << "1l2tau_SR = " << n_syncTree_1l2tau_SR << std::endl;
   std::cout << "1l2tau_Fake = " << n_syncTree_1l2tau_Fake << std::endl;
   std::cout << "2lSS1tau_SR = " << n_syncTree_2lSS1tau_SR << std::endl;
   std::cout << "2lSS1tau_Fake = " << n_syncTree_2lSS1tau_Fake << std::endl;
   std::cout << "2lSS1tau_Flip = " << n_syncTree_2lSS1tau_Flip << std::endl;
   std::cout << "2l2tau_SR = " << n_syncTree_2l2tau_SR << std::endl;
   std::cout << "2l2tau_Fake = " << n_syncTree_2l2tau_Fake << std::endl;
   std::cout << "3l1tau_SR = " << n_syncTree_3l1tau_SR << std::endl;
   std::cout << "3l1tau_Fake = " << n_syncTree_3l1tau_Fake << std::endl;
   std::cout << "ttWctrl_SR = " << n_syncTree_ttWctrl_SR << std::endl;
   std::cout << "ttWctrl_Fake = " << n_syncTree_ttWctrl_Fake << std::endl;
   std::cout << "ttWctrl_Flip = " << n_syncTree_ttWctrl_Flip << std::endl;
   std::cout << "ttZctrl_SR = " << n_syncTree_ttZctrl_SR << std::endl;
   std::cout << "ttZctrl_Fake = " << n_syncTree_ttZctrl_Fake << std::endl;
   
   gApplication->Terminate();
}
