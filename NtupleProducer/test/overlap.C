void overlap()
{
   gROOT->SetBatch();
   
   TFile *f = TFile::Open("output_sync.root");
   
   TH1F *h_overlap = (TH1F*)f->Get("overlap");
   
   gStyle->SetOptStat(0);
   
   TCanvas *c1 = new TCanvas();
   
   h_overlap->Scale(1./h_overlap->Integral());
   h_overlap->Draw("COLZ");
   
   c1->Print("overlap.eps");
   
   gApplication->Terminate();
}
