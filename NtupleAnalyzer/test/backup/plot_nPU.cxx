/* BASH COLORS */
#define RST   "[0m"
#define KRED  "[31m"
#define KGRN  "[32m"
#define KYEL  "[33m"
#define KBLU  "[34m"
#define KMAG  "[35m"
#define KCYN  "[36m"
#define KWHT  "[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "[1m" x RST
#define UNDL(x) "[4m" x RST

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TTree.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TString.h"
#include "TColor.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLine.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TObject.h"
#include "TRandom3.h"
#include "TRandom1.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TLegendEntry.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include <iostream>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Config.h"

#include <cassert>     //Can be used to terminate program if argument is not true.
//Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

using namespace std;

void plot_PU(bool corr)
{
    TFile* f = TFile::Open("./merged_ntuples/ttZ.root");
    TTree* t = (TTree*) f->Get("Tree");

    int nentries = t->GetEntries();

    Float_t weight;
    Float_t PU_SF;
    Int_t nPU;

    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("PU_SF", &PU_SF);
    t->SetBranchAddress("nPU", &nPU);

    TH1F* h = new TH1F("", "", 98, 0, 98);

    TCanvas* c = new TCanvas("", "", 1000, 800);

    for(int ientry=0; ientry<nentries; ientry++)
    {
        weight = 0; nPU = 0; PU_SF = 0;
        t->GetEntry(ientry);

        if(corr) {weight*= PU_SF;}

        h->Fill(nPU, weight);
    }

    h->Draw();

    if(corr) {c->SaveAs("nPU_withCorr.png");}
    else {c->SaveAs("nPU_noCorr.png");}

    delete c;
    delete h;

    return;
}

int main()
{
    plot_PU(true);
    plot_PU(false);
}
