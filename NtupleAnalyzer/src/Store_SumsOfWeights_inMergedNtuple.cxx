/* BASH COLORS */
#define RST   "[0m"
#define KRED  "[31m"
#define KGRN  "[32m"
#define KYEL  "[33m"
#define KBLU  "[34m"
#define KMAG  "[35m"
#define KCYN  "[36m"
#define KWHT  "[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "[1m" x RST
#define UNDL(x) "[4m" x RST

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TTree.h"
#include "TString.h"
#include "TColor.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLine.h"
#include "THStack.h"
#include "TString.h"
#include "TObject.h"


#include <iostream>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>


#include <cassert>     //Can be used to terminate program if argument is not true.
//Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

using namespace std;

//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
}

void Get_And_Save_Histogram(TString samplename)
{
    TString weightfilepath = "./weights_2017/merged_histograms/" + samplename + ".root";
    if(!Check_File_Existence(weightfilepath)) {cout<<BOLD(FRED("Error ! Weight file "<<weightfilepath<<" not found !"))<<endl; return;}

    TString outputfilepath = "./merged_ntuples/merged_" + samplename + ".root";
    if(!Check_File_Existence(outputfilepath)) {cout<<BOLD(FRED("Error ! Weight file "<<outputfilepath<<" not found !"))<<endl; return;}

    //Input file
    TFile* f_input = TFile::Open(weightfilepath);

    //Output file
    TFile *f_output = TFile::Open(outputfilepath, "UPDATE");

    TH1F* hSumWeights = (TH1F*) f_input->Get("hSumWeights");
    TH1F* hLHE = (TH1F*) f_input->Get("hLHE");
    // TH1F* hpileup = (TH1F*) f_input->Get("hpileup");

    f_output->cd();
    if(!f_input->GetListOfKeys()->Contains("hSumWeights") ) {cout<<"Histogram hSumWeights not found in file "<<weightfilepath<<" !"<<endl; return;}
    else {hSumWeights->Write("hSumWeights", TObject::kOverwrite);}
    if(f_input->GetListOfKeys()->Contains("hLHE") ) {hLHE->Write("hLHE", TObject::kOverwrite);} //Not necessarily saved in file

    f_input->Close();
    f_output->Close();

    return;
}

int main(int argc, char **argv)
{
    if(argc != 2) {return 1;}

    Get_And_Save_Histogram(argv[1]);

    return 0;
}
