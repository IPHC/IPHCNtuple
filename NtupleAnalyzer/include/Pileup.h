#ifndef PileUp_h
#define PileUp_h

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

#include <cassert>     //Can be used to terminate program if argument is not true //Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

using namespace std;

class Pileup
{
    public:

    //--- METHODS
    Pileup(); //Constructor
    ~Pileup(); //Destructor
    bool Init(TString); //Open files, histograms
    void Fill_MC_vector_from_HardCoded(); //Fill MC PU vector from hard-coded values
    bool Fill_MC_vector_from_Pileup_Profile();
    void Read_PU_Weights(); //Read PU weights from text file
    void Compute_PU_Weights();
    TString GetFilename(TString);
    void Print_Debug();
    // void Write_PU_Weight_File();

    //--- MEMBERS
    TString dir_parent;

    TString samplename;

    vector<double> v_MC_PU; //MC PU distribution //Take it from CMSSW config for 2017 campaign

    //Store final PU weights
    vector<double> v_PU_weights_Nom;
    vector<double> v_PU_weights_Up;
    vector<double> v_PU_weights_Down;

    //Paths of files and histograms
    TString fname_PU_Data_Nom;
    TString fname_PU_Data_Up;
    TString fname_PU_Data_Down;
    TString hname_PU_Data_Nom;
    TString hname_PU_Data_Up;
    TString hname_PU_Data_Down;

    //Files and histograms
    TFile* f_PU_Data_Nom;
    TFile* f_PU_Data_Up;
    TFile* f_PU_Data_Down;
    TH1F* h_PU_Data_Nom;
    TH1F* h_PU_Data_Up;
    TH1F* h_PU_Data_Down;

    //Vectors filled simultaneously --> Get correspondance b/w given sample name, its EOS path, and the nof associated files
    // vector<TString> v_samples; vector<TString> v_paths; vector<int> v_nof_files;
};


#endif
