#ifndef tHq_analysis_h
#define tHq_analysis_h


#include <TString.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TFile.h>
#include <TObject.h>
#include <iostream>
#include <fstream>

#include <sys/stat.h> //for mkdir, stat

class ScaleFactors
{
	public :

    ScaleFactors(); //Constructor
    ~ScaleFactors(); //Destructor

//-- Methods
    float Get_SF_RecoToLoose_Ele(float, float);
    float Get_SF_RecoToLoose_Mu(float, float);
    float Get_SF_LooseToTight_Leptons(int, int, float, float);
    float Get_Lepton_SF(int, int, float, float);
    float Get_Trigger_SF(int, int, float, int, float);

//-- Members
    TString path_dir;

    TFile* file_Eff_Reco_Ele_ptLt20;
    TFile* file_Eff_Reco_Ele_ptGt20;
    TFile* file_Eff_LooseID_Ele;
    TFile* file_Eff_Track_Mu_ptLt10;
    TFile* file_Eff_Track_Mu_ptGt10;
    TFile* file_Eff_LooseID_Mu_ptLt30;
    TFile* file_Eff_LooseID_Mu_ptGt30;
    TFile* file_Eff_Iso_Mu;
    TFile* file_Eff_LooseToTight_Ele_2lss;
    TFile* file_Eff_LooseToTight_Mu_2lss;
    TFile* file_Eff_LooseToTight_Ele_3l;
    TFile* file_Eff_LooseToTight_Mu_3l;

    TH2D* hEff_Reco_Ele_ptLt20;
    TH2D* hEff_Reco_Ele_ptGt20;
    TH2D* hEff_LooseID_Ele;
    TGraph* gEff_Track_Mu_ptLt10;
    TGraph* gEff_Track_Mu_ptGt10;
    TH2* hEff_LooseID_Mu_ptLt30;
    TH2* hEff_LooseID_Mu_ptGt30;
    TH2* hEff_Iso_Mu;
    TH2F* hEff_LooseToTight_Ele_2lss;
    TH2F* hEff_LooseToTight_Mu_2lss;
    TH2F* hEff_LooseToTight_Ele_3l;
    TH2F* hEff_LooseToTight_Mu_3l;
};

#endif
