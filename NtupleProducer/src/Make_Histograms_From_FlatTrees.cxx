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
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TTree.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TString.h"
#include "TString.h"
#include "TF1.h"

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

#include "../../NtupleAnalyzer/include/Pileup.h"

using namespace std;

//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
}


//Convert a double into a TString
// precision --> can choose if TString how many digits the TString should display
TString Convert_Number_To_TString(double number, int precision=5)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}

//Convert a TString into a float
double Convert_TString_To_Number(TString ts)
{
	double number = 0;
	string s = ts.Data();
	stringstream ss(s);
	ss >> number;
	return number;
}


void Make_PU_Histograms(vector<TString> v_process)
{
	cout<<"== WILL CREATE PILEUP PROFILES FOR ALL SAMPLES IN VECTOR =="<<endl<<endl;

	bool read_pileup_histograms_from_flatTrees = true; //true <-> read PU info from histogram ; //else produce histo from TTree (read all entries)

    TString list_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/lists/";
	TString output_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/PU_SF/profiles_MC/";

    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
		cout<<endl<<"-- Sample : "<<v_process[iproc]<<endl;

		TH1F* h_PU = 0;
		if(!read_pileup_histograms_from_flatTrees) {h_PU = new TH1F("", "", 100, 0, 100);}

		int nFiles_notFound = 0;
		int ifile = 0;
		while(nFiles_notFound < 5) //Stop if not finding 5 consecutive files (must have reached the end)
        {
            TString filename = list_dir + v_process[iproc] + "_ID" + Convert_Number_To_TString(ifile) + ".txt";
			ifile++;
            if(!Check_File_Existence(filename) ) {cout<<"File "<<filename<<" not found !"<<endl; nFiles_notFound++; continue;}

			ifstream file_in(filename.Data() );

			string line;
			while(getline(file_in, line) )
			{
				TString ts = line;
				cout<<endl<<"- Opening file "<<ts<<"..."<<endl;

				TFile* f = TFile::Open(ts);
				if(!f) {continue;}

                cout<<"-> Opened"<<endl<<endl;

				//--- OBSOLETE : DIRECTLY READ THE PILEUP VARIABLE IN TTree
				if(!read_pileup_histograms_from_flatTrees)
				{

					TTree* t = (TTree*) f->Get("FlatTree/tree");
					Int_t nPU;
					Float_t mc_weight;
					t->SetBranchAddress("mc_pu_trueNumInt", &nPU);
					t->SetBranchAddress("mc_weight", &mc_weight);
					int nentries = t->GetEntries();
					for(int ientry=0; ientry<nentries; ientry++)
					{
						// if(ientry>10000) {break;}

						if(ientry%20000==0) {cout<<"Event : "<<ientry<<" / "<<nentries<<endl;}

						t->GetEntry(ientry);

						h_PU->Fill(nPU, mc_weight);
					}
					delete t; t = NULL;
				}

				//--- NOW : store nPU histo at FlatTree level, so can just merge histos from all files
				if(read_pileup_histograms_from_flatTrees)
				{
                    TH1F* h_tmp = (TH1F*) f->Get("FlatTree/hpileup")->Clone();
					h_tmp->SetDirectory(0); //So that hist is not associated to file (& deleted when we close it)

                    if(!h_PU) {h_PU = (TH1F*) h_tmp;}
					else {h_PU->Add(h_tmp);}
				}

				f->Close();
			}
        }

		TString output_name = output_dir+v_process[iproc]+".root";
		TFile* f_output = 0;
		if(h_PU && h_PU->GetEntries() > 0)
		{
			f_output = new TFile(output_name, "RECREATE");
			f_output->cd();
			h_PU->Write("Pileup_MC");
			f_output->Close();
		}

		if(!read_pileup_histograms_from_flatTrees) {delete h_PU; h_PU = NULL;}

		cout<<"---> Saved file "<<output_name<<endl;
    }

	return;
}



void Make_SumWeightScales_Histograms(vector<TString> v_process)
{
	cout<<"== WILL CREATE SUM-OF-WEIGHT HISTOGRAMS OF SCALE VARIATIONS FOR ALL SAMPLES IN VECTOR =="<<endl<<endl;

    TString list_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/lists/";
	TString output_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/Scales/";

    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
		cout<<endl<<"-- Sample : "<<v_process[iproc]<<endl;

        TH1F* h_sumWeights_nominal = new TH1F("", "", 5, -2, 2);
        TH1F* h_sumWeightsScale_originalXWGTUP = new TH1F("", "", 5, -2, 2);
        TH1F* h_sumWeightsScale_muF0p5 = new TH1F("", "", 5, -2, 2);
        TH1F* h_sumWeightsScale_muF2 = new TH1F("", "", 5, -2, 2);
        TH1F* h_sumWeightsScale_muR0p5 = new TH1F("", "", 5, -2, 2);
        TH1F* h_sumWeightsScale_muR2 = new TH1F("", "", 5, -2, 2);
        TH1F* h_sumWeightsScale_muR2muF2 = new TH1F("", "", 5, -2, 2);
        TH1F* h_sumWeightsScale_muR0p5muF0p5 = new TH1F("", "", 5, -2, 2);

		int nFiles_notFound = 0;
		int ifile = 0;
		while(nFiles_notFound < 5) //Stop if not finding 5 consecutive files (must have reached the end)
        {
            // if(ifile == 1) {break;}

            TString filename = list_dir + v_process[iproc] + "_ID" + Convert_Number_To_TString(ifile) + ".txt";
			ifile++;
            if(!Check_File_Existence(filename) ) {cout<<"File "<<filename<<" not found !"<<endl; nFiles_notFound++; continue;}

			ifstream file_in(filename.Data() );

			string line;
			while(getline(file_in, line) )
			{
				TString ts = line;
				cout<<endl<<"- Opening file "<<ts<<"..."<<endl;

				TFile* f = TFile::Open(ts);
				if(!f) {continue;}
                cout<<"-> Opened"<<endl<<endl;

				TTree* t = (TTree*) f->Get("FlatTree/tree");

                Float_t weight_originalXWGTUP;
                Float_t weight_scale_muR0p5;
                Float_t weight_scale_muF0p5;
                Float_t weight_scale_muR2;
                Float_t weight_scale_muF2;
                Float_t weight_scale_muR0p5muF0p5;
                Float_t weight_scale_muR2muF2;
				Float_t mc_weight;

                t->SetBranchAddress("weight_originalXWGTUP", &weight_originalXWGTUP);
                t->SetBranchAddress("weight_scale_muR0p5", &weight_scale_muR0p5);
                t->SetBranchAddress("weight_scale_muF0p5", &weight_scale_muF0p5);
                t->SetBranchAddress("weight_scale_muR2", &weight_scale_muR2);
                t->SetBranchAddress("weight_scale_muF2", &weight_scale_muF2);
                t->SetBranchAddress("weight_scale_muR0p5muF0p5", &weight_scale_muR0p5muF0p5);
				t->SetBranchAddress("weight_scale_muR2muF2", &weight_scale_muR2muF2);
				t->SetBranchAddress("mc_weight", &mc_weight);

				int nentries = t->GetEntries();
				for(int ientry=0; ientry<nentries; ientry++)
				{
					// if(ientry>100) {break;}

					if(ientry%20000==0) {cout<<"Event : "<<ientry<<" / "<<nentries<<endl;}

					t->GetEntry(ientry);

                    h_sumWeights_nominal->Fill(mc_weight, 1);
                    h_sumWeightsScale_originalXWGTUP->Fill(mc_weight, weight_originalXWGTUP);
                    h_sumWeightsScale_muF0p5->Fill(mc_weight,weight_scale_muF0p5);
                    h_sumWeightsScale_muF2->Fill(mc_weight, weight_scale_muF2);
                    h_sumWeightsScale_muR0p5->Fill(mc_weight, weight_scale_muR0p5);
                    h_sumWeightsScale_muR2->Fill(mc_weight, weight_scale_muR2);
                    h_sumWeightsScale_muR2muF2->Fill(mc_weight, weight_scale_muR2muF2);
                    h_sumWeightsScale_muR0p5muF0p5->Fill(mc_weight, weight_scale_muR0p5muF0p5);
                }

				delete t; t = NULL;

				f->Close();
			} //loop on files
        } //loop on lists

		TString output_name = output_dir+v_process[iproc]+".root";
		TFile* f_output = 0;

		f_output = new TFile(output_name, "RECREATE");
		f_output->cd();

        h_sumWeights_nominal->Write("h_sumWeights_nominal");
        h_sumWeightsScale_originalXWGTUP->Write("h_sumWeightsScale_originalXWGTUP");
        h_sumWeightsScale_muF0p5->Write("h_sumWeightsScale_muF0p5");
        h_sumWeightsScale_muF2->Write("h_sumWeightsScale_muF2");
        h_sumWeightsScale_muR0p5->Write("h_sumWeightsScale_muR0p5");
        h_sumWeightsScale_muR2->Write("h_sumWeightsScale_muR2");
        h_sumWeightsScale_muR2muF2->Write("h_sumWeightsScale_muR2muF2");
        h_sumWeightsScale_muR0p5muF0p5->Write("h_sumWeightsScale_muR0p5muF0p5");

		f_output->Close();

        delete h_sumWeights_nominal;
        delete h_sumWeightsScale_originalXWGTUP;
        delete h_sumWeightsScale_muF0p5;
        delete h_sumWeightsScale_muF2;
        delete h_sumWeightsScale_muR0p5;
        delete h_sumWeightsScale_muR2;
        delete h_sumWeightsScale_muR2muF2;
        delete h_sumWeightsScale_muR0p5muF0p5;

		cout<<"---> Saved file "<<output_name<<endl;
    } //process loop

	return;
}


// void Create_PU_WeightFiles(vector<TString> v_process)
// {
// 	for(int iproc=0; iproc<v_process.size(); iproc++)
// 	{
// 		Pileup* pu = new Pileup();
// 		pu->Init(v_process[iproc]);
// 		pu->Fill_MC_vector_from_Pileup_Profile();
// 		pu->Write_PU_Weight_File();
//
// 		delete pu; pu = NULL;
// 	}
//
// 	return;
// }


int main(int argc, char *argv[])
{
	// if(argc != 2) {cout<<"Wrong arg ! Abort"<<endl; return 1;}
	// v_process.push_back(argv[1]);

//--------------------------------------------
    vector<TString> v_process;

	// v_process.push_back("THQ_4f_Hincl_13TeV_madgraph_pythia8");
	// v_process.push_back("THW_5f_Hincl_13TeV_madgraph_pythia8");

    //v_process.push_back("THQ_ctcvcp_4f_Hincl_13TeV_madgraph_pythia8");
    //v_process.push_back("THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8");
     //v_process.push_back("ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8");
    // v_process.push_back("ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8");
    // v_process.push_back("ttH_M125_TuneCP5_13TeV-powheg-pythia8");
    // v_process.push_back("ST_FCNC-TH_Tleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8");
    // v_process.push_back("ST_FCNC-TH_Tleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8");
    // v_process.push_back("TT_FCNC-TtoHJ_aTleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8");
    // v_process.push_back("TT_FCNC-TtoHJ_aTleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8");
    // v_process.push_back("TT_FCNC-aTtoHJ_Tleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8");
    // v_process.push_back("TT_FCNC-aTtoHJ_Tleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8");
    // v_process.push_back("TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8");
    // v_process.push_back("ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8");
    // v_process.push_back("TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8");
    // v_process.push_back("WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8");
	// v_process.push_back("tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8");
    // v_process.push_back("ZZTo2L2Nu_13TeV_powheg_pythia8");
    // v_process.push_back("ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v2_MINIAODSIM");
    // v_process.push_back("ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
    // v_process.push_back("TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8");
	// v_process.push_back("TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8");
    // v_process.push_back("WZZ_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("ZZZ_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8");
	// v_process.push_back("WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("TTTT_TuneCP5_13TeV-amcatnlo-pythia8");
     v_process.push_back("TTWW_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TTWH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TTTW_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8");
    // v_process.push_back("WpWpJJ_EWK-QCD_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("WW_DoubleScattering_13TeV-pythia8_TuneCP5");
    // v_process.push_back("WZG_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8_Fall17");
    // v_process.push_back("VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_Fall17");
    // v_process.push_back("DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_v1_MINIAODSIM");
    // v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
    // v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM");
    // v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
    // v_process.push_back("TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8");
    // v_process.push_back("TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8");
    // v_process.push_back("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8");
    // v_process.push_back("TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8");
    // v_process.push_back("TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8");
    // v_process.push_back("TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8");

	// v_process.push_back("");
//--------------------------------------------


//--------------------------------------------
    // Make_PU_Histograms(v_process);
//--------------------------------------------
	Make_SumWeightScales_Histograms(v_process);
//--------------------------------------------
	return 0;
}
