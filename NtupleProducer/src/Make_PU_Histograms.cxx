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
    TString list_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/lists/";
	TString output_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/PU_SF/profiles_MC/";

    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
        int nFiles_notFound = 0;

		TString output_name = output_dir+v_process[iproc]+".root";
		TFile* f_output = new TFile(output_name, "RECREATE");

		TH1F* h_PU = new TH1F("", "", 100, 0, 100);

	    for(int ifile=0; ifile<1000; ifile++)
        {
            TString filename = list_dir + v_process[iproc] + "_ID" + Convert_Number_To_TString(ifile) + ".txt";
            if(!Check_File_Existence(filename) ) {cout<<"File "<<filename<<"not found !"<<endl; nFiles_notFound++;}
            if(nFiles_notFound >= 5) {break;} //means that we reached the final file already

			ifstream file_in(filename.Data() );

			string line;
			while(getline(file_in, line) )
			{
				TString ts = line;
				cout<<endl<<"- Opening file "<<ts<<"..."<<endl;

				TFile* f = TFile::Open(ts);

				if(!f) {continue;}

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
				f->Close();
			}

        }

		f_output->cd();
		h_PU->Write("Pileup_MC");

		delete h_PU; h_PU = NULL;
		f_output->Close();

		cout<<"---> Saved file "<<output_name<<endl;
    }

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

    vector<TString> v_process;

	// v_process.push_back("THQ_4f_Hincl_13TeV_madgraph_pythia8"); //
	// v_process.push_back("THW_5f_Hincl_13TeV_madgraph_pythia8"); //
	// v_process.push_back("ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8"); //
	// v_process.push_back("tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8"); //
	// v_process.push_back("ZZZ_TuneCP5_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8"); //
	// v_process.push_back("ZZTo4L_13TeV_powheg_pythia8_ext");
	// v_process.push_back("ZZTo2L2Nu_13TeV_powheg_pythia8"); //
	// v_process.push_back("WZZ_TuneCP5_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8");//
	// v_process.push_back("WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8"); //
	// v_process.push_back("TTWW_TuneCP5_13TeV-madgraph-pythia8"); //
	// v_process.push_back("TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8");
	// v_process.push_back("TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8");
	// v_process.push_back("TTTT_TuneCP5_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8"); //
	 v_process.push_back("TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8"); //
	// v_process.push_back("TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"); //
	// v_process.push_back("TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8"); //
	// v_process.push_back("ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8"); //
	// v_process.push_back("ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"); //
	// v_process.push_back("ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"); //
	// v_process.push_back("ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8"); //
	// v_process.push_back("ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8"); //

    Make_PU_Histograms(v_process);

	return 0;
}
