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

#include <cassert> 	//Can be used to terminate program if argument is not true.
//Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

using namespace std;


// ------------------------------
//  ######   #######  ##     ## ##     ##  #######  ##    ##    ######## ##     ## ##    ##  ######
// ##    ## ##     ## ###   ### ###   ### ##     ## ###   ##    ##       ##     ## ###   ## ##    ##
// ##       ##     ## #### #### #### #### ##     ## ####  ##    ##       ##     ## ####  ## ##
// ##       ##     ## ## ### ## ## ### ## ##     ## ## ## ##    ######   ##     ## ## ## ## ##
// ##       ##     ## ##     ## ##     ## ##     ## ##  ####    ##       ##     ## ##  #### ##
// ##    ## ##     ## ##     ## ##     ## ##     ## ##   ###    ##       ##     ## ##   ### ##    ##
//  ######   #######  ##     ## ##     ##  #######  ##    ##    ##        #######  ##    ##  ######
// ------------------------------


//Convert a double into a TString
// precision --> can choose if TString how many digits the TString should display
TString Convert_Number_To_TString(double number, int precision=10)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}

//Convert a TString into a double
double Convert_TString_To_Number(TString ts)
{
	double number = 0;
	string s = ts.Data();
	stringstream ss(s);
	ss >> number;
	return number;
}


//Convert a Hexadecimal TString into an unsigned long
unsigned long Convert_Hexa_To_UnsignedLong(TString hex_TString)
{
	return std::strtoul(hex_TString.Data(), 0, 16);
}

//Convert an unsigned long into hexa TString
TString Convert_UnsignedLong_To_Hexa(unsigned long number)
{
	TString ts = "";
	stringstream ss;
	ss << std::hex << number;
	ts = ss.str();
	return ts;
}


//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
}







void Sum_Weights_FinalState(vector<TString> v_samples)
{
	cout<<FYEL("--- Will count the yields in Final State for each sample in the list ---")<<endl;
	cout<<"(Make sure you have properly merged the NTAnalyzer output files in 1 single file named after the sample !)"<<endl;

	TString dir_ntuples = "/opt/sbg/scratch1/cms/ntonon/Analyzer_ntuples_tHq/toy_tHqAnalysis";

	ofstream file_out("count_events_FS.txt");

	for(int isample=0; isample<v_samples.size(); isample++)
	{
		TString filename = v_samples[isample]+".root";

		TString filepath = dir_ntuples + "/" +  v_samples[isample] + "/" + filename;

		if(!Check_File_Existence(filepath) )
		{
			cout<<FRED("File "<<filepath<<" not found !")<<endl;
			continue;
		}

		TFile* f = new TFile(filepath);
		TTree* t = (TTree*) f->Get("Tree");

		double sum = 0;
		Float_t weight = 0;
		Char_t is_3l_THQ_SR = 0; //boolean foir category

		t->SetBranchAddress("is_3l_THQ_SR", &is_3l_THQ_SR);
		t->SetBranchAddress("weight", &weight);

		int nentries = t->GetEntries();

		for(int ientry=0; ientry<nentries; ientry++)
		{
			weight=0; is_3l_THQ_SR = false;

			t->GetEntry(ientry);

			if(is_3l_THQ_SR) sum+= weight;
		}

		file_out<<"** Sample : "<<v_samples[isample]<<endl;
		file_out<<" ----> nentries = "<<nentries<<" / SUM OF WEIGHTS = "<<sum<<endl<<endl;
	}
}

int main()
{
	vector<TString> v_samples;
	v_samples.push_back("THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1");
	// v_samples.push_back("THW_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1");
	// v_samples.push_back("TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM");
	// v_samples.push_back("TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM");
	// v_samples.push_back("ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix");
	// v_samples.push_back("TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8");
	// v_samples.push_back("TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8");
	// v_samples.push_back("WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM");
	// v_samples.push_back("WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8");
	// v_samples.push_back("WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8");
	// v_samples.push_back("ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8");
	// v_samples.push_back("ZZTo4L_13TeV_powheg_pythia8");
	// v_samples.push_back("tZq_ll_4f_13TeV-amcatnlo-herwigpp");
	// v_samples.push_back("ST_tWll_5f_LO_13TeV-MadGraph-pythia8");
	// v_samples.push_back("DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
	// v_samples.push_back("DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8");
	// v_samples.push_back("TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
	// v_samples.push_back("TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
	// v_samples.push_back("TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");

  	Sum_Weights_FinalState(v_samples);
}
