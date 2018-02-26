#ifndef tHq_analysis_h
#define tHq_analysis_h

/* BASH COLORS */
#define RST   "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST

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
#include "TMVA/DataLoader.h" //NEW -- TO BE USED INSTEAD OF FACTORY IN NEW ROOT VERSIONS ??
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Config.h"

#include <cassert> 	//Can be used to terminate program if argument is not true. Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

#include "Func_other.h" //Helper functions

using namespace std;

class tHq_analysis
{

	public :

	tHq_analysis(); //Default constructor
	tHq_analysis(vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<int>, TString, TString, TString, int, bool); //Overloaded constructor
	~tHq_analysis(); //Default destructor

//--- METHODS
	void Set_Luminosity(double);
	void Train_BDT(TString, TString, bool);
	void Produce_Templates(TString);
	void Produce_Control_Histograms();
	void Draw_Control_Plots(TString);

//--- MEMBERS
	bool stop_program;

	private :
//--- METHODS

//--- MEMBERS
	TMVA::Reader *reader;

	std::vector<TString> sample_list;
	std::vector<TString> syst_list;
	std::vector<TString> channel_list;
	std::vector<TString> stringv_list;

	std::vector<TString> var_list; std::vector<float> var_list_floats; //TMVA variables
	std::vector<TString> v_add_var_names; vector<float> v_add_var_floats; //Additional vars only for CR plots

	// std::vector<TString> var_list_BDTcut; std::vector<float> BDTcut_floats; //List of variables used in BDTfake, on which we cut
	// std::vector<TString> v_cut_name; std::vector<TString> v_cut_def; std::vector<float> v_cut_float; std::vector<bool> v_cut_IsUsedForBDT; //Variables for region cuts (e.g. NJets, ...)

	std::vector<int> color_list;
	std::vector<TColor*> v_custom_colors;

	TString dir_ntuples, t_name;
	TString plot_extension;
	double luminosity_rescale;
	int nbin;
	TString filename_suffix;
	bool use_custom_colorPalette;
};

#endif
