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


#include <iostream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>

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



/**
 * Executes bash command and returns the output as a Tstring
 * !! PROBLEMS WHEN PIPES INVOLVED (errno=ENOMEM, memory-related issue) --
 */
TString GetStdoutFromCommand(TString cmd_ts)
{
   string output = "";
   FILE* stream=0;
   const int max_buffer = 500; //Lack of buffer memory causes errors (errno=12)
   char buffer[max_buffer]; //Create buffer
   cmd_ts+= " 2>&1"; //Get both stdout and stderr outputs
   string cmd = cmd_ts.Data();

   stream = popen(cmd.c_str(), "r"); //Open read-only stream, run command

   if(stream)
   {
   	while(!feof(stream))
   	{
   		if(fgets(buffer, max_buffer, stream) != NULL) output.append(buffer); //Get the output
   	}

   	pclose(stream);
   }
   else //If stream was not opened <-> insufficient buffer memory. Retry with more !
   {
   	pclose(stream);

   	return "";
   }

   output.erase(std::remove(output.begin(), output.end(), '\n'), output.end()); //Remove the "end of line" character

   TString ts_output(output); //Convert string to TString

   return ts_output;
}

void Get_SumWeight_FlatTree(TString path_ft, int nfiles)
{
   cout<<FYEL("--- Will count Sum of Weighted Events (SWE) for FlatTree files in path : "<<path_ft<<" !")<<endl;
   // cout<<"(NB : make sure you did not do 'cmsenv', as it dis-activates the 'rfdir' command !)"<<endl;

   path_ft.Remove(TString::kTrailing, '/'); //Remove '/' at end of TString

   // TString command = "rfdir " + path_ft + " | wc -l";

   // ofstream file_out("list_tmp.txt");

   // int nfiles = Convert_TString_To_Number( GetStdoutFromCommand(command) );
	// if(nfiles == 0) {cout<<FRED("Error : no files found !")<<endl; return;}

   double sum_hweight = 0, sum_hcount = 0;

   for(int ifile=1; ifile<=nfiles; ifile++)
   {
      TString path_tmp = "root://sbgse1.in2p3.fr/" + path_ft + "/" + "output_" + Convert_Number_To_TString(ifile) + ".root";

      TFile* f =0;
	  f = TFile::Open(path_tmp);
	  if(!f) {continue;}

      TH1F* h = (TH1F*) f->Get("FlatTree/hweight");
	  sum_hweight+= h->GetBinContent(1);

	  h = (TH1F*) f->Get("FlatTree/hcount");
	  sum_hcount+= h->GetBinContent(1);

      f->Close();
   }

   // file_out<<path_ft<<" / nfiles = "<<nfiles<<" / sum_hweight = "<<setprecision(10)<<sum_hweight<<" / sum_hcount = "<<setprecision(10)<<sum_hcount<<endl<<endl;
	cout<<FMAG(" "<<path_ft<<" / nfiles = "<<nfiles<<" / sum_hweight = "<<setprecision(10)<<sum_hweight<<" ")<<endl<<endl;

   return;
}


int main()
{
   // TString path_ft = "/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/FlatTree/Walrus-patch2/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM/180122_094713/0000/";


	TString path_ft; int nfiles;

	cout<<FBLU("-- Enter path to Flat Tree files (e.g. /dpm/in2p3.fr/home/cms/phedex/store/....) :")<<endl;
	cin>>path_ft;
	cout<<FBLU("-- Enter total number of files to open :")<<endl;
	cin>>nfiles;


   Get_SumWeight_FlatTree(path_ft, nfiles);

  	return 0;
}
