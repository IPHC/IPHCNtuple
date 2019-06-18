/*
This code reads histograms from FlatTree files, merges and stores them
NB1 : more convenient to put it in the NtupleProducer/src dir., since we are reading the lists of FlatTree files stored in NTP/test
NB2 : output histograms are stored in NtupleAnalyzer/test/weights/... (more convenient, since they will be read by the NTA code)
*/

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
#include <TH1D.h>
#include <TH1D.h>
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






//--------------------------------------------
// ##     ## ######## ##       ########  ######## ########
// ##     ## ##       ##       ##     ## ##       ##     ##
// ##     ## ##       ##       ##     ## ##       ##     ##
// ######### ######   ##       ########  ######   ########
// ##     ## ##       ##       ##        ##       ##   ##
// ##     ## ##       ##       ##        ##       ##    ##
// ##     ## ######## ######## ##        ######## ##     ##
//--------------------------------------------



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






//--------------------------------------------
// ########  #### ##       ######## ##     ## ########
// ##     ##  ##  ##       ##       ##     ## ##     ##
// ##     ##  ##  ##       ##       ##     ## ##     ##
// ########   ##  ##       ######   ##     ## ########
// ##         ##  ##       ##       ##     ## ##
// ##         ##  ##       ##       ##     ## ##
// ##        #### ######## ########  #######  ##
//--------------------------------------------

/**
 * Readd all FlatTrees of process, and store PU histograms before any presel in output file
 */
void Get_Merged_PU_Histograms(vector<TString> v_process)
{
	cout<<"== WILL CREATE PILEUP PROFILES FOR ALL SAMPLES IN VECTOR (reading files lists) =="<<endl<<endl;

	//Should set to true, unless the PU histograms are missing from FlatTrees
	bool read_histograms_from_flatTrees = true; //true <-> read PU info from histogram ; //else produce histo from TTree (read all entries)

    TString list_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/lists/";
	TString output_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/PU_SF/profiles_MC/";

    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
		cout<<endl<<"-- Sample : "<<v_process[iproc]<<endl;

		TH1D* h_PU = 0;
		if(!read_histograms_from_flatTrees) {h_PU = new TH1D("", "", 100, 0, 100);}

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
				if(!read_histograms_from_flatTrees)
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
				if(read_histograms_from_flatTrees)
				{
                    TH1D* h_tmp = (TH1D*) f->Get("FlatTree/hpileup")->Clone();
					h_tmp->SetDirectory(0); //So that hist is not associated to file (& deleted when we close it)

                    if(!h_PU) {h_PU = (TH1D*) h_tmp;}
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

		if(!read_histograms_from_flatTrees) {delete h_PU; h_PU = NULL;}

		cout<<"---> Saved file "<<output_name<<endl;
    }

	return;
}





//--------------------------------------------
//  ######  ##     ## ##     ##    ##      ## ######## ####  ######   ##     ## ########  ######
// ##    ## ##     ## ###   ###    ##  ##  ## ##        ##  ##    ##  ##     ##    ##    ##    ##
// ##       ##     ## #### ####    ##  ##  ## ##        ##  ##        ##     ##    ##    ##
//  ######  ##     ## ## ### ##    ##  ##  ## ######    ##  ##   #### #########    ##     ######
//       ## ##     ## ##     ##    ##  ##  ## ##        ##  ##    ##  ##     ##    ##          ##
// ##    ## ##     ## ##     ##    ##  ##  ## ##        ##  ##    ##  ##     ##    ##    ##    ##
//  ######   #######  ##     ##     ###  ###  ######## ####  ######   ##     ##    ##     ######
//--------------------------------------------

/**
 *  * Readd all FlatTrees of process, and store Scale histograms before any presel in output file
 */
void Get_Merged_SumWeights_Histograms(vector<TString> v_process)
{
	//Should set to true, unless the SumWeights histograms are missing from FlatTrees
	bool read_histograms_from_flatTrees = true; //true <-> read info from histogram ; //else produce histo from TTree (read all entries)

    cout<<"== WILL CREATE SUM-OF-WEIGHT HISTOGRAMS FOR ALL SAMPLES IN VECTOR (reading files lists) =="<<endl<<endl;

    TString list_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/lists/";
	TString output_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/Scales/";

    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
		cout<<endl<<"-- Sample : "<<v_process[iproc]<<endl;

        TH1D* h_scales = 0;

        if(!read_histograms_from_flatTrees)
        {
            h_scales = new TH1D("", "", 10, 0, 10);
        }

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

                if(!read_histograms_from_flatTrees)
                {
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

                        h_scales->Fill(0.5, mc_weight);
                        h_scales->Fill(1+0.5, weight_originalXWGTUP);
                        h_scales->Fill(2+0.5, weight_scale_muR0p5);
                        h_scales->Fill(3+0.5,weight_scale_muF0p5);
                        h_scales->Fill(4+0.5, weight_scale_muR0p5muF0p5);
                        h_scales->Fill(5+0.5, weight_scale_muR2);
                        h_scales->Fill(6+0.5, weight_scale_muF2);
                        h_scales->Fill(7+0.5, weight_scale_muR2muF2);
                    }

    				delete t; t = NULL;
                }
                else
				{
                    TH1D* h_tmp = (TH1D*) f->Get("FlatTree/hscale")->Clone();
					h_tmp->SetDirectory(0); //So that hist is not associated to file (& deleted when we close it)

                    if(!h_scales) {h_scales = (TH1D*) h_tmp;}
					else {h_scales->Add(h_tmp);}
				}

				f->Close();

			} //loop on files
        } //loop on lists

		TString output_name = output_dir+v_process[iproc]+".root";
		TFile* f_output = 0;

		f_output = new TFile(output_name, "RECREATE");
		f_output->cd();

        h_scales->Write("h_sumWeightsScale");

		f_output->Close();

        if(!read_histograms_from_flatTrees) {delete h_scales;}

		cout<<"---> Saved file "<<output_name<<endl;
    } //process loop
	return;
}






//--------------------------------------------
// ##       ##     ## ########
// ##       ##     ## ##
// ##       ##     ## ##
// ##       ######### ######
// ##       ##     ## ##
// ##       ##     ## ##
// ######## ##     ## ########
//--------------------------------------------

/**
 * Read all FlatTrees of process, and store LHE sum weights histograms before any presel in output file
 * Only for THQ/THW samples ! (kT/kV reweighting)
 * NB : must have stored the 'hLHE' histogram in FlatTrees (didn't implement the re-reading of all FT entries for this func)
*/
void Get_Merged_SumWeightsLHE_Histograms(vector<TString> v_process)
{
    cout<<BOLD(FYEL("== WILL CREATE SUM-OF-WEIGHT HISTOGRAMS FOR ALL SAMPLES IN VECTOR (reading files lists) =="))<<endl<<endl;

    TString list_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/lists/";
	TString output_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/Scales/";

    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
		cout<<endl<<"-- Sample : "<<v_process[iproc]<<endl;

        TH1D* hLHE = 0;

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

                TH1D* h_tmp = (TH1D*) f->Get("FlatTree/hLHE")->Clone();
				h_tmp->SetDirectory(0); //So that hist is not associated to file (& deleted when we close it)

                if(!hLHE) {hLHE = (TH1D*) h_tmp;}
				else {hLHE->Add(h_tmp);}

				f->Close();
			} //loop on files
        } //loop on lists

		TString output_name = output_dir+v_process[iproc]+".root";
		TFile* f_output = 0;

		f_output = new TFile(output_name, "RECREATE");
		f_output->cd();

        hLHE->Write("hLHE");

		f_output->Close();

		cout<<"---> Saved file "<<output_name<<endl;
    } //process loop

	return;
}


















//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------



//--------------------------------------------
// ##     ## ######## ########   ######   ########       ###    ##       ##
// ###   ### ##       ##     ## ##    ##  ##            ## ##   ##       ##
// #### #### ##       ##     ## ##        ##           ##   ##  ##       ##
// ## ### ## ######   ########  ##   #### ######      ##     ## ##       ##
// ##     ## ##       ##   ##   ##    ##  ##          ######### ##       ##
// ##     ## ##       ##    ##  ##    ##  ##          ##     ## ##       ##
// ##     ## ######## ##     ##  ######   ########    ##     ## ######## ########

// ##     ## ####  ######  ########  #######   ######
// ##     ##  ##  ##    ##    ##    ##     ## ##    ##
// ##     ##  ##  ##          ##    ##     ## ##
// #########  ##   ######     ##    ##     ##  ######
// ##     ##  ##        ##    ##    ##     ##       ##
// ##     ##  ##  ##    ##    ##    ##     ## ##    ##
// ##     ## ####  ######     ##     #######   ######
//--------------------------------------------


/**
 * For each listed sample, merge the PU/SumWeights histograms & store them
 * NB : also merge sum of some LHE weights for THQ/THW samples
 */
void Get_Merged_Histograms_From_FlatTrees(vector<TString> v_process)
{
    cout<<endl<<endl<<endl<<BOLD(FYEL("== WILL MERGE HISTOGRAMS & STORE THEM =="))<<endl;

    TString list_dir = "./lists/";
	TString output_dir = "../../NtupleAnalyzer/test/weights_2017/merged_histograms/";
	cout<<"(NB : Sample names declared in main(), paths found in : "<<list_dir<<")"<<endl;
	mkdir(output_dir.Data(), 0777);

	//-- Loop on processes
    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
		cout<<endl<<UNDL("* Sample : "<<v_process[iproc]<<"")<<endl;

		//Output merged histograms
		TH1D* hpileup = 0;
		TH1D* hSumWeights = 0;
		TH1D* hLHE = 0;

		bool sample_found = false; //Verify that at least 1 file has been found for current sample
		int nFiles_notFound = 0;
		int ifile = 0;
		//--- Loop on lists of files
		while(nFiles_notFound < 5) //Stop if not finding 5 consecutive files (must have reached the end)
        {
            // if(ifile == 1) {break;} //debug

			//Read txt file containing list of FlatTree files to read
            TString filename = list_dir + v_process[iproc] + "_ID" + Convert_Number_To_TString(ifile) + ".txt";
			ifile++;
            if(!Check_File_Existence(filename) ) {cout<<"File "<<filename<<" not found !"<<endl; nFiles_notFound++; continue;}

			sample_found = true;

			ifstream file_in(filename.Data() );

			string line;
			TH1D *h_tmp1 = 0, *h_tmp2 = 0, *h_tmp3 = 0;
			//Loop on files to read
			while(getline(file_in, line) )
			{
				TString ts = line;
				cout<<endl<<"-- Opening file ["<<ts<<"] ..."<<endl;

				TFile* f = TFile::Open(ts);
				if(!f) {continue;}
                cout<<"... opened !"<<endl<<endl;

				//Get PU histo
				h_tmp1 = (TH1D*) f->Get("FlatTree/hpileup")->Clone();
				h_tmp1->SetDirectory(0); //So that hist is not associated to file (& deleted when we close it)
                if(!hpileup) {hpileup = (TH1D*) h_tmp1;}
				else {hpileup->Add(h_tmp1);}

				//Get sumWeights histo
				h_tmp2 = (TH1D*) f->Get("FlatTree/hSumWeights")->Clone();
				h_tmp2->SetDirectory(0); //So that hist is not associated to file (& deleted when we close it)
                if(!hSumWeights) {hSumWeights = (TH1D*) h_tmp2;}
				else {hSumWeights->Add(h_tmp2);}

				//Get LHE sum weights histo
				h_tmp3 = (TH1D*) f->Get("FlatTree/hLHE")->Clone();
				h_tmp3->SetDirectory(0); //So that hist is not associated to file (& deleted when we close it)
                if(!hLHE) {hLHE = (TH1D*) h_tmp3;}
				else {hLHE->Add(h_tmp3);}

				f->Close();
			} //end loop on files to read
        } //end loop on lists of files
		cout<<"(Must have reached the end of the loop, stop)"<<endl<<endl<<endl;

		if(!sample_found) {cout<<FRED("No file found for this sample ! Skip !")<<endl; continue;}

		TString output_name = output_dir+v_process[iproc]+".root";
		TFile* f_output = 0;

		f_output = new TFile(output_name, "RECREATE");
		f_output->cd();

		hpileup->Write("hpileup");
		hSumWeights->Write("hSumWeights");
		hLHE->Write("hLHE");

		f_output->Close();

		cout<<endl<<endl<<BOLD(FYEL("===> Saved file "<<output_name<<""))<<endl;
    } //end process loop

	return;
}




//--------------------------------------------
// ######## ########  ######  ######## #### ##    ##  ######
//    ##    ##       ##    ##    ##     ##  ###   ## ##    ##
//    ##    ##       ##          ##     ##  ####  ## ##
//    ##    ######    ######     ##     ##  ## ## ## ##   ####
//    ##    ##             ##    ##     ##  ##  #### ##    ##
//    ##    ##       ##    ##    ##     ##  ##   ### ##    ##
//    ##    ########  ######     ##    #### ##    ##  ######
//--------------------------------------------

/**
 * Can use this function to get the PU SF distribution, of PU distribution, before preselection, by running over all the FlatTree entries (~ no selection applied)
 * E.g. can store 1 standard PU distribution, and 1 PU distribution after the PU SF have been applied
 * Then, can use the output file to compare the Data/MC PU profiles before/after correction
 * NB : could automatically read the PU SF files... or, more simply, can just copy/hardcode here the PU SF vector for the process of interest
 */
void Get_PU_Distribution_before_presel(vector<TString> v_process)
{
    cout<<endl<<endl<<endl<<BOLD(FYEL("== WILL MERGE HISTOGRAMS & STORE THEM =="))<<endl;

    TString list_dir = "./lists/";
	TString output_dir = "../../NtupleAnalyzer/test/weights_2017/TEST/";
	cout<<"(NB : Sample names declared in main(), paths found in : "<<list_dir<<")"<<endl;
	mkdir(output_dir.Data(), 0777);

	bool debug = false; //read only 1 file
    int nrootfiles = 0; //counter to limit nof files

    vector<float> v_PU_SF; //HARDCODED PU_SF VALUES (AS FUNC OF NPU)

//ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8
    // v_PU_SF.push_back(0.000327303);
    // v_PU_SF.push_back(0.0412878);
    // v_PU_SF.push_back(0.0515754);
    // v_PU_SF.push_back(0.0640951);
    // v_PU_SF.push_back(0.10026);
    // v_PU_SF.push_back(0.126995);
    // v_PU_SF.push_back(0.148758);
    // v_PU_SF.push_back(0.198317);
    // v_PU_SF.push_back(0.107221);
    // v_PU_SF.push_back(0.355577);
    // v_PU_SF.push_back(0.461626);
    // v_PU_SF.push_back(0.654048);
    // v_PU_SF.push_back(0.710741);
    // v_PU_SF.push_back(0.757366);
    // v_PU_SF.push_back(0.753194);
    // v_PU_SF.push_back(0.834346);
    // v_PU_SF.push_back(0.927497);
    // v_PU_SF.push_back(0.996304);
    // v_PU_SF.push_back(1.10439);
    // v_PU_SF.push_back(1.12998);
    // v_PU_SF.push_back(1.22651);
    // v_PU_SF.push_back(1.27871);
    // v_PU_SF.push_back(1.32772);
    // v_PU_SF.push_back(1.32046);
    // v_PU_SF.push_back(1.37713);
    // v_PU_SF.push_back(1.36309);
    // v_PU_SF.push_back(1.34702);
    // v_PU_SF.push_back(1.35931);
    // v_PU_SF.push_back(1.35884);
    // v_PU_SF.push_back(1.35021);
    // v_PU_SF.push_back(1.3154);
    // v_PU_SF.push_back(1.24142);
    // v_PU_SF.push_back(1.19048);
    // v_PU_SF.push_back(1.11302);
    // v_PU_SF.push_back(1.03888);
    // v_PU_SF.push_back(0.994607);
    // v_PU_SF.push_back(0.935129);
    // v_PU_SF.push_back(0.886548);
    // v_PU_SF.push_back(0.849558);
    // v_PU_SF.push_back(0.813331);
    // v_PU_SF.push_back(0.811142);
    // v_PU_SF.push_back(0.87841);
    // v_PU_SF.push_back(0.962654);
    // v_PU_SF.push_back(1.041);
    // v_PU_SF.push_back(1.21077);
    // v_PU_SF.push_back(1.39067);
    // v_PU_SF.push_back(1.56733);
    // v_PU_SF.push_back(1.60314);
    // v_PU_SF.push_back(1.61175);
    // v_PU_SF.push_back(1.54243);
    // v_PU_SF.push_back(1.37284);
    // v_PU_SF.push_back(1.13766);
    // v_PU_SF.push_back(0.93614);
    // v_PU_SF.push_back(0.714744);
    // v_PU_SF.push_back(0.558522);
    // v_PU_SF.push_back(0.38709);
    // v_PU_SF.push_back(0.275052);
    // v_PU_SF.push_back(0.183913);
    // v_PU_SF.push_back(0.126931);
    // v_PU_SF.push_back(0.089078);
    // v_PU_SF.push_back(0.0618938);
    // v_PU_SF.push_back(0.0460299);
    // v_PU_SF.push_back(0.034228);
    // v_PU_SF.push_back(0.026902);
    // v_PU_SF.push_back(0.0208385);
    // v_PU_SF.push_back(0.0136676);
    // v_PU_SF.push_back(0.0111041);
    // v_PU_SF.push_back(0.00975514);
    // v_PU_SF.push_back(0.00820767);
    // v_PU_SF.push_back(0.00772725);
    // v_PU_SF.push_back(0.00668991);
    // v_PU_SF.push_back(0.00622559);
    // v_PU_SF.push_back(0.00380138);
    // v_PU_SF.push_back(0.00422235);
    // v_PU_SF.push_back(0.00253227);
    // v_PU_SF.push_back(0.00355532);
    // v_PU_SF.push_back(0.00303534);
    // v_PU_SF.push_back(0.00249638);
    // v_PU_SF.push_back(0.00254006);
    // v_PU_SF.push_back(0.000619584);
    // v_PU_SF.push_back(0.0022375);

//THQ_CTCVCP
	v_PU_SF.push_back(0.192241);
	v_PU_SF.push_back(3.67772);
	v_PU_SF.push_back(3.09023);
	v_PU_SF.push_back(2.22561);
	v_PU_SF.push_back(1.68218);
	v_PU_SF.push_back(1.55904);
	v_PU_SF.push_back(1.3062);
	v_PU_SF.push_back(1.28378);
	v_PU_SF.push_back(0.645943);
	v_PU_SF.push_back(1.52899);
	v_PU_SF.push_back(1.57005);
	v_PU_SF.push_back(1.55281);
	v_PU_SF.push_back(1.37511);
	v_PU_SF.push_back(1.21077);
	v_PU_SF.push_back(1.12449);
	v_PU_SF.push_back(1.10043);
	v_PU_SF.push_back(1.12409);
	v_PU_SF.push_back(1.1625);
	v_PU_SF.push_back(1.19364);
	v_PU_SF.push_back(1.20975);
	v_PU_SF.push_back(1.23154);
	v_PU_SF.push_back(1.25644);
	v_PU_SF.push_back(1.27511);
	v_PU_SF.push_back(1.28179);
	v_PU_SF.push_back(1.28754);
	v_PU_SF.push_back(1.27923);
	v_PU_SF.push_back(1.28099);
	v_PU_SF.push_back(1.27385);
	v_PU_SF.push_back(1.27918);
	v_PU_SF.push_back(1.25793);
	v_PU_SF.push_back(1.21859);
	v_PU_SF.push_back(1.16258);
	v_PU_SF.push_back(1.10516);
	v_PU_SF.push_back(1.03015);
	v_PU_SF.push_back(0.960988);
	v_PU_SF.push_back(0.905754);
	v_PU_SF.push_back(0.858262);
	v_PU_SF.push_back(0.827215);
	v_PU_SF.push_back(0.779085);
	v_PU_SF.push_back(0.746428);
	v_PU_SF.push_back(0.755271);
	v_PU_SF.push_back(0.794897);
	v_PU_SF.push_back(0.864377);
	v_PU_SF.push_back(0.973694);
	v_PU_SF.push_back(1.11036);
	v_PU_SF.push_back(1.28123);
	v_PU_SF.push_back(1.43428);
	v_PU_SF.push_back(1.50826);
	v_PU_SF.push_back(1.51495);
	v_PU_SF.push_back(1.43599);
	v_PU_SF.push_back(1.28821);
	v_PU_SF.push_back(1.1042);
	v_PU_SF.push_back(0.891568);
	v_PU_SF.push_back(0.692084);
	v_PU_SF.push_back(0.520517);
	v_PU_SF.push_back(0.36928);
	v_PU_SF.push_back(0.256988);
	v_PU_SF.push_back(0.174364);
	v_PU_SF.push_back(0.120457);
	v_PU_SF.push_back(0.0829975);
	v_PU_SF.push_back(0.0592459);
	v_PU_SF.push_back(0.0434596);
	v_PU_SF.push_back(0.0325716);
	v_PU_SF.push_back(0.0252627);
	v_PU_SF.push_back(0.0203882);
	v_PU_SF.push_back(0.014399);
	v_PU_SF.push_back(0.0104499);
	v_PU_SF.push_back(0.00908608);
	v_PU_SF.push_back(0.00803638);
	v_PU_SF.push_back(0.00742465);
	v_PU_SF.push_back(0.00702053);
	v_PU_SF.push_back(0.00674887);
	v_PU_SF.push_back(0.00673933);
	v_PU_SF.push_back(0.00680997);
	v_PU_SF.push_back(0.0053501);
	v_PU_SF.push_back(0.00442728);
	v_PU_SF.push_back(0.00468406);
	v_PU_SF.push_back(0.00469889);
	v_PU_SF.push_back(0.00479854);
	v_PU_SF.push_back(0.00477649);
	v_PU_SF.push_back(0.00532516);
	v_PU_SF.push_back(0.00603418);
	v_PU_SF.push_back(0.00512503);
	v_PU_SF.push_back(0.00360379);
	v_PU_SF.push_back(0.00426874);
	v_PU_SF.push_back(0.00345042);
	v_PU_SF.push_back(0.00344874);
	v_PU_SF.push_back(0.00257081);
	v_PU_SF.push_back(0.0031996);
	v_PU_SF.push_back(0.00225517);
	v_PU_SF.push_back(0.0032154);
	v_PU_SF.push_back(0.000973153);
	v_PU_SF.push_back(0.00175021);
	v_PU_SF.push_back(0.00103902);
	v_PU_SF.push_back(0.00122142);
	v_PU_SF.push_back(0.000355365);
	v_PU_SF.push_back(0);
	v_PU_SF.push_back(0);
	v_PU_SF.push_back(0);
	v_PU_SF.push_back(0);


	//-- Loop on processes
    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
		cout<<endl<<UNDL("* Sample : "<<v_process[iproc]<<"")<<endl;

		//Outputs
		TString output_name = output_dir + v_process[iproc] + ".root";
		TFile* f_output = new TFile(output_name, "RECREATE");
		TH1D* h_PU_noCorr = new TH1D("", "", 100, 0, 100);
		TH1D* h_PU_withCorr = new TH1D("", "", 100, 0, 100);

		//Output merged histograms
		TTree* t = 0;

		bool sample_found = false; //Verify that at least 1 file has been found for current sample
		int nFiles_notFound = 0;
		int ifile = 0;
		//--- Loop on lists of files
		while(nFiles_notFound < 5) //Stop if not finding 5 consecutive files (must have reached the end)
        {
			//Read txt file containing list of FlatTree files to read
            TString filename = list_dir + v_process[iproc] + "_ID" + Convert_Number_To_TString(ifile) + ".txt";
			ifile++;
            if(!Check_File_Existence(filename) ) {cout<<"File "<<filename<<" not found !"<<endl; nFiles_notFound++; continue;}

			sample_found = true;

			ifstream file_in(filename.Data() );

			string line;
			//Loop on files to read
			while(getline(file_in, line) )
			{
                // if(nrootfiles>40) {debug = true;}

				TString ts = line;
				cout<<endl<<"-- Opening file ["<<ts<<"] ..."<<endl;

                nrootfiles++;

				TFile* f = TFile::Open(ts);
				if(!f) {continue;}
                cout<<"... opened !"<<endl<<endl;

				//Get PU histo
                Int_t nPU = 0;
                Float_t mc_weight = 0;
				t = (TTree*) f->Get("FlatTree/tree");
				t->SetBranchStatus("*", 0);
				t->SetBranchStatus("mc_pu_trueNumInt", 1);
                t->SetBranchAddress("mc_pu_trueNumInt", &nPU);
				t->SetBranchStatus("mc_weight_originalValue", 1);
                t->SetBranchAddress("mc_weight_originalValue", &mc_weight);

				int nentries = t->GetEntries();
				// int nentries = 100;

                for(int ientry=0; ientry<nentries; ientry++)
                {
                    nPU = 0; mc_weight = 1;

                    t->GetEntry(ientry);

					// cout<<"mc_weight "<<mc_weight<<endl;
					// cout<<"nPU "<<nPU<<endl;
					// cout<<"v_PU_SF[nPU] "<<v_PU_SF[nPU]<<endl;

                    // if(nPU>=v_PU_SF.size()) {PU_SF = 0;}
                    // else {PU_SF = v_PU_SF[nPU];}

					h_PU_noCorr->Fill(nPU, mc_weight);
					h_PU_withCorr->Fill(nPU, mc_weight*v_PU_SF[nPU]);
				}

				f->Close();

				if(debug) {break;}
			} //end loop on files to read

			if(debug) {break;}
        } //end loop on lists of files
		cout<<"(Must have reached the end of the loop, stop)"<<endl<<endl<<endl;

		if(!sample_found) {cout<<FRED("No file found for this sample ! Skip !")<<endl; continue;}

		f_output->cd();

		h_PU_noCorr->Write("h_PU_noCorr");
		h_PU_withCorr->Write("h_PU_withCorr");

		// TCanvas* c = new TCanvas("", "", 1000, 800);
		// h_PU_noCorr->Draw("hist");
		// c->SaveAs("./test.png");

		delete h_PU_noCorr;
		delete h_PU_withCorr;

		f_output->Close();
		// delete c;

		cout<<endl<<endl<<BOLD(FYEL("===> Saved file "<<output_name<<""))<<endl;
    } //end process loop


	return;
}






/**
 * Read all FlatTree files for given process, and sum weights for each element of the LHE vector => Obtain sum of weight for each LHE weight, before preselection
 * Rather long (~20min for THQ_ctcvcp)
 * Idea to speed up process : use "Draw()"" + "TH1 *myh = (TH1*)gPad->GetPrimitive("htemp");" and then do mean*nentries for each index of the LHE vector? (NB : must also divide by generator weights). Problem : this is imperfect : depends on binning, etc. NB : actually, slower !!
 * Other possibilities : directly store these sums of weight in FlatTreeProducer.cc ; or read mapping LHE_ID <-> PDF_ID at beginning (must produce it beforehand), and then only store weights that we care about !
 */
void Get_SumOfWeights_LHEweights(vector<TString> v_process)
{
	// bool use_Draw_histograms = false; //Disactivate, using Draw() is too slow !

	cout<<"== WILL MERGE SUMS OF WEIGHTS, for all LHE weights in sample =="<<endl<<endl;

    TString list_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/lists/";
	TString output_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/LHE/";
	mkdir(output_dir.Data(), 0777);

    for(int iproc=0; iproc<v_process.size(); iproc++)
    {
		cout<<endl<<"-- Sample : "<<v_process[iproc]<<endl;

		TH1D* hLHE = 0;
		bool histo_already_created = false; //Set to true once we know the binning (= nof LHE weights in sample)
		int nbins = 0;

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

				TTree* t = (TTree*) f->Get("FlatTree/tree");
				t->SetBranchStatus("*", 0); //disable all branches, speed up
				t->SetBranchStatus("mc_pdfweights", 1);
				t->SetBranchStatus("mc_weight_originalValue", 1);
				t->SetBranchStatus("weight_originalXWGTUP", 1);

				vector<float>* mc_pdfweights;
				Float_t weight_originalXWGTUP = 0, mc_weight_originalValue = 0;
				// if(!use_Draw_histograms)
				{
					mc_pdfweights = new vector<float>;
					t->SetBranchAddress("mc_pdfweights", &mc_pdfweights);
					t->SetBranchAddress("mc_weight_originalValue", &mc_weight_originalValue);
					t->SetBranchAddress("weight_originalXWGTUP", &weight_originalXWGTUP);

					int nentries = t->GetEntries();
					for(int ientry=0; ientry<nentries; ientry++)
					{
						// if(ientry>10) {break;}

						mc_pdfweights->clear();
						weight_originalXWGTUP = 0;
						mc_weight_originalValue = 0;

						if(ientry%20000==0) {cout<<"Event : "<<ientry<<" / "<<nentries<<endl;}

						t->GetEntry(ientry);

						if(!histo_already_created)
						{
							hLHE = new TH1D("hLHE", "hLHE", mc_pdfweights->size(), 0, mc_pdfweights->size()); //bin width = 1
							hLHE->SetDirectory(0); //NECESSARY so that histo is not associated with TFile, and doesn't get deleted when file closed !
							histo_already_created = true;
						}

						// cout<<"mc_pdfweights->size() = "<<mc_pdfweights->size()<<endl;
						for(int iweight=0; iweight<mc_pdfweights->size(); iweight++)
						{
							// cout<<"iweight = "<<iweight<<" / mc_pdfweights->at(iweight) = "<<mc_pdfweights->at(iweight);
							// cout<<" / mc_weight_originalValue = "<<mc_weight_originalValue<<" / weight_originalXWGTUP = "<<weight_originalXWGTUP<<endl;
							// cout<<"--> weight = "<<mc_pdfweights->at(iweight)*mc_weight_originalValue/weight_originalXWGTUP<<endl;

							// hLHE->Fill(iweight+0.5, mc_pdfweights->at(iweight)); //bin width = 1
							hLHE->AddBinContent(iweight, mc_pdfweights->at(iweight)*mc_weight_originalValue/weight_originalXWGTUP); //bin width = 1 ; avoid to fill millions of entries ?
						}
					}

					delete t; t = NULL;
					delete mc_pdfweights;
				}
				/*
				else
				{
					TH1D* h_tmp = 0;

					if(!histo_already_created)
					{
						// t->Draw("mc_pdfweights");
						// h_tmp = (TH1D*) gPad->GetPrimitive("htemp"); //Retrieve histo in memory
						// // h_tmp->SetDirectory(0); //NECESSARY so that histo is not associated with TFile, and doesn't get deleted when file closed !
						// nbins = h_tmp->GetNbinsX();
						// cout<<"nbins = "<<nbins<<endl;
						// delete h_tmp; h_tmp = 0;

						nbins = 951; //FIXME -- hardcoded for testing, should get binning in some way

						hLHE = new TH1D("hLHE", "hLHE", nbins, 0, nbins); //bin width = 1
						hLHE->SetDirectory(0); //NECESSARY so that histo is not associated with TFile, and doesn't get deleted when file closed !
						histo_already_created = true;
					}

					for(int iweight=0; iweight<nbins; iweight++)
					{
						if(iweight>0) {break;} //FIXME

						TString var_tmp = "mc_pdfweights[" + Convert_Number_To_TString(iweight) + "]*mc_weight_originalValue/weight_originalXWGTUP";
						t->Draw(var_tmp); //Draw histo in memory
						h_tmp = (TH1D*) gPad->GetPrimitive("htemp"); //Retrieve histo in memory
						// h_tmp->SetDirectory(0); //NECESSARY so that histo is not associated with TFile, and doesn't get deleted when file closed !
						double sum_tmp = h_tmp->GetMean() * h_tmp->GetEntries(); //Approx. estimate of sum for given index
						cout<<"sum_tmp = "<<sum_tmp<<endl;

						hLHE->AddBinContent(iweight+0.5, sum_tmp); //bin width = 1 ; avoid to fill millions of entries ?
						delete h_tmp; h_tmp = 0;
					}
				}*/

				f->Close();
			}
        }

		TString output_name = output_dir+v_process[iproc]+".root";
		TFile* f_output = 0;
		f_output = new TFile(output_name, "RECREATE");
		f_output->cd();
		hLHE->Write("hLHE");
		f_output->Close();

		delete hLHE; hLHE = NULL;

		cout<<"---> Saved file "<<output_name<<endl;
    } //sample loop

	return;
}
















//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------


//--------------------------------------------
// ##     ##    ###    #### ##    ##
// ###   ###   ## ##    ##  ###   ##
// #### ####  ##   ##   ##  ####  ##
// ## ### ## ##     ##  ##  ## ## ##
// ##     ## #########  ##  ##  ####
// ##     ## ##     ##  ##  ##   ###
// ##     ## ##     ## #### ##    ##
//--------------------------------------------

int main(int argc, char *argv[])
{
	// if(argc != 2) {cout<<"Wrong arg ! Abort"<<endl; return 1;}
	// v_process.push_back(argv[1]);

//--------------------------------------------
//--------------------------------------------
//Un-comment all the samples for which you want to get/merge/store histograms
    vector<TString> v_process;
	// v_process.push_back("THQ_ctcvcp_4f_Hincl_13TeV_madgraph_pythia8");
	// v_process.push_back("THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8");
    // v_process.push_back("TTH_4f_ctcvcp_TuneCP5_13TeV_madgraph_pythia8");
	// v_process.push_back("ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8");
	// v_process.push_back("ttH_M125_TuneCP5_13TeV-powheg-pythia8");
    // v_process.push_back("TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8");
    // v_process.push_back("ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8");
    // v_process.push_back("TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8");
    // v_process.push_back("WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8");
	// v_process.push_back("tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8");
	// v_process.push_back("ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8");
    // v_process.push_back("ZZTo4L_13TeV_powheg_pythia8");
    // v_process.push_back("TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8");
	// v_process.push_back("TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8");
    // v_process.push_back("WZZ_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("ZZZ_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8");
	// v_process.push_back("WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("TTTT_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("TTWW_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TTWH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TTTW_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8");
    // v_process.push_back("WpWpJJ_EWK-QCD_TuneCP5_13TeV-madgraph-pythia8");
	// v_process.push_back("WZG_TuneCP5_13TeV-amcatnlo-pythia8");
    // v_process.push_back("VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8");
	// v_process.push_back("WW_DoubleScattering_13TeV-pythia8_TuneCP5");

    //GammaConv, DY, MC Fakes
	// v_process.push_back("W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8");
    // v_process.push_back("TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8");
    // v_process.push_back("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8");
    // v_process.push_back("TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8");
    // v_process.push_back("TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8");
    v_process.push_back("TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2_PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_v1_MINIAODSIM");
    // v_process.push_back("TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8");
	// v_process.push_back("DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8");
    // v_process.push_back("DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8");

    // v_process.push_back("ST_FCNC-TH_Tleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8");
    // v_process.push_back("ST_FCNC-TH_Tleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8");
    // v_process.push_back("TT_FCNC-TtoHJ_aTleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8");
    // v_process.push_back("TT_FCNC-TtoHJ_aTleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8");
    // v_process.push_back("TT_FCNC-aTtoHJ_Tleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8");
    // v_process.push_back("TT_FCNC-aTtoHJ_Tleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8");

    // v_process.push_back("ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8");
    // v_process.push_back("ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8");
    // v_process.push_back("ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8");

//--------------------------------------------
    //--- Obsolete ?
	// v_process.push_back("ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8"); --
    // v_process.push_back("TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"); --
    // v_process.push_back("WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8"); --
    // v_process.push_back("WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"); --
    // v_process.push_back("ST_tWll_5f_LO_TuneCP5_PSweights_13TeV_madgraph_pythia8_Fall17"); --
    // v_process.push_back("ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"); --
    // v_process.push_back("ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"); --
    // v_process.push_back("VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_Fall17"); --

	// v_process.push_back("DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_v1_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");

	// v_process.push_back("DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM");
	// v_process.push_back("DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
	// v_process.push_back("DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM");
	// v_process.push_back("DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
	// v_process.push_back("DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v2_MINIAODSIM");
	// v_process.push_back("DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v11_ext1_v1_MINIAODSIM");
	// v_process.push_back("DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8");
	// v_process.push_back("DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v2_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM");
	// v_process.push_back("DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8");

    // v_process.push_back("TprimeBToTH_M-1000_LH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TprimeBToTH_M-1100_LH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TprimeBToTH_M-1200_LH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TprimeBToTH_M-600_LH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TprimeBToTH_M-650_LH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TprimeBToTH_M-700_LH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TprimeBToTH_M-800_LH_TuneCP5_13TeV-madgraph-pythia8");
    // v_process.push_back("TprimeBToTH_M-900_LH_TuneCP5_13TeV-madgraph-pythia8");


//--------------------------------------------
//--------------------------------------------

//--------------------------------------------
//These function are only read/merge/store 1 kind of histogram each -- NOT NEEDED ANYMORE (use Get_Merged_Histograms_From_FlatTrees)
    // Get_Merged_PU_Histograms(v_process);
    // Get_Merged_SumWeights_Histograms(v_process);
	// Get_Merged_SumWeightsLHE_Histograms(v_process);

//Get PU_SF distri for ttH before presel -- for xchecking
    // Get_PU_Distribution_before_presel(v_process);

//This function reads/merges/stores several kinds of histos at once -- NEEDED TO COMPUTE EVENT WEIGHTS, ...
	Get_Merged_Histograms_From_FlatTrees(v_process);

//Get sums of weights for *all* LHE indices -- needed for PDF sets, scales, etc.
    // Get_SumOfWeights_LHEweights(v_process);

//--------------------------------------------
	return 0;
}
