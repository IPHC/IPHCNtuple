#include "../include/Pileup.h"

using namespace std;

bool use_my_own_PU_data_profile = false; //true <-> use my own PU for *data* ; false <-> use ttH file (different conventions)

//--------------------------------------------
// ##     ## ######## ##       ########  ######## ########
// ##     ## ##       ##       ##     ## ##       ##     ##
// ##     ## ##       ##       ##     ## ##       ##     ##
// ######### ######   ##       ########  ######   ########
// ##     ## ##       ##       ##        ##       ##   ##
// ##     ## ##       ##       ##        ##       ##    ##
// ##     ## ######## ######## ##        ######## ##     ##
//--------------------------------------------


//--------------------------------------------
// ########  ##     ##     ######  ##          ###     ######   ######
// ##     ## ##     ##    ##    ## ##         ## ##   ##    ## ##    ##
// ##     ## ##     ##    ##       ##        ##   ##  ##       ##
// ########  ##     ##    ##       ##       ##     ##  ######   ######
// ##        ##     ##    ##       ##       #########       ##       ##
// ##        ##     ##    ##    ## ##       ##     ## ##    ## ##    ##
// ##         #######      ######  ######## ##     ##  ######   ######
//--------------------------------------------


//Use stat function (from library sys/stat) to check if a file exists
bool Check_If_File_Exists(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
}


/**
 * Default constructor
 */
Pileup::Pileup()
{
    dir_parent = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/PU_SF/";

    output_dir = dir_parent + "profiles_MC/";

    if(!use_my_own_PU_data_profile)
    {
        fname_PU_Data_Nom  = dir_parent + "PileupData_ReRecoJSON_Full2017.root";
        fname_PU_Data_Up   = dir_parent + "PileupData_ReRecoJSON_Full2017.root";
        fname_PU_Data_Down = dir_parent + "PileupData_ReRecoJSON_Full2017.root";
        hname_PU_Data_Nom  = "pileup";
        hname_PU_Data_Up   = "pileup_plus";
        hname_PU_Data_Down = "pileup_minus";
    }
    else
    {
        fname_PU_Data_Nom  = dir_parent + "PileupNom.root";
        fname_PU_Data_Up   = dir_parent + "PileupUp.root";
        fname_PU_Data_Down = dir_parent + "PileupDown.root";
        hname_PU_Data_Nom  = "pileup";
        hname_PU_Data_Up   = "pileup";
        hname_PU_Data_Down = "pileup";
    }

    f_PU_Data_Nom  = NULL;
    f_PU_Data_Up   = NULL;
    f_PU_Data_Down = NULL;
    h_PU_Data_Nom  = NULL;
    h_PU_Data_Up   = NULL;
    h_PU_Data_Down = NULL;

	//Hard-code paths and number of FlatTree file, for each sample name
	// v_samples.push_back("tHq"); v_paths.push_back("THQ_4f_Hincl_13TeV_madgraph_pythia8/RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1_MINIAODSIM/180609_212538/"); v_nof_files.push_back(25);
}

/**
 * Default destructor
 */
Pileup::~Pileup()
{
//     delete h_PU_Data_Nom; h_PU_Data_Nom = NULL;
//     delete h_PU_Data_Up; h_PU_Data_Up = NULL;
//     delete h_PU_Data_Down; h_PU_Data_Down = NULL;
}


bool Pileup::Init(TString samplename)
{
    this->samplename = samplename;

    // Fill_MC_vector_from_HardCoded(); //Obsolete

    f_PU_Data_Nom = TFile::Open(fname_PU_Data_Nom);
    if(!f_PU_Data_Nom) {cout<<BOLD(FRED("Error : file "<<fname_PU_Data_Nom<<" not found !"))<<endl; return 0;}

    if(use_my_own_PU_data_profile) //Same histo name for all 3 variations (different files)
    {
        f_PU_Data_Up = TFile::Open(fname_PU_Data_Up);
        if(!f_PU_Data_Up) {cout<<BOLD(FRED("Error : file "<<fname_PU_Data_Up<<" not found !"))<<endl; return 0;}
        f_PU_Data_Down = TFile::Open(fname_PU_Data_Down);
        if(!f_PU_Data_Down) {cout<<BOLD(FRED("Error : file "<<fname_PU_Data_Down<<" not found !"))<<endl; return 0;}
    }
    else {f_PU_Data_Up = f_PU_Data_Nom; f_PU_Data_Down = f_PU_Data_Nom;} //Same file for all 3 variations (different histo names)

    h_PU_Data_Nom = (TH1F*) f_PU_Data_Nom->Get(hname_PU_Data_Nom);
    if(!h_PU_Data_Nom) {cout<<BOLD(FRED("Error : histo "<<hname_PU_Data_Nom<<" not found !"))<<endl; return 0;}
    h_PU_Data_Up = (TH1F*) f_PU_Data_Up->Get(hname_PU_Data_Up);
    if(!h_PU_Data_Up) {cout<<BOLD(FRED("Error : histo "<<hname_PU_Data_Up<<" not found !"))<<endl; return 0;}
    h_PU_Data_Down = (TH1F*) f_PU_Data_Down->Get(hname_PU_Data_Down);
    if(!h_PU_Data_Down) {cout<<BOLD(FRED("Error : histo "<<hname_PU_Data_Down<<" not found !"))<<endl; return 0;}

    return 1;
}


/**
 * Hard-code the PU MC distribution, taking the values from CMSSW, corresponding to production campaign
 * WARNING : in tHq2017 the samples we use have been produced with incorrect PU, can't use this simple method
 * Need to produce PU profile for each MC sample instead...
 * values taken here : https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py
*/
void Pileup::Fill_MC_vector_from_HardCoded()
{

    return;
}




 // ###### # #      #         #    #  ####     #    # ######  ####  #####  ####  #####
 // #      # #      #         ##  ## #    #    #    # #      #    #   #   #    # #    #
 // #####  # #      #         # ## # #         #    # #####  #        #   #    # #    #
 // #      # #      #         #    # #         #    # #      #        #   #    # #####
 // #      # #      #         #    # #    #     #  #  #      #    #   #   #    # #   #
 // #      # ###### ######    #    #  ####       ##   ######  ####    #    ####  #    #

bool Pileup::Fill_MC_vector_from_Pileup_Profile()
{
    v_MC_PU.resize(0);

    TString input_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/merged_histograms/"; //HARD-CODED
    TString filename = input_dir + samplename + ".root";

    if(!Check_If_File_Exists(filename) ) {cout<<BOLD(FRED("Error : MC PU profile file "<<filename<<" not found !"))<<endl; return 0;}

    TFile* f = TFile::Open(filename);
    if(!f->GetNkeys()) {cout<<BOLD(FRED("Error : file found, but has no keys !"))<<endl; return 0;}

    TH1F* h = (TH1F*) f->Get("hpileup");

    double integral = h->Integral();

    for(int ibin=1; ibin<100+1; ibin++) //start at bin 1
    {
      double content = 0;
      if(ibin == 1) {content = h->GetBinContent(0) + h->GetBinContent(1);} //I have underflow, put it in bin 1 (like data?)
      else {content = h->GetBinContent(ibin);}
      v_MC_PU.push_back(content / integral );
    }

    delete h; h = NULL;
    f->Close();

    cout<<endl<<"-> Pileup weightfile opened"<<endl<<endl<<endl;

    return 1;
}


 // #####  #    #    #    # ###### #  ####  #    # #####  ####
 // #    # #    #    #    # #      # #    # #    #   #   #
 // #    # #    #    #    # #####  # #      ######   #    ####
 // #####  #    #    # ## # #      # #  ### #    #   #        #
 // #      #    #    ##  ## #      # #    # #    #   #   #    #
 // #       ####     #    # ###### #  ####  #    #   #    ####


/**
 * Compute and store all PU weights in member vectors -- also write them in txt file
 */
void Pileup::Compute_PU_Weights()
{
    // cout<<BOLD(FYEL("--- Entering : Compute_PU_Weight()"))<<endl;

    //Write PU weights in output txt file, to cross check values (not used directly)
    ofstream textfile_PU_weights((output_dir + samplename + ".txt").Data());
    textfile_PU_weights<<"nPU \t weight \t weight_up \t weight_down"<<endl;

    if(!v_MC_PU.size()) {return;}

    double tot = 0, tot_up = 0, tot_down = 0; //For final normalisation
    tot = h_PU_Data_Nom->Integral();
    tot_up = h_PU_Data_Up->Integral();
    tot_down = h_PU_Data_Down->Integral();

    for(int i=0; i<v_MC_PU.size(); i++)
    {
        double pu_data_nom = h_PU_Data_Nom->GetBinContent(i+1); //start at bin 1
        // double pu_data_nom = h_PU_Data_Nom->GetBinContent(h_PU_Data_Nom->GetXaxis()->FindBin(i) );

        // cout<<endl<<"//--------------------------------------------"<<endl;
        // cout<<"Bin "<<i+1<<" :"<<endl;
        // cout<<"PU data = "<<pu_data_nom<<" / "<<tot<<" = "<<pu_data_nom/tot<<endl;
        // cout<<"v_MC_PU[i] = "<<v_MC_PU[i]<<endl;

        pu_data_nom/= tot;
        double weight = pu_data_nom / v_MC_PU[i]; //start at index 0
        if(isinf(weight) || isnan(weight) || weight < 0) {weight = 0;}

        // cout<<"=> weight = "<<weight<<endl<<endl;

        v_PU_weights_Nom.push_back(weight);

        double pu_data_up = h_PU_Data_Up->GetBinContent(h_PU_Data_Up->GetXaxis()->FindBin(i) );
        pu_data_up/= tot_up;
        double weight_up = pu_data_up / v_MC_PU[i];
        if(isinf(weight_up) || isnan(weight_up) || weight_up < 0) {weight_up = 0;}
        v_PU_weights_Up.push_back(weight_up);

        double pu_data_down = h_PU_Data_Down->GetBinContent(h_PU_Data_Down->GetXaxis()->FindBin(i) );
        pu_data_down/= tot_down;
        double weight_down = pu_data_down / v_MC_PU[i];
        if(isinf(weight_down) || isnan(weight_down) || weight_down < 0) {weight_down = 0;}
        v_PU_weights_Down.push_back(weight_down);

        textfile_PU_weights<<i<<"\t"<<weight<<"\t"<<weight_up<<"\t"<<weight_down<<endl;
    }

    return;
}






void Pileup::Print_Debug()
{

    //Print stored weights
    for(int i=0; i<v_PU_weights_Nom.size(); i++)
    {
        cout<<"//--------------------------------------------"<<endl;
        cout<<"v_PU_weights_Nom["<<i<<"] = "<<v_PU_weights_Nom[i]<<endl;
        cout<<"v_PU_weights_Up["<<i<<"] = "<<v_PU_weights_Up[i]<<endl;
        cout<<"v_PU_weights_Down["<<i<<"] = "<<v_PU_weights_Down[i]<<endl;
    }

}
