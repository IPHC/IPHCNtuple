#include "../include/Pileup.h"

using namespace std;

bool use_my_own_PU_data_profile = false; //true <-> use my own PU for data ; false <-> use ttH file (different conventions) //NB : using tth for now, deleted my own data profile

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
    // v_MC_PU.push_back(3.39597497605e-05);
    // v_MC_PU.push_back(6.63688402133e-06);
    // v_MC_PU.push_back(1.39533611284e-05);
    // v_MC_PU.push_back(3.64963078209e-05);
    // v_MC_PU.push_back(6.00872171664e-05);
    // v_MC_PU.push_back(9.33932578027e-05);
    // v_MC_PU.push_back(0.000120591524486);
    // v_MC_PU.push_back(0.000128694546198);
    // v_MC_PU.push_back(0.000361697233219);
    // v_MC_PU.push_back(0.000361796847553);
    // v_MC_PU.push_back(0.000702474896113);
    // v_MC_PU.push_back(0.00133766053707);
    // v_MC_PU.push_back(0.00237817050805);
    // v_MC_PU.push_back(0.00389825605651);
    // v_MC_PU.push_back(0.00594546732588);
    // v_MC_PU.push_back(0.00856825906255);
    // v_MC_PU.push_back(0.0116627396044);
    // v_MC_PU.push_back(0.0148793350787);
    // v_MC_PU.push_back(0.0179897368379);
    // v_MC_PU.push_back(0.0208723871946);
    // v_MC_PU.push_back(0.0232564170641);
    // v_MC_PU.push_back(0.0249826433945);
    // v_MC_PU.push_back(0.0262245860346);
    // v_MC_PU.push_back(0.0272704617569);
    // v_MC_PU.push_back(0.0283301107549);
    // v_MC_PU.push_back(0.0294006137386);
    // v_MC_PU.push_back(0.0303026836965);
    // v_MC_PU.push_back(0.0309692426278);
    // v_MC_PU.push_back(0.0308818046328);
    // v_MC_PU.push_back(0.0310566806228);
    // v_MC_PU.push_back(0.0309692426278);
    // v_MC_PU.push_back(0.0310566806228);
    // v_MC_PU.push_back(0.0310566806228);
    // v_MC_PU.push_back(0.0310566806228);
    // v_MC_PU.push_back(0.0307696426944);
    // v_MC_PU.push_back(0.0300103336052);
    // v_MC_PU.push_back(0.0288355370103);
    // v_MC_PU.push_back(0.0273233309106);
    // v_MC_PU.push_back(0.0264343533951);
    // v_MC_PU.push_back(0.0255453758796);
    // v_MC_PU.push_back(0.0235877272306);
    // v_MC_PU.push_back(0.0215627588047);
    // v_MC_PU.push_back(0.0195825559393);
    // v_MC_PU.push_back(0.0177296309658);
    // v_MC_PU.push_back(0.0160560731931);
    // v_MC_PU.push_back(0.0146022004183);
    // v_MC_PU.push_back(0.0134080690078);
    // v_MC_PU.push_back(0.0129586991411);
    // v_MC_PU.push_back(0.0125093292745);
    // v_MC_PU.push_back(0.0124360740539);
    // v_MC_PU.push_back(0.0123547104433);
    // v_MC_PU.push_back(0.0123953922486);
    // v_MC_PU.push_back(0.0124360740539);
    // v_MC_PU.push_back(0.0124360740539);
    // v_MC_PU.push_back(0.0123547104433);
    // v_MC_PU.push_back(0.0124360740539);
    // v_MC_PU.push_back(0.0123387597772);
    // v_MC_PU.push_back(0.0122414455005);
    // v_MC_PU.push_back(0.011705203844);
    // v_MC_PU.push_back(0.0108187105305);
    // v_MC_PU.push_back(0.00963985508986);
    // v_MC_PU.push_back(0.00827210065136);
    // v_MC_PU.push_back(0.00683770076341);
    // v_MC_PU.push_back(0.00545237697118);
    // v_MC_PU.push_back(0.00420456901556);
    // v_MC_PU.push_back(0.00367513566191);
    // v_MC_PU.push_back(0.00314570230825);
    // v_MC_PU.push_back(0.0022917978982);
    // v_MC_PU.push_back(0.00163221454973);
    // v_MC_PU.push_back(0.00114065309494);
    // v_MC_PU.push_back(0.000784838366118);
    // v_MC_PU.push_back(0.000533204105387);
    // v_MC_PU.push_back(0.000358474034915);
    // v_MC_PU.push_back(0.000238881117601);
    // v_MC_PU.push_back(0.0001984254989);
    // v_MC_PU.push_back(0.000157969880198);
    // v_MC_PU.push_back(0.00010375646169);
    // v_MC_PU.push_back(6.77366175538e-05);
    // v_MC_PU.push_back(4.39850477645e-05);
    // v_MC_PU.push_back(2.84298066026e-05);
    // v_MC_PU.push_back(1.83041729561e-05);
    // v_MC_PU.push_back(1.17473542058e-05);
    // v_MC_PU.push_back(7.51982735129e-06);
    // v_MC_PU.push_back(6.16160108867e-06);
    // v_MC_PU.push_back(4.80337482605e-06);
    // v_MC_PU.push_back(3.06235473369e-06);
    // v_MC_PU.push_back(1.94863396999e-06);
    // v_MC_PU.push_back(1.23726800704e-06);
    // v_MC_PU.push_back(7.83538083774e-07);
    // v_MC_PU.push_back(4.94602064224e-07);
    // v_MC_PU.push_back(3.10989480331e-07);
    // v_MC_PU.push_back(1.94628487765e-07);
    // v_MC_PU.push_back(1.57888581037e-07);
    // v_MC_PU.push_back(1.2114867431e-07);
    // v_MC_PU.push_back(7.49518929908e-08);
    // v_MC_PU.push_back(4.6060444984e-08);
    // v_MC_PU.push_back(2.81008884326e-08);
    // v_MC_PU.push_back(1.70121486128e-08);
    // v_MC_PU.push_back(1.02159894812e-08);
}


bool Pileup::Fill_MC_vector_from_Pileup_Profile()
{
    v_MC_PU.resize(0);

    TString output_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/PU_SF/profiles_MC/";
    TString filename = output_dir + samplename + ".root";

    if(!Check_If_File_Exists(filename) ) {cout<<BOLD(FRED("Error : MC PU profile file "<<filename<<" not found !"))<<endl; return 0;}

    TFile* f = TFile::Open(filename);
    if(!f->GetNkeys()) {cout<<BOLD(FRED("Error : file found, but has no keys !"))<<endl; return 0;}

    TH1F* h = (TH1F*) f->Get("Pileup_MC");

    double integral = h->Integral();

    for(int ibin=0; ibin<100; ibin++)
    {
        v_MC_PU.push_back(h->GetBinContent(ibin) / integral );
    }

    delete h; h = NULL;
    f->Close();

    return 1;
}


/**
 * Compute and store all PU weights in member vectors
 */
void Pileup::Compute_PU_Weights()
{
    // cout<<BOLD(FYEL("--- Entering : Compute_PU_Weight()"))<<endl;

    if(!v_MC_PU.size()) {return;}

    double tot = 0, tot_up = 0, tot_down = 0; //For final normalisation
    tot = h_PU_Data_Nom->Integral();
    tot_up = h_PU_Data_Up->Integral();
    tot_down = h_PU_Data_Down->Integral();

    for(int i=0; i<v_MC_PU.size(); i++)
    {
        double pu_data_nom = h_PU_Data_Nom->GetBinContent(h_PU_Data_Nom->GetXaxis()->FindBin(i) );
        // cout<<"pu_data_nom = "<<pu_data_nom<<" / "<<tot<<" = "<<pu_data_nom/tot<<endl;

        pu_data_nom/= tot;
        double weight = pu_data_nom / v_MC_PU[i];
        if(isinf(weight) || isnan(weight) || weight < 0) {weight = 0;} //FIXME

        // cout<<"v_MC_PU[i] = "<<v_MC_PU[i]<<endl;
        // cout<<"=> weight = "<<weight<<endl;

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
    }

    return;
}


/**
 * If PU weights have already been written in a text file (e.g. done by ttH2017), just need to read the file !
 */
void Pileup::Read_PU_Weights()
{
    cout<<"OBSOLETE ! Or update filepaths !"<<endl; return;


    TString filepath = dir_parent + "xxx/" + GetFilename(samplename); //obsolete
    if(!Check_If_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl; return;}
    else {cout<<endl<<"-- Opening PU weight file : "<<filepath<<endl<<endl<<endl;}

    ifstream file_in(filepath.Data());

    string line;
    while(getline(file_in, line))
    {
        stringstream ss(line);
        int npu; double weight;
        ss >> npu >> weight;

        if(npu > 99) {break;} //Only consider nPU up to 99

        v_PU_weights_Nom.push_back(weight);
    }

    return;
}


/**
 * Relate filename to samplename
 * //FIXME -- need to fix list !
 */
TString Pileup::GetFilename(TString name)
{
    TString _inputFile;

    if (name.Contains("WZTo3LNu") ) _inputFile = "WZTo3LNu.txt";
    else if (name.Contains("TTZ_M1to10") ) _inputFile = "TTZ_M1to10.txt";
    else if (name.Contains("TTZ") ) _inputFile = "TTZ.txt";

    else if (name.Contains("TTWJetsToLNu") ) _inputFile = "TTW_psw.txt";
    else if (name.Contains("tZq") ) _inputFile = "tZq.txt";
    else if (name.Contains("ttHJetToNonbb") ) _inputFile = "ttHJetToNonbb.txt";

    else if (name.Contains("DYJetsToLL_M-10to50") ) _inputFile = "DYJets_M10to50.txt";
    else if (name.Contains("DYJetsToLL_M-50_") ) _inputFile = "DYJets_M50.txt";
    else if (name.Contains("TTGJets") ) _inputFile = "TTGJets.txt";
    else if (name.Contains("TTTo2L2Nu_TuneCP5_PSweights") ) _inputFile = "TTToDiLep_psw.txt";
    else if (name.Contains("TTToSemiLeptonic_TuneCP5_PSweights") ) _inputFile = "TTToSemiLep_psw.txt";
    else if (name.Contains("TTTT") ) _inputFile = "TTTT.txt";
    else if (name.Contains("TTWW") ) _inputFile = "TTWW.txt";
    else if (name.Contains("WJets") ) _inputFile = "WJets.txt";
    else if (name.Contains("WWTo2L2Nu_NNPDF31") ) _inputFile = "WW.txt";
    else if (name.Contains("WWTo2L2Nu_DoubleScattering") ) _inputFile = "WWTo2L2Nuds.txt";
    else if (name.Contains("WWW") ) _inputFile = "WWW.txt";
    else if (name.Contains("WWZ") ) _inputFile = "WWZ.txt";
    else if (name.Contains("WZZ") ) _inputFile = "WZZ.txt";
    else if (name.Contains("ZZZ") ) _inputFile = "ZZZ.txt";
    else if (name.Contains("ZZTo") ) _inputFile = "ZZ.txt"; //FIXME : which ZZ sample is it ?

    else {return "xxx";} //Raise error

    return _inputFile;
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
