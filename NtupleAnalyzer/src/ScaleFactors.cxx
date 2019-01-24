// NB : using "max(1, min(Nbins, FindBin) )" ensures that we always read a bin between 1 and N (included)

#include "../include/ScaleFactors.h"

using namespace std;





//--------------------------------------------
// ##     ## ######## ##       ########  ######## ########     ######## ##     ## ##    ##  ######
// ##     ## ##       ##       ##     ## ##       ##     ##    ##       ##     ## ###   ## ##    ##
// ##     ## ##       ##       ##     ## ##       ##     ##    ##       ##     ## ####  ## ##
// ######### ######   ##       ########  ######   ########     ######   ##     ## ## ## ## ##
// ##     ## ##       ##       ##        ##       ##   ##      ##       ##     ## ##  #### ##
// ##     ## ##       ##       ##        ##       ##    ##     ##       ##     ## ##   ### ##    ## ###
// ##     ## ######## ######## ##        ######## ##     ##    ##        #######  ##    ##  ######  ###
//--------------------------------------------


//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Exists(const TString& name)
{
    struct stat buffer;
    return (stat (name.Data(), &buffer) == 0); //true if file exists
}










//--------------------------------------------
// #### ##    ## #### ########
//  ##  ###   ##  ##     ##
//  ##  ####  ##  ##     ##
//  ##  ## ## ##  ##     ##
//  ##  ##  ####  ##     ##
//  ##  ##   ###  ##     ##
// #### ##    ## ####    ##
//--------------------------------------------
/**
 * Constructor -- Open SF files
 */
ScaleFactors::ScaleFactors(TString samplename, bool debug_pileup)
{
    cout<<"-- Opening efficiency files..."<<endl;

    TString filepath;

//Open the lepton SF TFiles
//--------------------------------------------
    path_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/lepton_SF/";

    filepath = path_dir + "el_scaleFactors_gsf_ptLt20.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_Reco_Ele_ptLt20 = TFile::Open(filepath);}

    filepath = path_dir + "el_scaleFactors_gsf_ptGt20.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_Reco_Ele_ptGt20 = TFile::Open(filepath);}

    filepath = path_dir + "el_reco_loose_SF.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_LooseID_Ele = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_trkEff_ptLt10.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_Track_Mu_ptLt10 = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_trkEff_ptGt10.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_Track_Mu_ptGt10 = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_ptLt30.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_LooseID_Mu_ptLt30 = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_ptGt30.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_LooseID_Mu_ptGt30 = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_trkVtxCut_and_isoEff.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_Iso_Mu = TFile::Open(filepath);}

    filepath = path_dir + "lepMVAEffSF_e_2lss.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_LooseToTight_Ele_2lss = TFile::Open(filepath);}

    filepath = path_dir + "lepMVAEffSF_m_2lss.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_LooseToTight_Mu_2lss = TFile::Open(filepath);}

    filepath = path_dir + "lepMVAEffSF_e_3l.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_LooseToTight_Ele_3l = TFile::Open(filepath);}

    filepath = path_dir + "lepMVAEffSF_m_3l.root";
    if(!Check_File_Exists(filepath) ) {cout<<BOLD(FRED("Error : SF file "<<filepath<<" not found !"))<<endl;}
    else {file_Eff_LooseToTight_Mu_3l = TFile::Open(filepath);}

    //Read the histograms and tgraphs
    hEff_Reco_Ele_ptLt20       = (TH2D*) file_Eff_Reco_Ele_ptLt20->Get("EGamma_SF2D");
    hEff_Reco_Ele_ptGt20       = (TH2D*) file_Eff_Reco_Ele_ptGt20->Get("EGamma_SF2D");
    hEff_LooseID_Ele           = (TH2D*) file_Eff_LooseID_Ele->Get("EGamma_SF2D");

    gEff_Track_Mu_ptLt10       = (TGraph*) file_Eff_Track_Mu_ptLt10->Get("ratio_eff_eta3_tk0_dr030e030_corr");
    gEff_Track_Mu_ptGt10       = (TGraph*) file_Eff_Track_Mu_ptGt10->Get("ratio_eff_eta3_dr030e030_corr");
    hEff_LooseID_Mu_ptLt30     = (TH2*) file_Eff_LooseID_Mu_ptLt30->Get("NUM_LooseID_DEN_genTracks_pt_abseta");
    hEff_LooseID_Mu_ptGt30     = (TH2*) file_Eff_LooseID_Mu_ptGt30->Get("NUM_LooseID_DEN_genTracks_pt_abseta");
    hEff_Iso_Mu                = (TH2*) file_Eff_Iso_Mu->Get("NUM_ttHLoo_DEN_LooseID");

    hEff_LooseToTight_Ele_2lss = (TH2F*) file_Eff_LooseToTight_Ele_2lss->Get("sf");
    hEff_LooseToTight_Mu_2lss  = (TH2F*) file_Eff_LooseToTight_Mu_2lss->Get("sf");
    hEff_LooseToTight_Ele_3l   = (TH2F*) file_Eff_LooseToTight_Ele_3l->Get("sf");
    hEff_LooseToTight_Mu_3l    = (TH2F*) file_Eff_LooseToTight_Mu_3l->Get("sf");
//--------------------------------------------

//Open btagging efficiency file
//--------------------------------------------
    path_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/btag_SF/";

    TString btag_csv_filename = "DeepCSV_94XSF_V3_B_F.csv";

    const string& tagger_type = "deepcsv";
    const string& inputCSVfile = (const string&) path_dir.Data() + btag_csv_filename.Data();
    const string& measType = "iterativefit";

    if(!Check_File_Exists(inputCSVfile) ) {cout<<BOLD(FRED("Error : btag .csv file "<<inputCSVfile<<" not found !"))<<endl;}

    //Load the .csv file containing the SFs using BTagCalibration Standalone tool
    calib = BTagCalibrationX(tagger_type, inputCSVfile);

   //-- List the btagging systematics (+ central) //NB : "central" must always be in first position ! //NB2 : can not modify order, ... ?
    v_sysType.push_back("central");
    v_sysType.push_back("up_jes");
    v_sysType.push_back("down_jes");
    v_sysType.push_back("up_lf");
    v_sysType.push_back("down_lf");
    v_sysType.push_back("up_hf");
    v_sysType.push_back("down_hf");
    v_sysType.push_back("up_hfstats1");
    v_sysType.push_back("down_hfstats1");
    v_sysType.push_back("up_hfstats2");
    v_sysType.push_back("down_hfstats2");
    v_sysType.push_back("up_lfstats1");
    v_sysType.push_back("down_lfstats1");
    v_sysType.push_back("up_lfstats2");
    v_sysType.push_back("down_lfstats2");
    v_sysType.push_back("up_cferr1");
    v_sysType.push_back("down_cferr1");
    v_sysType.push_back("up_cferr2");
    v_sysType.push_back("down_cferr2");

    // * Systematics for bottom flavor jets:
   // Jet energy scale, "jes"
   // Light flavor contamination, "lf"
   // Linear and quadratic statistical fluctuations, "hfstats1" and "hfstats2"

   // * Systematics for light flavor jets:
   // Jet energy scale, "jes"
   // Heavy flavor contamination, "hf"
   // Linear and quadratic statistical fluctuations, "lfstats1" and "lfstats2"

   // * Systematics for charm flavor jets:
   // Linear and quadratic uncertainties, "cferr1" and "cferr2"

    //List the possible jet flavors
    v_jetFlavor.push_back(BTagEntryX::FLAV_B);
    v_jetFlavor.push_back(BTagEntryX::FLAV_C);
    v_jetFlavor.push_back(BTagEntryX::FLAV_UDSG);

    //To read & evaluate scale factors -- 1 per jet systematic
    v_btag_reader.resize(v_sysType.size());

    for(int isys=0; isys<v_sysType.size(); isys++) //all syst
    {
        const string& tmp = (const string&) v_sysType[isys].Data();

        v_btag_reader[isys] = BTagCalibrationXReader(BTagEntryX::OP_RESHAPING, tmp);

        for(int iflavor=0; iflavor<v_jetFlavor.size(); iflavor++) //all flavors
        {
            v_btag_reader[isys].load(calib, v_jetFlavor[iflavor], measType);
        }
    }

    cout<<FYEL("-> Efficiency files (lepton, btag) correctly opened")<<endl<<endl;

    //--- PILEUP
    my_pileup.Init(samplename);
    if(!my_pileup.Fill_MC_vector_from_Pileup_Profile()) {cout<<endl<<FRED("Pileup files not loaded...")<<endl<<endl<<endl;} //if files not found

    my_pileup.Compute_PU_Weights(); //Store PU weights in member vectors //obsolete
    if(debug_pileup) {my_pileup.Print_Debug();} //Print some debug infos
}

/**
 * Destructor -- close files
 */
ScaleFactors::~ScaleFactors()
{
    //Delete objects
    if(hEff_Reco_Ele_ptLt20) {delete hEff_Reco_Ele_ptLt20; hEff_Reco_Ele_ptLt20                   = NULL;}
    if(hEff_Reco_Ele_ptGt20) {delete hEff_Reco_Ele_ptGt20; hEff_Reco_Ele_ptGt20                   = NULL;}
    if(hEff_LooseID_Ele) {delete hEff_LooseID_Ele; hEff_LooseID_Ele                               = NULL;}
    if(gEff_Track_Mu_ptLt10) {delete gEff_Track_Mu_ptLt10; gEff_Track_Mu_ptLt10                   = NULL;}
    if(gEff_Track_Mu_ptGt10) {delete gEff_Track_Mu_ptGt10; gEff_Track_Mu_ptGt10                   = NULL;}
    if(hEff_LooseID_Mu_ptLt30) {delete hEff_LooseID_Mu_ptLt30; hEff_LooseID_Mu_ptLt30             = NULL;}
    if(hEff_LooseID_Mu_ptGt30) {delete hEff_LooseID_Mu_ptGt30; hEff_LooseID_Mu_ptGt30             = NULL;}
    if(hEff_Iso_Mu) {delete hEff_Iso_Mu; hEff_Iso_Mu                                              = NULL;}
    if(hEff_LooseToTight_Ele_2lss) {delete hEff_LooseToTight_Ele_2lss; hEff_LooseToTight_Ele_2lss = NULL;}
    if(hEff_LooseToTight_Mu_2lss) {delete hEff_LooseToTight_Mu_2lss; hEff_LooseToTight_Mu_2lss    = NULL;}
    if(hEff_LooseToTight_Ele_3l) {delete hEff_LooseToTight_Ele_3l; hEff_LooseToTight_Ele_3l       = NULL;}
    if(hEff_LooseToTight_Mu_3l) {delete hEff_LooseToTight_Mu_3l; hEff_LooseToTight_Mu_3l          = NULL;}

    //Close TFiles
    file_Eff_Reco_Ele_ptLt20->Close();
    file_Eff_Reco_Ele_ptGt20->Close();
    file_Eff_LooseID_Ele->Close();
    file_Eff_Track_Mu_ptLt10->Close();
    file_Eff_Track_Mu_ptGt10->Close();
    file_Eff_LooseID_Mu_ptLt30->Close();
    file_Eff_LooseID_Mu_ptGt30->Close();
    file_Eff_Iso_Mu->Close();
    file_Eff_LooseToTight_Ele_2lss->Close();
    file_Eff_LooseToTight_Mu_2lss->Close();
    file_Eff_LooseToTight_Ele_3l->Close();
    file_Eff_LooseToTight_Mu_3l->Close();
}











//--------------------------------------------
// ##       ######## ########  ########  #######  ##    ##     ######  ########
// ##       ##       ##     ##    ##    ##     ## ###   ##    ##    ## ##
// ##       ##       ##     ##    ##    ##     ## ####  ##    ##       ##
// ##       ######   ########     ##    ##     ## ## ## ##     ######  ######
// ##       ##       ##           ##    ##     ## ##  ####          ## ##
// ##       ##       ##           ##    ##     ## ##   ###    ##    ## ##
// ######## ######## ##           ##     #######  ##    ##     ######  ##
//--------------------------------------------




/**
 * Return efficiency SF for reconstruction of loose electrons
 * Var (to access uncert.) null by default
 */
float ScaleFactors::Get_SF_RecoToLoose_Ele(float pt, float eta, float var)
{
    TH2D* h = NULL; //Will point to different SF histos

    int ptbin = 0, etabin = 0;
    float eff_reco = 1, eff_looseID = 1;

    if(pt < 20) {h = hEff_Reco_Ele_ptLt20;}
    else {h = hEff_Reco_Ele_ptGt20;}

    etabin  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(eta) ) ); //X <-> eta
    ptbin = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(pt) ) ); //Y <-> pt
    eff_reco = h->GetBinContent(etabin,ptbin) + var*h->GetBinError(etabin,ptbin);

    h = hEff_LooseID_Ele;
    etabin = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(eta) ) );
    ptbin  = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(pt) ) );
    eff_looseID = h->GetBinContent(etabin,ptbin) + var*h->GetBinError(etabin,ptbin);

    // cout<<"Get_SF_RecoToLoose_Ele = "<<eff_reco*eff_looseID<<endl;

    return eff_reco*eff_looseID;
}

/**
 * Return efficiency SF for reconstruction of loose muons
 * Var (to access uncert.) null by default
 */
float ScaleFactors::Get_SF_RecoToLoose_Mu(float pt, float eta, float var)
{
    TGraph* g = NULL; //Will point to different SF graphs
    TH2* h = NULL; //Will point to different SF histos

    int ptbin = 0, etabin = 0;
    float eff_track = 1, eff_looseID = 1, eff_iso = 1;

    if(pt < 10) {g = gEff_Track_Mu_ptLt10;}
    else {g = gEff_Track_Mu_ptGt10;}
    eff_track = g->Eval(eta);

    if(pt < 30) {h = hEff_LooseID_Mu_ptLt30;}
    else {h = hEff_LooseID_Mu_ptGt30;}
    ptbin  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(pt))); //X <-> pt
    etabin = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(fabs(eta)))); //Y <-> abs(eta)
    eff_looseID = h->GetBinContent(ptbin,etabin) + var*h->GetBinError(ptbin,etabin);
    if(pt >= 15 && pt < 30 && fabs(eta) >= 2.1) {eff_looseID = 1.;} //Hard-coded fix : POG put some null bins in eff file !

    h = hEff_Iso_Mu;
    ptbin  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(pt))); //X <-> pt
    etabin = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(fabs(eta)))); //Y <-> abs(eta)
    eff_iso = h->GetBinContent(ptbin,etabin) + var*h->GetBinError(ptbin,etabin);

    // cout<<"Get_SF_RecoToLoose_Mu = "<<eff_track*eff_looseID*eff_iso<<endl;

    return eff_track*eff_looseID*eff_iso;
}


/**
 * Return efficiency SF for identifying loose leptons as "tight-ID"
 */
float ScaleFactors::Get_SF_LooseToTight_Leptons(int pdgid, int nlep, float pt, float eta, float var)
{
    if(abs(pdgid) != 11 && abs(pdgid) != 13) {return 1;}
    if(nlep != 2 && nlep != 3) {return 1;}

    TH2* h = NULL; //Will point to different SF histos

    int ptbin = 0, etabin = 0;
    // float eff_looseToTight = 1;

    if(abs(pdgid) == 11)
    {
        if(nlep == 2) {h = hEff_LooseToTight_Ele_2lss;}
        else {h = hEff_LooseToTight_Ele_3l;}
    }
    else
    {
        if(nlep == 2) {h = hEff_LooseToTight_Mu_2lss;}
        else {h = hEff_LooseToTight_Mu_3l;}
    }

    ptbin  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(pt))); //X <-> pt
    etabin = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(fabs(eta)))); //Y <-> abs(eta)

    // cout<<"Get_SF_LooseToTight_Leptons = "<<h->GetBinContent(ptbin,etabin)<<endl;

    float result = h->GetBinContent(ptbin,etabin);
    if(abs(pdgid) == 11) {result+= var * h->GetBinError(ptbin,etabin);}

    return result;
}


/**
 * Get global lepton SF, by combining other SFs
 * Var (to access uncert.) null by default
*/
float ScaleFactors::Get_Lepton_SF(int nlep, int pdgid, float pt, float eta, TString var_type)
{
    if(var_type != "" && var_type != "looseUp" && var_type != "looseDown" && var_type != "tightUp" && var_type != "tightDown") {return 1;}
    if(abs(pdgid) != 11 && abs(pdgid) != 13) {return 1;}
    if(nlep != 2 && nlep != 3) {return 1;}

    float var=0;
    float eff_recoToLoose = 1;

    if(var_type == "looseUp") {var = +1;}
    else if(var_type == "looseDown") {var = -1;}
    else {var = 0;}
    if(abs(pdgid) == 11) {eff_recoToLoose = Get_SF_RecoToLoose_Ele(pt, eta, var);}
    else {eff_recoToLoose = Get_SF_RecoToLoose_Mu(pt, eta, var);}

    if(var_type == "looseUp") {var = +1;}
    else if(var_type == "looseDown") {var = -1;}
    else {var = 0;}
    float eff_looseToTight = Get_SF_LooseToTight_Leptons(pdgid, nlep, pt, eta, var); // var is ignored in all cases for the tight part (systematics handled as nuisance parameter)

    float sf = eff_recoToLoose * eff_looseToTight;

    if(!sf) {cout<<FRED("Problem in ScaleFactors.cxx : Get_Lepton_SF() = 0")<<endl;}

    // cout<<"Get_Lepton_SF = "<<sf<<endl;

    return sf;
}





//--------------------------------------------
// ######## ########  ####  ######    ######   ######## ########      ######  ########
//    ##    ##     ##  ##  ##    ##  ##    ##  ##       ##     ##    ##    ## ##
//    ##    ##     ##  ##  ##        ##        ##       ##     ##    ##       ##
//    ##    ########   ##  ##   #### ##   #### ######   ########      ######  ######
//    ##    ##   ##    ##  ##    ##  ##    ##  ##       ##   ##            ## ##
//    ##    ##    ##   ##  ##    ##  ##    ##  ##       ##    ##     ##    ## ##
//    ##    ##     ## ####  ######    ######   ######## ##     ##     ######  ##
//--------------------------------------------

/**
 * Get trigger SF, based on nof leptons, ele or muons, and pT of leading lepton
 * Depends on nof leptons (2 or 3), lepton flavors, and pT oh hardest lepton.
 * Shift variable = +/- 1 ---> get values shift by +/-1 sigma. Null by default
 */
float ScaleFactors::Get_Trigger_SF(int nlep, int pdgid1, float pt1, int pdgid2, float shift/*=0*/)
{
    if(nlep == 3) {return 1. + shift*0.05;}
    else if(nlep == 2)
    {
        if(abs(pdgid1) == 13 && abs(pdgid2) == 13) //uu
        {
            if(pt1 < 35) {return 0.972 + shift*0.006;}
            else {return 0.994 + shift*0.001;}
        }
        else if( (abs(pdgid1) == 13 && abs(pdgid2) == 11) || (abs(pdgid1) == 11 && abs(pdgid2) == 13) ) //ue, eu
        {
            if(pt1 < 35) {return 0.952 + shift*0.008;}
            else if(pt1 < 50) {return 0.983 + shift*0.003;}
            else {return 1. + shift*0.001;}
        }
        else //ee
        {
            if(pt1 < 30) {return 0.937 + shift*0.027;}
            else {return 0.991 + shift*0.002;}
        }
    }

    cout<<BOLD(FRED("Error : wrong nlep value !"))<<endl; return 1.;
}







//--------------------------------------------
// ########  ########    ###     ######       ######  ########
// ##     ##    ##      ## ##   ##    ##     ##    ## ##
// ##     ##    ##     ##   ##  ##           ##       ##
// ########     ##    ##     ## ##   ####     ######  ######
// ##     ##    ##    ######### ##    ##           ## ##
// ##     ##    ##    ##     ## ##    ##     ##    ## ##
// ########     ##    ##     ##  ######       ######  ##
//--------------------------------------------


/**
 * Get Btag SF for given jet
 * Inspired from : https://github.com/cmarper/ttH/blob/a5a62199b624ea2cde07d8d79d204b8bdc5de499/macros/bTagSF_CSVshape.cc
 */
float ScaleFactors::Get_Btag_SF(TString syst_tmp, int jetFlavor, float pt, float eta, float deepcsv)
{
    const string& syst = (const string&) syst_tmp.Data(); //Convert TString to required type

    float SF = 1.;

    BTagEntryX::JetFlavor flav;
    bool isBFlav = false;
    bool isCFlav = false;
    bool isLFlav = false;
    // cout<<"isBFlav "<<isBFlav<<" / isCFlav "<<isCFlav<<" / isLFlav "<<isLFlav<<endl;

    if(abs(jetFlavor) == 5)
    {
        flav = BTagEntryX::FLAV_B;
        isBFlav = true;
    }
    else if(abs(jetFlavor) == 4)
    {
        flav = BTagEntryX::FLAV_C;
        isCFlav = true;
    }
    else
    {
        flav = BTagEntryX::FLAV_UDSG;
        isLFlav = true;
    }

    for(int isys=0; isys<v_sysType.size(); isys++)
    {
        int index_reader; //To get the correct BTagCalibrationXReader element from vector

        const string& tmp = (const string&) v_sysType[isys].Data();

        if(syst == tmp)
        {
            if((syst == "up_jes" && isCFlav)
            || (syst == "down_jes" && isCFlav)
            || (syst == "up_lf" && !isBFlav)
            || (syst == "down_lf" && !isBFlav)
            || (syst == "up_hf" && !isLFlav)
            || (syst == "down_hf" && !isLFlav)
            || (syst == "up_hfstats1" && !isBFlav)
            || (syst == "down_hfstats1" && !isBFlav)
            || (syst == "up_hfstats2" && !isBFlav)
            || (syst == "down_hfstats2" && !isBFlav)
            || (syst == "up_lfstats1" && !isLFlav)
            || (syst == "down_lfstats1" && !isLFlav)
            || (syst == "up_lfstats2" && !isLFlav)
            || (syst == "down_lfstats2" && !isLFlav)
            || (syst == "up_cferr1" && !isCFlav)
            || (syst == "down_cferr1" && !isCFlav)
            || (syst == "up_cferr2" && !isCFlav)
            || (syst == "down_cferr2" && !isCFlav) ) {index_reader = 0;} //If incorrect sysType/flav association, use "central" (first element)

            else {index_reader = isys;} //If correct sysType/flav association, use associated reader
            // cout<<"index_reader = "<<index_reader<<endl;
            // SF = v_btag_reader[index_reader].eval(flav, eta, pt, deepcsv);

            //CHANGED -- use eval_auto_bounds, as recommended
            //NB1 : Uncert. increased by factor of two if pt value outside provided range
            //NB2 : Must use fabs(eta) ; returns 1 above eta > 2.4 ?
            if(!index_reader) {SF = v_btag_reader[index_reader].eval_auto_bounds("central", flav, fabs(eta), pt, deepcsv);} //If calling 'central' reader, must call it with 'central' as systName
            else {SF = v_btag_reader[index_reader].eval_auto_bounds(syst, flav, fabs(eta), pt, deepcsv);}

            break; //syst found
        }
    }

    if(!SF) {SF = 1.;}

    return SF;
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
 * Get PU weight for given MC event
 * "var" can be 'nom' / 'up' / 'down'
 */
float ScaleFactors::Get_Pileup_SF(int nPU, TString var)
{
    // cout<<BOLD(FYEL("--- Entering : Compute_PU_Weight()"))<<endl;

    if(nPU==99) {return 1;} //Fix binning issue (I stopped at 98)

    if(nPU < 0 || nPU >= my_pileup.v_PU_weights_Nom.size())
    {
        // cout<<BOLD(FRED("ERROR ! Wrong nPU argument in : Pileup::Get_PU_Weight() ! ==> nPU = "))<<nPU<<endl;
        return 1;
    }

    if(var=="nom")
    {
        return my_pileup.v_PU_weights_Nom[nPU];
    }
    else if(var=="up")
    {
        return my_pileup.v_PU_weights_Up[nPU];
    }
    else if(var=="down")
    {
        return my_pileup.v_PU_weights_Down[nPU];
    }
    else {cout<<BOLD(FRED("Error : wrong arg 'var= "<<var<<" ' in Pileup::Get_PU_Weight() !"))<<endl; return 1;}

    return 1;
}




//--------------------------------------------
//  ######   ######     ###    ##       ########  ######
// ##    ## ##    ##   ## ##   ##       ##       ##    ##
// ##       ##        ##   ##  ##       ##       ##
//  ######  ##       ##     ## ##       ######    ######
//       ## ##       ######### ##       ##             ##
// ##    ## ##    ## ##     ## ##       ##       ##    ##
//  ######   ######  ##     ## ######## ########  ######
//--------------------------------------------

void ScaleFactors::Read_Scale_SumWeights(TString samplename, float& sumWeights_nominal, float& sumWeights_scale_originalXWGTUP, float& sumWeights_scale_muF0p5, float& sumWeights_scale_muF2, float& sumWeights_scale_muR0p5, float& sumWeights_scale_muR2, float& sumWeights_scale_muR2muF2, float& sumWeights_scale_muR0p5muF0p5)
{
    TString filename = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/Scales/" + samplename + ".root";
    if(!Check_File_Exists(filename))
    {
        cout<<endl<<endl<<FRED("File "<<filename<<" containing the Sum of Weights for the scale variations was not found. The scale variations are set to 1.")<<endl<<endl;
        return;
    }
    else {cout<<FYEL("-> File containing sum of weights for scale variations correctly opened")<<endl<<endl;}


    TFile* f = TFile::Open(filename);
    TH1F* h_tmp = 0;

    //Nominal
    h_tmp = (TH1F*) f->Get("h_sumWeights_nominal");
    if(!h_tmp) {cout<<FRED("Histogram not found (sum weights scale variation) ! Abort")<<endl; return;}
    sumWeights_nominal = h_tmp->Integral();

    delete h_tmp; h_tmp = 0;

    //originalXWGTUP -- useless ?
    h_tmp = (TH1F*) f->Get("h_sumWeightsScale_originalXWGTUP");
    if(!h_tmp) {cout<<FRED("Histogram not found (sum weights scale variation) ! Abort")<<endl; return;}
    sumWeights_scale_originalXWGTUP = h_tmp->Integral();
    delete h_tmp; h_tmp = 0;

    //muR = 0.5 / muF = 1
    h_tmp = (TH1F*) f->Get("h_sumWeightsScale_muR0p5");
    if(!h_tmp) {cout<<FRED("Histogram not found (sum weights scale variation) ! Abort")<<endl; return;}
    sumWeights_scale_muR0p5 = h_tmp->Integral();
    delete h_tmp; h_tmp = 0;

    //muR = 2 / muF = 1
    h_tmp = (TH1F*) f->Get("h_sumWeightsScale_muR2");
    if(!h_tmp) {cout<<FRED("Histogram not found (sum weights scale variation) ! Abort")<<endl; return;}
    sumWeights_scale_muR2 = h_tmp->Integral();
    delete h_tmp; h_tmp = 0;

    //muR = 1 / muF = 0.5
    h_tmp = (TH1F*) f->Get("h_sumWeightsScale_muF0p5");
    if(!h_tmp) {cout<<FRED("Histogram not found (sum weights scale variation) ! Abort")<<endl; return;}
    sumWeights_scale_muF0p5 = h_tmp->Integral();
    delete h_tmp; h_tmp = 0;

    //muR = 1 / muF = 2
    h_tmp = (TH1F*) f->Get("h_sumWeightsScale_muF2");
    if(!h_tmp) {cout<<FRED("Histogram not found (sum weights scale variation) ! Abort")<<endl; return;}
    sumWeights_scale_muF2 = h_tmp->Integral();
    delete h_tmp; h_tmp = 0;

    //muR = 0.5 / muF = 0.5
    h_tmp = (TH1F*) f->Get("h_sumWeightsScale_muR0p5muF0p5");
    if(!h_tmp) {cout<<FRED("Histogram not found (sum weights scale variation) ! Abort")<<endl; return;}
    sumWeights_scale_muR0p5muF0p5 = h_tmp->Integral();
    delete h_tmp; h_tmp = 0;

    //muR = 2 / muF = 2
    h_tmp = (TH1F*) f->Get("h_sumWeightsScale_muR2muF2");
    if(!h_tmp) {cout<<FRED("Histogram not found (sum weights scale variation) ! Abort")<<endl; return;}
    sumWeights_scale_muR2muF2 = h_tmp->Integral();
    delete h_tmp; h_tmp = 0;

    f->Close();

    return;
}
