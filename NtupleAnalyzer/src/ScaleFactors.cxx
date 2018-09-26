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
//  ######  ########     ######  ##          ###     ######   ######
// ##    ## ##          ##    ## ##         ## ##   ##    ## ##    ##
// ##       ##          ##       ##        ##   ##  ##       ##
//  ######  ######      ##       ##       ##     ##  ######   ######
//       ## ##          ##       ##       #########       ##       ##
// ##    ## ##          ##    ## ##       ##     ## ##    ## ##    ##
//  ######  ##           ######  ######## ##     ##  ######   ######
//--------------------------------------------

/**
 * Constructor -- Open SF files
 */
ScaleFactors::ScaleFactors()
{
    path_dir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/SF_files/";

    //Open the SF TFiles
    TString filepath;

    filepath = path_dir + "el_scaleFactors_gsf_ptLt20.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_Reco_Ele_ptLt20 = TFile::Open(filepath);}

    filepath = path_dir + "el_scaleFactors_gsf_ptGt20.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_Reco_Ele_ptGt20 = TFile::Open(filepath);}

    filepath = path_dir + "el_reco_loose_SF.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_LooseID_Ele = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_trkEff_ptLt10.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_Track_Mu_ptLt10 = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_trkEff_ptGt10.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_Track_Mu_ptGt10 = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_ptLt30.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_LooseID_Mu_ptLt30 = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_ptGt30.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_LooseID_Mu_ptGt30 = TFile::Open(filepath);}

    filepath = path_dir + "mu_scaleFactors_trkVtxCut_and_isoEff.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_Iso_Mu = TFile::Open(filepath);}

    filepath = path_dir + "lepMVAEffSF_e_2lss.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_LooseToTight_Ele_2lss = TFile::Open(filepath);}

    filepath = path_dir + "lepMVAEffSF_m_2lss.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_LooseToTight_Mu_2lss = TFile::Open(filepath);}

    filepath = path_dir + "lepMVAEffSF_e_3l.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
    else {file_Eff_LooseToTight_Ele_3l = TFile::Open(filepath);}

    filepath = path_dir + "lepMVAEffSF_m_3l.root";
    if(!Check_File_Exists(filepath) ) {cout<<"Error : SF file "<<filepath<<" not found !"<<endl;}
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

float ScaleFactors::Get_SF_RecoToLoose_Ele(float pt, float eta)
{
    TH2D* h = NULL; //Will point to different SF histos

    int ptbin = 0, etabin = 0;
    float eff_reco = 1, eff_looseID = 1;

    if(pt < 20) {h = hEff_Reco_Ele_ptLt20;}
    else {h = hEff_Reco_Ele_ptGt20;}

    etabin  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(eta) ) ); //X <-> eta
    ptbin = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(pt) ) ); //Y <-> pt
    eff_reco = h->GetBinContent(etabin,ptbin);

    h = hEff_LooseID_Ele;
    etabin = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(eta) ) );
    ptbin  = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(pt) ) );
    eff_looseID = h->GetBinContent(etabin,ptbin);

    // cout<<"Get_SF_RecoToLoose_Ele = "<<eff_reco*eff_looseID<<endl;

    return eff_reco*eff_looseID;
}


float ScaleFactors::Get_SF_RecoToLoose_Mu(float pt, float eta)
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
    eff_looseID = h->GetBinContent(ptbin,etabin);

    h = hEff_Iso_Mu;

    ptbin  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(pt))); //X <-> pt
    etabin = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(fabs(eta)))); //Y <-> abs(eta)
    eff_iso = h->GetBinContent(ptbin,etabin);

    // cout<<"Get_SF_RecoToLoose_Mu = "<<eff_track*eff_looseID*eff_iso<<endl;

    return eff_track*eff_looseID*eff_iso;
}



float ScaleFactors::Get_SF_LooseToTight_Leptons(int pdgid, int nlep, float pt, float eta)
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

    return h->GetBinContent(ptbin,etabin);
}



float ScaleFactors::Get_Lepton_SF(int nlep, int pdgid, float pt, float eta)
{
    if(abs(pdgid) != 11 && abs(pdgid) != 13) {return 1;}
    if(nlep != 2 && nlep != 3) {return 1;}

    float eff_recoToLoose = 1;

    if(abs(pdgid) == 11) {eff_recoToLoose = Get_SF_RecoToLoose_Ele(pt, eta);}
    else {eff_recoToLoose = Get_SF_RecoToLoose_Mu(pt, eta);}

    float eff_looseToTight = Get_SF_LooseToTight_Leptons(pdgid, nlep, pt, eta);

    // cout<<"Get_Lepton_SF = "<<eff_recoToLoose * eff_looseToTight<<endl;

    return eff_recoToLoose * eff_looseToTight;
}


float ScaleFactors::Get_Trigger_SF(int nlep, int pdgid1, float pt1, int pdgid2, float pt2)
{
    if(nlep == 3) {return 1;}
    else if(nlep == 2)
    {
        if(abs(pdgid1) == 13 && abs(pdgid2) == 13) //uu
        {
            if(pt1 < 35 || pt2 < 35) {return 0.972;}
            else {return 0.994;}
        }
        else if( (abs(pdgid1) == 13 && abs(pdgid2) == 11) || (abs(pdgid1) == 11 && abs(pdgid2) == 13) ) //ue, eu
        {
            if(pt1 < 35 || pt2 < 35) {return 0.952;}
            else if(pt1 < 50 || pt2 < 50) {return 0.983;}
            else {return 1.;}
        }
        else //ee
        {
            if(pt1 < 30 || pt2 < 30) {return 0.937;}
            else {return 0.991;}
        }
    }

    cout<<"Error : wrong nlep value !"<<endl; return 1;
}


// float triggerSF_ttH(int pdgid1, float pt1, int pdgid2, float pt2, int nlep, float var_ee=0){
//     if (var_ee!=0) assert(0); // NOT IMPLEMENTED
//     if (nlep>2) return 1;
//     if (abs(pdgid1)==11 && abs(pdgid2)==11){
//         if (std::max(pt1,pt2)<40) return 0.95;
//         else return 0.99;
//     }
//     else if (abs(pdgid1)==13 && abs(pdgid2)==13) {
//         return 1.;
//     }
//     else return 0.98;
// }
