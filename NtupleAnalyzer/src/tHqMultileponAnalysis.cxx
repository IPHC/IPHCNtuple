//ThreeLeptonSelection_THQ3l
// tHqMultileponAnalysis::fillOutputTree()


//----------  FIXMES  -------------------------
// -- mHT value for MEM ?
//--------------------------------------------

//--------------------------------------------
// -- changed SR selections from ttH to tHq
// -- added boolean to choose tHq or ttH selections
//--------------------------------------------

#include "../include/tHqMultileponAnalysis.h"

#include "TSystem.h"
#include "SignalExtractionMVA.cxx"
#include "Helper.cxx"
#include "FakeRate.cxx"
#include "ChargeFlip.cxx"
// #include "BTagging.cxx"

#define kCat_3l_2b_2j   0
#define kCat_3l_1b_2j   1
#define kCat_3l_2b_1j   2
#define kCat_3l_1b_1j   3
#define kCat_3l_2b_0j   4
#define kCat_4l_2b      5
#define kCat_4l_1b      6
#define kCat_2lss_2b_4j 7
#define kCat_2lss_1b_4j 8
#define kCat_2lss_2b_3j 9
#define kCat_2lss_1b_3j 10
#define kCat_2lss_2b_2j 11

#define kCat_3l_1b_0j 12
#define kCat_3l_0b_1j 13
#define kCat_3l_0b_0j 14

#define kCat_2lss_2b_1j 15
#define kCat_2lss_1b_2j 16
#define kCat_2lss_1b_1j 17

#define DEBUG false

using namespace std;

//Cutflow vector --> Count events for THQ_3l_SR selection, after each cut
vector<double> v_cutflow(15);

bool use_ttH_selections = true; //FIXME
//--------------------------------------------









//--------------------------------------------
// #### ##    ## #### ######## ####    ###    ##       #### ########    ###    ######## ####  #######  ##    ##
//  ##  ###   ##  ##     ##     ##    ## ##   ##        ##       ##    ## ##      ##     ##  ##     ## ###   ##
//  ##  ####  ##  ##     ##     ##   ##   ##  ##        ##      ##    ##   ##     ##     ##  ##     ## ####  ##
//  ##  ## ## ##  ##     ##     ##  ##     ## ##        ##     ##    ##     ##    ##     ##  ##     ## ## ## ##
//  ##  ##  ####  ##     ##     ##  ######### ##        ##    ##     #########    ##     ##  ##     ## ##  ####
//  ##  ##   ###  ##     ##     ##  ##     ## ##        ##   ##      ##     ##    ##     ##  ##     ## ##   ###
// #### ##    ## ####    ##    #### ##     ## ######## #### ######## ##     ##    ##    ####  #######  ##    ##
//--------------------------------------------



/**
 * Default constructor
 */
// tHqMultileponAnalysis::tHqMultileponAnalysis()
// {
//     cout<<BOLD(FRED("DO NOT USE THE DEFAULT CONSTRUCTOR !"))<<endl<<endl;
// }


/**
 * Destructor
 */
tHqMultileponAnalysis::~tHqMultileponAnalysis()
{
    delete vEvent;
    delete vElectronLoose ; delete vElectronFakeable ; delete vElectronTight;
    delete vMuonLoose; delete vMuonFakeable; delete vMuonTight;
    delete vTauLoose; delete vTauMedium; delete vTauFakeable;
    delete vJetLoose;
    delete vTruth;

    delete tOutput;
    delete outputfile;

    // delete f_CSVwgt_HF; delete f_CSVwgt_LF; delete f_BTag_eff;
	delete f_QFwgt; delete f_FRwgt;

    sf->~ScaleFactors(); //Destroy SF object
}

/**
 * Overloaded constructor -- to be used
 * @param inputFileName  File containing paths to input root files
 * @param sampleName     name of sample
 * @param treeName       name of tree to be read in file
 * @param outputFileName name of output file
 * @param isdata         data or mc
 * @param doSystCombine  ?
 * @param xsec           xsec of process
 * @param lumi           lumi
 * @param nowe           number of weighted events in Flatree (needed to compute correct event weights)
 * @param nmax           max number of processed evnts
 */
tHqMultileponAnalysis::tHqMultileponAnalysis(TString inputFileName, TString sampleName, TString treeName, TString outputFileName, bool isdata, bool doSystCombine, float xsec, float lumi, int nowe, int nmax)
{
    //
    _isdata = isdata;
    _xsec = xsec;
    _lumi = lumi;
    _nowe = nowe;
    _nmax = nmax;
    _outputFileName = outputFileName;
    _sampleName = sampleName;
    _process = "toto";

    fChain = new TChain(treeName);

    std::ifstream infile;
    infile.open(inputFileName);
    std::string ifile = "";
    while( getline(infile, ifile) )
    {
        std::string fnameStr = std::string(ifile);
        fChain->Add(fnameStr.c_str());
    }
    infile.close();
    InitTree();

    TString outputfileNameRoot = _outputFileName+".root";
    outputfile = new TFile(outputfileNameRoot, "RECREATE");

    sf = new ScaleFactors(); //Opens all the SF files, fills SF histograms

    InitFiles();
    initializeOutputTree();
}


/**
 * Initialize vectors contained in input TTree
 */
void tHqMultileponAnalysis::InitTree()
{
    // Set branch addresses and branch pointers
    if (!fChain) return;

	//Object collections defined in NtupleProducer
    vEvent            = new std::vector<Event>();
    vElectronLoose    = new std::vector<Electron>();
    vElectronFakeable = new std::vector<Electron>();
    vElectronTight    = new std::vector<Electron>();
    vMuonLoose        = new std::vector<Muon>();
    vMuonFakeable     = new std::vector<Muon>();
    vMuonTight        = new std::vector<Muon>();
    vTauLoose         = new std::vector<Tau>();
    vTauFakeable      = new std::vector<Tau>();
    vTauMedium      = new std::vector<Tau>();
    vJetLoose         = new std::vector<Jet>();
    vTruth            = new std::vector<Truth>();


    fChain->SetBranchAddress("Event", &vEvent);
    fChain->SetBranchAddress("ElectronLoose", &vElectronLoose);
    fChain->SetBranchAddress("ElectronFakeable", &vElectronFakeable);
    fChain->SetBranchAddress("ElectronTight", &vElectronTight);
    fChain->SetBranchAddress("MuonLoose", &vMuonLoose);
    fChain->SetBranchAddress("MuonFakeable", &vMuonFakeable);
    fChain->SetBranchAddress("MuonTight", &vMuonTight);
    fChain->SetBranchAddress("TauLoose", &vTauLoose);
    fChain->SetBranchAddress("TauFakeable", &vTauFakeable);
    fChain->SetBranchAddress("TauMedium", &vTauMedium);
    fChain->SetBranchAddress("JetLoose", &vJetLoose);
    fChain->SetBranchAddress("Truth", &vTruth);
    // fChain->SetBranchAddress("GenJet", &vGenJet);

    return;
}

/**
 * Open the files containing weights and corrections
 */
void tHqMultileponAnalysis::InitFiles()
{
    //--- Load MVA weight files
    Load_MVA();

    //--- Loading weight files and creating corresponding histograms

/*
    // b-tagging
    std::string inputFileHF = "/opt/sbg/scratch1/cms/TTH/weight/csv_rwt_fit_hf_76x_2016_02_08.root"; //-- needs update ??
    std::string inputFileLF = "/opt/sbg/scratch1/cms/TTH/weight/csv_rwt_fit_lf_76x_2016_02_08.root";
    f_CSVwgt_HF = new TFile ((inputFileHF).c_str());
    f_CSVwgt_LF = new TFile ((inputFileLF).c_str());
    fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF); //Fill histos with btagging weights

    // new b-tagging (using BTagCalibrationXStandaloneWhatever)
    // setup calibration + reader
    // BTagCalibrationX *
    calib = BTagCalibrationX("csvv2", "/opt/sbg/scratch1/cms/TTH/weight/CSVv2_Moriond17_B_H.csv");

    // BTagCalibrationXReader *
    reader = BTagCalibrationXReader(  BTagEntryX::OP_LOOSE,  // operating point
                                                            "central",            // central sys type
                                                            {"up", "down"});      // other sys types

    reader.load(   calib,                // calibration instance
                    BTagEntryX::FLAV_B,    // btag flavour
                    "comb");               // measurement type

    reader.load(    calib,                // calibration instance
                    BTagEntryX::FLAV_C,    // btag flavour
                    "comb");              // measurement type

    reader.load(    calib,                // calibration instance
                    BTagEntryX::FLAV_UDSG, // btag flavour
                    "comb");              // measurement type


    // BTag Efficiencies
    TString inputFileBTagEff = "/opt/sbg/scratch1/cms/TTH/weight/output_tt_effBtag.root"; //HARD-CODED
    f_BTag_eff = new TFile (inputFileBTagEff);
    fill_eff_btagging_histos(f_BTag_eff);
*/

    // charge flip
    // TString inputFileQF = "/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleAnalyzer/test/FR_weights/QF_data_el.root"; //HARD-CODED
    TString inputFileQF = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/ElectronChargeMisIdRates_2017_2018Jun22.root"; //HARD-CODED
    f_QFwgt    = new TFile (inputFileQF);
    fillQFhistos(f_QFwgt);

    // fake rate
    TString inputFileFR = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/FR_data_ttH_mva.root"; //HARD-CODED
    f_FRwgt    = new TFile (inputFileFR);
    fillFRhistos(f_FRwgt);
}


/**
 * Re-initialize vectors of objects
 */
void tHqMultileponAnalysis::InitCollections()
{
    //--- Re-init all vectors
	vLeptonLoose.clear();
	vLeptonFakeable.clear();
	vLeptonTight.clear();

	vLooseBTagJets.clear();
	vLightJets.clear();
	vLightJets_FwdPtCut.clear();

	return;
}

/**
 * Initialize the values of the ttH-predefined categories (computed at NTProd level)
 */
void tHqMultileponAnalysis::InitPredefinedCategories()
{
    //ttH2017 categories (implemented at NTP level)
    is_ttH_2lSS              = 0;
    is_ttH_2lSS_SR           = 0;
    is_ttH_2lSS_SR_Data      = 0;
    is_ttH_2lSS_Fake         = 0;
    is_ttH_2lSS_Flip         = 0;
    is_ttH_2lSS_Flip_Data    = 0;
    is_ttH_3l                = 0;
    is_ttH_3l_SR             = 0;
    is_ttH_3l_SR_Data        = 0;
    is_ttH_3l_Fake           = 0;
    is_ttH_ttWctrl           = 0;
    is_ttH_ttWctrl_SR        = 0;
    is_ttH_ttWctrl_SR_Data   = 0;
    is_ttH_ttWctrl_Fake      = 0;
    is_ttH_ttWctrl_Flip      = 0;
    is_ttH_ttWctrl_Flip_Data = 0;
    is_ttH_ttZctrl           = 0;
    is_ttH_ttZctrl_SR        = 0;
    is_ttH_ttZctrl_SR_Data   = 0;
    is_ttH_ttZctrl_Fake      = 0;
    is_ttH_WZctrl            = 0;
    is_ttH_WZctrl_SR         = 0;
    is_ttH_WZctrl_SR_Data    = 0;
    is_ttH_WZctrl_Fake       = 0;
}


/**
 * Initialize the values of my custom categories
 */
void tHqMultileponAnalysis::InitCustomCategories()
{
    //Custom categories
    is_tHq_2lSS              = 0;
    is_tHq_2lSS_SR           = 0;
    is_tHq_2lSS_Training     = 0;
    is_tHq_2lSS_Fake         = 0;
    is_tHq_2lSS_Flip         = 0;
    is_tHq_2lSS_GammaConv    = 0;
    is_tHq_3l                = 0;
    is_tHq_3l_SR             = 0;
    is_tHq_3l_Training       = 0;
    is_tHq_3l_Fake           = 0;
    is_tHq_ttWctrl           = 0;
    is_tHq_ttWctrl_SR        = 0;
    is_tHq_ttWctrl_Fake      = 0;
    is_tHq_ttWctrl_Flip      = 0;
    is_tHq_ttWctrl_GammaConv = 0;
    is_tHq_ttZctrl           = 0;
    is_tHq_ttZctrl_SR        = 0;
    is_tHq_ttZctrl_Fake      = 0;
    is_tHq_WZctrl            = 0;
    is_tHq_WZctrl_SR         = 0;
    is_tHq_WZctrl_Fake       = 0;

    return;
}



/**
 * Re-initialize TLorentzVectors
 * To be called between each sel function (use different objects, so different output)
 */
void tHqMultileponAnalysis::InitTLorentzVectors()
{
    //---- Re-init all TLorentzVectors (call to constructor)
    multilepton_Lepton1_P4 = TLorentzVector();
    multilepton_Lepton2_P4 = TLorentzVector();
    multilepton_Lepton3_P4 = TLorentzVector();
    multilepton_Lepton4_P4 = TLorentzVector();
    multilepton_Lepton1_P4_Matched = TLorentzVector();
    multilepton_Lepton2_P4_Matched = TLorentzVector();
    multilepton_Lepton3_P4_Matched = TLorentzVector();
    multilepton_Lepton4_P4_Matched = TLorentzVector();
    multilepton_Bjet1_P4 = TLorentzVector();
    multilepton_Bjet2_P4 = TLorentzVector();
    multilepton_Bjet1_P4_Matched = TLorentzVector();
    multilepton_Bjet2_P4_Matched = TLorentzVector();
    multilepton_JetHighestPt1_P4 = TLorentzVector();
    multilepton_JetHighestPt2_P4 = TLorentzVector();
    multilepton_JetClosestMw1_P4 = TLorentzVector();
    multilepton_JetClosestMw2_P4 = TLorentzVector();
    multilepton_JetLowestMjj1_P4 = TLorentzVector();
    multilepton_JetLowestMjj2_P4 = TLorentzVector();
    multilepton_JetHighestPt1_2ndPair_P4 = TLorentzVector();
    multilepton_JetHighestPt2_2ndPair_P4 = TLorentzVector();
    multilepton_JetClosestMw1_2ndPair_P4 = TLorentzVector();
    multilepton_JetClosestMw2_2ndPair_P4 = TLorentzVector();
    multilepton_JetLowestMjj1_2ndPair_P4 = TLorentzVector();
    multilepton_JetLowestMjj2_2ndPair_P4 = TLorentzVector();
    multilepton_JetHighestEta1_P4 = TLorentzVector();
    multilepton_JetHighestEta2_P4 = TLorentzVector();
    multilepton_h0_P4 = TLorentzVector();
    multilepton_t1_P4 = TLorentzVector();
    multilepton_t2_P4 = TLorentzVector();
    multilepton_mET = TLorentzVector();
    multilepton_Ptot = TLorentzVector();

	//--- Re-init all jet CSV values
    multilepton_Bjet1_CSV = -999;
    multilepton_Bjet2_CSV = -999;
    multilepton_JetHighestPt1_CSV = -999;
    multilepton_JetHighestPt2_CSV = -999;
    multilepton_JetClosestMw1_CSV = -999;
    multilepton_JetClosestMw2_CSV = -999;
    multilepton_JetLowestMjj1_CSV = -999;
    multilepton_JetLowestMjj2_CSV = -999;
    multilepton_JetHighestPt1_2ndPair_CSV = -999;
    multilepton_JetHighestPt2_2ndPair_CSV = -999;
    multilepton_JetClosestMw1_2ndPair_CSV = -999;
    multilepton_JetClosestMw2_2ndPair_CSV = -999;
    multilepton_JetLowestMjj1_2ndPair_CSV = -999;
    multilepton_JetLowestMjj2_2ndPair_CSV = -999;
    multilepton_JetHighestEta1_CSV = -999;
    multilepton_JetHighestEta2_CSV = -999;
}


/**
 * Initialize all variables and category : to be called between before each selection function (else event can activate 2 categories at once)
 */
void tHqMultileponAnalysis::InitVariables()
{
	weightfake = 0;
	weightflip = 0;

	// signal_3l_TT_MVA = -9; signal_3l_TTV_MVA = -9;
	// signal_2lss_TT_MVA = -9; signal_2lss_TTV_MVA = -9;

	channel = -1;

    //tHq2016 input vars
    nJet25               = 0.;
    nJetLoose            = 0.;
    maxEtaJet25          = 0.;
    lepCharge            = 0.;
    nJetEta1             = 0.;
    dEtaFwdJetBJet       = 0.;
    dEtaFwdJet2BJet      = 0.;
    dEtaFwdJetClosestLep = 0.;
    dPhiHighestPtSSPair  = 0.;
    minDRll              = 0.;
    Lep3Pt               = 0.;

    //ttH2017 input vars
    lep1_conePt    = 0;
    lep2_conePt    = 0;
    mindr_lep1_jet = 0;
    mindr_lep2_jet = 0;
    mT_lep1        = 0;
    mT_lep2        = 0;
    max_lep_eta    = 0;

    //New input vars, to be tested
    minv_FwdJetBJet = 0;
    etaFwdJet       = 0;
    ptFwdJet        = 0;
    LeadJetEta      = 0;
    LeadJetPt       = 0;
    dRjj            = 0;
    deepCSV_max     = 0;
    deepCSV_2nd     = 0;
    Mjj_max         = 0;
    dPhiLepBJet_max = 0;
    dPhijj_max      = 0;

    //Additional vars
    lep1Pt = 0; lep2Pt = 0; lep3Pt = 0; hardestBjetPt = 0; hardestBjetEta = 0; fwdJetPt = 0;

	return;
}

























//--------------------------------------------
// ##        #######   #######  ########     ######## ##     ## ######## ##    ## ########  ######
// ##       ##     ## ##     ## ##     ##    ##       ##     ## ##       ###   ##    ##    ##    ##
// ##       ##     ## ##     ## ##     ##    ##       ##     ## ##       ####  ##    ##    ##
// ##       ##     ## ##     ## ########     ######   ##     ## ######   ## ## ##    ##     ######
// ##       ##     ## ##     ## ##           ##        ##   ##  ##       ##  ####    ##          ##
// ##       ##     ## ##     ## ##           ##         ## ##   ##       ##   ###    ##    ##    ##
// ########  #######   #######  ##           ########    ###    ######## ##    ##    ##     ######
//--------------------------------------------


/**
 * Main function : Loop on all input events, store objects, call selection functions
 */
void tHqMultileponAnalysis::Loop()
{
    if (!fChain) return;

    for(int i=0; i<15; i++)
    {
        v_cutflow[i] = 0;
    }

    Long64_t nentries = fChain->GetEntries();
    int nentries_max = nentries;
    if ( _nmax != -1 && _nmax < nentries ) nentries_max = _nmax;

    cout<<endl<<FYEL("--------------------------------------------")<<endl;
    cout<<BOLD(FYEL(""<< "Will process "<<nentries_max<<"/"<<nentries<<" events" <<""))<<endl;
    cout<<FYEL("--------------------------------------------")<<endl<<endl;;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries_max;jentry++)
    {
        // cout<<endl<<jentry<< " / "<<nentries_max<<endl;

        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%200000 == 0) std::cout << "--- "<<jentry<< " / "<<nentries_max<<endl;

        //if(jentry > 100000) break;

        nb = fChain->GetEntry(jentry);   nbytes += nb;

        InitCollections(); //Re-init object collections for each event
        InitPredefinedCategories(); //Re-init ttH categories
        InitCustomCategories(); //Re-init custom categories

        int nlep_tth = Read_Predefined_Categories(); //Read booleans for ttH2017 categories, implemented at NTP level

        event_id = vEvent->at(0).id;
        MET = vEvent->at(0).metpt;
        // metLD = Compute_METLD(MET);
        metLD = vEvent->at(0).metLD;
        event_run = vEvent->at(0).run;
        // cout<<"event id "<<event_id<<endl;

        weights_pdf.clear();
        ids_pdf.clear();

        if(!_isdata )
        {
            weight = _lumi*_xsec/_nowe;

            mc_weight = vEvent->at(0).mc_weight;
            weight = weight * mc_weight;
        }
        else
        {
            weight    = 1.;
            mc_weight = 1.;
            weight_PV = 1.;
			weightfake = 1.;
			weightflip = 1.;
        }

        // for(int itruth=0; itruth<vTruth->size(); itruth++)
        // {
        //     if(abs(vTruth->at(0).mc_truth_id.at(itruth)) == 12) ) {cout<<vTruth->at(i).mc_truth_id<<" !!!"<<endl;
        //     if(abs(vTruth->at(0).mc_truth_id.at(itruth)) == 14) ) {cout<<vTruth->at(i).mc_truth_id<<" !!!"<<endl;
        //     if(abs(vTruth->at(0).mc_truth_id.at(itruth)) == 16) ) {cout<<vTruth->at(i).mc_truth_id<<" !!!"<<endl;
        //     }
        // }





        //--------------------------------------------
        // #    # #    #  ####  #    #  ####
        // ##  ## #    # #    # ##   # #
        // # ## # #    # #    # # #  #  ####
        // #    # #    # #    # #  # #      #
        // #    # #    # #    # #   ## #    #
        // #    #  ####   ####  #    #  ####
        //--------------------------------------------


        //-- Can use here custom definition for FO/tight
        // for(unsigned int imuon=0; imuon < vMuonLoose->size() ; imuon++)
        // {
        //     Lepton l; l.setLepton(&vMuonLoose->at(imuon),imuon,0,1);
        //
        //     vLeptonLoose.push_back(l);
        //
        //     if(vMuonLoose->at(imuon).isFakeableTTH)
        //     {
        //         vLeptonFakeable.push_back(l);
        //     }
        //
        //     if(vMuonLoose->at(imuon).isTightTTH)
        //     {
        //         vLeptonTight.push_back(l);
        //     }
        // }

        //-- Or can directly use definitions from ttH2017 (coded in NtupleProducer)
        for(unsigned int imuon=0; imuon < vMuonLoose->size() ; imuon++)
        {
            Lepton l; l.setLepton(&vMuonLoose->at(imuon),imuon,0,1);
            vLeptonLoose.push_back(l);
        }
        for(unsigned int imuon=0; imuon < vMuonFakeable->size() ; imuon++)
        {
            Lepton l; l.setLepton(&vMuonFakeable->at(imuon),imuon,0,1);
            vLeptonFakeable.push_back(l);
        }
        for(unsigned int imuon=0; imuon < vMuonTight->size() ; imuon++)
        {
            Lepton l; l.setLepton(&vMuonTight->at(imuon),imuon,0,1);
            vLeptonTight.push_back(l);
        }

        //--------------------------------------------
        // ###### #      ######  ####  ##### #####   ####  #    #  ####
        // #      #      #      #    #   #   #    # #    # ##   # #
        // #####  #      #####  #        #   #    # #    # # #  #  ####
        // #      #      #      #        #   #####  #    # #  # #      #
        // #      #      #      #    #   #   #   #  #    # #   ## #    #
        // ###### ###### ######  ####    #   #    #  ####  #    #  ####
        //--------------------------------------------

        //-- Can use here custom definition for FO/tight
        // for(unsigned int ielectron=0; ielectron < vElectronLoose->size() ; ielectron++)
        // {
        //     Lepton l; l.setLepton(&vElectronLoose->at(ielectron),ielectron,1,0);
        //
        //     vLeptonLoose.push_back(l);
        //
        //     if(vElectronLoose->at(ielectron).isFakeableTTH)
        //     {
        //         vLeptonFakeable.push_back(l);
        //     }
        //
        //     if(vElectronLoose->at(ielectron).isTightTTH)
        //     {
        //         vLeptonTight.push_back(l);
        //     }
        // }

        //-- Or can directly use definitions from ttH2017 (coded in NtupleProducer)
        for(unsigned int ielectron=0; ielectron < vElectronLoose->size() ; ielectron++)
        {
            Lepton l; l.setLepton(&vElectronLoose->at(ielectron),ielectron,1,0);
            vLeptonLoose.push_back(l);
        }
        for(unsigned int ielectron=0; ielectron < vElectronFakeable->size() ; ielectron++)
        {
            Lepton l; l.setLepton(&vElectronFakeable->at(ielectron),ielectron,1,0);
            vLeptonFakeable.push_back(l);
        }
        for(unsigned int ielectron=0; ielectron < vElectronTight->size() ; ielectron++)
        {
            Lepton l; l.setLepton(&vElectronTight->at(ielectron),ielectron,1,0);
            vLeptonTight.push_back(l);
        }

        //--------------------------------------------
        // #####   ##   #    #  ####
        //   #    #  #  #    # #
        //   #   #    # #    #  ####
        //   #   ###### #    #      #
        //   #   #    # #    # #    #
        //   #   #    #  ####   ####
        //--------------------------------------------

        //-- Can use here custom definition for FO/tight
        // for(unsigned int itau=0; itau < vTauLoose->size() ; itau++)
        // {
        //     Lepton l; l.setLepton(&vTauLoose->at(itau),itau,0,0);
        //
        //     vLeptonLoose.push_back(l);
        //
        //     if(vTauLoose->at(itau).isFakeableTTH)
        //     {
        //         vTauFakeable->push_back(vTauLoose->at(itau));
        //     }
        //     if(vTauLoose->at(itau).isMediumTTH)
        //     {
        //         vTauMedium->push_back(vTauLoose->at(itau));
        //     }
        // }

        //-- Or can directly use definitions from ttH2017 (coded in NtupleProducer)
        for(unsigned int itau=0; itau < vTauLoose->size() ; itau++)
        {
            Lepton l; l.setLepton(&vTauLoose->at(itau),itau,0,0);
            vLeptonLoose.push_back(l);
        }
        for(unsigned int itau=0; itau < vTauFakeable->size() ; itau++)
        {
            Lepton l; l.setLepton(&vTauFakeable->at(itau),itau,0,0);
            vLeptonFakeable.push_back(l);
        }




        //--------------------------------------------
        // #####             # ###### #####  ####
        // #    #            # #        #   #
        // #####  #####      # #####    #    ####
        // #    #            # #        #        #
        // #    #       #    # #        #   #    #
        // #####         ####  ######   #    ####
        //--------------------------------------------

        nLooseBJets  = 0;
        nMediumBJets = 0;
        nLightJets = 0;
        nLightJets_Fwd40 = 0;

        is_hasJetTransitionRegion = 0;

        //-- ttH btag definitions
        // isLooseBTag  =  (deepCSVb+deepCSVbb > 0.1522);
        // isMediumBTag =  (deepCSVb+deepCSVbb > 0.4941);
        // isTightBTag = (deepCSVb+deepCSVbb > 0.8001);

        //NB : btag def must be consistent w/ SelectBjets()
        for(unsigned int ijet=0; ijet < vJetLoose->size() ; ijet++)
        {
            if(use_ttH_selections)
            {
                if(fabs(vJetLoose->at(ijet).eta) > 2.4) {continue;} // -- don't relax eta cut (ttH2017)
            }
            else
            {
                if(fabs(vJetLoose->at(ijet).eta) > 4.7) {continue;} // -- relax eta cut (tHq2017)

                if(fabs(vJetLoose->at(ijet).eta) > 2.7 && fabs(vJetLoose->at(ijet).eta) < 3.0) {is_hasJetTransitionRegion = 1;} //Noisy region, might want to veto these events
            }

            //--- B-tagged jets
            if(vJetLoose->at(ijet).isLooseBTag)
            {
                vLooseBTagJets.push_back(vJetLoose->at(ijet)); nLooseBJets++;

                if(vJetLoose->at(ijet).isMediumBTag) {nMediumBJets++;}
            }
            else
            {
                vLightJets.push_back(vJetLoose->at(ijet)); nLightJets++;

                if(fabs(vJetLoose->at(ijet).eta ) < 2.4)
                {
                    vLightJets_FwdPtCut.push_back(vJetLoose->at(ijet)); nLightJets_Fwd40++;
                }
                else if(vJetLoose->at(ijet).pt > 40 && fabs(vJetLoose->at(ijet).eta ) > 2.4)
                {
                    vLightJets_FwdPtCut.push_back(vJetLoose->at(ijet)); nLightJets_Fwd40++;
                }
            }
        }

        //--------------------------------------------
        //  ####  #####  #####  ###### #####  # #    #  ####
        // #    # #    # #    # #      #    # # ##   # #    #
        // #    # #    # #    # #####  #    # # # #  # #
        // #    # #####  #    # #      #####  # #  # # #  ###
        // #    # #   #  #    # #      #   #  # #   ## #    #
        //  ####  #    # #####  ###### #    # # #    #  ####
        //--------------------------------------------

        //Sort object vectors by Pt
        //-- moved from 'SortingLeptonPt' to 'SortingLeptonConePt'
        std::sort(vLeptonLoose.begin(), vLeptonLoose.end(), SortingLeptonConePt);
        std::sort(vLeptonLoose.begin(), vLeptonLoose.end(), SortingLeptonConePt);
        std::sort(vLeptonFakeable.begin(), vLeptonFakeable.end(), SortingLeptonConePt);
        std::sort(vLeptonLoose.begin(), vLeptonLoose.end(), SortingLeptonConePt);
        std::sort(vLeptonTight.begin(), vLeptonTight.end(), SortingLeptonConePt);

        nTightLep = vLeptonTight.size();
        nFakeableLep = vLeptonFakeable.size();



        //--------------------------------------------
        // ##### #####  #  ####   ####  ###### #####
        //   #   #    # # #    # #    # #      #    #
        //   #   #    # # #      #      #####  #    #
        //   #   #####  # #  ### #  ### #      #####
        //   #   #   #  # #    # #    # #      #   #
        //   #   #    # #  ####   ####  ###### #    #
        //--------------------------------------------


        bool TRIGm   = vEvent->at(0).trig_m  ;
        bool TRIGe   = vEvent->at(0).trig_e  ;
        bool TRIGee  = vEvent->at(0).trig_ee ;
        bool TRIGmm  = vEvent->at(0).trig_mm ;
        bool TRIGem  = vEvent->at(0).trig_em ;
        bool TRIGeee = vEvent->at(0).trig_eee;
        bool TRIGmme = vEvent->at(0).trig_emm;
        bool TRIGeem = vEvent->at(0).trig_eem;
        bool TRIGmmm = vEvent->at(0).trig_mmm;

        //Previous
        // bool E = false, M = false, EE = false, MM = false, EM = false;
        // if(TRIGme || TRIGeem || TRIGmme ) {EM = true;}
        // if(TRIGmm || TRIGmmm)                       {MM = true;}
        // if(TRIGee || TRIGeee)	                    {EE = true;}
        // if(TRIGm)                                   {M  = true;}
        // if(TRIGe)                                   {E  = true;}

        //Updated
        bool E = false, M = false, EE = false, MM = false, EM = false, MME = false, EEM = false, MMM = false, EEE = false;
        // if(TRIGmmm)                                 {MMM = true;}
        // if(TRIGeee)                                 {EEE = true;}
        // if(TRIGmme)                                 {MME = true;}
        // if(TRIGeem)                                 {EEM = true;}
        // if(TRIGmm)                                  {MM = true;}
        // if(TRIGee)          	                    {EE = true;}
        // if(TRIGem)          	                    {EM = true;}
        // if(TRIGm)                                   {M  = true;}
        // if(TRIGe)                                   {E  = true;}

        if(TRIGeee || TRIGee)               {EE = true;}
        if(TRIGmmm || TRIGmm)               {MM = true;}
        if(TRIGmme || TRIGeem || TRIGem)    {EM = true;}
        if(TRIGm)                           {M = true;}
        if(TRIGe)                           {E = true;}

        bool emdataset = _sampleName.Contains("MuonE");
        bool mmdataset = _sampleName.Contains("DoubleM");
        bool eedataset = _sampleName.Contains("DoubleE");
        bool mdataset  = _sampleName.Contains("SingleM");
        bool edataset  = _sampleName.Contains("SingleE");

        is_trigger     = false;
        is_trigger_ttH = false;

        //Hypothesis : events satisfying EE will end up in eedataset, etc.
        if(_isdata) //For data
        {
            if      ( EE  &&                               (eedataset) ) {is_trigger = true;}
            else if ( !EE && MM  &&                        (mmdataset) ) {is_trigger = true;}
            else if ( !EE && !MM && EM  &&                 (emdataset) ) {is_trigger = true;}
            else if ( !EE && !MM && !EM && E  &&           (edataset ) ) {is_trigger = true;}
            else if ( !EE && !MM && !EM && !E && M &&      (mdataset ) ) {is_trigger = true;}
        }
        else //For MC
        {
            if(EM || MM || EE || M || E) {is_trigger = true;}
        }


        //Kirill's 2lss trigger logic
        bool passTrigMC = (vEvent->at(0).trig_e || vEvent->at(0).trig_m || vEvent->at(0).trig_ee || vEvent->at(0).trig_mm || vEvent->at(0).trig_em);
        bool passTrigData =
        (edataset && vEvent->at(0).trig_e && !vEvent->at(0).trig_ee && !vEvent->at(0).trig_m && !vEvent->at(0).trig_mm && !vEvent->at(0).trig_em) ||
        (eedataset && vEvent->at(0).trig_ee && !vEvent->at(0).trig_mm && !vEvent->at(0).trig_em) ||
        (mdataset && vEvent->at(0).trig_m && !vEvent->at(0).trig_ee && !vEvent->at(0).trig_e && !vEvent->at(0).trig_mm && !vEvent->at(0).trig_em) ||
        (mmdataset && vEvent->at(0).trig_mm && !vEvent->at(0).trig_ee && !vEvent->at(0).trig_em) ||
        (emdataset && vEvent->at(0).trig_em && !vEvent->at(0).trig_ee && !vEvent->at(0).trig_mm);

        is_trigger_ttH = (!_isdata && passTrigMC) || (_isdata && passTrigData);

        //The event will be rejected in sel. func., but counted
        if(!is_trigger)
        {
            weight = 0;
            continue; //Skip data events failing trigger
        }


        // if(!is_trigger && !is_trigger_ttH)
        // {
        //     weight = 0;
        //     continue; //Skip data events failing trigger
        // }




        //--------------------------------------------
        // #    #   ##   #####  #   ##   #####  #      ######  ####
        // #    #  #  #  #    # #  #  #  #    # #      #      #
        // #    # #    # #    # # #    # #####  #      #####   ####
        // #    # ###### #####  # ###### #    # #      #           #
        //  #  #  #    # #   #  # #    # #    # #      #      #    #
        //   ##   #    # #    # # #    # #####  ###### ######  ####
        //--------------------------------------------

        if(vLeptonLoose.size() < 2) {continue;} //Else, infinite loop when checking Zveto

        // ##########
        // # OSSF pair - Z veto #
        // ##########
        pass_Zveto = false;
        pass_Zveto_ee = false;
        pass_cleanup = false;

        nSFOS = 0;
        mllll = -1;
        inv_mll = 0;

        //-- Taken from Kirill (NTP code, Sync.cxx)
        float mll_min = 10E+10;
        float mll_z_min = 10E+10;
        for(int il=0; il<vLeptonLoose.size(); il++)
        {
            TLorentzVector l1;
            l1.SetPtEtaPhiE(vLeptonLoose.at(il).pt, vLeptonLoose.at(il).eta, vLeptonLoose.at(il).phi, vLeptonLoose.at(il).E);

            for(int ill=il+1; ill<vLeptonLoose.size(); ill++)
            {
                TLorentzVector l2;
                l2.SetPtEtaPhiE(vLeptonLoose.at(ill).pt, vLeptonLoose.at(ill).eta, vLeptonLoose.at(ill).phi, vLeptonLoose.at(ill).E);

                float mll = (l1+l2).M();
                float mll_z = fabs((l1+l2).M()-91.2);

                bool SFOS = (vLeptonLoose.at(il).id == -vLeptonLoose.at(ill).id);

                if(SFOS) {nSFOS++;}

                if(mll < mll_min) {mll_min = mll;}
                if(mll_z < mll_z_min && SFOS)
                {
                    mll_z_min = mll_z;
                    inv_mll = (l1+l2).M();
                }

                for(int illl=ill+1; illl<vLeptonLoose.size(); illl++)
                {
                    TLorentzVector l3;
                    l3.SetPtEtaPhiE(vLeptonLoose.at(illl).pt, vLeptonLoose.at(illl).eta, vLeptonLoose.at(illl).phi, vLeptonLoose.at(illl).E);

                    for(int illll=illl+1; illll<vLeptonLoose.size(); illll++)
                    {
                        TLorentzVector l4;
                        l4.SetPtEtaPhiE(vLeptonLoose.at(illll).pt, vLeptonLoose.at(illll).eta, vLeptonLoose.at(illll).phi, vLeptonLoose.at(illll).E);

                        bool SFOS34 = (vLeptonLoose.at(illl).id == -vLeptonLoose.at(illll).id);

                        if(!SFOS34) {continue;}

                        float mllll = (l1+l2+l3+l4).M();
                    }
                }
            }
        }
        if(mll_min > 12) {pass_cleanup = true;}
        if(mll_z_min > 10) {pass_Zveto = true;}


        //--- Taken from Kirill (NTP code, Sync.cxx)
        float mee_z_min = 10E+10;
        for(int ie=0; ie<vElectronLoose->size(); ie++)
        {
            TLorentzVector e1;
            e1.SetPtEtaPhiE(vElectronLoose->at(ie).pt, vElectronLoose->at(ie).eta, vElectronLoose->at(ie).phi, vElectronLoose->at(ie).E);

            for(int iee=ie+1; iee<vElectronLoose->size(); iee++)
            {
                TLorentzVector e2;
                e2.SetPtEtaPhiE(vElectronLoose->at(iee).pt, vElectronLoose->at(iee).eta, vElectronLoose->at(iee).phi, vElectronLoose->at(iee).E);

                float mee = fabs((e1+e2).M()-91.2);
                if(mee < mee_z_min) {mee_z_min = mee;}
            }
        }
        if(mee_z_min > 10) {pass_Zveto_ee = true;}



        //--------------------------------------------
        // #####  ###### #####  #    #  ####     #####  #####  # #    # #####  ####  #    # #####
        // #    # #      #    # #    # #    #    #    # #    # # ##   #   #   #    # #    #   #
        // #    # #####  #####  #    # #         #    # #    # # # #  #   #   #    # #    #   #
        // #    # #      #    # #    # #  ###    #####  #####  # #  # #   #   #    # #    #   #
        // #    # #      #    # #    # #    #    #      #   #  # #   ##   #   #    # #    #   #
        // #####  ###### #####   ####   ####     #      #    # # #    #   #    ####   ####    #
        //--------------------------------------------

        // bool debug_printout = false;
        if(DEBUG)
        {
            cout<<endl<<endl<<BOLD(FBLU("------------ EVENT "<<event_id<<" -----------"))<<endl;

            cout<<BOLD(FBLU("--- Loose Btag Jets : "))<<endl;
            for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
            {
                cout<<FBLU(""<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt<<" / eta = "<<vLooseBTagJets.at(ijet).eta<<" / phi = "<<vLooseBTagJets.at(ijet).phi<<" / CSV = "<<vLooseBTagJets.at(ijet).DeepCSVbtag<<"")<<endl;
            }
            cout<<endl;
            cout<<BOLD(FCYN("### Light Jets : "))<<endl;
            for(int ijet=0; ijet<vLightJets.size(); ijet++)
            {
                cout<<FCYN("* "<<ijet<<" pT = "<<vLightJets.at(ijet).pt<<" / eta = "<<vLightJets.at(ijet).eta<<" / phi = "<<vLightJets.at(ijet).phi<<" / CSV = "<<vLightJets.at(ijet).DeepCSVbtag<<"")<<endl;
            }
            cout<<endl;
            cout<<BOLD(FYEL("### Fakeable Leptons : "))<<endl<<endl;
            for(int ilep=0; ilep<vLeptonFakeable.size(); ilep++)
            {
                cout<<"-- vLeptonFakeable["<<ilep<<"].id = "<<vLeptonFakeable[ilep].id<<endl;


                TString type = "";
                if(vLeptonFakeable[ilep].id < 0) type = "anti";
                if( abs(vLeptonFakeable[ilep].id) == 11) type+= "Ele";
                else if( abs(vLeptonFakeable[ilep].id) == 13) type+= "Mu";

                cout<<FYEL(<<" "<<ilep<<" "<<type<<": FO ? "<<vLeptonFakeable[ilep].isFakeableTTH<<", Tight ? "<<vLeptonFakeable[ilep].isTightTTH<<" / pT = "<<vLeptonFakeable.at(ilep).pt<<" / eta = "<<vLeptonFakeable.at(ilep).eta<<" / phi = "<<vLeptonFakeable.at(ilep).phi<<"")<<endl;
            }

            // cout<<endl;
            // cout<<BOLD(FYEL("### ttH predefined categories : "))<<endl<<endl;
            // cout<<"is_ttH_2lSS = "<<(bool) is_ttH_2lSS<<endl;
            cout<<BOLD(FBLU("------------------------------------"))<<endl<<endl;
        }


        //--------------------------------------------
        //  ####  ###### #      ######  ####  ##### #  ####  #    #  ####
        // #      #      #      #      #    #   #   # #    # ##   # #
        //  ####  #####  #      #####  #        #   # #    # # #  #  ####
        //      # #      #      #      #        #   # #    # #  # #      #
        // #    # #      #      #      #    #   #   # #    # #   ## #    #
        //  ####  ###### ###### ######  ####    #   #  ####  #    #  ####
        //--------------------------------------------



        //Regions
        //NB : all regions not entirely orthogonal, since an event can be e.g. both 3l_SR_Fake and 2lSS_ttWctrl !

        bool is_3l_custom_cat = false; //Check if 1 of my custom categories is true
        bool is_2l_custom_cat = false; //Check if 1 of my custom categories is true

//--------------------------------------------
        if(ThreeLeptonSelection_THQ3l(jentry) ) {is_3l_custom_cat = true;}
        if(TwoLeptonSelection_THQ2lSSl(jentry) ) {is_2l_custom_cat = true;}

        if(Sample_isUsed_forTraining() ) //Looser training selection (training samples only)
        {
            if(ThreeLeptonSelection_THQ3l_TrainingSelection(jentry) ) {is_3l_custom_cat = true;}
            if(TwoLeptonSelection_THQ2lSS_TrainingSelection(jentry) ) {is_2l_custom_cat = true;}
        }

//--------------------------------------------

        //If event satisfies any 2l-cat, will be saved as such
        //Creates a bias if also saved as a 3l-cat event (e.g. lepCharge will be wrong)... but need avoid double-count
        if(is_2l_custom_cat)
        {
            Compute_Variables("2l");
            Apply_ScaleFactors(2);
            fillOutputTree();
        }
        else if(is_3l_custom_cat)
        {
            Compute_Variables("3l");
            Apply_ScaleFactors(3);
            fillOutputTree();
        }
        else if(nlep_tth != 0)
        {
            if(nlep_tth == 3) {Compute_Variables("3l"); Apply_ScaleFactors(3);}
            else if(nlep_tth == 2) {Compute_Variables("2l"); Apply_ScaleFactors(2);}
            fillOutputTree();
        }
    }




	// ofstream file_out("cutflow.txt");
    // for(int i=0; i<v_cutflow.size(); i++)
    // {
	// 	file_out<<"v_cutflow "<<i<<" = "<<v_cutflow[i]<<endl;
    // }


	outputfile->cd();
    tOutput->Write();
}





























//--------------------------------------------
//  ######  ######## ##              ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
// ##    ## ##       ##              ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
// ##       ##       ##              ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
//  ######  ######   ##              ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
//       ## ##       ##              ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
// ##    ## ##       ##       ###    ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
//  ######  ######## ######## ###    ##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######
//--------------------------------------------






//--------------------------------------------
//  #####
// #     # #
//       # #
//  #####  #
//       # #
// #     # #
//  #####  ######
//--------------------------------------------

/**
 * 3l event selection
 */
bool tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l(int evt)
{
	if(DEBUG) {cout<<FYEL("-- ThreeLeptonSelection_THQ3l")<<endl;}

    InitVariables();
	v_cutflow[0]++;

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return 0;}

	//At least 3 FO leptons
    if(vLeptonFakeable.size() < 3) {return 0;}
    v_cutflow[1]++;

    //pT cuts
    if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15 || vLeptonFakeable.at(2).conept < 10) {return 0;}
    v_cutflow[2]++;

    //No pair of loose leptons with mll < 12
    if(!pass_cleanup) {return 0;}
    v_cutflow[3]++;

    //metLD cuts
    if(vJetLoose->size() < 4)
    {
        // met LD : nJet25 >= 4 || met_pt*0.6 + mhtJet25*0.4 > 30 + 15*(mZ1 > 0) //from marco's plots...
        if(nSFOS > 0 && metLD < 0.3) {return 0;}
        else if(!nSFOS && metLD < 0.2) {return 0;}
    }
    v_cutflow[4]++;

    if( fabs(vLeptonFakeable.at(0).charge + vLeptonFakeable.at(1).charge + vLeptonFakeable.at(2).charge) != 1) {return 0;}
    v_cutflow[5]++;

    if(vJetLoose->size() < 2) {return 0;}

    v_cutflow[6]++;

//--------------------------------------------

    if(!ThreeLeptonSelection_THQ3l_Regions(evt) ) {return 0;}

    return 1;
}



/**
 * 3l categories definitions
 */
bool tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l_Regions(int evt)
{
    bool pass_isTight = false;
    bool pass_MCMatching = false;
    bool pass_m4l = false;
    bool has_lessThan4Tight = false;
    bool has_noIsoTau = false;
    bool pass_nBJet_ttH = false;
    bool pass_njet_tHq = false;
    bool pass_tighter_ptCuts = false;
    if(vLeptonFakeable.at(0).isTightTTH && vLeptonFakeable.at(1).isTightTTH && vLeptonFakeable.at(2).isTightTTH) {pass_isTight = true;}
    if(_isdata || (vLeptonFakeable.at(0).hasMCMatch && vLeptonFakeable.at(1).hasMCMatch && vLeptonFakeable.at(2).hasMCMatch) ) {pass_MCMatching = true;}
    if(mllll > 140 || mllll < 0) {pass_m4l = true;}
    if(vLeptonTight.size() < 4) {has_lessThan4Tight = true;}
    if(vTauLoose->size() == 0) {has_noIsoTau = true;} //NB : tight or loose taus ?? //FIXME
    bool pass_basic_SR_cuts(pass_m4l && has_lessThan4Tight && has_noIsoTau);
    if(nLooseBJets >= 2 || nMediumBJets >= 1) {pass_nBJet_ttH = true;} //at least 2 b-loose or 1 b-medium
    if(nMediumBJets > 0 && nLightJets > 0) {pass_njet_tHq = true;} //at least 1 b-medium & 1 fwd (light) jet
    if(vLeptonFakeable.at(0).conept > 25 && vLeptonFakeable.at(1).conept > 15 && vLeptonFakeable.at(2).conept > 15) {pass_tighter_ptCuts = true;}

    //Determine category
    if(use_ttH_selections && pass_Zveto && pass_tighter_ptCuts && pass_basic_SR_cuts && pass_nBJet_ttH) {is_tHq_3l = 1;} //ttH2017 selections
    else if(!use_ttH_selections && pass_Zveto && pass_basic_SR_cuts && pass_njet_tHq) {is_tHq_3l = 1;}
    if(is_tHq_3l && pass_isTight && pass_MCMatching) {is_tHq_3l_SR = 1;}
    if(is_tHq_3l && !pass_isTight) {is_tHq_3l_Fake = 1; Apply_FakeRate_Weight("3l");}

    if(!pass_Zveto && pass_basic_SR_cuts && pass_nBJet_ttH && nLooseBJets>=2 && nMediumBJets>=1) {is_tHq_ttZctrl = 1;}
    if(is_tHq_ttZctrl && pass_isTight && pass_MCMatching) {is_tHq_ttZctrl_SR = 1;}
    if(is_tHq_ttZctrl && !pass_isTight) {is_tHq_ttZctrl_Fake = 1; Apply_FakeRate_Weight("3l");}

    if(!pass_Zveto && pass_tighter_ptCuts && pass_basic_SR_cuts && vJetLoose->size() >= 2 && !pass_nBJet_ttH) {is_tHq_WZctrl = 1;}
    if(is_tHq_WZctrl && pass_isTight && pass_MCMatching) {is_tHq_WZctrl_SR = 1;} //removed basic cuts for now
    if(is_tHq_WZctrl && !pass_isTight) {is_tHq_WZctrl_Fake = 1; Apply_FakeRate_Weight("3l");}

//--------------------------------------------
    if(is_tHq_3l || is_tHq_ttZctrl || is_tHq_WZctrl) {return true;}

    return false;
}















//--------------------------------------------
//  #####
// #     # #
//       # #
//  #####  #
// #       #
// #       #
// ####### ######
//--------------------------------------------

/**
 * 2l selection
 */
bool tHqMultileponAnalysis::TwoLeptonSelection_THQ2lSSl(int evt)
{
	if(DEBUG) {cout<<FYEL("-- TwoLeptonSelection_THQ2lSSl")<<endl;}

    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return 0;}

    //At least 2 leptons
    if(vLeptonFakeable.size() < 2) {return 0;}
    if(vLeptonTight.size() > 2) {return 0;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    if(!vLeptonFakeable.at(0).passTightCharge || !vLeptonFakeable.at(1).passTightCharge ) {return 0;}

    //pT cuts
    if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15) {return 0;}

    //No pair of loose leptons with mll < 12
    if(!pass_cleanup) {return 0;}

    //No loose tau
    if(vTauLoose->size() > 0) {return 0;} //FIXME or tau tight ??

    //No ee pair with mll-mZ < 10
    if(!pass_Zveto_ee) {return 0;}

    //metLD cuts
    if(fabs(vLeptonFakeable.at(0).id) == 11 && fabs(vLeptonFakeable.at(1).id) == 11)
    {
        // metLDee : met_pt*0.6 + mhtJet25*0.4 > 30 //from marco's plots
        if(metLD < 0.2) {return 0;}
    }

    if(vJetLoose->size() < 2) {return 0;}

//--------------------------------------------

    if(!TwoLeptonSelection_THQ2lSS_Regions(evt) ) {return 0;} //Check categories

    return 1;
}





/**
 * 2l categories selection
 */
bool tHqMultileponAnalysis::TwoLeptonSelection_THQ2lSS_Regions(int evt)
{
    bool pass_isTight = false;
    bool pass_MCMatching = false;
    bool is_SS_pair = false;
    bool contains_ele = false;
    bool pass_nBJet_ttH = false;
    bool pass_njet_tHq = false;
    if(vLeptonFakeable.at(0).isTightTTH && vLeptonFakeable.at(1).isTightTTH) {pass_isTight = true;}
    if(_isdata || (vLeptonFakeable.at(0).hasMCMatch && vLeptonFakeable.at(1).hasMCMatch) ) {pass_MCMatching = true;}
    if(vLeptonFakeable.at(0).charge * vLeptonFakeable.at(1).charge > 0) {is_SS_pair = true;}
    if( fabs(vLeptonFakeable.at(0).id) == 11 || fabs(vLeptonFakeable.at(1).id) == 11 ) {contains_ele = true;}
    if(vJetLoose->size() >= 3 && (nLooseBJets >= 2 || nMediumBJets >= 1) ) {pass_nBJet_ttH = true;}
    if(vJetLoose->size() >= 2 && nMediumBJets > 0 && nLightJets > 0)  {pass_njet_tHq = true;}


    //Determine category
    if(use_ttH_selections && vJetLoose->size() >= 4 && pass_nBJet_ttH) {is_tHq_2lSS = 1;} //ttH2017 cuts
    else if(!use_ttH_selections && pass_njet_tHq) {is_tHq_2lSS = 1;} //tHq2017 cuts
    if(is_tHq_2lSS && is_SS_pair && pass_isTight && pass_MCMatching) {is_tHq_2lSS_SR = 1;}
    if(is_tHq_2lSS && is_SS_pair && !pass_isTight) {is_tHq_2lSS_Fake = 1; Apply_FakeRate_Weight("2l");}
    if(is_tHq_2lSS && !is_SS_pair && contains_ele && pass_isTight && pass_MCMatching) {is_tHq_2lSS_Flip = 1; Apply_FlipRate_Weight();}
    if(is_tHq_2lSS && !is_tHq_2lSS_Fake && !is_tHq_2lSS_Flip && contains_ele && is_GammaConv_Event() ) {is_tHq_2lSS_GammaConv = 1;} //FIXME

    if(pass_nBJet_ttH && vJetLoose->size() == 3) {is_tHq_ttWctrl = 1;}
	if(is_tHq_ttWctrl && is_SS_pair && pass_isTight && pass_MCMatching) {is_tHq_ttWctrl_SR = 1;}
    if(is_tHq_ttWctrl && is_SS_pair && !pass_isTight) {is_tHq_ttWctrl_Fake = 1; Apply_FakeRate_Weight("2l");}
    if(is_tHq_ttWctrl && !is_SS_pair && contains_ele && pass_isTight && pass_MCMatching) {is_tHq_ttWctrl_Flip = 1; Apply_FlipRate_Weight();}
    if(is_tHq_ttWctrl && !is_tHq_ttWctrl_Fake && !is_tHq_ttWctrl_Flip && contains_ele && is_GammaConv_Event() ) {is_tHq_ttWctrl_GammaConv = 1;} //FIXME

//--------------------------------------------
	if(is_tHq_2lSS || is_tHq_ttWctrl)
    {
        return true;
    }

    return false;
}











//--------------------------------------------
// ######## ########     ###    #### ##    ## #### ##    ##  ######       ######  ######## ##
//    ##    ##     ##   ## ##    ##  ###   ##  ##  ###   ## ##    ##     ##    ## ##       ##
//    ##    ##     ##  ##   ##   ##  ####  ##  ##  ####  ## ##           ##       ##       ##
//    ##    ########  ##     ##  ##  ## ## ##  ##  ## ## ## ##   ####     ######  ######   ##
//    ##    ##   ##   #########  ##  ##  ####  ##  ##  #### ##    ##           ## ##       ##
//    ##    ##    ##  ##     ##  ##  ##   ###  ##  ##   ### ##    ##     ##    ## ##       ##       ###
//    ##    ##     ## ##     ## #### ##    ## #### ##    ##  ######       ######  ######## ######## ###
//--------------------------------------------

/**
 * Selection function implementing the ttH2017 training selection, for 3l cat.
 */
bool tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l_TrainingSelection(int evt)
{
    if(DEBUG) {cout<<FYEL("-- ThreeLeptonSelection_THQ3l_TrainingSelection")<<endl;}

    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return 0;}

	//At least 3 FO leptons
    if(vLeptonFakeable.size() < 3) {return 0;}

    if(!pass_Zveto) {return 0;}

    if(use_ttH_selections) //ttH training sel
    {
        if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15 || vLeptonFakeable.at(2).conept < 15) {return 0;}

        //metLD cuts //maybe a bit different !
        // allcuts += "(nJet25_Recl >= 4 || (met_pt*0.00397 + mhtJet25_Recl*0.00265 - 0.184 > 0.0 + 0.1*(mZ1_Recl > 0)))"
        if(vJetLoose->size() < 4)
        {
            if(nSFOS > 0 && metLD < 0.3) {return 0;}
            else if(!nSFOS && metLD < 0.2) {return 0;}
        }

        if(nLooseBJets < 2) {return 0;}
    }
    else //Try even looser cuts for tHq
    {
        if(vLeptonFakeable.at(0).conept < 20 || vLeptonFakeable.at(1).conept < 10 || vLeptonFakeable.at(2).conept < 10) {return 0;}

        if(!nLooseBJets || !nLightJets) {return 0;}
    }


    //--------------------------------------------
    is_tHq_3l_Training = 1;

    return 1;
}


/**
 *  * Selection function implementing the ttH2017 training selection, for 2lSS cat.
 */
bool tHqMultileponAnalysis::TwoLeptonSelection_THQ2lSS_TrainingSelection(int evt)
{
    if(DEBUG) {cout<<FYEL("-- TwoLeptonSelection_THQ2lSS_TrainingSelection")<<endl;}

    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return 0;}

    //At least 3 FO leptons
    if(vLeptonFakeable.size() < 2) {return 0;}

    //At most 2 tight leptons
    if(vLeptonTight.size() > 2) {return 0;}

    if(vLeptonFakeable.at(0).charge != vLeptonFakeable.at(1).charge) {return 0;}

    //No ee pair with mll-mZ < 10
    if(!pass_Zveto_ee) {return 0;}

    if(use_ttH_selections) //ttH training sel
    {
        if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15) {return 0;}

        if(vJetLoose->size() < 4) {return 0;}

        if(nLooseBJets < 2 && nMediumBJets == 0) {return 0;}
    }
    else //Try even looser cuts for tHq
    {
        if(vLeptonFakeable.at(0).conept < 20 || vLeptonFakeable.at(1).conept < 10) {return 0;}

        if(!nLooseBJets || !nLightJets) {return 0;}
    }


    //--------------------------------------------
    is_tHq_2lSS_Training = 1;

    return 1;
}












//--------------------------------------------
//  ######      ###    ##     ## ##     ##    ###        ######   #######  ##    ## ##     ##
// ##    ##    ## ##   ###   ### ###   ###   ## ##      ##    ## ##     ## ###   ## ##     ##
// ##         ##   ##  #### #### #### ####  ##   ##     ##       ##     ## ####  ## ##     ##
// ##   #### ##     ## ## ### ## ## ### ## ##     ##    ##       ##     ## ## ## ## ##     ##
// ##    ##  ######### ##     ## ##     ## #########    ##       ##     ## ##  ####  ##   ##
// ##    ##  ##     ## ##     ## ##     ## ##     ##    ##    ## ##     ## ##   ###   ## ##   ###
//  ######   ##     ## ##     ## ##     ## ##     ##     ######   #######  ##    ##    ###    ###
//--------------------------------------------

bool tHqMultileponAnalysis::is_GammaConv_Event()
{
    if(_isdata) {return false;}
    if(is_tHq_2lSS_Fake || is_tHq_2lSS_Flip) {return false;}
    if(_isdata) {return false;}
    if(vLeptonFakeable.size() < 2) {return false;}

    if(abs(vLeptonFakeable[0].id) == 11 && is_Electron_Matched_To_Gamma(0) ) {return true;}
    if(abs(vLeptonFakeable[1].id) == 11 && is_Electron_Matched_To_Gamma(1) ) {return true;}

    return false;
}


bool tHqMultileponAnalysis::is_Electron_Matched_To_Gamma(int lepRec_idx)
{
    if(_isdata) {return false;}
    if( abs(vLeptonFakeable[lepRec_idx].id) != 11) {return -1;} //Check ele only

    int genLep_idx = Ele_To_Lepton_Matching(lepRec_idx);
    int genPho_idx = Ele_To_Gamma_Matching(lepRec_idx);

    if(genLep_idx < 0 && genPho_idx > 0) {return true;}
    else if(genLep_idx > 0 && genPho_idx > 0)
    {
        double dlep = Get_Distance(vLeptonFakeable[lepRec_idx].pt, vLeptonFakeable[lepRec_idx].eta, vLeptonFakeable[lepRec_idx].phi, vTruth->at(0).gen_pt.at(genLep_idx), vTruth->at(0).gen_eta.at(genLep_idx), vTruth->at(0).gen_phi.at(genLep_idx) );
        double dpho = Get_Distance(vLeptonFakeable[lepRec_idx].pt, vLeptonFakeable[lepRec_idx].eta, vLeptonFakeable[lepRec_idx].phi, vTruth->at(0).gen_pt.at(genPho_idx), vTruth->at(0).gen_eta.at(genPho_idx), vTruth->at(0).gen_phi.at(genPho_idx) );

        if(dpho < dlep) {return true;}
    }

    return false;
}



int tHqMultileponAnalysis::Ele_To_Gamma_Matching(int lepRec_idx)
{
    if(_isdata) {return -1;}
    if( abs(vLeptonFakeable[lepRec_idx].id) != 11) {return -1;} //Check ele only

    double dR_min = 999;
    int idx_bestMatch = -1;

    //Loop on all truth objects
    for(int itruth = 0; itruth < vTruth->at(0).gen_id.size(); itruth++)
    {
        if(abs(vTruth->at(0).gen_id.at(itruth)) != 22) {continue;} //Check only gen photons
        if(vTruth->at(0).gen_status.at(itruth) != 1) {continue;} //Check that gen photon belongs to final state
        if(vTruth->at(0).gen_pt.at(itruth) < 1) {continue;} //Check that gen photon has pT>1

        if(vLeptonFakeable[lepRec_idx].conept < 0.3 * vTruth->at(0).gen_pt.at(itruth) || vLeptonFakeable[lepRec_idx].conept > 1.5 * vTruth->at(0).gen_pt.at(itruth) ) {continue;}

        double dR = GetDeltaR(vLeptonFakeable[lepRec_idx].eta, vLeptonFakeable[lepRec_idx].phi, vTruth->at(0).gen_eta.at(itruth), vTruth->at(0).gen_phi.at(itruth) );

        if(dR > 0.3) {continue;}

        if(dR < dR_min)
        {
            dR_min = dR;
            idx_bestMatch = itruth;
        }
    }

    return idx_bestMatch;
}

int tHqMultileponAnalysis::Ele_To_Lepton_Matching(int lepRec_idx)
{
    if(_isdata) {return -1;}
    if( abs(vLeptonFakeable[lepRec_idx].id) != 11) {return -1;} //Check ele only

    double dR_min = 999;
    int idx_bestMatch = -1;

    //Loop on all truth objects
    for(int itruth = 0; itruth < vTruth->at(0).gen_id.size(); itruth++)
    {
        if(abs(vTruth->at(0).gen_id.at(itruth)) != abs(vLeptonFakeable[lepRec_idx].id) ) {continue;} //Check gen leptons with same ID

        double dR = GetDeltaR(vLeptonFakeable[lepRec_idx].eta, vLeptonFakeable[lepRec_idx].phi, vTruth->at(0).gen_eta.at(itruth), vTruth->at(0).gen_phi.at(itruth) );

        if(dR > 0.3) {continue;}

        if(dR < dR_min)
        {
            dR_min = dR;
            idx_bestMatch = itruth;
        }

        // if(dR < 0.3) {return itruth;}
        // else if(vLeptonFakeable[lepRec_idx].pt < 10 && abs(vLeptonFakeable[lepRec_idx].id) == 13 && vLeptonFakeable[lepRec_idx].id != vTruth->at(0).gen_id.at(itruth) ) {return -1;}
        // else if(dR < 0.7) {return itruth;}
        // else if(std::min(vLeptonFakeable[lepRec_idx].pt, vTruth->at(0).gen_pt.at(itruth)) / std::max(vLeptonFakeable[lepRec_idx].pt, vTruth->at(0).gen_pt.at(itruth)) < 0.3) {return -1;}
        // else if(dR < 1.2) {return itruth;}
    }

    return idx_bestMatch;
}









































//--------------------------------------------
// ######## #### ##       ##          ######## ######## ########  ######## ########
// ##        ##  ##       ##             ##       ##    ##     ## ##       ##
// ##        ##  ##       ##             ##       ##    ##     ## ##       ##
// ######    ##  ##       ##             ##       ##    ########  ######   ######
// ##        ##  ##       ##             ##       ##    ##   ##   ##       ##
// ##        ##  ##       ##             ##       ##    ##    ##  ##       ##
// ##       #### ######## ########       ##       ##    ##     ## ######## ########

//   ##
//   ##
// ######
//   ##
//   ##

// ##     ## ######## ##     ##    ##     ##    ###    ########   ######
// ###   ### ##       ###   ###    ##     ##   ## ##   ##     ## ##    ##
// #### #### ##       #### ####    ##     ##  ##   ##  ##     ## ##
// ## ### ## ######   ## ### ##    ##     ## ##     ## ########   ######
// ##     ## ##       ##     ##     ##   ##  ######### ##   ##         ##
// ##     ## ##       ##     ##      ## ##   ##     ## ##    ##  ##    ##
// ##     ## ######## ##     ##       ###    ##     ## ##     ##  ######
//--------------------------------------------





/**
 * Fill output tree. Compute lot of variables, needed either for MVA analysis of for MEM computation
 */
void tHqMultileponAnalysis::fillOutputTree()
{
	InitTLorentzVectors(); //Re-init the MEM inputs at each call of the function

    //--------------------------------------------
    // #####       # ###### #####  ####
    // #    #      # #        #   #
    // #####       # #####    #    ####
    // #    #      # #        #        #
    // #    # #    # #        #   #    #
    // #####   ####  ######   #    ####
    //--------------------------------------------

    //Select 2 b-jets (highest CSV)
    bool doSelectOnlyBjets = true;
    // TLorentzVector Bjet1, Bjet2;
    multilepton_Bjet1_Id = -999; multilepton_Bjet2_Id = -999;
    int ib1=-1, ib2=-1;
    SelectBjets(ib1, ib2, doSelectOnlyBjets);
    if (ib1!=-1)
    {
        // Bjet1.SetPtEtaPhiE(vJetLoose->at(ib1).pt, vJetLoose->at(ib1).eta, vJetLoose->at(ib1).phi, vJetLoose->at(ib1).E);
        // FillJetInfoOutputTree(&multilepton_Bjet1_Id, 5, &multilepton_Bjet1_P4, Bjet1, &multilepton_Bjet1_CSV, vJetLoose->at(ib1).CSVv2, &multilepton_Bjet1_JEC_Up, &multilepton_Bjet1_JEC_Down, vJetLoose->at(ib1).JES_uncert(), &multilepton_Bjet1_JER_Up, &multilepton_Bjet1_JER_Down, vJetLoose->at(ib1).pt_JER(), vJetLoose->at(ib1).pt_JER_up(), vJetLoose->at(ib1).pt_JER_down());
        multilepton_Bjet1_Id = 5;
        multilepton_Bjet1_CSV = vJetLoose->at(ib1).DeepCSVbtag;
        multilepton_Bjet1_P4.SetPtEtaPhiE(vJetLoose->at(ib1).pt, vJetLoose->at(ib1).eta, vJetLoose->at(ib1).phi, vJetLoose->at(ib1).E);

    }
    if (ib2!=-1)
    {
        // Bjet2.SetPtEtaPhiE(vJetLoose->at(ib2).pt, vJetLoose->at(ib2).eta, vJetLoose->at(ib2).phi, vJetLoose->at(ib2).E);
        // FillJetInfoOutputTree(&multilepton_Bjet2_Id, 5, &multilepton_Bjet2_P4, Bjet2, &multilepton_Bjet2_CSV, vJetLoose->at(ib2).CSVv2, &multilepton_Bjet2_JEC_Up, &multilepton_Bjet2_JEC_Down, vJetLoose->at(ib2).JES_uncert(), &multilepton_Bjet2_JER_Up, &multilepton_Bjet2_JER_Down, vJetLoose->at(ib2).pt_JER(), vJetLoose->at(ib2).pt_JER_up(), vJetLoose->at(ib2).pt_JER_down());
        multilepton_Bjet2_Id = 5;
        multilepton_Bjet2_CSV = vJetLoose->at(ib2).DeepCSVbtag;
        multilepton_Bjet2_P4.SetPtEtaPhiE(vJetLoose->at(ib2).pt, vJetLoose->at(ib2).eta, vJetLoose->at(ib2).phi, vJetLoose->at(ib2).E);
    }


    //Define jet category
	//2lss
    if (is_tHq_2lSS && ib1!=-1 && ib2!=-1 && vJetLoose->size()-2>=4) catJets = kCat_2lss_2b_4j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2==-1 && vJetLoose->size()-1>=4) catJets = kCat_2lss_1b_4j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2!=-1 && vJetLoose->size()-2==3) catJets = kCat_2lss_2b_3j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2==-1 && vJetLoose->size()-1==3) catJets = kCat_2lss_1b_3j;
	else if (is_tHq_2lSS && ib1!=-1 && ib2!=-1 && vJetLoose->size()-2==2) catJets = kCat_2lss_2b_2j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2!=-1 && vJetLoose->size()-2==1) catJets = kCat_2lss_2b_1j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2==-1 && vJetLoose->size()-1==2) catJets = kCat_2lss_1b_2j;
	else if (is_tHq_2lSS && ib1!=-1 && ib2==-1 && vJetLoose->size()-1==1) catJets = kCat_2lss_1b_1j;
    //4l
    // else if (is4l && ib1!=-1 && ib2!=-1) catJets = kCat_4l_2b;
    // else if (is4l && ib1!=-1 && ib2==-1) catJets = kCat_4l_1b;
    //3l
    else if (is_tHq_3l && ib1!=-1 && ib2!=-1 && vJetLoose->size()-2>=2) catJets = kCat_3l_2b_2j;
    else if (is_tHq_3l && ib1!=-1 && ib2==-1 && vJetLoose->size()-1>=2) catJets = kCat_3l_1b_2j;
    else if (is_tHq_3l && ib1!=-1 && ib2!=-1 && vJetLoose->size()-2==1) catJets = kCat_3l_2b_1j;
    else if (is_tHq_3l && ib1!=-1 && ib2==-1 && vJetLoose->size()-1==1) catJets = kCat_3l_1b_1j;
	else if (is_tHq_3l && ib1!=-1 && ib2!=-1 && vJetLoose->size()-2==0) catJets = kCat_3l_2b_0j;
	else if (is_tHq_3l && ib1!=-1 && ib2==-1 && vJetLoose->size()-1==0) catJets = kCat_3l_1b_0j;
	else if (is_tHq_3l && ib1==-1 && ib2==-1 && vJetLoose->size()==1) catJets = kCat_3l_0b_1j;
	else if (is_tHq_3l && ib1==-1 && ib2==-1 && vJetLoose->size()==0) catJets = kCat_3l_0b_0j;

	else catJets = -1;


    //--------------------------------------------
    // #      ###### #####  #####  ####  #    #  ####
    // #      #      #    #   #   #    # ##   # #
    // #      #####  #    #   #   #    # # #  #  ####
    // #      #      #####    #   #    # #  # #      #
    // #      #      #        #   #    # #   ## #    #
    // ###### ###### #        #    ####  #    #  ####
    //--------------------------------------------

    //-- Select 2 or 3 hardest leptons
    multilepton_Lepton1_Id = -999;
    multilepton_Lepton2_Id = -999;
    multilepton_Lepton3_Id = -999;
    multilepton_Lepton4_Id = -999;


    //Always fill at least 2 leptons
    multilepton_Lepton1_P4 = vLeptonFakeable.at(0).p4;
    multilepton_Lepton1_Id = vLeptonFakeable.at(0).id;
    multilepton_Lepton2_P4 = vLeptonFakeable.at(1).p4;
    multilepton_Lepton2_Id = vLeptonFakeable.at(1).id;


    if(is_tHq_3l) //If corresponds to a 3l event
    {
        multilepton_Lepton3_P4 = vLeptonFakeable.at(2).p4;
        multilepton_Lepton3_Id = vLeptonFakeable.at(2).id;
    }




    // ###############################################################################
    // #                  _       _     _               _                            #
    // #  _ __ ___   __ _| |_ ___| |__ (_)_ __   __ _  | |_ ___     __ _  ___ _ __   #
    // # | '_ ` _ \ / _` | __/ __| '_ \| | '_ \ / _` | | __/ _ \   / _` |/ _ \ '_ \  #
    // # | | | | | | (_| | || (__| | | | | | | | (_| | | || (_) | | (_| |  __/ | | | #
    // # |_| |_| |_|\__,_|\__\___|_| |_|_|_| |_|\__, |  \__\___/   \__, |\___|_| |_| #
    // #                                        |___/              |___/             #
    // #                                                                             #
    // ###############################################################################

    multilepton_Lepton1_Id_Matched = -999; multilepton_Lepton1_Label_Matched = -999; multilepton_Lepton1_DeltaR_Matched = -999;
    multilepton_Lepton2_Id_Matched = -999; multilepton_Lepton2_Label_Matched = -999; multilepton_Lepton2_DeltaR_Matched = -999;
    multilepton_Lepton3_Id_Matched = -999; multilepton_Lepton3_Label_Matched = -999; multilepton_Lepton3_DeltaR_Matched = -999;
    multilepton_Bjet1_Id_Matched = -999; multilepton_Bjet1_Label_Matched = -999; multilepton_Bjet1_DeltaR_Matched = -999;
    multilepton_Bjet2_Id_Matched = -999; multilepton_Bjet2_Label_Matched = -999; multilepton_Bjet2_DeltaR_Matched = -999;

    if ( !_isdata )
    {
        float lep1_dr_gen       = 100.,     lep2_dr_gen     = 100.,     lep3_dr_gen     = 100.,     lep4_dr_gen     = 100. ;
        float jet1_dr_gen       = 100.,     jet2_dr_gen     = 100.;
        float lep1_dr_gen_min   = 100.,     lep2_dr_gen_min = 100.,     lep3_dr_gen_min = 100.,     lep4_dr_gen_min = 100. ;
        float jet1_dr_gen_min   = 100.,     jet2_dr_gen_min = 100.;
        int   lep1_matched      = -1,       lep2_matched    = -1,       lep3_matched   = -1,       lep4_matched    = -1;
        int   jet1_matched      = -1,       jet2_matched    = -1;

        TLorentzVector LeptonX;

        for(unsigned int itruth = 0; itruth < vTruth->at(0).mc_truth_label.size() ; itruth++)
        {
            if( abs(vTruth->at(0).mc_truth_id.at(itruth)) < 18 )
            {
                lep1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vLeptonFakeable.at(0).eta, vLeptonFakeable.at(0).phi );
                if( lep1_dr_gen < lep1_dr_gen_min)
                {
                    lep1_dr_gen_min = lep1_dr_gen;  lep1_matched = itruth;
                }

                lep2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vLeptonFakeable.at(1).eta, vLeptonFakeable.at(1).phi );
                if( lep2_dr_gen < lep2_dr_gen_min)
                {   lep2_dr_gen_min = lep2_dr_gen;  lep2_matched = itruth;  }

                if(vLeptonFakeable.size()>=3)
                {
                    lep3_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vLeptonFakeable.at(2).eta, vLeptonFakeable.at(2).phi );
                    if( lep3_dr_gen < lep3_dr_gen_min)
                    {   lep3_dr_gen_min = lep3_dr_gen;  lep3_matched = itruth;  }
                }

                if(ib1!=-1)
                {
                    jet1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vJetLoose->at(ib1).eta, vJetLoose->at(ib1).phi );
                    if( jet1_dr_gen < jet1_dr_gen_min) {jet1_dr_gen_min = jet1_dr_gen;  jet1_matched = itruth;}
                }

                if(ib2!=-1)
                {
                    jet2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vJetLoose->at(ib2).eta, vJetLoose->at(ib2).phi );
                    if( jet2_dr_gen < jet2_dr_gen_min)
                    {   jet2_dr_gen_min = jet2_dr_gen;  jet2_matched = itruth;  }
                }

            }
        }

        if(lep1_matched >= 0)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt.at(lep1_matched),       vTruth->at(0).mc_truth_eta.at(lep1_matched),
                                    vTruth->at(0).mc_truth_phi.at(lep1_matched),      vTruth->at(0).mc_truth_E.at(lep1_matched)     );
            multilepton_Lepton1_P4_Matched      = LeptonX;
            multilepton_Lepton1_Id_Matched      = vTruth->at(0).mc_truth_id.at(lep1_matched);
            multilepton_Lepton1_Label_Matched   = vTruth->at(0).mc_truth_label.at(lep1_matched);
            multilepton_Lepton1_DeltaR_Matched  = lep1_dr_gen_min;
        }
        if(lep2_matched >= 0)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt.at(lep2_matched),       vTruth->at(0).mc_truth_eta.at(lep2_matched),
                                    vTruth->at(0).mc_truth_phi.at(lep2_matched),      vTruth->at(0).mc_truth_E.at(lep2_matched)     );
            multilepton_Lepton2_P4_Matched      = LeptonX;
            multilepton_Lepton2_Id_Matched      = vTruth->at(0).mc_truth_id.at(lep2_matched);
            multilepton_Lepton2_Label_Matched   = vTruth->at(0).mc_truth_label.at(lep2_matched);
            multilepton_Lepton2_DeltaR_Matched  = lep2_dr_gen_min;
        }
        if(lep3_matched >= 0)
        {
            if(vLeptonFakeable.size()>=3)
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt.at(lep3_matched),       vTruth->at(0).mc_truth_eta.at(lep3_matched),
                                        vTruth->at(0).mc_truth_phi.at(lep3_matched),      vTruth->at(0).mc_truth_E.at(lep3_matched)     );
                multilepton_Lepton3_P4_Matched      = LeptonX;
                multilepton_Lepton3_Id_Matched      = vTruth->at(0).mc_truth_id.at(lep3_matched);
                multilepton_Lepton3_Label_Matched   = vTruth->at(0).mc_truth_label.at(lep3_matched);
                multilepton_Lepton3_DeltaR_Matched  = lep3_dr_gen_min;
            }
        }
        if(jet1_matched >= 0)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt.at(jet1_matched),       vTruth->at(0).mc_truth_eta.at(jet1_matched),
                                    vTruth->at(0).mc_truth_phi.at(jet1_matched),      vTruth->at(0).mc_truth_E.at(jet1_matched)     );
            multilepton_Bjet1_P4_Matched        = LeptonX;
            multilepton_Bjet1_Id_Matched        = vTruth->at(0).mc_truth_id.at(jet1_matched);
            multilepton_Bjet1_Label_Matched     = vTruth->at(0).mc_truth_label.at(jet1_matched);
            multilepton_Bjet1_DeltaR_Matched    = jet1_dr_gen_min;
        }
        if(jet2_matched >= 0)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt.at(jet2_matched),       vTruth->at(0).mc_truth_eta.at(jet2_matched),
                                    vTruth->at(0).mc_truth_phi.at(jet2_matched),      vTruth->at(0).mc_truth_E.at(jet2_matched)     );
            multilepton_Bjet2_P4_Matched        = LeptonX;
            multilepton_Bjet2_Id_Matched        = vTruth->at(0).mc_truth_id.at(jet2_matched);
            multilepton_Bjet2_Label_Matched     = vTruth->at(0).mc_truth_label.at(jet2_matched);
            multilepton_Bjet2_DeltaR_Matched    = jet2_dr_gen_min;
        }

    }

    // ========================

    //NB : jet IDs do not have meaning
    int ij1=-1, ij2=-1, ik1=-1, ik2=-1, ie1=-1, ie2=-1, il1=-1, il2=-1, im1=-1, im2=-1, io1=-1, io2=-1, ip1=-1, ip2=-1;
    SelectOtherJets(ib1, ib2, ij1, ij2, ik1, ik2, ie1, ie2, il1, il2, im1, im2, io1, io2, ip1, ip2);

    multilepton_JetHighestPt1_Id = -999;
    multilepton_JetHighestPt2_Id = -999;
    multilepton_JetClosestMw1_Id = -999;
    multilepton_JetClosestMw2_Id = -999;
    multilepton_JetLowestMjj1_Id = -999;
    multilepton_JetLowestMjj2_Id = -999;
    multilepton_JetHighestEta1_Id = -999;
    multilepton_JetHighestEta2_Id = -999;
    multilepton_JetHighestPt1_2ndPair_Id = -999;
    multilepton_JetHighestPt2_2ndPair_Id = -999;
    multilepton_JetClosestMw1_2ndPair_Id = -999;
    multilepton_JetClosestMw2_2ndPair_Id = -999;
    multilepton_JetLowestMjj1_2ndPair_Id = -999;
    multilepton_JetLowestMjj2_2ndPair_Id = -999;

    // TLorentzVector Jet1, Jet2;
    if (ij1!=-1)
    {
        // Jet1.SetPtEtaPhiE(vJetLoose->at(ij1).pt, vJetLoose->at(ij1).eta, vJetLoose->at(ij1).phi, vJetLoose->at(ij1).E);
        // FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vJetLoose->at(ij1).CSVv2, &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vJetLoose->at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vJetLoose->at(ij1).pt_JER(), vJetLoose->at(ij1).pt_JER_up(), vJetLoose->at(ij1).pt_JER_down());
        multilepton_JetHighestPt1_Id = 1;
        multilepton_JetHighestPt1_CSV = vJetLoose->at(ij1).DeepCSVbtag;
        multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vJetLoose->at(ij1).pt, vJetLoose->at(ij1).eta, vJetLoose->at(ij1).phi, vJetLoose->at(ij1).E);
    }
    if (ij2!=-1)
    {
        // Jet2.SetPtEtaPhiE(vJetLoose->at(ij2).pt, vJetLoose->at(ij2).eta, vJetLoose->at(ij2).phi, vJetLoose->at(ij2).E);
        // FillJetInfoOutputTree(&multilepton_JetHighestPt2_Id, 1, &multilepton_JetHighestPt2_P4, Jet2, &multilepton_JetHighestPt2_CSV, vJetLoose->at(ij2).CSVv2, &multilepton_JetHighestPt2_JEC_Up, &multilepton_JetHighestPt2_JEC_Down, vJetLoose->at(ij2).JES_uncert(), &multilepton_JetHighestPt2_JER_Up, &multilepton_JetHighestPt2_JER_Down, vJetLoose->at(ij2).pt_JER(), vJetLoose->at(ij2).pt_JER_up(), vJetLoose->at(ij2).pt_JER_down());
        multilepton_JetHighestPt2_Id = 1;
        multilepton_JetHighestPt2_CSV = vJetLoose->at(ij2).DeepCSVbtag;
        multilepton_JetHighestPt2_P4.SetPtEtaPhiE(vJetLoose->at(ij2).pt, vJetLoose->at(ij2).eta, vJetLoose->at(ij2).phi, vJetLoose->at(ij2).E);
    }

    if (ik1!=-1 && ik2!=-1){
        // Jet1.SetPtEtaPhiE(vJetLoose->at(ik1).pt, vJetLoose->at(ik1).eta, vJetLoose->at(ik1).phi, vJetLoose->at(ik1).E);
        // Jet2.SetPtEtaPhiE(vJetLoose->at(ik2).pt, vJetLoose->at(ik2).eta, vJetLoose->at(ik2).phi, vJetLoose->at(ik2).E);
        // FillJetInfoOutputTree(&multilepton_JetClosestMw1_Id, 2, &multilepton_JetClosestMw1_P4, Jet1, &multilepton_JetClosestMw1_CSV, vJetLoose->at(ik1).CSVv2, &multilepton_JetClosestMw1_JEC_Up, &multilepton_JetClosestMw1_JEC_Down, vJetLoose->at(ik1).JES_uncert(), &multilepton_JetClosestMw1_JER_Up, &multilepton_JetClosestMw1_JER_Down, vJetLoose->at(ik1).pt_JER(), vJetLoose->at(ik1).pt_JER_up(), vJetLoose->at(ik1).pt_JER_down());
        // FillJetInfoOutputTree(&multilepton_JetClosestMw2_Id, 2, &multilepton_JetClosestMw2_P4, Jet2, &multilepton_JetClosestMw2_CSV, vJetLoose->at(ik2).CSVv2, &multilepton_JetClosestMw2_JEC_Up, &multilepton_JetClosestMw2_JEC_Down, vJetLoose->at(ik2).JES_uncert(), &multilepton_JetClosestMw2_JER_Up, &multilepton_JetClosestMw2_JER_Down, vJetLoose->at(ik2).pt_JER(), vJetLoose->at(ik2).pt_JER_up(), vJetLoose->at(ik2).pt_JER_down());
        multilepton_JetClosestMw1_Id = 2;
        multilepton_JetClosestMw2_Id = 2;
        multilepton_JetClosestMw1_CSV = vJetLoose->at(ik1).DeepCSVbtag;
        multilepton_JetClosestMw2_CSV = vJetLoose->at(ik2).DeepCSVbtag;
        multilepton_JetClosestMw1_P4.SetPtEtaPhiE(vJetLoose->at(ik1).pt, vJetLoose->at(ik1).eta, vJetLoose->at(ik1).phi, vJetLoose->at(ik1).E);
        multilepton_JetClosestMw2_P4.SetPtEtaPhiE(vJetLoose->at(ik2).pt, vJetLoose->at(ik2).eta, vJetLoose->at(ik2).phi, vJetLoose->at(ik2).E);
    }
    if (il1!=-1 && il2!=-1){
        // Jet1.SetPtEtaPhiE(vJetLoose->at(il1).pt, vJetLoose->at(il1).eta, vJetLoose->at(il1).phi, vJetLoose->at(il1).E);
        // Jet2.SetPtEtaPhiE(vJetLoose->at(il2).pt, vJetLoose->at(il2).eta, vJetLoose->at(il2).phi, vJetLoose->at(il2).E);
        // FillJetInfoOutputTree(&multilepton_JetLowestMjj1_Id, 3, &multilepton_JetLowestMjj1_P4, Jet1, &multilepton_JetLowestMjj1_CSV, vJetLoose->at(il1).CSVv2, &multilepton_JetLowestMjj1_JEC_Up, &multilepton_JetLowestMjj1_JEC_Down, vJetLoose->at(il1).JES_uncert(), &multilepton_JetLowestMjj1_JER_Up, &multilepton_JetLowestMjj1_JER_Down, vJetLoose->at(il1).pt_JER(), vJetLoose->at(il1).pt_JER_up(), vJetLoose->at(il1).pt_JER_down());
        // FillJetInfoOutputTree(&multilepton_JetLowestMjj2_Id, 3, &multilepton_JetLowestMjj2_P4, Jet2, &multilepton_JetLowestMjj2_CSV, vJetLoose->at(il2).CSVv2, &multilepton_JetLowestMjj2_JEC_Up, &multilepton_JetLowestMjj2_JEC_Down, vJetLoose->at(il2).JES_uncert(), &multilepton_JetLowestMjj2_JER_Up, &multilepton_JetLowestMjj2_JER_Down, vJetLoose->at(il2).pt_JER(), vJetLoose->at(il2).pt_JER_up(), vJetLoose->at(il2).pt_JER_down());
        multilepton_JetLowestMjj1_Id = 3;
        multilepton_JetLowestMjj2_Id = 3;
        multilepton_JetLowestMjj1_CSV = vJetLoose->at(il1).DeepCSVbtag;
        multilepton_JetLowestMjj2_CSV = vJetLoose->at(il2).DeepCSVbtag;
        multilepton_JetLowestMjj1_P4.SetPtEtaPhiE(vJetLoose->at(il1).pt, vJetLoose->at(il1).eta, vJetLoose->at(il1).phi, vJetLoose->at(il1).E);
        multilepton_JetLowestMjj2_P4.SetPtEtaPhiE(vJetLoose->at(il2).pt, vJetLoose->at(il2).eta, vJetLoose->at(il2).phi, vJetLoose->at(il2).E);
    }
    if(ie1!=-1)
    {
        multilepton_JetHighestEta1_Id = 4;
        multilepton_JetHighestEta1_CSV = vJetLoose->at(ie1).DeepCSVbtag;
        multilepton_JetHighestEta1_P4.SetPtEtaPhiE(vJetLoose->at(ie1).pt, vJetLoose->at(ie1).eta, vJetLoose->at(ie1).phi, vJetLoose->at(ie1).E );
    }
    if(ie2!=-1) //2jets
    {
        multilepton_JetHighestEta2_Id = 4;
        multilepton_JetHighestEta2_CSV = vJetLoose->at(ie2).DeepCSVbtag;
        multilepton_JetHighestEta2_P4.SetPtEtaPhiE(vJetLoose->at(ie2).pt, vJetLoose->at(ie2).eta, vJetLoose->at(ie2).phi, vJetLoose->at(ie2).E );
    }


    //--- Fill 2nd pairs (first one is closest to mW) -- needed for 2l only (more jets)
    if(is_tHq_2lSS && ij1!=-1 && ij2!=-1)
    {
        if (im1!=-1)
        {
            // Jet1.SetPtEtaPhiE(vJetLoose->at(im1).pt, vJetLoose->at(im1).eta, vJetLoose->at(im1).phi, vJetLoose->at(im1).E);
            // FillJetInfoOutputTree(&multilepton_JetHighestPt1_2ndPair_Id, 1, &multilepton_JetHighestPt1_2ndPair_P4, Jet1, &multilepton_JetHighestPt1_2ndPair_CSV, vJetLoose->at(im1).CSVv2, &multilepton_JetHighestPt1_2ndPair_JEC_Up, &multilepton_JetHighestPt1_2ndPair_JEC_Down, vJetLoose->at(im1).JES_uncert(), &multilepton_JetHighestPt1_2ndPair_JER_Up, &multilepton_JetHighestPt1_2ndPair_JER_Down, vJetLoose->at(im1).pt_JER(), vJetLoose->at(im1).pt_JER_up(), vJetLoose->at(im1).pt_JER_down());
            multilepton_JetHighestPt1_2ndPair_Id = 1;
            multilepton_JetHighestPt1_2ndPair_CSV = vJetLoose->at(im1).DeepCSVbtag;
            multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vJetLoose->at(im1).pt, vJetLoose->at(im1).eta, vJetLoose->at(im1).phi, vJetLoose->at(im1).E);
        }
        if(im2!=-1){
            // Jet2.SetPtEtaPhiE(vJetLoose->at(im2).pt, vJetLoose->at(im2).eta, vJetLoose->at(im2).phi, vJetLoose->at(im2).E);
            // FillJetInfoOutputTree(&multilepton_JetHighestPt2_2ndPair_Id, 1, &multilepton_JetHighestPt2_2ndPair_P4, Jet2, &multilepton_JetHighestPt2_2ndPair_CSV, vJetLoose->at(im2).CSVv2, &multilepton_JetHighestPt2_2ndPair_JEC_Up, &multilepton_JetHighestPt2_2ndPair_JEC_Down, vJetLoose->at(im2).JES_uncert(), &multilepton_JetHighestPt2_2ndPair_JER_Up, &multilepton_JetHighestPt2_2ndPair_JER_Down, vJetLoose->at(im2).pt_JER(), vJetLoose->at(im2).pt_JER_up(), vJetLoose->at(im2).pt_JER_down());
            multilepton_JetHighestPt2_2ndPair_Id = 1;
            multilepton_JetHighestPt1_2ndPair_CSV = vJetLoose->at(im2).DeepCSVbtag;
            multilepton_JetHighestPt2_2ndPair_P4.SetPtEtaPhiE(vJetLoose->at(im2).pt, vJetLoose->at(im2).eta, vJetLoose->at(im2).phi, vJetLoose->at(im2).E);
        }
        if (io1!=-1 && io2!=-1){
            // Jet1.SetPtEtaPhiE(vJetLoose->at(ip1).pt, vJetLoose->at(ip1).eta, vJetLoose->at(ip1).phi, vJetLoose->at(ip1).E);
            // Jet2.SetPtEtaPhiE(vJetLoose->at(io2).pt, vJetLoose->at(io2).eta, vJetLoose->at(io2).phi, vJetLoose->at(io2).E);
            // FillJetInfoOutputTree(&multilepton_JetClosestMw1_2ndPair_Id, 2, &multilepton_JetClosestMw1_2ndPair_P4, Jet1, &multilepton_JetClosestMw1_2ndPair_CSV, vJetLoose->at(io1).CSVv2, &multilepton_JetClosestMw1_2ndPair_JEC_Up, &multilepton_JetClosestMw1_2ndPair_JEC_Down, vJetLoose->at(io1).JES_uncert(), &multilepton_JetClosestMw1_2ndPair_JER_Up, &multilepton_JetClosestMw1_2ndPair_JER_Down, vJetLoose->at(io1).pt_JER(), vJetLoose->at(io1).pt_JER_up(), vJetLoose->at(io1).pt_JER_down());
            // FillJetInfoOutputTree(&multilepton_JetClosestMw2_2ndPair_Id, 2, &multilepton_JetClosestMw2_2ndPair_P4, Jet2, &multilepton_JetClosestMw2_2ndPair_CSV, vJetLoose->at(io2).CSVv2, &multilepton_JetClosestMw2_2ndPair_JEC_Up, &multilepton_JetClosestMw2_2ndPair_JEC_Down, vJetLoose->at(io2).JES_uncert(), &multilepton_JetClosestMw2_2ndPair_JER_Up, &multilepton_JetClosestMw2_2ndPair_JER_Down, vJetLoose->at(io2).pt_JER(), vJetLoose->at(io2).pt_JER_up(), vJetLoose->at(io2).pt_JER_down());
            multilepton_JetClosestMw1_2ndPair_Id = 2;
            multilepton_JetClosestMw2_2ndPair_Id = 2;
            multilepton_JetClosestMw1_2ndPair_CSV = vJetLoose->at(io1).DeepCSVbtag;
            multilepton_JetClosestMw2_2ndPair_CSV = vJetLoose->at(io2).DeepCSVbtag;
            multilepton_JetClosestMw1_2ndPair_P4.SetPtEtaPhiE(vJetLoose->at(io1).pt, vJetLoose->at(io1).eta, vJetLoose->at(io1).phi, vJetLoose->at(io1).E);
            multilepton_JetClosestMw2_2ndPair_P4.SetPtEtaPhiE(vJetLoose->at(io2).pt, vJetLoose->at(io2).eta, vJetLoose->at(io2).phi, vJetLoose->at(io2).E);
        }
        if (ip1!=-1 && ip2!=-1){
            // Jet1.SetPtEtaPhiE(vJetLoose->at(ip1).pt, vJetLoose->at(ip1).eta, vJetLoose->at(ip1).phi, vJetLoose->at(ip1).E);
            // Jet2.SetPtEtaPhiE(vJetLoose->at(ip2).pt, vJetLoose->at(ip2).eta, vJetLoose->at(ip2).phi, vJetLoose->at(ip2).E);
            // FillJetInfoOutputTree(&multilepton_JetLowestMjj1_2ndPair_Id, 3, &multilepton_JetLowestMjj1_2ndPair_P4, Jet1, &multilepton_JetLowestMjj1_2ndPair_CSV, vJetLoose->at(ip1).CSVv2, &multilepton_JetLowestMjj1_2ndPair_JEC_Up, &multilepton_JetLowestMjj1_2ndPair_JEC_Down, vJetLoose->at(ip1).JES_uncert(), &multilepton_JetLowestMjj1_2ndPair_JER_Up, &multilepton_JetLowestMjj1_2ndPair_JER_Down, vJetLoose->at(ip1).pt_JER(), vJetLoose->at(ip1).pt_JER_up(), vJetLoose->at(ip1).pt_JER_down());
            // FillJetInfoOutputTree(&multilepton_JetLowestMjj2_2ndPair_Id, 3, &multilepton_JetLowestMjj2_2ndPair_P4, Jet2, &multilepton_JetLowestMjj2_2ndPair_CSV, vJetLoose->at(ip2).CSVv2, &multilepton_JetLowestMjj2_2ndPair_JEC_Up, &multilepton_JetLowestMjj2_2ndPair_JEC_Down, vJetLoose->at(ip2).JES_uncert(), &multilepton_JetLowestMjj2_2ndPair_JER_Up, &multilepton_JetLowestMjj2_2ndPair_JER_Down, vJetLoose->at(ip2).pt_JER(), vJetLoose->at(ip2).pt_JER_up(), vJetLoose->at(ip2).pt_JER_down());
            multilepton_JetLowestMjj1_2ndPair_Id = 3;
            multilepton_JetLowestMjj2_2ndPair_Id = 3;
            multilepton_JetLowestMjj1_2ndPair_CSV = vJetLoose->at(ip1).DeepCSVbtag;
            multilepton_JetLowestMjj2_2ndPair_CSV = vJetLoose->at(ip2).DeepCSVbtag;
            multilepton_JetLowestMjj1_2ndPair_P4.SetPtEtaPhiE(vJetLoose->at(ip1).pt, vJetLoose->at(ip1).eta, vJetLoose->at(ip1).phi, vJetLoose->at(ip1).E);
            multilepton_JetLowestMjj2_2ndPair_P4.SetPtEtaPhiE(vJetLoose->at(ip2).pt, vJetLoose->at(ip2).eta, vJetLoose->at(ip2).phi, vJetLoose->at(ip2).E);
        }
    }


    // ##########################################################################
    // #      _                  _       _                                      #
    // #  ___| |_ __ _ _ __   __| | __ _| | ___  _ __   ___    __ _  ___ _ __   #
    // # / __| __/ _` | '_ \ / _` |/ _` | |/ _ \| '_ \ / _ \  / _` |/ _ \ '_ \  #
    // # \__ \ || (_| | | | | (_| | (_| | | (_) | | | |  __/ | (_| |  __/ | | | #
    // # |___/\__\__,_|_| |_|\__,_|\__,_|_|\___/|_| |_|\___|  \__, |\___|_| |_| #
    // #                                                      |___/             #
    // #                                                                        #
    // ##########################################################################

    if( !_isdata )
    {
        for(unsigned int itruth = 0; itruth < vTruth->at(0).mc_truth_label.size() ; itruth++)
        {
            TLorentzVector LeptonX;

            if( vTruth->at(0).mc_truth_label.at(itruth) == 1 )
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt.at(itruth),
                                        vTruth->at(0).mc_truth_eta.at(itruth),
                                        vTruth->at(0).mc_truth_phi.at(itruth),
                                        vTruth->at(0).mc_truth_E.at(itruth) );

                multilepton_h0_P4 = LeptonX;
                multilepton_h0_Id = vTruth->at(0).mc_truth_id.at(itruth);
            }

            if( vTruth->at(0).mc_truth_label.at(itruth) == 2 )
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt.at(itruth),
                                        vTruth->at(0).mc_truth_eta.at(itruth),
                                        vTruth->at(0).mc_truth_phi.at(itruth),
                                        vTruth->at(0).mc_truth_E.at(itruth) );

                multilepton_t1_P4 = LeptonX;
                multilepton_t1_Id = vTruth->at(0).mc_truth_id.at(itruth);
            }

            if( vTruth->at(0).mc_truth_label.at(itruth) == 3 )
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt.at(itruth),
                                        vTruth->at(0).mc_truth_eta.at(itruth),
                                        vTruth->at(0).mc_truth_phi.at(itruth),
                                        vTruth->at(0).mc_truth_E.at(itruth) );

                multilepton_t2_P4 = LeptonX;
                multilepton_t2_Id = vTruth->at(0).mc_truth_id.at(itruth);
            }
        }
    }

    multilepton_mET.SetPtEtaPhiE(vEvent->at(0).metpt, 0, vEvent->at(0).metphi, vEvent->at(0).metpt);
    multilepton_mETcov00 = vEvent->at(0).metcov00;
    multilepton_mETcov01 = vEvent->at(0).metcov01;
    multilepton_mETcov10 = vEvent->at(0).metcov10;
    multilepton_mETcov11 = vEvent->at(0).metcov11;
    multilepton_mHT = vEvent->at(0).metsumet; //FIXME -- correct ?

    mc_ttZhypAllowed = 0;
    if (multilepton_Lepton1_Id!=-999 && multilepton_Lepton2_Id!=-999 && multilepton_Lepton3_Id!=-999)
    {
        //+++ && --- cases
        if (multilepton_Lepton1_Id*multilepton_Lepton2_Id>0 && multilepton_Lepton2_Id*multilepton_Lepton3_Id>0) {mc_ttZhypAllowed =-1;}
        //OSSF pair case
        else if ( (multilepton_Lepton1_Id==-multilepton_Lepton2_Id)
                || (multilepton_Lepton1_Id==-multilepton_Lepton3_Id)
                || (multilepton_Lepton2_Id==-multilepton_Lepton3_Id)) {mc_ttZhypAllowed = 1;}
    }

    tOutput->Fill();

    return;
}



























//--------------------------------------------
// #### ##    ## ########  ##     ## ########    ##     ##    ###    ########
//  ##  ###   ## ##     ## ##     ##    ##       ##     ##   ## ##   ##     ##
//  ##  ####  ## ##     ## ##     ##    ##       ##     ##  ##   ##  ##     ##
//  ##  ## ## ## ########  ##     ##    ##       ##     ## ##     ## ########
//  ##  ##  #### ##        ##     ##    ##        ##   ##  ######### ##   ##
//  ##  ##   ### ##        ##     ##    ##         ## ##   ##     ## ##    ##
// #### ##    ## ##         #######     ##          ###    ##     ## ##     ##
//--------------------------------------------


/**
 * Compute input variables and other useful variables (channel, ...)
 * Compute both tHq2016 & ttH2017 input variables, for tests
 * @param region (3l or 2l) ==> To compute properly "sum_id" & "Lep3Pt/lep3Pt" & lepCharge
 */
void tHqMultileponAnalysis::Compute_Variables(TString region)
{
	if(region != "3l" && region != "2l")
	{
		cout<<FRED("Error ! Wrong region name at Compute_Variables() call !")<<endl;
	}

    int nleptons = 0;
    if(region == "3l") {nleptons = 3;}
    else {nleptons = 2;}

    //--- Determine leptonic channel of event
	int sum_id = 0;
	if(region == "3l")
    {
        sum_id = fabs(vLeptonFakeable.at(0).id)+fabs(vLeptonFakeable.at(1).id)+fabs(vLeptonFakeable.at(2).id);
        if(sum_id == 39) {channel = 0;} //uuu
        else if(sum_id == 37) {channel = 1;} //uue
        else if(sum_id == 35) {channel = 2;} //eeu
        else if(sum_id == 33) {channel = 3;} //eee
    }
    else
    {
        sum_id = fabs(vLeptonFakeable.at(0).id)+fabs(vLeptonFakeable.at(1).id);
        if(sum_id == 26) {channel = 0;} //uu
        else if(sum_id == 24) {channel = 1;} //ue + eu
        else if(sum_id == 22) {channel = 2;} //ee -- not used in tHq2016
    }



    int ijet_forward=-1, ijet_hardest_btag=-1, ijet_2nd_hardest_btag=-1;
    double tmp = -999;

    //--- Find "forward jet"
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta) > tmp)
        {
            tmp = fabs(vLightJets.at(ijet).eta);
            ijet_forward = ijet;
        }
    }

    //--- Find "hardest" and "second hardest" tagged jets
    double pt1=0, pt2=0;
    for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
    {
        if(fabs(vLooseBTagJets.at(ijet).pt) > pt1)
        {
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt);
            ijet_hardest_btag = ijet;
        }
        else if(fabs(vLooseBTagJets.at(ijet).pt) > pt2)
        {
            pt2 = fabs(vLooseBTagJets.at(ijet).pt);
            ijet_2nd_hardest_btag = ijet;
        }
    }

    TLorentzVector fwdJet, BJet;
    if(ijet_forward >= 0) {fwdJet.SetPtEtaPhiE(vLightJets.at(ijet_forward).pt, vLightJets.at(ijet_forward).eta, vLightJets.at(ijet_forward).phi, vLightJets.at(ijet_forward).E);}
    if(ijet_hardest_btag >= 0) {BJet.SetPtEtaPhiE(vLooseBTagJets.at(ijet_hardest_btag).pt, vLooseBTagJets.at(ijet_hardest_btag).eta, vLooseBTagJets.at(ijet_hardest_btag).phi, vLooseBTagJets.at(ijet_hardest_btag).E);}

    // cout<<"ijet_forward "<<ijet_forward<<endl;
    // cout<<"ijet_hardest_btag "<<ijet_hardest_btag<<endl;
    // cout<<"ijet_2nd_hardest_btag "<<ijet_2nd_hardest_btag<<endl;


//--------------------------------------------
//       #     #         #####    ###     #    #####
// ##### #     #  ####  #     #  #   #   ##   #     #
//   #   #     # #    #       # #     # # #   #
//   #   ####### #    #  #####  #     #   #   ######
//   #   #     # #  # # #       #     #   #   #     #
//   #   #     # #   #  #        #   #    #   #     #
//   #   #     #  ### # #######   ###   #####  #####
//--------------------------------------------

    //--- Compute input variables (tHq 2016)
    //NB : "forward jet" = most forward non-CSV loose jet

	//--- Var1 : nof jets with pT>25 and |eta|<2.4
    nJet25 = 0;
    for(int ijet=0; ijet<vJetLoose->size(); ijet++)
    {
        if(vJetLoose->at(ijet).pt > 25 && fabs(vJetLoose->at(ijet).eta) < 2.4) {nJet25++;}
    }

    //NB : nJet25 != vJetLoose->size(), because vJetLoose accounts for JES variations (pT could be slightly < 25 !)
    nJetLoose = vJetLoose->size();


    //--- Var2: max eta of any 'non-CSV-loose' jet
    tmp = -999;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta) > tmp) {tmp = fabs(vLightJets.at(ijet).eta);}
    }
    maxEtaJet25 = tmp;

	//--- Var3 : sum of leptons charges
    for(int ilep=0; ilep<nleptons; ilep++)
    {
        lepCharge+= vLeptonFakeable.at(ilep).charge;
    }

	//--- Var4 : nof 'non-csv-loose' jets with eta>1.0
    nJetEta1 = 0;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if( fabs(vLightJets.at(ijet).eta ) > 1.0) {nJetEta1++;}
    }

	//--- Var5 : dEta between forward light jet and hardest tagged jet
    dEtaFwdJetBJet = -999;
    if(ijet_hardest_btag >= 0 && ijet_forward >= 0) {dEtaFwdJetBJet = fabs( vLightJets.at(ijet_forward).eta - vLooseBTagJets.at(ijet_hardest_btag).eta );}

	//--- Var6 : dEta between forward and 2nd hardest tagged jet
    if(ijet_2nd_hardest_btag < 0 || ijet_forward < 0) {dEtaFwdJet2BJet = -1;}
    else {dEtaFwdJet2BJet = fabs( vLightJets.at(ijet_forward).eta - vLooseBTagJets.at(ijet_2nd_hardest_btag).eta );}


	//--- Var7 : dEta between forward light jet and closet lepton (angular dist.)
    dEtaFwdJetClosestLep = -999; tmp = 999;
    for(int ilep=0; ilep<vLeptonFakeable.size(); ilep++)
    {
        if(ijet_forward >= 0)
        {
            if( fabs(vLightJets.at(ijet_forward).eta - vLeptonFakeable.at(ilep).eta ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta - vLeptonFakeable.at(ilep).eta ); tmp = dEtaFwdJetClosestLep;}
        }
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = -999; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vLeptonFakeable.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vLeptonFakeable.at(i).pt, vLeptonFakeable.at(i).eta, vLeptonFakeable.at(i).phi, vLeptonFakeable.at(i).E );
        for(int j=i+1; j<vLeptonFakeable.size(); j++)
        {
            lepj.SetPtEtaPhiE(vLeptonFakeable.at(j).pt, vLeptonFakeable.at(j).eta, vLeptonFakeable.at(j).phi, vLeptonFakeable.at(j).E );
            if(vLeptonFakeable.at(i).charge==vLeptonFakeable.at(j).charge && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vLeptonFakeable.at(i).phi - vLeptonFakeable.at(j).phi ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    lepi = TLorentzVector();
    lepj = TLorentzVector();
    for(int i=0; i<vLeptonFakeable.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vLeptonFakeable.at(i).pt, vLeptonFakeable.at(i).eta, vLeptonFakeable.at(i).phi, vLeptonFakeable.at(i).E );
        for(int j=i+1; j<vLeptonFakeable.size(); j++)
        {
            lepj.SetPtEtaPhiE(vLeptonFakeable.at(j).pt, vLeptonFakeable.at(j).eta, vLeptonFakeable.at(j).phi, vLeptonFakeable.at(j).E );
            if(lepi.DeltaR(lepj) < minDRll)
            {

                minDRll = lepi.DeltaR(lepj);
            }
        }
    }

    //--- Var10 : pT of 3rd hardest lepton
    if(region == "3l") {Lep3Pt = vLeptonFakeable.at(2).pt;}
    else {Lep3Pt = vLeptonFakeable.at(1).pt;}

    //--- Fill additionnal variables, used for control only
    lep1Pt = vLeptonFakeable.at(0).pt; lep2Pt = vLeptonFakeable.at(1).pt; lep3Pt = vLeptonFakeable.at(1).pt;
    if(region == "3l") {lep3Pt = vLeptonFakeable.at(2).pt;}
    if(ijet_hardest_btag >= 0)
    {
        hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt;
        hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta;
    }
    if(ijet_forward >= 0) {fwdJetPt = vLightJets.at(ijet_forward).pt;}


    //--------------------------------------------
    //             #     #  #####    ###     #   #######
    // ##### ##### #     # #     #  #   #   ##   #    #
    //   #     #   #     #       # #     # # #       #
    //   #     #   #######  #####  #     #   #      #
    //   #     #   #     # #       #     #   #     #
    //   #     #   #     # #        #   #    #     #
    //   #     #   #     # #######   ###   #####   #
    //--------------------------------------------

    lep1_conePt = vLeptonFakeable.at(0).conept;
    if(nleptons == 3) {lep2_conePt = vLeptonFakeable.at(2).conept;} //'Trailing' lepton
    else {lep2_conePt = vLeptonFakeable.at(1).conept;}

    mindr_lep1_jet = 999;
    mindr_lep2_jet = 999;
    for(int i=0; i<vJetLoose->size(); i++)
    {
        float dr1 = GetDeltaR(vLeptonFakeable.at(0).eta,vLeptonFakeable.at(0).phi,vJetLoose->at(i).eta,vJetLoose->at(i).phi);
        if(dr1 < mindr_lep1_jet) {mindr_lep1_jet = dr1;}

        float dr2 = GetDeltaR(vLeptonFakeable.at(1).eta,vLeptonFakeable.at(1).phi,vJetLoose->at(i).eta,vJetLoose->at(i).phi);
        if(dr2 < mindr_lep2_jet) {mindr_lep2_jet = dr2;}
    }

    // mT_lep1 = comp_MT_met_lep(lep_tmp, vEvent->at(0).metpt, vEvent->at(0).metphi);
    mT_lep1 = sqrt( 2*vLeptonFakeable.at(0).conept * MET * (1 - cos(vLeptonFakeable.at(0).phi - vEvent->at(0).metphi) ) );
    // mT_lep2 = comp_MT_met_lep(lep_tmp, vEvent->at(0).metpt, vEvent->at(0).metphi);
    mT_lep2 = sqrt( 2*vLeptonFakeable.at(1).conept * MET * (1 - cos(vLeptonFakeable.at(1).phi - vEvent->at(0).metphi) ) );

    max_lep_eta = std::max( fabs(vLeptonFakeable.at(0).eta), fabs(vLeptonFakeable.at(1).eta) );
    if(is_tHq_3l) {max_lep_eta = std::max(std::max(fabs(vLeptonFakeable.at(0).eta), fabs(vLeptonFakeable.at(1).eta)), fabs(vLeptonFakeable.at(2).eta));}

    //--------------------------------------------
    //                               #
    // #    # ###### #    #         #     ##### ######  ####  ##### # #    #  ####
    // ##   # #      #    #        #        #   #      #        #   # ##   # #    #
    // # #  # #####  #    #       #         #   #####   ####    #   # # #  # #
    // #  # # #      # ## #      #          #   #           #   #   # #  # # #  ###
    // #   ## #      ##  ##     #           #   #      #    #   #   # #   ## #    #
    // #    # ###### #    #    #            #   ######  ####    #   # #    #  ####
    //--------------------------------------------

    minv_FwdJetBJet = 0;
    if(ijet_forward >= 0 && ijet_hardest_btag >= 0) {minv_FwdJetBJet = (fwdJet + BJet).M();}

    etaFwdJet = -10;
    if(ijet_forward >= 0) {etaFwdJet = fwdJet.Eta();}

    ptFwdJet = 0;
    if(ijet_forward >= 0) {ptFwdJet = fwdJet.Pt();}

    LeadJetEta = -999;
    LeadJetPt = 0;
    for(int j=0; j<vJetLoose->size(); j++)
    {
        if(vJetLoose->at(j).pt > LeadJetPt)
        {
            LeadJetPt = vJetLoose->at(j).pt;
            LeadJetEta = vJetLoose->at(j).eta;
        }
    }

    dRjj = -999;
    Mjj_max = 0;
    dPhijj_max = -999;
    TLorentzVector jet_tmp;
    for(int i=0; i<vJetLoose->size(); i++)
    {
        jet_tmp.SetPtEtaPhiE(vJetLoose->at(i).pt, vJetLoose->at(i).eta, vJetLoose->at(i).phi, vJetLoose->at(i).E );
        for(int j=i+1; j<vJetLoose->size(); j++)
        {
            TLorentzVector jet_tmp_2;
            jet_tmp_2.SetPtEtaPhiE(vJetLoose->at(j).pt, vJetLoose->at(j).eta, vJetLoose->at(j).phi, vJetLoose->at(j).E );

            if( GetDeltaR(vJetLoose->at(i).eta, vJetLoose->at(i).phi, vJetLoose->at(j).eta, vJetLoose->at(j).phi) > dRjj)
            {
                dRjj = GetDeltaR(vJetLoose->at(i).eta, vJetLoose->at(i).phi, vJetLoose->at(j).eta, vJetLoose->at(j).phi);
            }
            if( (jet_tmp + jet_tmp_2).M() >  Mjj_max)
            {
                Mjj_max = (jet_tmp + jet_tmp_2).M();
            }
            if(Phi_MPi_Pi(vJetLoose->at(i).phi - vJetLoose->at(j).phi ) > dPhijj_max)
            {
                dPhijj_max = Phi_MPi_Pi(vJetLoose->at(i).phi - vJetLoose->at(j).phi );
            }
        }
    }

    deepCSV_max = 0;
    deepCSV_2nd = 0;
    for(int i=0; i<vLooseBTagJets.size(); i++)
    {
        if(vLooseBTagJets.at(i).DeepCSVbtag > deepCSV_max) {deepCSV_max = vLooseBTagJets.at(i).DeepCSVbtag;}
        else if(vLooseBTagJets.at(i).DeepCSVbtag > deepCSV_2nd) {deepCSV_2nd = vLooseBTagJets.at(i).DeepCSVbtag;}
    }

    dPhiLepBJet_max = -999;
    if(ijet_hardest_btag >= 0)
    {
        for(int i=0; i<nleptons; i++)
        {
            if(Phi_MPi_Pi(vLooseBTagJets.at(ijet_hardest_btag).phi - vLeptonFakeable.at(i).phi ) > dPhiLepBJet_max)
            {
                dPhiLepBJet_max = Phi_MPi_Pi(vLooseBTagJets.at(ijet_hardest_btag).phi - vLeptonFakeable.at(i).phi );
            }
        }
    }

    //--------------------------------------------
    // #####  #####  # #    # #####  ####  #    # #####
    // #    # #    # # ##   #   #   #    # #    #   #
    // #    # #    # # # #  #   #   #    # #    #   #
    // #####  #####  # #  # #   #   #    # #    #   #
    // #      #   #  # #   ##   #   #    # #    #   #
    // #      #    # # #    #   #    ####   ####    #
    //--------------------------------------------

    bool do_printout = false;

    if(do_printout)
    {
        cout<<endl<<endl<<BOLD(FBLU("------------ EVENT -----------"))<<endl;

        cout<<FYEL("--- Tagged Jets : ")<<endl;
        for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt<<" / eta = "<<vLooseBTagJets.at(ijet).eta<<" / phi = "<<vLooseBTagJets.at(ijet).phi<<" / CSV = "<<vLooseBTagJets.at(ijet).DeepCSVbtag<<endl;
        }
        cout<<"Hardest & 2nd hardest jets are "<<ijet_hardest_btag<<", "<<ijet_2nd_hardest_btag<<endl<<endl;

        cout<<FYEL("--- Forward Jets : ")<<endl;
        for(int ijet=0; ijet<vLightJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLightJets.at(ijet).pt<<" / eta = "<<vLightJets.at(ijet).eta<<" / phi = "<<vLightJets.at(ijet).phi<<" / CSV = "<<vLightJets.at(ijet).DeepCSVbtag<<endl;
        }
        cout<<"Forwardest jet is "<<ijet_forward<<endl<<endl;

        cout<<FYEL("--- Selected Leptons : ")<<endl;
        for(int ilep=0; ilep<vLeptonFakeable.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vLeptonFakeable.at(ilep).pt<<" / eta = "<<vLeptonFakeable.at(ilep).eta<<" / phi = "<<vLeptonFakeable.at(ilep).phi<<" : Charge = "<<vLeptonFakeable.at(ilep).charge<<endl;
        }

        cout<<FYEL("--- Input variables : ")<<endl;

        cout<<"nJet25 = "<<nJet25<<endl;
        cout<<"maxEtaJet25 = "<<maxEtaJet25<<endl;
        cout<<"lepCharge = "<<lepCharge<<endl;
        cout<<"nJetEta1 = "<<nJetEta1 <<endl;
        cout<<"dEtaFwdJetBJet = "<<dEtaFwdJetBJet<<endl;
        cout<<"dEtaFwdJet2BJet = "<<dEtaFwdJet2BJet<<endl;
        cout<<"dEtaFwdJetClosestLep = "<<dEtaFwdJetClosestLep<<endl;
        cout<<"dPhiHighestPtSSPair = "<<dPhiHighestPtSSPair<<endl;
        cout<<"minDRll = "<<minDRll<<endl;
        cout<<"Lep3Pt = "<<Lep3Pt<<endl;

        cout<<"------------------"<<endl<<endl;
    }


    //--- COMPUTE BDT OUTPUT from weight files -- tHq2016 weights
    // if(region == "3l")
    // {
        // signal_3l_TT_MVA   = mva_3l_tt->EvaluateMVA(""BDT"G method");
        // signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");
    // }
    // else
    // {
        // signal_2lss_TT_MVA   = mva_2lss_tt->EvaluateMVA("BDTG method");
        // signal_2lss_TTV_MVA   = mva_2lss_ttV->EvaluateMVA("BDTG method");
    // }
    // ======================================================================================================

    return;
}






















//--------------------------------------------
// ######## ######## ########  ######## ########       ########  ########     ###    ##    ##  ######  ##     ##
//    ##       ##    ##     ## ##       ##             ##     ## ##     ##   ## ##   ###   ## ##    ## ##     ##
//    ##       ##    ##     ## ##       ##             ##     ## ##     ##  ##   ##  ####  ## ##       ##     ##
//    ##       ##    ########  ######   ######         ########  ########  ##     ## ## ## ## ##       #########
//    ##       ##    ##   ##   ##       ##             ##     ## ##   ##   ######### ##  #### ##       ##     ##
//    ##       ##    ##    ##  ##       ##             ##     ## ##    ##  ##     ## ##   ### ##    ## ##     ##
//    ##       ##    ##     ## ######## ########       ########  ##     ## ##     ## ##    ##  ######  ##     ##
//--------------------------------------------

/**
 * Create all output branches
 */
void tHqMultileponAnalysis::initializeOutputTree()
{
    outputfile->cd();
    tOutput = new TTree("Tree", "Tree");

	//-- Event main infos
	tOutput->Branch("channel",&channel,"channel/F");
    tOutput->Branch("weight",&weight,"weight/F");
	tOutput->Branch("weightfake",&weightfake,"weightfake/F");
	tOutput->Branch("weightflip",&weightflip,"weightflip/F");
    tOutput->Branch("event_id",&event_id,"event_id/F");
    tOutput->Branch("event_run",&event_run,"event_run/F");
	tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");
    tOutput->Branch("is_trigger",&is_trigger,"is_trigger/O");
    tOutput->Branch("is_trigger_ttH",&is_trigger_ttH,"is_trigger_ttH/O");
    tOutput->Branch("is_hasJetTransitionRegion",&is_hasJetTransitionRegion,"is_hasJetTransitionRegion/O");


	//--- Categories & MVA

    //tHq 2017 categories (implemented in NTA)
    tOutput->Branch("is_tHq_2lSS",&is_tHq_2lSS,"is_tHq_2lSS/B");
	tOutput->Branch("is_tHq_2lSS_SR",&is_tHq_2lSS_SR,"is_tHq_2lSS_SR/B");
	tOutput->Branch("is_tHq_2lSS_Training",&is_tHq_2lSS_Training,"is_tHq_2lSS_Training/B");
    tOutput->Branch("is_tHq_2lSS_Fake",&is_tHq_2lSS_Fake,"is_tHq_2lSS_Fake/B");
    tOutput->Branch("is_tHq_2lSS_Flip",&is_tHq_2lSS_Flip,"is_tHq_2lSS_Flip/B");
    tOutput->Branch("is_tHq_2lSS_GammaConv",&is_tHq_2lSS_GammaConv,"is_tHq_2lSS_GammaConv/B");
    tOutput->Branch("is_tHq_3l",&is_tHq_3l,"is_tHq_3l/B");
	tOutput->Branch("is_tHq_3l_SR",&is_tHq_3l_SR,"is_tHq_3l_SR/B");
	tOutput->Branch("is_tHq_3l_Training",&is_tHq_3l_Training,"is_tHq_3l_Training/B");
    tOutput->Branch("is_tHq_3l_Fake",&is_tHq_3l_Fake,"is_tHq_3l_Fake/B");
    tOutput->Branch("is_tHq_ttWctrl",&is_tHq_ttWctrl,"is_tHq_ttWctrl/B");
    tOutput->Branch("is_tHq_ttWctrl_SR",&is_tHq_ttWctrl_SR,"is_tHq_ttWctrl_SR/B");
	tOutput->Branch("is_tHq_ttWctrl_Fake",&is_tHq_ttWctrl_Fake,"is_tHq_ttWctrl_Fake/B");
    tOutput->Branch("is_tHq_ttWctrl_Flip",&is_tHq_ttWctrl_Flip,"is_tHq_ttWctrl_Flip/B");
    tOutput->Branch("is_tHq_ttWctrl_GammaConv",&is_tHq_ttWctrl_GammaConv,"is_tHq_ttWctrl_GammaConv/B");
    tOutput->Branch("is_tHq_ttZctrl",&is_tHq_ttZctrl,"is_tHq_ttZctrl/B");
    tOutput->Branch("is_tHq_ttZctrl_SR",&is_tHq_ttZctrl_SR,"is_tHq_ttZctrl_SR/B");
    tOutput->Branch("is_tHq_ttZctrl_Fake",&is_tHq_ttZctrl_Fake,"is_tHq_ttZctrl_Fake/B");
    tOutput->Branch("is_tHq_WZctrl",&is_tHq_WZctrl,"is_tHq_WZctrl/B");
    tOutput->Branch("is_tHq_WZctrl_SR",&is_tHq_WZctrl_SR,"is_tHq_WZctrl_SR/B");
    tOutput->Branch("is_tHq_WZctrl_Fake",&is_tHq_WZctrl_Fake,"is_tHq_WZctrl_Fake/B");

    //ttH2017 predefined categories (implemented in NTP)
    tOutput->Branch("is_ttH_2lSS",&is_ttH_2lSS,"is_ttH_2lSS/B");
    tOutput->Branch("is_ttH_2lSS_SR",&is_ttH_2lSS_SR,"is_ttH_2lSS_SR/B");
    tOutput->Branch("is_ttH_2lSS_SR_Data",&is_ttH_2lSS_SR_Data,"is_ttH_2lSS_SR_Data/B");
    tOutput->Branch("is_ttH_2lSS_Fake",&is_ttH_2lSS_Fake,"is_ttH_2lSS_Fake/B");
    tOutput->Branch("is_ttH_2lSS_Flip",&is_ttH_2lSS_Flip,"is_ttH_2lSS_Flip/B");
    tOutput->Branch("is_ttH_2lSS_Flip_Data",&is_ttH_2lSS_Flip_Data,"is_ttH_2lSS_Flip_Data/B");
    tOutput->Branch("is_ttH_3l",&is_ttH_3l,"is_ttH_3l/B");
    tOutput->Branch("is_ttH_3l_SR",&is_ttH_3l_SR,"is_ttH_3l_SR/B");
    tOutput->Branch("is_ttH_3l_SR_Data",&is_ttH_3l_SR_Data,"is_ttH_3l_SR_Data/B");
    tOutput->Branch("is_ttH_3l_Fake",&is_ttH_3l_Fake,"is_ttH_3l_Fake/B");
    tOutput->Branch("is_ttH_ttWctrl",&is_ttH_ttWctrl,"is_ttH_ttWctrl/B");
    tOutput->Branch("is_ttH_ttWctrl_SR",&is_ttH_ttWctrl_SR,"is_ttH_ttWctrl_SR/B");
    tOutput->Branch("is_ttH_ttWctrl_SR_Data",&is_ttH_ttWctrl_SR_Data,"is_ttH_ttWctrl_SR_Data/B");
    tOutput->Branch("is_ttH_ttWctrl_Fake",&is_ttH_ttWctrl_Fake,"is_ttH_ttWctrl_Fake/B");
    tOutput->Branch("is_ttH_ttWctrl_Flip",&is_ttH_ttWctrl_Flip,"is_ttH_ttWctrl_Flip/B");
    tOutput->Branch("is_ttH_ttWctrl_Flip_Data",&is_ttH_ttWctrl_Flip_Data,"is_ttH_ttWctrl_Flip_Data/B");
    tOutput->Branch("is_ttH_ttZctrl",&is_ttH_ttZctrl,"is_ttH_ttZctrl/B");
    tOutput->Branch("is_ttH_ttZctrl_SR",&is_ttH_ttZctrl_SR,"is_ttH_ttZctrl_SR/B");
    tOutput->Branch("is_ttH_ttZctrl_SR_Data",&is_ttH_ttZctrl_SR_Data,"is_ttH_ttZctrl_SR_Data/B");
    tOutput->Branch("is_ttH_ttZctrl_Fake",&is_ttH_ttZctrl_Fake,"is_ttH_ttZctrl_Fake/B");
    tOutput->Branch("is_ttH_WZctrl",&is_ttH_WZctrl,"is_ttH_WZctrl/B");
    tOutput->Branch("is_ttH_WZctrl_SR",&is_ttH_WZctrl_SR,"is_ttH_WZctrl_SR/B");
    tOutput->Branch("is_ttH_WZctrl_SR_Data",&is_ttH_WZctrl_SR_Data,"is_ttH_WZctrl_SR_Data/B");
    tOutput->Branch("is_ttH_WZctrl_Fake",&is_ttH_WZctrl_Fake,"is_ttH_WZctrl_Fake/B");

	// tOutput->Branch("signal_3l_TT_MVA",&signal_3l_TT_MVA,"signal_3l_TT_MVA/F");
	// tOutput->Branch("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,"signal_3l_TTV_MVA/F");
    // tOutput->Branch("signal_2lss_TT_MVA",&signal_2lss_TT_MVA,"signal_2lss_TT_MVA/F");
    // tOutput->Branch("signal_2lss_TTV_MVA",&signal_2lss_TTV_MVA,"signal_2lss_TTV_MVA/F");


    //-- Input variables from tHq2016 analysis
    tOutput->Branch("nJet25",&nJet25,"nJet25/F");
    tOutput->Branch("nJetLoose",&nJetLoose,"nJetLoose/F");
    tOutput->Branch("maxEtaJet25",&maxEtaJet25,"maxEtaJet25/F");
    tOutput->Branch("lepCharge",&lepCharge,"lepCharge/F");
    tOutput->Branch("nJetEta1",&nJetEta1,"nJetEta1/F");
    tOutput->Branch("dEtaFwdJetBJet",&dEtaFwdJetBJet,"dEtaFwdJetBJet/F");
    tOutput->Branch("dEtaFwdJet2BJet",&dEtaFwdJet2BJet,"dEtaFwdJet2BJet/F");
    tOutput->Branch("dEtaFwdJetClosestLep",&dEtaFwdJetClosestLep,"dEtaFwdJetClosestLep/F");
    tOutput->Branch("dPhiHighestPtSSPair",&dPhiHighestPtSSPair,"dPhiHighestPtSSPair/F");
    tOutput->Branch("minDRll",&minDRll,"minDRll/F");
    tOutput->Branch("Lep3Pt",&Lep3Pt,"Lep3Pt/F");

    //-- Input variables from ttH2017 analysis
    tOutput->Branch("lep1_conePt",&lep1_conePt,"lep1_conePt/F");
    tOutput->Branch("lep2_conePt",&lep2_conePt,"lep2_conePt/F");
    tOutput->Branch("mindr_lep1_jet",&mindr_lep1_jet,"mindr_lep1_jet/F");
    tOutput->Branch("mindr_lep2_jet",&mindr_lep2_jet,"mindr_lep2_jet/F");
    tOutput->Branch("mT_lep1",&mT_lep1,"mT_lep1/F");
    tOutput->Branch("mT_lep2",&mT_lep2,"mT_lep2/F");
    tOutput->Branch("max_lep_eta",&max_lep_eta,"max_lep_eta/F");

    //-- new variables
    tOutput->Branch("minv_FwdJetBJet",&minv_FwdJetBJet,"minv_FwdJetBJet/F");
    tOutput->Branch("etaFwdJet",&etaFwdJet,"etaFwdJet/F");
    tOutput->Branch("ptFwdJet",&ptFwdJet,"ptFwdJet/F");
    tOutput->Branch("LeadJetEta",&LeadJetEta,"LeadJetEta/F");
    tOutput->Branch("LeadJetPt",&LeadJetPt,"LeadJetPt/F");
    tOutput->Branch("dRjj",&dRjj,"dRjj/F");
    tOutput->Branch("deepCSV_max",&deepCSV_max,"deepCSV_max/F");
    tOutput->Branch("deepCSV_2nd",&deepCSV_2nd,"deepCSV_2nd/F");
    tOutput->Branch("Mjj_max",&Mjj_max,"Mjj_max/F");
    tOutput->Branch("dPhiLepBJet_max",&dPhiLepBJet_max,"dPhiLepBJet_max/F");
    tOutput->Branch("dPhijj_max",&dPhijj_max,"dPhijj_max/F");


	//-- More control vars
    tOutput->Branch("inv_mll",&inv_mll,"inv_mll/F");
    tOutput->Branch("hardestBjetPt",&hardestBjetPt,"hardestBjetPt/F");
    tOutput->Branch("hardestBjetEta",&hardestBjetEta,"hardestBjetEta/F");
    tOutput->Branch("fwdJetPt",&fwdJetPt,"fwdJetPt/F");
    tOutput->Branch("lep1Pt",&lep1Pt,"lep1Pt/F");
    tOutput->Branch("lep2Pt",&lep2Pt,"lep2Pt/F");
    tOutput->Branch("lep3Pt",&lep3Pt,"lep3Pt/F");
    tOutput->Branch("MET",&MET,"MET/F");

	//-- Nof objects variables
	tOutput->Branch("nLooseBJets",&nLooseBJets,"nLooseBJets/F");
	tOutput->Branch("nMediumBJets",&nMediumBJets,"nMediumBJets/F");
	tOutput->Branch("nTightLep",&nTightLep,"nTightLep/F");
	tOutput->Branch("nFakeableLep",&nFakeableLep,"nFakeableLep/F");
    tOutput->Branch("nLightJets",&nLightJets,"nLightJets/F");
	tOutput->Branch("nLightJets_Fwd40",&nLightJets_Fwd40,"nLightJets_Fwd40/F");


	//-- Other, MEM necessary vars, ...
    tOutput->Branch("PV_weight",&weight_PV,"PV_weight/F");
    tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/B");

    tOutput->Branch("catJets",&catJets,"catJets/I");

    tOutput->Branch("multilepton_Lepton1_Id",               &multilepton_Lepton1_Id,                "multilepton_Lepton1_Id/I");
    tOutput->Branch("multilepton_Lepton1_P4",               "TLorentzVector",                       &multilepton_Lepton1_P4);
    tOutput->Branch("multilepton_Lepton1_DeltaR_Matched",   &multilepton_Lepton1_DeltaR_Matched,    "multilepton_Lepton1_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton1_Label_Matched",    &multilepton_Lepton1_Label_Matched,     "multilepton_Lepton1_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton1_Id_Matched",       &multilepton_Lepton1_Id_Matched,        "multilepton_Lepton1_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton1_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton1_P4_Matched);
    tOutput->Branch("multilepton_Lepton2_Id",               &multilepton_Lepton2_Id,                "multilepton_Lepton2_Id/I");
    tOutput->Branch("multilepton_Lepton2_P4",               "TLorentzVector",                       &multilepton_Lepton2_P4);
    tOutput->Branch("multilepton_Lepton2_DeltaR_Matched",   &multilepton_Lepton2_DeltaR_Matched,    "multilepton_Lepton2_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton2_Label_Matched",    &multilepton_Lepton2_Label_Matched,     "multilepton_Lepton2_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton2_Id_Matched",       &multilepton_Lepton2_Id_Matched,        "multilepton_Lepton2_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton2_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton2_P4_Matched);
    tOutput->Branch("multilepton_Lepton3_Id",               &multilepton_Lepton3_Id,                "multilepton_Lepton3_Id/I");
    tOutput->Branch("multilepton_Lepton3_P4",               "TLorentzVector",                       &multilepton_Lepton3_P4);
    tOutput->Branch("multilepton_Lepton3_DeltaR_Matched",   &multilepton_Lepton3_DeltaR_Matched,    "multilepton_Lepton3_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton3_Label_Matched",    &multilepton_Lepton3_Label_Matched,     "multilepton_Lepton3_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton3_Id_Matched",       &multilepton_Lepton3_Id_Matched,        "multilepton_Lepton3_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton3_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton3_P4_Matched);
    tOutput->Branch("multilepton_Lepton4_Id",               &multilepton_Lepton4_Id,                "multilepton_Lepton4_Id/I");
    tOutput->Branch("multilepton_Lepton4_P4",               "TLorentzVector",                       &multilepton_Lepton4_P4);
    tOutput->Branch("multilepton_Lepton4_DeltaR_Matched",   &multilepton_Lepton4_DeltaR_Matched,    "multilepton_Lepton4_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton4_Label_Matched",    &multilepton_Lepton4_Label_Matched,     "multilepton_Lepton4_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton4_Id_Matched",       &multilepton_Lepton4_Id_Matched,        "multilepton_Lepton4_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton4_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton4_P4_Matched);

    tOutput->Branch("multilepton_Bjet1_Id",                 &multilepton_Bjet1_Id,                  "multilepton_Bjet1_Id/I");
    tOutput->Branch("multilepton_Bjet1_P4",                 "TLorentzVector",                       &multilepton_Bjet1_P4);
    tOutput->Branch("multilepton_Bjet1_CSV",                &multilepton_Bjet1_CSV,                 "multilepton_Bjet1_CSV/F");
    tOutput->Branch("multilepton_Bjet1_JEC_Up",             &multilepton_Bjet1_JEC_Up,              "multilepton_Bjet1_JEC_Up/F");
    tOutput->Branch("multilepton_Bjet1_JEC_Down",           &multilepton_Bjet1_JEC_Down,            "multilepton_Bjet1_JEC_Down/F");
    tOutput->Branch("multilepton_Bjet1_JER_Up",             &multilepton_Bjet1_JER_Up,              "multilepton_Bjet1_JER_Up/F");
    tOutput->Branch("multilepton_Bjet1_JER_Down",           &multilepton_Bjet1_JER_Down,            "multilepton_Bjet1_JER_Down/F");
    tOutput->Branch("multilepton_Bjet1_DeltaR_Matched",     &multilepton_Bjet1_DeltaR_Matched,      "multilepton_Bjet1_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Bjet1_Label_Matched",      &multilepton_Bjet1_Label_Matched,       "multilepton_Bjet1_Label_Matched/I");
    tOutput->Branch("multilepton_Bjet1_Id_Matched",         &multilepton_Bjet1_Id_Matched,          "multilepton_Bjet1_Id_Matched/I");
    tOutput->Branch("multilepton_Bjet1_P4_Matched",         "TLorentzVector",                       &multilepton_Bjet1_P4_Matched);

    tOutput->Branch("multilepton_Bjet2_Id",                 &multilepton_Bjet2_Id,                  "multilepton_Bjet2_Id/I");
    tOutput->Branch("multilepton_Bjet2_P4",                 "TLorentzVector",                       &multilepton_Bjet2_P4);
    tOutput->Branch("multilepton_Bjet2_CSV",                &multilepton_Bjet2_CSV,                 "multilepton_Bjet2_CSV/F");
    tOutput->Branch("multilepton_Bjet2_JEC_Up",             &multilepton_Bjet2_JEC_Up,              "multilepton_Bjet2_JEC_Up/F");
    tOutput->Branch("multilepton_Bjet2_JEC_Down",           &multilepton_Bjet2_JEC_Down,            "multilepton_Bjet2_JEC_Down/F");
    tOutput->Branch("multilepton_Bjet2_JER_Up",             &multilepton_Bjet2_JER_Up,              "multilepton_Bjet2_JER_Up/F");
    tOutput->Branch("multilepton_Bjet2_JER_Down",           &multilepton_Bjet2_JER_Down,            "multilepton_Bjet2_JER_Down/F");
    tOutput->Branch("multilepton_Bjet2_DeltaR_Matched",     &multilepton_Bjet2_DeltaR_Matched,      "multilepton_Bjet2_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Bjet2_Label_Matched",      &multilepton_Bjet2_Label_Matched,       "multilepton_Bjet2_Label_Matched/I");
    tOutput->Branch("multilepton_Bjet2_Id_Matched",         &multilepton_Bjet2_Id_Matched,          "multilepton_Bjet2_Id_Matched/I");
    tOutput->Branch("multilepton_Bjet2_P4_Matched",         "TLorentzVector",                       &multilepton_Bjet2_P4_Matched);

    tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
    tOutput->Branch("multilepton_JetHighestPt1_CSV",&multilepton_JetHighestPt1_CSV,"multilepton_JetHighestPt1_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt1_JEC_Up",&multilepton_JetHighestPt1_JEC_Up,"multilepton_JetHighestPt1_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_JEC_Down",&multilepton_JetHighestPt1_JEC_Down,"multilepton_JetHighestPt1_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt1_JER_Up",&multilepton_JetHighestPt1_JER_Up,"multilepton_JetHighestPt1_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_JER_Down",&multilepton_JetHighestPt1_JER_Down,"multilepton_JetHighestPt1_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
    tOutput->Branch("multilepton_JetHighestPt2_CSV",&multilepton_JetHighestPt2_CSV,"multilepton_JetHighestPt2_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt2_JEC_Up",&multilepton_JetHighestPt2_JEC_Up,"multilepton_JetHighestPt2_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_JEC_Down",&multilepton_JetHighestPt2_JEC_Down,"multilepton_JetHighestPt2_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt2_JER_Up",&multilepton_JetHighestPt2_JER_Up,"multilepton_JetHighestPt2_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_JER_Down",&multilepton_JetHighestPt2_JER_Down,"multilepton_JetHighestPt2_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
    tOutput->Branch("multilepton_JetClosestMw1_CSV",&multilepton_JetClosestMw1_CSV,"multilepton_JetClosestMw1_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw1_JEC_Up",&multilepton_JetClosestMw1_JEC_Up,"multilepton_JetClosestMw1_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_JEC_Down",&multilepton_JetClosestMw1_JEC_Down,"multilepton_JetClosestMw1_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw1_JER_Up",&multilepton_JetClosestMw1_JER_Up,"multilepton_JetClosestMw1_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_JER_Down",&multilepton_JetClosestMw1_JER_Down,"multilepton_JetClosestMw1_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
    tOutput->Branch("multilepton_JetClosestMw2_CSV",&multilepton_JetClosestMw2_CSV,"multilepton_JetClosestMw2_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw2_JEC_Up",&multilepton_JetClosestMw2_JEC_Up,"multilepton_JetClosestMw2_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_JEC_Down",&multilepton_JetClosestMw2_JEC_Down,"multilepton_JetClosestMw2_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw2_JER_Up",&multilepton_JetClosestMw2_JER_Up,"multilepton_JetClosestMw2_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_JER_Down",&multilepton_JetClosestMw2_JER_Down,"multilepton_JetClosestMw2_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_CSV",&multilepton_JetLowestMjj1_CSV,"multilepton_JetLowestMjj1_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JEC_Up",&multilepton_JetLowestMjj1_JEC_Up,"multilepton_JetLowestMjj1_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JEC_Down",&multilepton_JetLowestMjj1_JEC_Down,"multilepton_JetLowestMjj1_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JER_Up",&multilepton_JetLowestMjj1_JER_Up,"multilepton_JetLowestMjj1_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JER_Down",&multilepton_JetLowestMjj1_JER_Down,"multilepton_JetLowestMjj1_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_CSV",&multilepton_JetLowestMjj2_CSV,"multilepton_JetLowestMjj2_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JEC_Up",&multilepton_JetLowestMjj2_JEC_Up,"multilepton_JetLowestMjj2_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JEC_Down",&multilepton_JetLowestMjj2_JEC_Down,"multilepton_JetLowestMjj2_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JER_Up",&multilepton_JetLowestMjj2_JER_Up,"multilepton_JetLowestMjj2_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JER_Down",&multilepton_JetLowestMjj2_JER_Down,"multilepton_JetLowestMjj2_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestEta1_Id",&multilepton_JetHighestEta1_Id,"multilepton_JetHighestEta1_Id/I");
    tOutput->Branch("multilepton_JetHighestEta1_P4","TLorentzVector",&multilepton_JetHighestEta1_P4);
    tOutput->Branch("multilepton_JetHighestEta1_CSV",&multilepton_JetHighestEta1_CSV,"multilepton_JetHighestEta1_CSV/F");
    tOutput->Branch("multilepton_JetHighestEta2_Id",&multilepton_JetHighestEta2_Id,"multilepton_JetHighestEta2_Id/I");
    tOutput->Branch("multilepton_JetHighestEta2_P4","TLorentzVector",&multilepton_JetHighestEta2_P4);
    tOutput->Branch("multilepton_JetHighestEta2_CSV",&multilepton_JetHighestEta2_CSV,"multilepton_JetHighestEta2_CSV/F");
    tOutput->Branch("multilepton_JetHighestEta1_JEC_Up",&multilepton_JetHighestEta1_JEC_Up,"multilepton_JetHighestEta1_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestEta1_JEC_Down",&multilepton_JetHighestEta1_JEC_Down,"multilepton_JetHighestEta1_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestEta1_JER_Up",&multilepton_JetHighestEta1_JER_Up,"multilepton_JetHighestEta1_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestEta1_JER_Down",&multilepton_JetHighestEta1_JER_Down,"multilepton_JetHighestEta1_JER_Down/F");
    tOutput->Branch("multilepton_JetHighestEta2_JEC_Up",&multilepton_JetHighestEta2_JEC_Up,"multilepton_JetHighestEta2_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestEta2_JEC_Down",&multilepton_JetHighestEta2_JEC_Down,"multilepton_JetHighestEta2_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestEta2_JER_Up",&multilepton_JetHighestEta2_JER_Up,"multilepton_JetHighestEta2_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestEta2_JER_Down",&multilepton_JetHighestEta2_JER_Down,"multilepton_JetHighestEta2_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_Id",&multilepton_JetHighestPt1_2ndPair_Id,"multilepton_JetHighestPt1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt1_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_CSV",&multilepton_JetHighestPt1_2ndPair_CSV,"multilepton_JetHighestPt1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Up",&multilepton_JetHighestPt1_2ndPair_JEC_Up,"multilepton_JetHighestPt1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Down",&multilepton_JetHighestPt1_2ndPair_JEC_Down,"multilepton_JetHighestPt1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JER_Up",&multilepton_JetHighestPt1_2ndPair_JER_Up,"multilepton_JetHighestPt1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JER_Down",&multilepton_JetHighestPt1_2ndPair_JER_Down,"multilepton_JetHighestPt1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_Id",&multilepton_JetHighestPt2_2ndPair_Id,"multilepton_JetHighestPt2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt2_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_CSV",&multilepton_JetHighestPt2_2ndPair_CSV,"multilepton_JetHighestPt2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Up",&multilepton_JetHighestPt2_2ndPair_JEC_Up,"multilepton_JetHighestPt2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Down",&multilepton_JetHighestPt2_2ndPair_JEC_Down,"multilepton_JetHighestPt2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JER_Up",&multilepton_JetHighestPt2_2ndPair_JER_Up,"multilepton_JetHighestPt2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JER_Down",&multilepton_JetHighestPt2_2ndPair_JER_Down,"multilepton_JetHighestPt2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_Id",&multilepton_JetClosestMw1_2ndPair_Id,"multilepton_JetClosestMw1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw1_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_CSV",&multilepton_JetClosestMw1_2ndPair_CSV,"multilepton_JetClosestMw1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Up",&multilepton_JetClosestMw1_2ndPair_JEC_Up,"multilepton_JetClosestMw1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Down",&multilepton_JetClosestMw1_2ndPair_JEC_Down,"multilepton_JetClosestMw1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JER_Up",&multilepton_JetClosestMw1_2ndPair_JER_Up,"multilepton_JetClosestMw1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JER_Down",&multilepton_JetClosestMw1_2ndPair_JER_Down,"multilepton_JetClosestMw1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_Id",&multilepton_JetClosestMw2_2ndPair_Id,"multilepton_JetClosestMw2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw2_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_CSV",&multilepton_JetClosestMw2_2ndPair_CSV,"multilepton_JetClosestMw2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Up",&multilepton_JetClosestMw2_2ndPair_JEC_Up,"multilepton_JetClosestMw2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Down",&multilepton_JetClosestMw2_2ndPair_JEC_Down,"multilepton_JetClosestMw2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JER_Up",&multilepton_JetClosestMw2_2ndPair_JER_Up,"multilepton_JetClosestMw2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JER_Down",&multilepton_JetClosestMw2_2ndPair_JER_Down,"multilepton_JetClosestMw2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_Id",&multilepton_JetLowestMjj1_2ndPair_Id,"multilepton_JetLowestMjj1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj1_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_CSV",&multilepton_JetLowestMjj1_2ndPair_CSV,"multilepton_JetLowestMjj1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Up",&multilepton_JetLowestMjj1_2ndPair_JEC_Up,"multilepton_JetLowestMjj1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Down",&multilepton_JetLowestMjj1_2ndPair_JEC_Down,"multilepton_JetLowestMjj1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Up",&multilepton_JetLowestMjj1_2ndPair_JER_Up,"multilepton_JetLowestMjj1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Down",&multilepton_JetLowestMjj1_2ndPair_JER_Down,"multilepton_JetLowestMjj1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_Id",&multilepton_JetLowestMjj2_2ndPair_Id,"multilepton_JetLowestMjj2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj2_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_CSV",&multilepton_JetLowestMjj2_2ndPair_CSV,"multilepton_JetLowestMjj2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Up",&multilepton_JetLowestMjj2_2ndPair_JEC_Up,"multilepton_JetLowestMjj2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Down",&multilepton_JetLowestMjj2_2ndPair_JEC_Down,"multilepton_JetLowestMjj2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Up",&multilepton_JetLowestMjj2_2ndPair_JER_Up,"multilepton_JetLowestMjj2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Down",&multilepton_JetLowestMjj2_2ndPair_JER_Down,"multilepton_JetLowestMjj2_2ndPair_JER_Down/F");

    // Test adding truth information

    tOutput->Branch("multilepton_h0_Id",                    &multilepton_h0_Id,                 "multilepton_h0_Id/I");
    tOutput->Branch("multilepton_h0_P4",                    "TLorentzVector",                   &multilepton_h0_P4);
    tOutput->Branch("multilepton_h0_Label",                    &multilepton_h0_Label,                 "multilepton_h0_Label/I");
    tOutput->Branch("multilepton_t1_Id",                    &multilepton_t1_Id,                 "multilepton_t1_Id/I");
    tOutput->Branch("multilepton_t1_P4",                    "TLorentzVector",                   &multilepton_t1_P4);
    tOutput->Branch("multilepton_t1_Label",                    &multilepton_t1_Label,                 "multilepton_h0_Label/I");
    tOutput->Branch("multilepton_t2_Id",                    &multilepton_t2_Id,                 "multilepton_t2_Id/I");
    tOutput->Branch("multilepton_t2_P4",                    "TLorentzVector",                   &multilepton_t2_P4);
    tOutput->Branch("multilepton_t2_Label",                    &multilepton_t2_Label,                 "multilepton_h0_Label/I");

    // End test adding truth information

    tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
    tOutput->Branch("multilepton_mETcov00",&multilepton_mETcov00,"multilepton_mETcov00/D");
    tOutput->Branch("multilepton_mETcov01",&multilepton_mETcov01,"multilepton_mETcov01/D");
    tOutput->Branch("multilepton_mETcov10",&multilepton_mETcov10,"multilepton_mETcov10/D");
    tOutput->Branch("multilepton_mETcov11",&multilepton_mETcov11,"multilepton_mETcov11/D");
    tOutput->Branch("multilepton_mHT",&multilepton_mHT,"multilepton_mHT/F");
    tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

    return;
}
































//--------------------------------------------
// ##     ## ######## ##       ########  ######## ########     ######## ##     ## ##    ##  ######   ######
// ##     ## ##       ##       ##     ## ##       ##     ##    ##       ##     ## ###   ## ##    ## ##    ##
// ##     ## ##       ##       ##     ## ##       ##     ##    ##       ##     ## ####  ## ##       ##
// ######### ######   ##       ########  ######   ########     ######   ##     ## ## ## ## ##        ######
// ##     ## ##       ##       ##        ##       ##   ##      ##       ##     ## ##  #### ##             ##
// ##     ## ##       ##       ##        ##       ##    ##     ##       ##     ## ##   ### ##    ## ##    ##
// ##     ## ######## ######## ##        ######## ##     ##    ##        #######  ##    ##  ######   ######
//--------------------------------------------


bool Pass_FakeApp_Requirements(vector<Lepton> vLeptonFakeable)
{

}




void tHqMultileponAnalysis::FillJetInfoOutputTree(int* tree_Id, int Id, TLorentzVector* tree_P4, TLorentzVector P4, float* tree_CSV, float CSV, float* tree_JEC_Up, float* tree_JEC_Down, float JEC_value, float* tree_JER_Up, float* tree_JER_Down, float JER, float JER_Up, float JER_Down)
{
    *tree_Id = Id;
    *tree_P4 = P4;

    *tree_CSV = CSV;

    *tree_JEC_Up = P4.E()*(1.+JEC_value);
    *tree_JEC_Down = P4.E()*(1.-JEC_value);

    *tree_JER_Up = P4.E()*JER_Up/JER;
    *tree_JER_Down = P4.E()*JER_Down/JER;

    return;
}



 //  ####  ###### #             #####       # ###### #####  ####
 // #      #      #             #    #      # #        #   #
 //  ####  #####  #             #####       # #####    #    ####
 //      # #      #      ###    #    #      # #        #        #
 // #    # #      #      ###    #    # #    # #        #   #    #
 //  ####  ###### ###### ###    #####   ####  ######   #    ####


/**
//Selects the two highest b-tag jets. If only one b-tag, selects just this one.
 * @param vJetLoose     Contains all selected jets (b & light )
 * @param ibsel1            Return bjet index 1
 * @param ibsel2            Return bjet index 2
 * @param doSelectOnlyBjets True if consider only bjets in vector
 */
void tHqMultileponAnalysis::SelectBjets(int &ibsel1, int &ibsel2, bool doSelectOnlyBjets=true)
{
    int ib1=-1, ib2=-1;

    // double CSV_threshold = 0.5426; //Loose WP -- old // Now using (DeepCSVbtag = deepCSVb+deepCSVbb)

    Float_t btag_max=-9999, btag_max2=-9999;
    for (int ib=0; ib<vJetLoose->size(); ib++)

    {
        if(doSelectOnlyBjets && !vJetLoose->at(ib).isLooseBTag) {continue;}

        if(vJetLoose->at(ib).DeepCSVbtag>btag_max)
        {
          btag_max2 = btag_max;
          ib2 = ib1;
          btag_max = vJetLoose->at(ib).DeepCSVbtag;
          ib1 = ib;
        }
        else if (vJetLoose->at(ib).DeepCSVbtag<btag_max && vJetLoose->at(ib).DeepCSVbtag>btag_max2)
        {
          btag_max2 = vJetLoose->at(ib).DeepCSVbtag;
          ib2 = ib;
        }
    }


    ibsel1 = ib1;
    ibsel2 = ib2;

    return;
}




//-------------------------------------------
 //  ####  ###### #             #      #  ####  #    # #####         # ###### #####  ####
 // #      #      #             #      # #    # #    #   #           # #        #   #
 //  ####  #####  #             #      # #      ######   #           # #####    #    ####
 //      # #      #      ###    #      # #  ### #    #   #           # #        #        #
 // #    # #      #      ###    #      # #    # #    #   #      #    # #        #   #    #
 //  ####  ###### ###### ###    ###### #  ####  #    #   #       ####  ######   #    ####
//-------------------------------------------
/**
 * Determine highest pt, 2nd highest pt, highest eta jets, dijet with inv. mass closest to W mass, dijet with lowest inv mass
 * @param vJetLoose [Contains all jets (btagged or not) passing the object selection]
 */
void tHqMultileponAnalysis::SelectOtherJets(const int ib1, const int ib2, int &ij1, int &ij2, int &ik1, int &ik2, int &ie1, int &ie2, int &il1, int &il2, int &im1, int &im2, int &io1, int &io2, int &ip1, int &ip2)
{
    TLorentzVector Pjet1, Pjet2;
    float pt_max=0, pt_max2=0;
    float diffmass_min = 10000, mass_min = 10000;
    float eta_max=0, eta_max2=0;

    for(int ij=0; ij<vJetLoose->size(); ij++)
    {
        if (ij==ib1 || ij==ib2) {continue;} //Don't take bjets into account

        if(vJetLoose->at(ij).pt > pt_max ) //Highest pT
        {
            pt_max2 = pt_max;
            ij2 = ij1;
            pt_max = vJetLoose->at(ij).pt;
            ij1 = ij;
        }
        else if(vJetLoose->at(ij).pt > pt_max2) //2nd Highest pT
        {
            pt_max2 = vJetLoose->at(ij).pt;
            ij2 = ij;
        }

        if(fabs(vJetLoose->at(ij).eta) > eta_max ) //Highest eta
        {
            eta_max2 = eta_max;
            ie2 = ie1;
            eta_max = fabs(vJetLoose->at(ij).eta);
            ie1 = ij;
        }
        else if(fabs(vJetLoose->at(ij).eta ) > eta_max2) //2nd Highest eta
        {
            eta_max2 = fabs(vJetLoose->at(ij).eta);
            ie2 = ij;
        }

        for(int ik=ij+1; ik<vJetLoose->size(); ik++) //Dijet w/ lowest (m - mW)
        {
            if (ik==ib1 || ik==ib2) {continue;} //Don't take bjets into account

            Pjet1.SetPtEtaPhiE(vJetLoose->at(ij).pt, vJetLoose->at(ij).eta, vJetLoose->at(ij).phi, vJetLoose->at(ij).E);
            Pjet2.SetPtEtaPhiE(vJetLoose->at(ik).pt, vJetLoose->at(ik).eta, vJetLoose->at(ik).phi, vJetLoose->at(ik).E);

            if(TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min) //Lowest (m-mW) difference
            {
                ik1=ij;
                ik2=ik;
                diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
            }
            if((Pjet1+Pjet2).M()<mass_min)
            {
                il1=ij;
                il2=ik;
                mass_min = (Pjet1+Pjet2).M();
            }
        }
    }


    //Choose 2 more jets (1rst pair is closest to mW) -- needed for 2l (more jets)
    pt_max=0; pt_max2=0;
    diffmass_min = 10000; mass_min = 10000;
    for(int ij=0; ij<vJetLoose->size(); ij++)
    {
        if (ij==ib1 || ij==ib2 || ij==ik1 || ij==ik2) {continue;} //Don't take bjets and 1rst mW pair into account

        if(vJetLoose->at(ij).pt > pt_max ) //Highest pT
        {
            pt_max2 = pt_max;
            im2 = im1;
            pt_max = vJetLoose->at(ij).pt;
            im1 = ij;
        }
        else if(vJetLoose->at(ij).pt > pt_max2) //2nd Highest pT
        {
            pt_max2 = vJetLoose->at(ij).pt;
            im2 = ij;
        }

        for(int ik=ij+1; ik<vJetLoose->size(); ik++) //Dijet w/ lowest (m - mW)
        {
            if (ik==ib1 || ik==ib2 || ik==ik1 || ik==ik2) {continue;} //Don't take bjets and 1rst pair into account

            Pjet1.SetPtEtaPhiE(vJetLoose->at(ij).pt, vJetLoose->at(ij).eta, vJetLoose->at(ij).phi, vJetLoose->at(ij).E);
            Pjet2.SetPtEtaPhiE(vJetLoose->at(ik).pt, vJetLoose->at(ik).eta, vJetLoose->at(ik).phi, vJetLoose->at(ik).E);

            if(TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min) //Lowest (m-mW) difference
            {
                io1=ij;
                io2=ik;
                diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
            }
            if((Pjet1+Pjet2).M()<mass_min)
            {
                ip1=ij;
                ip2=ik;
                mass_min = (Pjet1+Pjet2).M();
            }
        }
    }

    //cout<<ij1<<", "<<ij2<<", "<<ik1<<", "<<ik2<<", "<<ie1<<", "<<ie2<<endl;

    return;
}




float tHqMultileponAnalysis::Phi_0_2Pi(float phi)
{
    float phi_0_2pi = phi;
    if (phi>= TMath::Pi()) phi_0_2pi  -= 2.*TMath::Pi();
    if (phi<0.)            phi_0_2pi  += 2.*TMath::Pi();
    return phi_0_2pi;
}

float tHqMultileponAnalysis::Phi_MPi_Pi(float phi)
{
    float phi_MPi_Pi = phi;
    while(phi_MPi_Pi > TMath::Pi()) {phi_MPi_Pi -= 2.*TMath::Pi();}
    while(phi_MPi_Pi < - TMath::Pi()) {phi_MPi_Pi += 2.*TMath::Pi();}

	return phi_MPi_Pi;
}

float tHqMultileponAnalysis::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
    if(fabs(phi2 - phi1) > pow(10, 3) ) {return -999;} //Some phi values are huge in vTruth (bad init ?)

    float DeltaPhi = Phi_MPi_Pi( phi2 - phi1 );
    return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

//Decides if given sample is to be used for training or not (else, don't call training selection function)
bool tHqMultileponAnalysis::Sample_isUsed_forTraining()
{
	// if(_sampleName.Contains("ttZJets_13TeV_madgraphMLM") || _sampleName.Contains("ttWJets_13TeV_madgraphMLM") || _sampleName.Contains("THQ") || _sampleName.Contains("TTJets") ) //2016 samples
    {
        if(_sampleName.Contains("ttZ", TString::kIgnoreCase) || _sampleName.Contains("ttW", TString::kIgnoreCase) || _sampleName.Contains("THQ", TString::kIgnoreCase) || _sampleName.Contains("THW", TString::kIgnoreCase) || _sampleName.Contains("ttH", TString::kIgnoreCase) || _sampleName.Contains("TTJet", TString::kIgnoreCase) || _sampleName.Contains("TTTo2", TString::kIgnoreCase) || _sampleName.Contains("TTToSemi", TString::kIgnoreCase) ) //2017 samples
        return true;
    }

    return false;
}


//Func. taken from ttH2017
double tHqMultileponAnalysis::comp_MT_met_lep(TLorentzVector leptonP4, double met_pt, double met_phi)
{
    double met_px = met_pt*std::cos(met_phi);
    double met_py = met_pt*std::sin(met_phi);
    double mT = std::sqrt(std::max(0., pow(leptonP4.Et() + met_pt, 2) - (pow(leptonP4.Px() + met_px, 2) + pow(leptonP4.Py() + met_py, 2))));
    return mT;
}

//Compute met_LD variable, as it's done in Sync.cxx (should include it in NTP code)
double tHqMultileponAnalysis::Compute_METLD(double MET)
{
    TLorentzVector jet;
    float jet_px = 0;
    float jet_py = 0;
    for(int ij=0;ij<vJetLoose->size();ij++)
    {
        jet.SetPtEtaPhiE(vJetLoose->at(ij).pt, vJetLoose->at(ij).eta, vJetLoose->at(ij).phi, vJetLoose->at(ij).E);
        jet_px += jet.Px();
        jet_py += jet.Py();
    }
    for(int ij=0;ij<vLeptonFakeable.size();ij++)
    {
        jet.SetPtEtaPhiE(vLeptonFakeable.at(ij).pt, vLeptonFakeable.at(ij).eta, vLeptonFakeable.at(ij).phi, vLeptonFakeable.at(ij).E);
        jet_px += jet.Px();
        jet_py += jet.Py();
    }
    for(int ij=0;ij<vTauLoose->size();ij++)
    {
        jet.SetPtEtaPhiE(vTauLoose->at(ij).pt, vTauLoose->at(ij).eta, vTauLoose->at(ij).phi, vTauLoose->at(ij).E);
        jet_px += jet.Px();
        jet_py += jet.Py();
    }

    double MHT = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );
    return MET*0.00397 + MHT*0.00265;
}


/**
 * Read the ttH2017 predefined categories values (implemented at NTP level)
 */
int tHqMultileponAnalysis::Read_Predefined_Categories()
{
    int nlep_tth = 0;

    //Can't convert directly boolean into Char_t
    //NB : using Char_t because booleans should not be used in vectors (later)
    if(vEvent->at(0).is_2lSS) {is_ttH_2lSS                           = 1;}
    if(vEvent->at(0).is_2lSS_SR) {is_ttH_2lSS_SR                     = 1;}
    if(vEvent->at(0).is_2lSS_SR_Data) {is_ttH_2lSS_SR_Data           = 1;}
    if(vEvent->at(0).is_2lSS_Fake) {is_ttH_2lSS_Fake                 = 1;}
    if(vEvent->at(0).is_2lSS_Flip) {is_ttH_2lSS_Flip                 = 1;}
    if(vEvent->at(0).is_2lSS_Flip_Data) {is_ttH_2lSS_Flip_Data       = 1;}
    if(vEvent->at(0).is_3l) {is_ttH_3l                               = 1;}
    if(vEvent->at(0).is_3l_SR) {is_ttH_3l_SR                         = 1;}
    if(vEvent->at(0).is_3l_SR_Data) {is_ttH_3l_SR_Data               = 1;}
    if(vEvent->at(0).is_3l_Fake) {is_ttH_3l_Fake                     = 1;}
    if(vEvent->at(0).is_ttWctrl) {is_ttH_ttWctrl                     = 1;}
    if(vEvent->at(0).is_ttWctrl_SR) {is_ttH_ttWctrl_SR               = 1;}
    if(vEvent->at(0).is_ttWctrl_SR_Data) {is_ttH_ttWctrl_SR_Data     = 1;}
    if(vEvent->at(0).is_ttWctrl_Fake) {is_ttH_ttWctrl_Fake           = 1;}
    if(vEvent->at(0).is_ttWctrl_Flip) {is_ttH_ttWctrl_Flip           = 1;}
    if(vEvent->at(0).is_ttWctrl_Flip_Data) {is_ttH_ttWctrl_Flip_Data = 1;}
    if(vEvent->at(0).is_ttZctrl) {is_ttH_ttZctrl                     = 1;}
    if(vEvent->at(0).is_ttZctrl_SR) {is_ttH_ttZctrl_SR               = 1;}
    if(vEvent->at(0).is_ttZctrl_SR_Data) {is_ttH_ttZctrl_SR_Data     = 1;}
    if(vEvent->at(0).is_ttZctrl_Fake) {is_ttH_ttZctrl_Fake           = 1;}
    if(vEvent->at(0).is_WZctrl) {is_ttH_WZctrl                       = 1;}
    if(vEvent->at(0).is_WZctrl_SR) {is_ttH_WZctrl_SR                 = 1;}
    if(vEvent->at(0).is_WZctrl_SR_Data) {is_ttH_WZctrl_SR_Data       = 1;}
    if(vEvent->at(0).is_WZctrl_Fake) {is_ttH_WZctrl_Fake             = 1;}

    //Determine if event is to be saved as a 2lep or 3lep category (needed for variable computation)
    if(is_ttH_2lSS || is_ttH_2lSS_Fake || is_ttH_2lSS_Flip || (is_ttH_2lSS_Flip_Data && _isdata) || is_ttH_ttWctrl || is_ttH_ttWctrl_Fake || is_ttH_ttWctrl_Flip || (is_ttH_ttWctrl_Flip_Data && _isdata) )
    {
        nlep_tth = 2;
    }
    else if(is_ttH_3l || is_ttH_3l_Fake || is_ttH_ttZctrl || is_ttH_ttZctrl_Fake) //WZ CR not implemented
    {
        nlep_tth = 3;
    }

    return nlep_tth;
}


/**
 * Set the value of "weightfake" variable. Weight read from ttH FR file
 * -- nLep_cat : "3l" or "2l" event
 */
void tHqMultileponAnalysis::Apply_FakeRate_Weight(TString nLep_cat)
{
    // if(!_isdata) {return;}

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;

    int nlep = 2;
    if(nLep_cat == "3l") {nlep = 3;}
    else if(nLep_cat != "2l") {cout<<"Wrong FakeRate category !"<<endl; return;}

    //Consider only 2 or 3 hardest FO leptons
    for(int i=0; i<nlep; i++)
    {
        if( !vLeptonFakeable.at(i).isTightTTH ) //For each lepton failing the tight requirements, multiply event by a weight
        {
            leptonsPts.push_back(vLeptonFakeable.at(i).pt);
            leptonsEtas.push_back(vLeptonFakeable.at(i).eta);
            leptonsIds.push_back(vLeptonFakeable.at(i).id);
        }
    }

    weightfake = get_FR_weight(leptonsPts, leptonsEtas, leptonsIds);

    // cout<<"weightfake = "<<weightfake<<endl;

    return;
}

/**
 * Set the value of "weightflip" variable. Weight read from ttH FR file
 * -- nLep_cat : "3l" or "2l" event
 */
void tHqMultileponAnalysis::Apply_FlipRate_Weight()
{
    // if(!_isdata) {return;}

    weightflip = get_QF_weight(vLeptonFakeable.at(0).pt, fabs(vLeptonFakeable.at(0).eta), vLeptonFakeable.at(0).id, vLeptonFakeable.at(1).pt, fabs(vLeptonFakeable.at(1).eta), vLeptonFakeable.at(1).id);

    // cout<<"weightflip = "<<weightflip<<endl;

    return;
}



double tHqMultileponAnalysis::Get_Distance(double rec_pt, double rec_eta, double rec_phi, double gen_pt, double gen_eta, double gen_phi)
{
    double dr = GetDeltaR(rec_eta, rec_phi, gen_eta, gen_phi);
    double dptRel = fabs(rec_pt - gen_pt) / gen_pt;

    return dr + 0.2 * dptRel;
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
 * Apply scale factors (lepton, trigger, ...) to weight, using the ScaleFactors object
 * @param nlep : read different files depending if 2lss or 3l
 */
void tHqMultileponAnalysis::Apply_ScaleFactors(int nlep)
{
    if(nlep != 2 && nlep != 3) {cout<<"Error : wrong nlep value !"<<endl; return;}

    //Lepton SF
    //--------------------------------------------
    float lepton_SF = 1;

    lepton_SF*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[0].id, vLeptonFakeable[0].pt, vLeptonFakeable[0].eta);
    lepton_SF*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[1].id, vLeptonFakeable[1].pt, vLeptonFakeable[1].eta);

    if(nlep == 3) {lepton_SF*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[2].id, vLeptonFakeable[2].pt, vLeptonFakeable[2].eta);}

    // cout<<"lepton_SF = "<<lepton_SF<<endl;
    // cout<<"Weight : "<<weight<<" --> "<<weight * lepton_SF<<endl;

    // weight*= lepton_SF;
    //--------------------------------------------


    //Trigger SF
    //--------------------------------------------
    float trigger_SF = sf->Get_Trigger_SF(nlep, vLeptonFakeable[0].id, vLeptonFakeable[0].pt, vLeptonFakeable[1].id, vLeptonFakeable[1].pt);

    // cout<<"trigger_SF = "<<trigger_SF<<endl;
    // cout<<"Weight : "<<weight<<" --> "<<weight * trigger_SF<<endl;

    // weight*= trigger_SF;
    //--------------------------------------------



    return;
}
