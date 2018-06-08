//ThreeLeptonSelection_THQ3l
//TwoLeptonSelection_THQ2l
//Selection_GammaConv
// tHqMultileponAnalysis::fillOutputTree()

#include "../include/tHqMultileponAnalysis.h"

#include "TSystem.h"
#include "SignalExtractionMVA.cxx"
#include "Helper.cxx"
#include "BTagging.cxx"
#include "FakeRate.cxx"
#include "ChargeFlip.cxx"

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

#define DEBUG           false

using namespace std;

//Cutflow vector --> Count events for THQ_3l_SR selection, after each cut
vector<double> v_cutflow(15);
// bool is_debug;
// vector<double> list_ids;
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
tHqMultileponAnalysis::tHqMultileponAnalysis()
{
    cout<<BOLD(FRED("DO NOT USE THE DEFAULT CONSTRUCTOR !"))<<endl<<endl;
}

/**
 * Destructor
 */
tHqMultileponAnalysis::~tHqMultileponAnalysis()
{
    delete vEvent; delete vElectron; delete vMuon; delete vTau; delete vJet; delete vTruth;

    delete tOutput;
    delete outputfile;

    delete f_CSVwgt_HF; delete f_CSVwgt_LF; delete f_BTag_eff; delete f_QFwgt; delete f_FRwgt;
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

    vEvent    = new std::vector<Event>();
    vElectron = new std::vector<Electron>();
    vMuon     = new std::vector<Muon>();
    vTau      = new std::vector<Tau>();
    vJet      = new std::vector<Jet>();
    vTruth    = new std::vector<Truth>();

    fChain->SetBranchAddress("Event",    &vEvent   );
    fChain->SetBranchAddress("Electron", &vElectron);
    fChain->SetBranchAddress("Muon",     &vMuon    );
    fChain->SetBranchAddress("Tau",      &vTau     );
    fChain->SetBranchAddress("Jet",      &vJet     );
    fChain->SetBranchAddress("Truth",    &vTruth   );
}

/**
 * Open the files containing weights and corrections
 */
void tHqMultileponAnalysis::InitFiles()
{
    //--- Load MVA weight files
    Load_MVA();

    //--- Loading weight files and creating corresponding histograms

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
    std::string inputFileBTagEff = "/opt/sbg/scratch1/cms/TTH/weight/output_tt_effBtag.root";
    f_BTag_eff = new TFile ((inputFileBTagEff).c_str());
    fill_eff_btagging_histos(f_BTag_eff);

    // charge flip
    // std::string inputFileQF = "/opt/sbg/scratch1/cms/TTH/weightMoriond2017/FakeFlip/QF_data_el.root";
    std::string inputFileQF = "/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleAnalyzer/test/FR_weights/QF_data_el.root";
    f_QFwgt    = new TFile ((inputFileQF).c_str());
    fillQFhistos(f_QFwgt);

    // fake rate
    // std::string inputFileFR = "/opt/sbg/scratch1/cms/TTH/weightMoriond2017/FakeFlip/FR_data_ttH_mva.root";
    std::string inputFileFR = "/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleAnalyzer/test/FR_weights/FR_data_ttH_mva.root"; //FIXME
    f_FRwgt    = new TFile ((inputFileFR).c_str());
    fillFRhistos(f_FRwgt);
}


/**
 * Re-initialize vectors of objects
 */
void tHqMultileponAnalysis::InitCollections()
{
    //--- Re-init all vectors
    vLeptons.clear();
    vSelectedMuons.clear();
    vSelectedElectrons.clear();
    vSelectedTaus.clear();
    vSelectedLeptons.clear();
    vFakeMuons.clear();
    vFakeElectrons.clear();
    vFakeLeptons.clear();
    vFakeableLeptons.clear();
    vLooseLeptons.clear();
    vTightLeptons.clear();

    vSelectedJets.clear();
    vLooseBTagJets.clear();
    vLightJets.clear();
    vLightJets_FwdPtCut.clear();
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

    nJet25 = 0.;
    maxEtaJet25 = 0.;
    lepCharge = 0.;
    nJetEta1  = 0.;
    dEtaFwdJetBJet = 0.;
    dEtaFwdJet2BJet = 0.;
    dEtaFwdJetClosestLep = 0.;
    dPhiHighestPtSSPair = 0.;
    minDRll = 0.;
    Lep3Pt = 0.;

    lep1Pt = 0; lep2Pt = 0; lep3Pt = 0; inv_mll = 0; hardestBjetPt = 0; hardestBjetEta = 0; fwdJetPt = 0;

    //Some categories (e.g. 3l SR & Training) overlap, whereas all the others are orthogonal
    //---> Make sure that when we fill the TTree, only 1 category is TRUE !
    is_3l_THQ_SR = 0;
    is_3l_THQ_Training = 0;
    is_3l_Z_CR = 0;
    is_3l_AppFakes_SR = 0;
    is_3l_GammaConv_SR = 0;
    is_2l_THQ_SR = 0;
    is_2l_THQ_Training = 0;
    is_2l_EMU_CR = 0;
    is_2l_AppFakes_SR = 0;
    is_2l_GammaConv_SR = 0;
    is_2l_QFlip_SR = 0;

    // is_debug = false;

	return;
}




















//--------------------------------------------
// ########  ########    ###    ########     ######## ##     ## ######## ##    ## ########  ######
// ##     ## ##         ## ##   ##     ##    ##       ##     ## ##       ###   ##    ##    ##    ##
// ##     ## ##        ##   ##  ##     ##    ##       ##     ## ##       ####  ##    ##    ##
// ########  ######   ##     ## ##     ##    ######   ##     ## ######   ## ## ##    ##     ######
// ##   ##   ##       ######### ##     ##    ##        ##   ##  ##       ##  ####    ##          ##
// ##    ##  ##       ##     ## ##     ##    ##         ## ##   ##       ##   ###    ##    ##    ##
// ##     ## ######## ##     ## ########     ########    ###    ######## ##    ##    ##     ######
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

    std::cout << "Number of input events = " << nentries << std::endl;
    std::cout << "Will process " << nentries_max << " events" << std::endl << std::endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries_max;jentry++)
    {
        // cout<<endl<<jentry<< " / "<<nentries_max<<endl;

        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%200000 == 0) std::cout << "number of processed events " << jentry << std::endl;

        if(DEBUG) std::cout << "New event ===" << std::endl;

        //if(jentry > 100000) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        // int pvn = vEvent->at(0).pv_n();

        event_id = vEvent->at(0).id;
        // if(event_id < 0) cout<<"event id "<<event_id<<endl;

        MET = vEvent->at(0).metpt;

        event_run = vEvent->at(0).run;

        weights_pdf.clear();
        ids_pdf.clear();

        if(!_isdata )
        {
            weight = _lumi*_xsec/_nowe;

            mc_weight = vEvent->at(0).mc_weight;
            weight = weight * mc_weight;

            // weight_scale_muF0p5 = vEvent->at(0).weight_scale_muF0p5();
            // weight_scale_muF2 = vEvent->at(0).weight_scale_muF2();
            // weight_scale_muR0p5 = vEvent->at(0).weight_scale_muR0p5();
            // weight_scale_muR2 = vEvent->at(0).weight_scale_muR2();
            // weights_pdf = vEvent->at(0).pdf_weights();
            // ids_pdf = vEvent->at(0).pdf_ids();
        }
        else
        {
            weight    = 1.;
            mc_weight = 1.;
            weight_PV = 1.;
			weightfake = 1.;
			weightflip = 1.;
        }





        // ###########################################################
        // #  _       _ _   _       _ _           _   _              #
        // # (_)_ __ (_) |_(_) __ _| (_)___  __ _| |_(_) ___  _ __   #
        // # | | '_ \| | __| |/ _` | | / __|/ _` | __| |/ _ \| '_ \  #
        // # | | | | | | |_| | (_| | | \__ \ (_| | |_| | (_) | | | | #
        // # |_|_| |_|_|\__|_|\__,_|_|_|___/\__,_|\__|_|\___/|_| |_| #
        // #                                                         #
        // ###########################################################

        InitCollections(); //Re-init object collections for each event


        // ######################################
        // #  _        _                        #
        // # | |_ _ __(_) __ _  __ _  ___ _ __  #
        // # | __| '__| |/ _` |/ _` |/ _ \ '__| #
        // # | |_| |  | | (_| | (_| |  __/ |    #
        // #  \__|_|  |_|\__, |\__, |\___|_|    #
        // #             |___/ |___/            #
        // #                                    #
        // ######################################

        {
            bool TRIGm   = vEvent->at(0).trig_m  ;
            bool TRIGe   = vEvent->at(0).trig_e  ;
            bool TRIGee  = vEvent->at(0).trig_ee ;
            bool TRIGmm  = vEvent->at(0).trig_mm ;
            bool TRIGme  = vEvent->at(0).trig_em ;
            bool TRIGeee = vEvent->at(0).trig_eee;
            bool TRIGmme = vEvent->at(0).trig_emm;
            bool TRIGeem = vEvent->at(0).trig_eem;
            bool TRIGmmm = vEvent->at(0).trig_mmm;

            bool E = false, M = false, EE = false, MM = false, EM = false;
            if(TRIGme || TRIGeem || TRIGmme ) {EM = true;}
            if(TRIGmm || TRIGmmm)                       {MM = true;}
            if(TRIGee || TRIGeee)	                    {EE = true;}
            if(TRIGm)                                   {M  = true;}
            if(TRIGe)                                   {E  = true;}

            bool emdataset = _sampleName.Contains("MuonE");
            bool mmdataset = _sampleName.Contains("DoubleM");
            bool eedataset = _sampleName.Contains("DoubleE");
            bool mdataset  = _sampleName.Contains("SingleM");
            bool edataset  = _sampleName.Contains("SingleE");

            is_trigger = false;

            // if(_isdata)
            // {
            //     if      ( EM  &&                               (emdataset) ) {is_trigger = true;}
            //     else if ( !EM && MM  &&                        (mmdataset) ) {is_trigger = true;}
            //     else if ( !EM && !MM && EE  &&                 (eedataset) ) {is_trigger = true;}
            //     else if ( !EM && !MM && !EE && M  &&           (mdataset ) ) {is_trigger = true;}
            //     else if ( !EM && !MM && !EE && !M && E &&      (edataset ) ) {is_trigger = true;}
            // }
            if(_isdata) //Modified datasets order
            {
                if      ( EE  &&                               (eedataset) ) {is_trigger = true;}
                else if ( !EE && MM  &&                        (mmdataset) ) {is_trigger = true;}
                else if ( !EE && !MM && EM  &&                 (emdataset) ) {is_trigger = true;}
                else if ( !EE && !MM && !EM && E  &&           (edataset ) ) {is_trigger = true;}
                else if ( !EE && !MM && !EM && !E && M &&      (mdataset ) ) {is_trigger = true;}
            }
            else
            {
                if(EM || MM || EE || M || E) {is_trigger = true;}
            }

            //The event will be rejected in sel. func., but counted
            if(!is_trigger)
            {
                weight = 0;
                continue; //Skip data events failing trigger
            }
        }

        // #####################################
        // #  _ __ ___  _   _  ___  _ __  ___  #
        // # | '_ ` _ \| | | |/ _ \| '_ \/ __| #
        // # | | | | | | |_| | (_) | | | \__ \ #
        // # |_| |_| |_|\__,_|\___/|_| |_|___/ #
        // #                                   #
        // #####################################

        for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {
            Lepton l; l.setLepton(&vMuon->at(imuon),imuon,0,1);

            vLeptons.push_back(l);

            if(vMuon->at(imuon).isLooseTTH)
            {
                vLooseLeptons.push_back(l);
            }

            if(vMuon->at(imuon).isFakeableTTH)
            {
                vSelectedLeptons.push_back(l);
                vFakeableLeptons.push_back(l);
            }

            if(vMuon->at(imuon).isTightTTH)
            {
                vTightLeptons.push_back(l);
            }
        }


        // ##############################################
        // #       _           _                        #
        // #   ___| | ___  ___| |_ _ __ ___  _ __  ___  #
        // #  / _ \ |/ _ \/ __| __| '__/ _ \| '_ \/ __| #
        // # |  __/ |  __/ (__| |_| | | (_) | | | \__ \ #
        // #  \___|_|\___|\___|\__|_|  \___/|_| |_|___/ #
        // #                                            #
        // ##############################################

        for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++)
        {
            Lepton l; l.setLepton(&vElectron->at(ielectron),ielectron,1,0);

            vLeptons.push_back(l);

            if(vElectron->at(ielectron).isLooseTTH)
            {
                vLooseLeptons.push_back(l);
            }


            if(vElectron->at(ielectron).isFakeableTTH)
            {
                vSelectedLeptons.push_back(l);
                vFakeableLeptons.push_back(l);
            }

            if(vElectron->at(ielectron).isTightTTH)
            {
                // if(!vElectron->at(ielectron).cutEventSel() || !vElectron->at(ielectron).noLostHits() ) {continue;} //conversion veto, always true for muons
                vTightLeptons.push_back(l);
            }
        }

        // ########################
        // #  _                   #
        // # | |_ __ _ _   _ ___  #
        // # | __/ _` | | | / __| #
        // # | || (_| | |_| \__ \ #
        // #  \__\__,_|\__,_|___/ #
        // #                      #
        // ########################

        //NB : not needed for tHq analysis

        for(unsigned int itau=0; itau < vTau->size() ; itau++)
        {
            Lepton l; l.setLepton(&vTau->at(itau),itau,0,0);

            if( vTau->at(itau).isTightTTH)
            {
                vSelectedTaus.push_back(vTau->at(itau));
            }
        }

        // #############################################
        // #                _           _              #
        // #   ___  _ __ __| | ___ _ __(_)_ __   __ _  #
        // #  / _ \| '__/ _` |/ _ \ '__| | '_ \ / _` | #
        // # | (_) | | | (_| |  __/ |  | | | | | (_| | #
        // #  \___/|_|  \__,_|\___|_|  |_|_| |_|\__, | #
        // #                                    |___/  #
        // #                                           #
        // #############################################

        //Sort object vectors by Pt
        std::sort(vLeptons.begin(), vLeptons.end(), SortingLeptonPt);
        std::sort(vSelectedLeptons.begin(), vSelectedLeptons.end(), SortingLeptonPt);
        std::sort(vFakeableLeptons.begin(), vFakeableLeptons.end(), SortingLeptonPt);
        std::sort(vLooseLeptons.begin(), vLooseLeptons.end(), SortingLeptonPt);
        std::sort(vTightLeptons.begin(), vTightLeptons.end(), SortingLeptonPt);

        nTightLep = vTightLeptons.size();
        nFakeableLep = vFakeableLeptons.size();

        // ################################
        // #                              #
        // #  _           _      _        #
        // # | |__       (_) ___| |_ ___  #
        // # | '_ \ _____| |/ _ \ __/ __| #
        // # | |_) |_____| |  __/ |_\__ \ #
        // # |_.__/     _/ |\___|\__|___/ #
        // #           |__/               #
        // #                              #
        // ################################

        nLooseBJets  = 0;
        nMediumBJets = 0;
        nLightJets = 0;
        nLightJets_Fwd40 = 0;

        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {
            if(vJet->at(ijet).pt < 25) {continue;}
            // if( fabs(vJet->at(ijet).eta ) > 4.7) {continue;} //Need to clean jets with eta>4.7 ?

			//---------
			//All jets with pT>25 in vJets are "selected jets"
			vSelectedJets.push_back(vJet->at(ijet));

            //--- B-tagged jets
            if(vJet->at(ijet).CSVv2 > 0.5426 && fabs(vJet->at(ijet).eta) < 2.4)
            {
                vLooseBTagJets.push_back(vJet->at(ijet)); nLooseBJets++;

                if(vJet->at(ijet).CSVv2 > 0.8484) {nMediumBJets++;}
            }
            //--- Light jets
            else if(vJet->at(ijet).CSVv2 < 0.5426)
            {
                vLightJets.push_back(vJet->at(ijet)); nLightJets++;

                if(fabs(vJet->at(ijet).eta ) < 2.4)
                {
                    vLightJets_FwdPtCut.push_back(vJet->at(ijet)); nLightJets_Fwd40++;
                }
                else if(vJet->at(ijet).pt > 40 && fabs(vJet->at(ijet).eta ) > 2.4)
                {
                    vLightJets_FwdPtCut.push_back(vJet->at(ijet)); nLightJets_Fwd40++;
                }
            }
        }



        //--------------------------------------------
        // #####  ###### #####  #    #  ####     #####  #####  # #    # #####  ####  #    # #####
        // #    # #      #    # #    # #    #    #    # #    # # ##   #   #   #    # #    #   #
        // #    # #####  #####  #    # #         #    # #    # # # #  #   #   #    # #    #   #
        // #    # #      #    # #    # #  ###    #####  #####  # #  # #   #   #    # #    #   #
        // #    # #      #    # #    # #    #    #      #   #  # #   ##   #   #    # #    #   #
        // #####  ###### #####   ####   ####     #      #    # # #    #   #    ####   ####    #
        //--------------------------------------------

        bool debug_printout = false;
        if(debug_printout)
        {
            cout<<endl<<endl<<BOLD(FBLU("------------ EVENT -----------"))<<endl;

            cout<<FYEL("--- Loose Tagged Jets : ")<<endl;
            for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
            {
                cout<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt<<" / eta = "<<vLooseBTagJets.at(ijet).eta<<" / phi = "<<vLooseBTagJets.at(ijet).phi<<" / CSV = "<<vLooseBTagJets.at(ijet).CSVv2<<endl;
            }

            cout<<FYEL("--- Light Jets : ")<<endl;
            for(int ijet=0; ijet<vLightJets.size(); ijet++)
            {
                cout<<ijet<<" pT = "<<vLightJets.at(ijet).pt<<" / eta = "<<vLightJets.at(ijet).eta<<" / phi = "<<vLightJets.at(ijet).phi<<" / CSV = "<<vLightJets.at(ijet).CSVv2<<endl;
            }

            cout<<FYEL("--- Loose Leptons : ")<<endl;
            for(int ilep=0; ilep<vLooseLeptons.size(); ilep++)
            {
                TString type = "";
                if( abs(vLooseLeptons[ilep].id) == 11) type = "Ele";
                else if( abs(vLooseLeptons[ilep].id) == 13) type = "Mu";

                cout<<type<<" "<<ilep<<" FO ? "<<vLooseLeptons[ilep].isFakeableTTH<<", Tight ? "<<vLooseLeptons[ilep].isTightTTH<<" / pT = "<<vLooseLeptons.at(ilep).pt<<" / eta = "<<vLooseLeptons.at(ilep).eta<<" / phi = "<<vLooseLeptons.at(ilep).phi<<" : Charge = "<<vLooseLeptons.at(ilep).charge<<endl;
            }
        }

        /*
        if( _isdata )
        {
            for(int j=0; j<list_ids.size(); j++)
            {
                // if( (int) event_id == (int) list_ids[j])
                if(event_id == list_ids[j])
                {
                    Debug_Selection(jentry);
                    break;
                }
            }

            ThreeLeptonSelection_THQ3l_SR(jentry);
            TwoLeptonSelection_THQ2l_SR(jentry);
            continue;
        }*/


        // ############################################
        // #           _           _   _              #
        // #  ___  ___| | ___  ___| |_(_) ___  _ __   #
        // # / __|/ _ \ |/ _ \/ __| __| |/ _ \| '_ \  #
        // # \__ \  __/ |  __/ (__| |_| | (_) | | | | #
        // # |___/\___|_|\___|\___|\__|_|\___/|_| |_| #
        // #                                          #
        // ############################################

        if(Sample_isUsed_forGammaConv() == true)
        {
            ThreeLeptonSelection_GammaConv(jentry);
            TwoLeptonSelection_GammaConv(jentry);
        }

        //Training
        // if(Sample_isUsed_forTraining() == true) //Only few samples are used for training //FIXME
        {
            ThreeLeptonSelection_THQ3l_Training(jentry);
            TwoLeptonSelection_THQ2l_Training(jentry);
        }
        if(_sampleName.Contains("ttWJets_13TeV_madgraphMLM") || _sampleName.Contains("ttZJets_13TeV_madgraphMLM") ) {continue;} //Used *only* for training

        //Data-driven Application regions
        if( _isdata ) //Apply to datasets only
        {
            ThreeLeptonSelection_ApplicationFakes(jentry);
            TwoLeptonSelection_ApplicationFakes(jentry);

            TwoLeptonSelection_ApplicationChargeFlip(jentry);
        }

        //SR & CR
        ThreeLeptonSelection_THQ3l_SR(jentry);
        ThreeLeptonSelection_THQ3l_Z_CR(jentry);

        TwoLeptonSelection_THQ2l_SR(jentry);
        // TwoLeptonSelection_THQ2l_EMU_OS_CR(jentry); // -- disactivate it for now : only TTbar + ST, >1M events !
    }


	ofstream file_out("cutflow.txt");
    for(int i=0; i<v_cutflow.size(); i++)
    {
		file_out<<"v_cutflow "<<i<<" = "<<v_cutflow[i]<<endl;
    }

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




 //  #####             #####  ######
 // #     # #         #     # #     #
 //       # #         #       #     #
 //  #####  #          #####  ######
 //       # #               # #   #
 // #     # #         #     # #    #
 //  #####  ######     #####  #     #

//---------------------------------------------------------------------------

/**
 * 3l signal region selection
 */
void tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l_SR(int evt)
{
    InitVariables();

	v_cutflow[0]++;

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}
    v_cutflow[1]++;

    // ####################
    // # Common selection #
    // ####################
	//At least 3 FO leptons

    if(vFakeableLeptons.size() < 3) {return;}
    v_cutflow[2]++;

    //Ask exactly 3 tight leptons, which must be the 3 hardest FO leptons
	if(vTightLeptons.size() != 3) {return;}
    if( !vFakeableLeptons.at(0).isTightTTH || !vFakeableLeptons.at(1).isTightTTH || !vFakeableLeptons.at(2).isTightTTH ) {return;}
    // if( !(vFakeableLeptons.at(0).isTightTTH && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) || !(vFakeableLeptons.at(2).isTightTTH && vFakeableLeptons.at(2).cutEventSel() && vFakeableLeptons.at(2).noLostHits()) ) {return;}
    v_cutflow[3]++;


    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 15 || vFakeableLeptons.at(2).pt < 15) {return;}

    if(vSelectedJets.size() < 2) {return;}
    v_cutflow[4]++;

    if(nMediumBJets == 0) {return;}
    v_cutflow[5]++;

    if(vLightJets_FwdPtCut.size() == 0)  {return;}
    // if(vLightJets.size() == 0)  {return;}
    v_cutflow[6]++;

    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}
    v_cutflow[7]++;


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && vFakeableLeptons.at(i).charge == -vFakeableLeptons.at(j).charge )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 15) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}
    v_cutflow[8]++;


//--------------------------------------------
    is_3l_THQ_SR = 1; //3l SR category

    Compute_Variables("3l");

    fillOutputTree();

    return;
}






 //  #####
 // #     # #         ##### #####    ##   # #    # # #    #  ####
 //       # #           #   #    #  #  #  # ##   # # ##   # #    #
 //  #####  #           #   #    # #    # # # #  # # # #  # #
 //       # #           #   #####  ###### # #  # # # #  # # #  ###
 // #     # #           #   #   #  #    # # #   ## # #   ## #    #
 //  #####  ######      #   #    # #    # # #    # # #    #  ####

/**
 * 3l training selection
 */
void tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l_Training(int evt)
{
    InitVariables();

    if(weight==0) {return;}

    if(vFakeableLeptons.size() < 3) {return;}

    //-- Only consider the 3 hardest 'FakeableObject' leptons
    if(vFakeableLeptons.at(0).pt < 20 || vFakeableLeptons.at(1).pt < 10 || vFakeableLeptons.at(2).pt < 10)
    {
        return;
    }

	int jets_tmp = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
    	if( fabs(vSelectedJets.at(ijet).eta) < 2.4 ) jets_tmp++;
    }
    if(jets_tmp < 2) {return;} //At least 2 jets with pT>25 & eta<2.4

    if(nLooseBJets == 0) {return;}

    if(vLightJets.size() == 0) {return;}

    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && vFakeableLeptons.at(i).charge == -vFakeableLeptons.at(j).charge )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


//--------------------------------------------
    is_3l_THQ_Training = 1;

    Compute_Variables("3l");

    fillOutputTree();

    return;
}




 //  #####            #######     #####  ######
 // #     # #              #     #     # #     #
 //       # #             #      #       #     #
 //  #####  #            #       #       ######
 //       # #           #        #       #   #
 // #     # #          #         #     # #    #
 //  #####  ######    #######     #####  #     #

/**
 * 3l Z control region selection
 */
void tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l_Z_CR(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}

    // ####################
    // # Common selection #
    // ####################

    //At least 3 leptons
    if(vFakeableLeptons.size() < 3) {return;}

    //Ask exactly 3 tight leptons, which must be the 3 hardest FO leptons
	if(vTightLeptons.size() != 3) {return;}
    if( !vFakeableLeptons.at(0).isTightTTH || !vFakeableLeptons.at(1).isTightTTH || !vFakeableLeptons.at(2).isTightTTH ) {return;}
	// if( !(vFakeableLeptons.at(0).isTightTTH && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) || !(vFakeableLeptons.at(2).isTightTTH && vFakeableLeptons.at(2).cutEventSel() && vFakeableLeptons.at(2).noLostHits()) ) {return;}

    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 15 || vFakeableLeptons.at(2).pt < 15) {return;}

    if( fabs(vFakeableLeptons.at(0).charge + vFakeableLeptons.at(1).charge + vFakeableLeptons.at(2).charge ) != 1) {return;}

    if(vSelectedJets.size() < 2) {return;}

    if(nLooseBJets == 0) {return;}

    if(vLightJets_FwdPtCut.size() == 0)  {return;}
    // if(vLightJets.size() == 0)  {return;}


    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}


    // ##########
    // # OSSF pair - Z INVERTED veto (selection enriched in ttZ/tZ, WZ events)
    // ##########
    bool pass_inverted_Zveto = false;
    float mll_tmp = 0; float best_mll=999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && vFakeableLeptons.at(i).charge == -vFakeableLeptons.at(j).charge )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if(fabs(mll_tmp - 91.188) < 15)
                {
                    pass_inverted_Zveto = true ;
                    if(fabs(mll_tmp - 91.188) < fabs(best_mll - 91.188)) {best_mll = mll_tmp;}
                }
            }

        }
    }
    if(!pass_inverted_Zveto) {return;}

//--------------------------------------------
    is_3l_Z_CR = 1; //3l SR category

    Compute_Variables("3l");

    fillOutputTree();
}




 //  #####             #####  ######
 // #     # #         #     # #     #
 //       # #         #       #     #
 //  #####  #          #####  ######
 // #       #               # #   #
 // #       #         #     # #    #
 // ####### ######     #####  #     #

/**
 * 2l Signal region selection
 */
void tHqMultileponAnalysis::TwoLeptonSelection_THQ2l_SR(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 leptons
    if(vFakeableLeptons.size() < 2) {return;}

    //==2 tight leptons, which must be the 2 hardest FO leptons
	if(vTightLeptons.size() != 2) {return;}
    if( !vFakeableLeptons.at(0).isTightTTH || !vFakeableLeptons.at(1).isTightTTH ) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    // if(!vFakeableLeptons.at(0).tightCharge || !vFakeableLeptons.at(1).tightCharge ) {return;}

    //SS lepton pair
    if(vFakeableLeptons.at(0).charge * vFakeableLeptons.at(1).charge < 0) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 15) {return;}

    //Jet cuts
    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    if(vLightJets_FwdPtCut.size() == 0)  {return;}
    // if(vLightJets.size() == 0)  {return;}


    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}


    // ##########
    // # OSSF pair - Z veto ee #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && fabs(vFakeableLeptons.at(i).id)==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}

//--------------------------------------------
    is_2l_THQ_SR = 1; //2l SR category

    Compute_Variables("2l");

    fillOutputTree();

    return;
}




 //  #####
 // #     # #         ##### #####    ##   # #    # # #    #  ####
 //       # #           #   #    #  #  #  # ##   # # ##   # #    #
 //  #####  #           #   #    # #    # # # #  # # # #  # #
 // #       #           #   #####  ###### # #  # # # #  # # #  ###
 // #       #           #   #   #  #    # # #   ## # #   ## #    #
 // ####### ######      #   #    # #    # # #    # # #    #  ####

/**
 * 2l training selection
 */
void tHqMultileponAnalysis::TwoLeptonSelection_THQ2l_Training(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 FO leptons
    if(vFakeableLeptons.size() < 2) {return;}


    //At most 2 tight leptons
	if(vTightLeptons.size() > 2) {return;}

    //SS lepton pair
    if(vFakeableLeptons.at(0).charge * vFakeableLeptons.at(1).charge < 0) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    // if(!vFakeableLeptons.at(0).tightCharge ) {return;}
    // if(!vFakeableLeptons.at(1).tightCharge ) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt < 20 || vFakeableLeptons.at(1).pt < 10) {return;}

    //Jet cuts
    int jets_tmp = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if( fabs(vSelectedJets.at(ijet).eta) < 2.4 ) jets_tmp++;
    }
    if(jets_tmp < 2) {return;} //At least 2 jets with pT>25 & eta<2.4

    if(nLooseBJets == 0) {return;}

    if(vLightJets.size() == 0)  {return;}

//--------------------------------------------
    is_2l_THQ_Training = 1; //2l Training category

    Compute_Variables("2l");

    fillOutputTree();

    return;
}






 //  #####                                     #####  ######
 // #     # #         ###### #    # #    #    #     # #     #
 //       # #         #      ##  ## #    #    #       #     #
 //  #####  #         #####  # ## # #    #    #       ######
 // #       #         #      #    # #    #    #       #   #
 // #       #         #      #    # #    #    #     # #    #
 // ####### ######    ###### #    #  ####      #####  #     #

/**
 * 2l EMU OS CR : enriched in ttbar events
 */
void tHqMultileponAnalysis::TwoLeptonSelection_THQ2l_EMU_OS_CR(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 leptons
    if(vFakeableLeptons.size() < 2) {return;}


    //==2 tight leptons, which must be the 2 hardest leptons
	if(vTightLeptons.size() != 2) {return;}
    if( !vFakeableLeptons.at(0).isTightTTH || !vFakeableLeptons.at(1).isTightTTH ) {return;}
    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    // if(!vFakeableLeptons.at(0).tightCharge || !vFakeableLeptons.at(1).tightCharge) {return;}
    if(vFakeableLeptons.at(0).charge * vFakeableLeptons.at(1).charge > 0) {return;}

    //e-mu pair
    if( fabs(vFakeableLeptons.at(0).id * vFakeableLeptons.at(1).id) != 143 ) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 10) {return;}
    if( fabs(vFakeableLeptons.at(1).id) == 11 && vFakeableLeptons.at(1).pt < 15) {return;}

    //Jet cuts
    int jets_tmp = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if( fabs(vSelectedJets.at(ijet).eta) < 2.4 ) jets_tmp++;
    }
    if(jets_tmp < 2) {return;} //At least 2 jets with pT>25 & eta<2.4

    if(nMediumBJets == 0) {return;}

	if(vLightJets_FwdPtCut.size() == 0) {return;}
    // if(vLightJets.size() == 0) {return;}


    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}

    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && fabs(vFakeableLeptons.at(i).id)==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}

//--------------------------------------------
    is_2l_EMU_CR = 1; //2l SR category

    Compute_Variables("2l");

    fillOutputTree();

    return;
}
















//--------------------------------------------
// ########    ###    ##    ## ########  ######        ###    ########  ########  ##       ####
// ##         ## ##   ##   ##  ##       ##    ##      ## ##   ##     ## ##     ## ##        ##
// ##        ##   ##  ##  ##   ##       ##           ##   ##  ##     ## ##     ## ##        ##
// ######   ##     ## #####    ######    ######     ##     ## ########  ########  ##        ##
// ##       ######### ##  ##   ##             ##    ######### ##        ##        ##        ##
// ##       ##     ## ##   ##  ##       ##    ##    ##     ## ##        ##        ##        ##  ###
// ##       ##     ## ##    ## ########  ######     ##     ## ##        ##        ######## #### ###
//--------------------------------------------



 //  #####
 // #     # #         ######   ##   #    # ######  ####
 //       # #         #       #  #  #   #  #      #
 //  #####  #         #####  #    # ####   #####   ####
 //       # #         #      ###### #  #   #           #
 // #     # #         #      #    # #   #  #      #    #
 //  #####  ######    #      #    # #    # ######  ####

/**
 * 3l fakes application region selection
 */
void tHqMultileponAnalysis::ThreeLeptonSelection_ApplicationFakes(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}

    // ####################
    // # Common selection #
    // ####################
	//At least 3 FO leptons

    if(vFakeableLeptons.size() < 3) {return;}

    if(vTightLeptons.size() > 2) {return;} //At most 2 tight leptons

	// if( !(vFakeableLeptons.at(0).isTightTTH && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) || !(vFakeableLeptons.at(2).isTightTTH && vFakeableLeptons.at(2).cutEventSel() && vFakeableLeptons.at(2).noLostHits()) ) {return;}

    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 15 || vFakeableLeptons.at(2).pt < 15) {return;}

    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    if(vLightJets_FwdPtCut.size() == 0)  {return;}
    // if(vLightJets.size() == 0)  {return;}

    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && vFakeableLeptons.at(i).charge == -vFakeableLeptons.at(j).charge )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 15) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    // #########################
    // # fake rate reweighting #
    // #########################

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;

    //only considering the 3 hardest FO leptons
    for(int i=0; i<3; i++)
    {
        // if( !vFakeableLeptons.at(i).isTightTTH ) //For each lepton failing the tight requirements, multiply event by a weight
        {
            leptonsPts.push_back(     vFakeableLeptons.at(i).pt );
            leptonsEtas.push_back(    vFakeableLeptons.at(i).eta );
            leptonsIds.push_back(     vFakeableLeptons.at(i).id );
        }
    }

    weightfake = get_FR_weight(leptonsPts, leptonsEtas, leptonsIds);

    // cout<<"weightfake = "<<weightfake<<endl;

//--------------------------------------------
    is_3l_AppFakes_SR = 1; //3l AppFakes SR category

    Compute_Variables("3l");

    fillOutputTree();

    return;
}





 //  #####            #######
 // #     # #         #         ##   #    # ######  ####
 //       # #         #        #  #  #   #  #      #
 //  #####  #         #####   #    # ####   #####   ####
 // #       #         #       ###### #  #   #           #
 // #       #         #       #    # #   #  #      #    #
 // ####### ######    #       #    # #    # ######  ####

/**
 * 2l fakes application region selection
 */
void tHqMultileponAnalysis::TwoLeptonSelection_ApplicationFakes(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 FO leptons
    if(vFakeableLeptons.size() < 2) {return;}

    if(vTightLeptons.size() > 1) {return;} //at most 1 tight lepton
    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    // if(!vFakeableLeptons.at(0).tightCharge || !vFakeableLeptons.at(1).tightCharge ) {return;}

    //SS lepton pair
    if(vFakeableLeptons.at(0).charge * vFakeableLeptons.at(1).charge < 0) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 15) {return;}

    //Jet cuts
    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    if(vLightJets_FwdPtCut.size() == 0)  {return;}
    // if(vLightJets.size() == 0)  {return;}


    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && fabs(vFakeableLeptons.at(i).id)==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    // #########################
    // # fake rate reweighting #
    // #########################

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;

    //Only consider the 2 hardest FO leptons
    for(int i=0; i<2; i++)
    {
        if( !vFakeableLeptons.at(i).isTightTTH ) //For each lepton failing the tight requirements, multiply event by a weight
        {
            leptonsPts.push_back(     vFakeableLeptons.at(i).pt                );
            leptonsEtas.push_back(    vFakeableLeptons.at(i).eta               );
            leptonsIds.push_back(     vFakeableLeptons.at(i).id                );
        }
    }

    weightfake = get_FR_weight(leptonsPts, leptonsEtas, leptonsIds);

    // cout<<"weightfake = "<<weightfake<<endl;

//--------------------------------------------
    is_2l_AppFakes_SR = 1; //2l SR category

    Compute_Variables("2l");

    fillOutputTree();

    return;
}





















//--------------------------------------------
//  #######  ######## ##       #### ########
// ##     ## ##       ##        ##  ##     ##
// ##     ## ##       ##        ##  ##     ##
// ##     ## ######   ##        ##  ########
// ##  ## ## ##       ##        ##  ##
// ##    ##  ##       ##        ##  ##
//  ##### ## ##       ######## #### ##
//--------------------------------------------

/**
 * 2l charge flip application region selection (data-driven)
 */
void tHqMultileponAnalysis::TwoLeptonSelection_ApplicationChargeFlip(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 FO leptons
    if(vFakeableLeptons.size() < 2) {return;}
    if(vTightLeptons.size() > 2) {return;}

    if( !vFakeableLeptons.at(0).isTightTTH || !vFakeableLeptons.at(1).isTightTTH ) {return;}
    // if( !(vFakeableLeptons.at(0).isTightTTH && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) ) {return;}

    //Neglect uu channel
    if( abs(vFakeableLeptons.at(0).id) == 13 && abs(vFakeableLeptons.at(1).id) == 13 ) {return;}

    //OS lepton pair
    if(vFakeableLeptons.at(0).charge * vFakeableLeptons.at(1).charge > 0) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele) //DISACTIVATE ??
    // if(!vFakeableLeptons.at(0).tightCharge ) {return;}
    // if(!vFakeableLeptons.at(1).tightCharge ) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 15) {return;}

    //Jet cuts
    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    if(vLightJets_FwdPtCut.size() == 0)  {return;}
    // if(vLightJets.size() == 0)  {return;}


    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && fabs(vFakeableLeptons.at(i).id)==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    // #########################
    // # QFlip reweighting     #
    // #########################

    weightflip = get_QF_weight(vFakeableLeptons.at(0).pt, fabs(vFakeableLeptons.at(0).eta), vFakeableLeptons.at(0).id, vFakeableLeptons.at(1).pt, fabs(vFakeableLeptons.at(1).eta), vFakeableLeptons.at(1).id);

    // cout<<"weightflip = "<<weightflip<<endl;
	// weight*= weightflip;


//--------------------------------------------
    is_2l_QFlip_SR = 1; //2l SR category

    Compute_Variables("2l");

    fillOutputTree();

    return;
}

















//--------------------------------------------
//  ######      ###    ##     ## ##     ##    ###             ######   #######  ##    ## ##     ##
// ##    ##    ## ##   ###   ### ###   ###   ## ##           ##    ## ##     ## ###   ## ##     ##
// ##         ##   ##  #### #### #### ####  ##   ##          ##       ##     ## ####  ## ##     ##
// ##   #### ##     ## ## ### ## ## ### ## ##     ## ####### ##       ##     ## ## ## ## ##     ##
// ##    ##  ######### ##     ## ##     ## #########         ##       ##     ## ##  ####  ##   ##
// ##    ##  ##     ## ##     ## ##     ## ##     ##         ##    ## ##     ## ##   ###   ## ##   ###
//  ######   ##     ## ##     ## ##     ## ##     ##          ######   #######  ##    ##    ###    ###
//--------------------------------------------




 //  #####
 // #     # #          ####    ##   #    # #    #   ##
 //       # #         #    #  #  #  ##  ## ##  ##  #  #
 //  #####  #         #      #    # # ## # # ## # #    #
 //       # #         #  ### ###### #    # #    # ######
 // #     # #         #    # #    # #    # #    # #    #
 //  #####  ######     ####  #    # #    # #    # #    #

/**
 * 3l gamma-conv selection
 */
void tHqMultileponAnalysis::ThreeLeptonSelection_GammaConv(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}

    // ####################
    // # Common selection #
    // ####################
	//At least 3 FO leptons

    if(vFakeableLeptons.size() < 3) {return;}

    //Ask exactly 3 tight leptons, which must be the 3 hardest FO leptons
	if(vTightLeptons.size() != 3) {return;} //FIXME
    if( !vFakeableLeptons.at(0).isTightTTH || !vFakeableLeptons.at(1).isTightTTH || !vFakeableLeptons.at(2).isTightTTH ) {return;}

    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 15 || vFakeableLeptons.at(2).pt < 15) {return;}

    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    if(vLightJets_FwdPtCut.size() == 0)  {return;}
    // if(vLightJets.size() == 0)  {return;}

    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && vFakeableLeptons.at(i).charge == -vFakeableLeptons.at(j).charge )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 15) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


//---------------------------------------------------------------------------

// ------------------------------------------------
// #    #   ##   #####  ####  #    # # #    #  ####
// ##  ##  #  #    #   #    # #    # # ##   # #    #
// # ## # #    #   #   #      ###### # # #  # #
// #    # ######   #   #      #    # # #  # # #  ###
// #    # #    #   #   #    # #    # # #   ## #    #
// #    # #    #   #    ####  #    # # #    #  ####
// ------------------------------------------------

    bool is_GammaConv = false;

    if(!_isdata)
    {
        float lep1_dr_gen = 100, lep2_dr_gen = 100 , lep3_dr_gen = 100;
        float lep1_dr_gen_min = 0.3, lep2_dr_gen_min = 0.3, lep3_dr_gen_min = 0.3;
        int lep1_matched = -1, lep2_matched = -1, lep3_matched = -1;

        //Loop on all truth objects
        for(int itruth = 0; itruth < vTruth->at(0).gen_id.size(); itruth++)
        {
            // if( abs(vTruth->at(0).gen_id.at(itruth)) != 22 && abs(vTruth->at(0).gen_id.at(itruth)) != 11 && abs(vTruth->at(0).gen_id.at(itruth)) != 13) {continue;} //Check only gen photons & leptons
            if( abs(vTruth->at(0).gen_id.at(itruth)) != 22) {continue;} //Check only gen photons

            if(vTruth->at(0).gen_pt.at(itruth) < 1) {continue;}
            // if(vTruth->at(0).gen_status().at(itruth) != 1) {continue;}

            if( abs(vFakeableLeptons[0].id) == 11) //If electron
            {
                lep1_dr_gen = GetDeltaR(vTruth->at(0).gen_eta.at(itruth),  vTruth->at(0).gen_phi.at(itruth), vFakeableLeptons.at(0).eta, vFakeableLeptons.at(0).phi );

                if(lep1_dr_gen < lep1_dr_gen_min)
                {
                    if(vFakeableLeptons.at(0).pt > 0.3*vTruth->at(0).gen_pt.at(itruth) && vFakeableLeptons.at(0).pt < 1.5*vTruth->at(0).gen_pt.at(itruth) )
                    {
                        lep1_dr_gen_min = lep1_dr_gen;  lep1_matched = itruth;
                        // cout<<"-- Matched FO 0 to itruth "<<itruth<<endl;
                    }
                }
            }
            if( abs(vFakeableLeptons[1].id) == 11) //If electron
            {
                lep2_dr_gen = GetDeltaR(vTruth->at(0).gen_eta.at(itruth),  vTruth->at(0).gen_phi.at(itruth), vFakeableLeptons.at(1).eta, vFakeableLeptons.at(1).phi );

                if( lep2_dr_gen < lep2_dr_gen_min)
                {
                    if(vFakeableLeptons.at(1).pt > 0.3*vTruth->at(0).gen_pt.at(itruth) && vFakeableLeptons.at(1).pt < 1.5*vTruth->at(0).gen_pt.at(itruth) )
                    {
                        lep2_dr_gen_min = lep2_dr_gen;  lep2_matched = itruth;
                        // cout<<"-- Matched FO 1 to itruth "<<itruth<<endl;
                    }
                }
            }
            if( abs(vFakeableLeptons[2].id) == 11) //If electron
            {
                lep3_dr_gen = GetDeltaR(vTruth->at(0).gen_eta.at(itruth),  vTruth->at(0).gen_phi.at(itruth), vFakeableLeptons.at(2).eta, vFakeableLeptons.at(2).phi );

                if( lep3_dr_gen < lep3_dr_gen_min)
                {
                    if(vFakeableLeptons.at(2).pt > 0.3*vTruth->at(0).gen_pt.at(itruth) && vFakeableLeptons.at(2).pt < 1.5*vTruth->at(0).gen_pt.at(itruth) )
                    {
                        lep3_dr_gen_min = lep3_dr_gen;  lep3_matched = itruth;
                        // cout<<"-- Matched FO 1 to itruth "<<itruth<<endl;
                    }
                }
            }
        }

        // if(lep1_matched == lep2_matched || lep1_matched == lep3_matched || lep3_matched == lep2_matched) {return;}
        // if(lep1_matched >= 0 || lep2_matched >= 0 || lep3_matched >= 0) {is_GammaConv = true;}

        if(lep1_matched >= 0 && abs(vTruth->at(0).gen_id.at(lep1_matched)) == 22) {is_GammaConv = true;}
        else  if(lep2_matched >= 0 && abs(vTruth->at(0).gen_id.at(lep2_matched)) == 22) {is_GammaConv = true;}
        else  if(lep3_matched >= 0 && abs(vTruth->at(0).gen_id.at(lep3_matched)) == 22) {is_GammaConv = true;}
    }

    if(!is_GammaConv) {return;}

//--------------------------------------------
    is_3l_GammaConv_SR = 1; //3l SR category

    Compute_Variables("3l");

    fillOutputTree();

    return;
}




 //  #####
 // #     # #          ####    ##   #    # #    #   ##
 //       # #         #    #  #  #  ##  ## ##  ##  #  #
 //  #####  #         #      #    # # ## # # ## # #    #
 // #       #         #  ### ###### #    # #    # ######
 // #       #         #    # #    # #    # #    # #    #
 // ####### ######     ####  #    # #    # #    # #    #


/**
 * 2l gamma-conv selection
 */
void tHqMultileponAnalysis::TwoLeptonSelection_GammaConv(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}

    // ####################
    // # Common selection #
    // ####################

    //At least 2 leptons
    if(vFakeableLeptons.size() < 2) {return;}

    //==2 tight leptons, which must be the 2 hardest FO leptons
    if(vTightLeptons.size() != 2) {return;} //FIXME
    if( !vFakeableLeptons.at(0).isTightTTH || !vFakeableLeptons.at(1).isTightTTH ) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    // if(!vFakeableLeptons.at(0).tightCharge || !vFakeableLeptons.at(1).tightCharge ) {return;}

    //SS lepton pair
    if(vFakeableLeptons.at(0).charge * vFakeableLeptons.at(1).charge < 0) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt < 25 || vFakeableLeptons.at(1).pt < 15) {return;}

    //Jet cuts
    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    if(vLightJets_FwdPtCut.size() == 0)  {return;}
    // if(vLightJets.size() == 0)  {return;}


    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && fabs(vFakeableLeptons.at(i).id)==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}

//---------------------------------------------------------------------------

// ------------------------------------------------
// #    #   ##   #####  ####  #    # # #    #  ####
// ##  ##  #  #    #   #    # #    # # ##   # #    #
// # ## # #    #   #   #      ###### # # #  # #
// #    # ######   #   #      #    # # #  # # #  ###
// #    # #    #   #   #    # #    # # #   ## #    #
// #    # #    #   #    ####  #    # # #    #  ####
// ------------------------------------------------

    bool is_GammaConv = false;

    if(!_isdata)
    {
        float lep1_dr_gen       = 100.,     lep2_dr_gen     = 100.;
        float lep1_dr_gen_min   = 0.3,     lep2_dr_gen_min = 0.3;
        int   lep1_matched      = -1,       lep2_matched    = -1;

        //Loop on all truth objects
        for(int itruth = 0; itruth < vTruth->at(0).gen_id.size(); itruth++)
        {
            // if( abs(vTruth->at(0).gen_id.at(itruth)) != 22 && abs(vTruth->at(0).gen_id.at(itruth)) != 11 && abs(vTruth->at(0).gen_id.at(itruth)) != 13) {continue;} //Check only gen photons
            if( abs(vTruth->at(0).gen_id.at(itruth)) != 22) {continue;} //Check only gen photons

            if(vTruth->at(0).gen_pt.at(itruth) < 1) {continue;}
            // if(vTruth->at(0).gen_status().at(itruth) != 1) {continue;}

            if( abs(vFakeableLeptons[0].id) == 11) //If electron
            {
                lep1_dr_gen = GetDeltaR(vTruth->at(0).gen_eta.at(itruth),  vTruth->at(0).gen_phi.at(itruth), vFakeableLeptons.at(0).eta, vFakeableLeptons.at(0).phi );

                if(lep1_dr_gen < lep1_dr_gen_min)
                {
                    if(vFakeableLeptons.at(0).pt > 0.3*vTruth->at(0).gen_pt.at(itruth) && vFakeableLeptons.at(0).pt < 1.5*vTruth->at(0).gen_pt.at(itruth) )
                    {
                        lep1_dr_gen_min = lep1_dr_gen;  lep1_matched = itruth;
                        // cout<<"-- Matched FO 0 to itruth "<<itruth<<endl;
                    }
                }
            }
            if( abs(vFakeableLeptons[1].id) == 11) //If electron
            {
                lep2_dr_gen = GetDeltaR(vTruth->at(0).gen_eta.at(itruth),  vTruth->at(0).gen_phi.at(itruth), vFakeableLeptons.at(1).eta, vFakeableLeptons.at(1).phi );
                if( lep2_dr_gen < lep2_dr_gen_min)
                {
                    if(vFakeableLeptons.at(1).pt > 0.3*vTruth->at(0).gen_pt.at(itruth) && vFakeableLeptons.at(1).pt < 1.5*vTruth->at(0).gen_pt.at(itruth) )
                    {
                        lep2_dr_gen_min = lep2_dr_gen;  lep2_matched = itruth;
                        // cout<<"-- Matched FO 1 to itruth "<<itruth<<endl;
                    }
                }
            }
        }

        // if(lep1_matched == lep2_matched) {return;}
        // if(lep1_matched >= 0 || lep2_matched >= 0) {is_GammaConv = true;}

        if(lep1_matched >= 0 && abs(vTruth->at(0).gen_id.at(lep1_matched)) == 22) {is_GammaConv = true;}
        else  if(lep2_matched >= 0 && abs(vTruth->at(0).gen_id.at(lep2_matched)) == 22) {is_GammaConv = true;}
    }

    if(!is_GammaConv) {return;}

//--------------------------------------------
    is_2l_GammaConv_SR = 1; //2l SR category

    Compute_Variables("2l");

    fillOutputTree();

    return;
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

	bool is2lss=false, is3l=false, is4l=false;

	if(is_3l_THQ_SR || is_3l_THQ_Training || is_3l_Z_CR || is_3l_AppFakes_SR || is_3l_GammaConv_SR) {is3l = true;}
	else if(is_2l_THQ_SR || is_2l_THQ_Training || is_2l_EMU_CR || is_2l_AppFakes_SR || is_2l_GammaConv_SR || is_2l_QFlip_SR) {is2lss = true;}

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
        // Bjet1.SetPtEtaPhiE(vSelectedJets.at(ib1).pt, vSelectedJets.at(ib1).eta, vSelectedJets.at(ib1).phi, vSelectedJets.at(ib1).E);
        // FillJetInfoOutputTree(&multilepton_Bjet1_Id, 5, &multilepton_Bjet1_P4, Bjet1, &multilepton_Bjet1_CSV, vSelectedJets.at(ib1).CSVv2, &multilepton_Bjet1_JEC_Up, &multilepton_Bjet1_JEC_Down, vSelectedJets.at(ib1).JES_uncert(), &multilepton_Bjet1_JER_Up, &multilepton_Bjet1_JER_Down, vSelectedJets.at(ib1).pt_JER(), vSelectedJets.at(ib1).pt_JER_up(), vSelectedJets.at(ib1).pt_JER_down());
        multilepton_Bjet1_Id = 5;
        multilepton_Bjet1_CSV = vSelectedJets.at(ib1).CSVv2;
        multilepton_Bjet1_P4.SetPtEtaPhiE(vSelectedJets.at(ib1).pt, vSelectedJets.at(ib1).eta, vSelectedJets.at(ib1).phi, vSelectedJets.at(ib1).E);

    }
    if (ib2!=-1)
    {
        // Bjet2.SetPtEtaPhiE(vSelectedJets.at(ib2).pt, vSelectedJets.at(ib2).eta, vSelectedJets.at(ib2).phi, vSelectedJets.at(ib2).E);
        // FillJetInfoOutputTree(&multilepton_Bjet2_Id, 5, &multilepton_Bjet2_P4, Bjet2, &multilepton_Bjet2_CSV, vSelectedJets.at(ib2).CSVv2, &multilepton_Bjet2_JEC_Up, &multilepton_Bjet2_JEC_Down, vSelectedJets.at(ib2).JES_uncert(), &multilepton_Bjet2_JER_Up, &multilepton_Bjet2_JER_Down, vSelectedJets.at(ib2).pt_JER(), vSelectedJets.at(ib2).pt_JER_up(), vSelectedJets.at(ib2).pt_JER_down());
        multilepton_Bjet2_Id = 5;
        multilepton_Bjet2_CSV = vSelectedJets.at(ib2).CSVv2;
        multilepton_Bjet2_P4.SetPtEtaPhiE(vSelectedJets.at(ib2).pt, vSelectedJets.at(ib2).eta, vSelectedJets.at(ib2).phi, vSelectedJets.at(ib2).E);
    }


    //Define jet category
	//2lss
    if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2>=4) catJets = kCat_2lss_2b_4j;
    else if (is2lss && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1>=4) catJets = kCat_2lss_1b_4j;
    else if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==3) catJets = kCat_2lss_2b_3j;
    else if (is2lss && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==3) catJets = kCat_2lss_1b_3j;
	else if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==2) catJets = kCat_2lss_2b_2j;
    else if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==1) catJets = kCat_2lss_2b_1j;
    else if (is2lss && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==2) catJets = kCat_2lss_1b_2j;
	else if (is2lss && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==1) catJets = kCat_2lss_1b_1j;
    //4l
    else if (is4l && ib1!=-1 && ib2!=-1) catJets = kCat_4l_2b;
    else if (is4l && ib1!=-1 && ib2==-1) catJets = kCat_4l_1b;
    //3l
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2>=2) catJets = kCat_3l_2b_2j;
    else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1>=2) catJets = kCat_3l_1b_2j;
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==1) catJets = kCat_3l_2b_1j;
    else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==1) catJets = kCat_3l_1b_1j;
	else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==0) catJets = kCat_3l_2b_0j;
	else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==0) catJets = kCat_3l_1b_0j;
	else if (is3l && ib1==-1 && ib2==-1 && vSelectedJets.size()==1) catJets = kCat_3l_0b_1j;
	else if (is3l && ib1==-1 && ib2==-1 && vSelectedJets.size()==0) catJets = kCat_3l_0b_0j;

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
    multilepton_Lepton1_P4 = vFakeableLeptons.at(0).p4;
    multilepton_Lepton1_Id = vFakeableLeptons.at(0).id;
    multilepton_Lepton2_P4 = vFakeableLeptons.at(1).p4;
    multilepton_Lepton2_Id = vFakeableLeptons.at(1).id;

    if(is3l) //If corresponds to a 3l event
    {
        multilepton_Lepton3_P4 = vFakeableLeptons.at(2).p4;
        multilepton_Lepton3_Id = vFakeableLeptons.at(2).id;
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
                lep1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vFakeableLeptons.at(0).eta, vFakeableLeptons.at(0).phi );
                if( lep1_dr_gen < lep1_dr_gen_min)
                {
                    lep1_dr_gen_min = lep1_dr_gen;  lep1_matched = itruth;
                }

                lep2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vFakeableLeptons.at(1).eta, vFakeableLeptons.at(1).phi );
                if( lep2_dr_gen < lep2_dr_gen_min)
                {   lep2_dr_gen_min = lep2_dr_gen;  lep2_matched = itruth;  }

                if(vFakeableLeptons.size()>=3)
                {
                    lep3_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vFakeableLeptons.at(2).eta, vFakeableLeptons.at(2).phi );
                    if( lep3_dr_gen < lep3_dr_gen_min)
                    {   lep3_dr_gen_min = lep3_dr_gen;  lep3_matched = itruth;  }
                }

                jet1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vSelectedJets.at(ib1).eta, vSelectedJets.at(ib1).phi );
                if( jet1_dr_gen < jet1_dr_gen_min)
                {   jet1_dr_gen_min = jet1_dr_gen;  jet1_matched = itruth;  }


                if(ib2!=-1)
                {
                    jet2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vSelectedJets.at(ib2).eta, vSelectedJets.at(ib2).phi );
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
            if(vFakeableLeptons.size()>=3)
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
        // Jet1.SetPtEtaPhiE(vSelectedJets.at(ij1).pt, vSelectedJets.at(ij1).eta, vSelectedJets.at(ij1).phi, vSelectedJets.at(ij1).E);
        // FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vSelectedJets.at(ij1).CSVv2, &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vSelectedJets.at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vSelectedJets.at(ij1).pt_JER(), vSelectedJets.at(ij1).pt_JER_up(), vSelectedJets.at(ij1).pt_JER_down());
        multilepton_JetHighestPt1_Id = 1;
        multilepton_JetHighestPt1_CSV = vSelectedJets.at(ij1).CSVv2;
        multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt, vSelectedJets.at(ij1).eta, vSelectedJets.at(ij1).phi, vSelectedJets.at(ij1).E);
    }
    if (ij2!=-1)
    {
        // Jet2.SetPtEtaPhiE(vSelectedJets.at(ij2).pt, vSelectedJets.at(ij2).eta, vSelectedJets.at(ij2).phi, vSelectedJets.at(ij2).E);
        // FillJetInfoOutputTree(&multilepton_JetHighestPt2_Id, 1, &multilepton_JetHighestPt2_P4, Jet2, &multilepton_JetHighestPt2_CSV, vSelectedJets.at(ij2).CSVv2, &multilepton_JetHighestPt2_JEC_Up, &multilepton_JetHighestPt2_JEC_Down, vSelectedJets.at(ij2).JES_uncert(), &multilepton_JetHighestPt2_JER_Up, &multilepton_JetHighestPt2_JER_Down, vSelectedJets.at(ij2).pt_JER(), vSelectedJets.at(ij2).pt_JER_up(), vSelectedJets.at(ij2).pt_JER_down());
        multilepton_JetHighestPt2_Id = 1;
        multilepton_JetHighestPt2_CSV = vSelectedJets.at(ij2).CSVv2;
        multilepton_JetHighestPt2_P4.SetPtEtaPhiE(vSelectedJets.at(ij2).pt, vSelectedJets.at(ij2).eta, vSelectedJets.at(ij2).phi, vSelectedJets.at(ij2).E);
    }

    if (ik1!=-1 && ik2!=-1){
        // Jet1.SetPtEtaPhiE(vSelectedJets.at(ik1).pt, vSelectedJets.at(ik1).eta, vSelectedJets.at(ik1).phi, vSelectedJets.at(ik1).E);
        // Jet2.SetPtEtaPhiE(vSelectedJets.at(ik2).pt, vSelectedJets.at(ik2).eta, vSelectedJets.at(ik2).phi, vSelectedJets.at(ik2).E);
        // FillJetInfoOutputTree(&multilepton_JetClosestMw1_Id, 2, &multilepton_JetClosestMw1_P4, Jet1, &multilepton_JetClosestMw1_CSV, vSelectedJets.at(ik1).CSVv2, &multilepton_JetClosestMw1_JEC_Up, &multilepton_JetClosestMw1_JEC_Down, vSelectedJets.at(ik1).JES_uncert(), &multilepton_JetClosestMw1_JER_Up, &multilepton_JetClosestMw1_JER_Down, vSelectedJets.at(ik1).pt_JER(), vSelectedJets.at(ik1).pt_JER_up(), vSelectedJets.at(ik1).pt_JER_down());
        // FillJetInfoOutputTree(&multilepton_JetClosestMw2_Id, 2, &multilepton_JetClosestMw2_P4, Jet2, &multilepton_JetClosestMw2_CSV, vSelectedJets.at(ik2).CSVv2, &multilepton_JetClosestMw2_JEC_Up, &multilepton_JetClosestMw2_JEC_Down, vSelectedJets.at(ik2).JES_uncert(), &multilepton_JetClosestMw2_JER_Up, &multilepton_JetClosestMw2_JER_Down, vSelectedJets.at(ik2).pt_JER(), vSelectedJets.at(ik2).pt_JER_up(), vSelectedJets.at(ik2).pt_JER_down());
        multilepton_JetClosestMw1_Id = 2;
        multilepton_JetClosestMw2_Id = 2;
        multilepton_JetClosestMw1_CSV = vSelectedJets.at(ik1).CSVv2;
        multilepton_JetClosestMw2_CSV = vSelectedJets.at(ik2).CSVv2;
        multilepton_JetClosestMw1_P4.SetPtEtaPhiE(vSelectedJets.at(ik1).pt, vSelectedJets.at(ik1).eta, vSelectedJets.at(ik1).phi, vSelectedJets.at(ik1).E);
        multilepton_JetClosestMw2_P4.SetPtEtaPhiE(vSelectedJets.at(ik2).pt, vSelectedJets.at(ik2).eta, vSelectedJets.at(ik2).phi, vSelectedJets.at(ik2).E);
    }
    if (il1!=-1 && il2!=-1){
        // Jet1.SetPtEtaPhiE(vSelectedJets.at(il1).pt, vSelectedJets.at(il1).eta, vSelectedJets.at(il1).phi, vSelectedJets.at(il1).E);
        // Jet2.SetPtEtaPhiE(vSelectedJets.at(il2).pt, vSelectedJets.at(il2).eta, vSelectedJets.at(il2).phi, vSelectedJets.at(il2).E);
        // FillJetInfoOutputTree(&multilepton_JetLowestMjj1_Id, 3, &multilepton_JetLowestMjj1_P4, Jet1, &multilepton_JetLowestMjj1_CSV, vSelectedJets.at(il1).CSVv2, &multilepton_JetLowestMjj1_JEC_Up, &multilepton_JetLowestMjj1_JEC_Down, vSelectedJets.at(il1).JES_uncert(), &multilepton_JetLowestMjj1_JER_Up, &multilepton_JetLowestMjj1_JER_Down, vSelectedJets.at(il1).pt_JER(), vSelectedJets.at(il1).pt_JER_up(), vSelectedJets.at(il1).pt_JER_down());
        // FillJetInfoOutputTree(&multilepton_JetLowestMjj2_Id, 3, &multilepton_JetLowestMjj2_P4, Jet2, &multilepton_JetLowestMjj2_CSV, vSelectedJets.at(il2).CSVv2, &multilepton_JetLowestMjj2_JEC_Up, &multilepton_JetLowestMjj2_JEC_Down, vSelectedJets.at(il2).JES_uncert(), &multilepton_JetLowestMjj2_JER_Up, &multilepton_JetLowestMjj2_JER_Down, vSelectedJets.at(il2).pt_JER(), vSelectedJets.at(il2).pt_JER_up(), vSelectedJets.at(il2).pt_JER_down());
        multilepton_JetLowestMjj1_Id = 3;
        multilepton_JetLowestMjj2_Id = 3;
        multilepton_JetLowestMjj1_CSV = vSelectedJets.at(il1).CSVv2;
        multilepton_JetLowestMjj2_CSV = vSelectedJets.at(il2).CSVv2;
        multilepton_JetLowestMjj1_P4.SetPtEtaPhiE(vSelectedJets.at(il1).pt, vSelectedJets.at(il1).eta, vSelectedJets.at(il1).phi, vSelectedJets.at(il1).E);
        multilepton_JetLowestMjj2_P4.SetPtEtaPhiE(vSelectedJets.at(il2).pt, vSelectedJets.at(il2).eta, vSelectedJets.at(il2).phi, vSelectedJets.at(il2).E);
    }
    if(ie1!=-1)
    {
        multilepton_JetHighestEta1_Id = 4;
        multilepton_JetHighestEta1_CSV = vSelectedJets.at(ie1).CSVv2;
        multilepton_JetHighestEta1_P4.SetPtEtaPhiE(vSelectedJets.at(ie1).pt, vSelectedJets.at(ie1).eta, vSelectedJets.at(ie1).phi, vSelectedJets.at(ie1).E );
    }
    if(ie2!=-1) //2jets
    {
        multilepton_JetHighestEta2_Id = 4;
        multilepton_JetHighestEta2_CSV = vSelectedJets.at(ie2).CSVv2;
        multilepton_JetHighestEta2_P4.SetPtEtaPhiE(vSelectedJets.at(ie2).pt, vSelectedJets.at(ie2).eta, vSelectedJets.at(ie2).phi, vSelectedJets.at(ie2).E );
    }


    //--- Fill 2nd pairs (first one is closest to mW) -- needed for 2l only (more jets)
    if(is2lss && ij1!=-1 && ij2!=-1)
    {
        if (im1!=-1)
        {
            // Jet1.SetPtEtaPhiE(vSelectedJets.at(im1).pt, vSelectedJets.at(im1).eta, vSelectedJets.at(im1).phi, vSelectedJets.at(im1).E);
            // FillJetInfoOutputTree(&multilepton_JetHighestPt1_2ndPair_Id, 1, &multilepton_JetHighestPt1_2ndPair_P4, Jet1, &multilepton_JetHighestPt1_2ndPair_CSV, vSelectedJets.at(im1).CSVv2, &multilepton_JetHighestPt1_2ndPair_JEC_Up, &multilepton_JetHighestPt1_2ndPair_JEC_Down, vSelectedJets.at(im1).JES_uncert(), &multilepton_JetHighestPt1_2ndPair_JER_Up, &multilepton_JetHighestPt1_2ndPair_JER_Down, vSelectedJets.at(im1).pt_JER(), vSelectedJets.at(im1).pt_JER_up(), vSelectedJets.at(im1).pt_JER_down());
            multilepton_JetHighestPt1_2ndPair_Id = 1;
            multilepton_JetHighestPt1_2ndPair_CSV = vSelectedJets.at(im1).CSVv2;
            multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im1).pt, vSelectedJets.at(im1).eta, vSelectedJets.at(im1).phi, vSelectedJets.at(im1).E);
        }
        if(im2!=-1){
            // Jet2.SetPtEtaPhiE(vSelectedJets.at(im2).pt, vSelectedJets.at(im2).eta, vSelectedJets.at(im2).phi, vSelectedJets.at(im2).E);
            // FillJetInfoOutputTree(&multilepton_JetHighestPt2_2ndPair_Id, 1, &multilepton_JetHighestPt2_2ndPair_P4, Jet2, &multilepton_JetHighestPt2_2ndPair_CSV, vSelectedJets.at(im2).CSVv2, &multilepton_JetHighestPt2_2ndPair_JEC_Up, &multilepton_JetHighestPt2_2ndPair_JEC_Down, vSelectedJets.at(im2).JES_uncert(), &multilepton_JetHighestPt2_2ndPair_JER_Up, &multilepton_JetHighestPt2_2ndPair_JER_Down, vSelectedJets.at(im2).pt_JER(), vSelectedJets.at(im2).pt_JER_up(), vSelectedJets.at(im2).pt_JER_down());
            multilepton_JetHighestPt2_2ndPair_Id = 1;
            multilepton_JetHighestPt1_2ndPair_CSV = vSelectedJets.at(im2).CSVv2;
            multilepton_JetHighestPt2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im2).pt, vSelectedJets.at(im2).eta, vSelectedJets.at(im2).phi, vSelectedJets.at(im2).E);
        }
        if (io1!=-1 && io2!=-1){
            // Jet1.SetPtEtaPhiE(vSelectedJets.at(ip1).pt, vSelectedJets.at(ip1).eta, vSelectedJets.at(ip1).phi, vSelectedJets.at(ip1).E);
            // Jet2.SetPtEtaPhiE(vSelectedJets.at(io2).pt, vSelectedJets.at(io2).eta, vSelectedJets.at(io2).phi, vSelectedJets.at(io2).E);
            // FillJetInfoOutputTree(&multilepton_JetClosestMw1_2ndPair_Id, 2, &multilepton_JetClosestMw1_2ndPair_P4, Jet1, &multilepton_JetClosestMw1_2ndPair_CSV, vSelectedJets.at(io1).CSVv2, &multilepton_JetClosestMw1_2ndPair_JEC_Up, &multilepton_JetClosestMw1_2ndPair_JEC_Down, vSelectedJets.at(io1).JES_uncert(), &multilepton_JetClosestMw1_2ndPair_JER_Up, &multilepton_JetClosestMw1_2ndPair_JER_Down, vSelectedJets.at(io1).pt_JER(), vSelectedJets.at(io1).pt_JER_up(), vSelectedJets.at(io1).pt_JER_down());
            // FillJetInfoOutputTree(&multilepton_JetClosestMw2_2ndPair_Id, 2, &multilepton_JetClosestMw2_2ndPair_P4, Jet2, &multilepton_JetClosestMw2_2ndPair_CSV, vSelectedJets.at(io2).CSVv2, &multilepton_JetClosestMw2_2ndPair_JEC_Up, &multilepton_JetClosestMw2_2ndPair_JEC_Down, vSelectedJets.at(io2).JES_uncert(), &multilepton_JetClosestMw2_2ndPair_JER_Up, &multilepton_JetClosestMw2_2ndPair_JER_Down, vSelectedJets.at(io2).pt_JER(), vSelectedJets.at(io2).pt_JER_up(), vSelectedJets.at(io2).pt_JER_down());
            multilepton_JetClosestMw1_2ndPair_Id = 2;
            multilepton_JetClosestMw2_2ndPair_Id = 2;
            multilepton_JetClosestMw1_2ndPair_CSV = vSelectedJets.at(io1).CSVv2;
            multilepton_JetClosestMw2_2ndPair_CSV = vSelectedJets.at(io2).CSVv2;
            multilepton_JetClosestMw1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(io1).pt, vSelectedJets.at(io1).eta, vSelectedJets.at(io1).phi, vSelectedJets.at(io1).E);
            multilepton_JetClosestMw2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(io2).pt, vSelectedJets.at(io2).eta, vSelectedJets.at(io2).phi, vSelectedJets.at(io2).E);
        }
        if (ip1!=-1 && ip2!=-1){
            // Jet1.SetPtEtaPhiE(vSelectedJets.at(ip1).pt, vSelectedJets.at(ip1).eta, vSelectedJets.at(ip1).phi, vSelectedJets.at(ip1).E);
            // Jet2.SetPtEtaPhiE(vSelectedJets.at(ip2).pt, vSelectedJets.at(ip2).eta, vSelectedJets.at(ip2).phi, vSelectedJets.at(ip2).E);
            // FillJetInfoOutputTree(&multilepton_JetLowestMjj1_2ndPair_Id, 3, &multilepton_JetLowestMjj1_2ndPair_P4, Jet1, &multilepton_JetLowestMjj1_2ndPair_CSV, vSelectedJets.at(ip1).CSVv2, &multilepton_JetLowestMjj1_2ndPair_JEC_Up, &multilepton_JetLowestMjj1_2ndPair_JEC_Down, vSelectedJets.at(ip1).JES_uncert(), &multilepton_JetLowestMjj1_2ndPair_JER_Up, &multilepton_JetLowestMjj1_2ndPair_JER_Down, vSelectedJets.at(ip1).pt_JER(), vSelectedJets.at(ip1).pt_JER_up(), vSelectedJets.at(ip1).pt_JER_down());
            // FillJetInfoOutputTree(&multilepton_JetLowestMjj2_2ndPair_Id, 3, &multilepton_JetLowestMjj2_2ndPair_P4, Jet2, &multilepton_JetLowestMjj2_2ndPair_CSV, vSelectedJets.at(ip2).CSVv2, &multilepton_JetLowestMjj2_2ndPair_JEC_Up, &multilepton_JetLowestMjj2_2ndPair_JEC_Down, vSelectedJets.at(ip2).JES_uncert(), &multilepton_JetLowestMjj2_2ndPair_JER_Up, &multilepton_JetLowestMjj2_2ndPair_JER_Down, vSelectedJets.at(ip2).pt_JER(), vSelectedJets.at(ip2).pt_JER_up(), vSelectedJets.at(ip2).pt_JER_down());
            multilepton_JetLowestMjj1_2ndPair_Id = 3;
            multilepton_JetLowestMjj2_2ndPair_Id = 3;
            multilepton_JetLowestMjj1_2ndPair_CSV = vSelectedJets.at(ip1).CSVv2;
            multilepton_JetLowestMjj2_2ndPair_CSV = vSelectedJets.at(ip2).CSVv2;
            multilepton_JetLowestMjj1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(ip1).pt, vSelectedJets.at(ip1).eta, vSelectedJets.at(ip1).phi, vSelectedJets.at(ip1).E);
            multilepton_JetLowestMjj2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(ip2).pt, vSelectedJets.at(ip2).eta, vSelectedJets.at(ip2).phi, vSelectedJets.at(ip2).E);
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
    multilepton_mHT = vEvent->at(0).metsumet;

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
 * @param region (3l or 2l)
 */
void tHqMultileponAnalysis::Compute_Variables(TString region)
{
	if(region != "3l" && region != "2l")
	{
		cout<<FRED("Error ! Wrong region name at Compute_Variables() call !")<<endl;
	}

    int nleptons = 0;
    if(region == "3l") nleptons = 3;
    else nleptons = 2;

    //--- Determine leptonic channel of event
	int sum_id = 0;
	if(region == "3l")
    {
        sum_id = fabs(vFakeableLeptons.at(0).id)+fabs(vFakeableLeptons.at(1).id)+fabs(vFakeableLeptons.at(2).id);
        if(sum_id == 39) {channel = 0;} //uuu
        else if(sum_id == 37) {channel = 1;} //uue
        else if(sum_id == 35) {channel = 2;} //eeu
        else if(sum_id == 33) {channel = 3;} //eee
    }
    else
    {
        sum_id = fabs(vFakeableLeptons.at(0).id)+fabs(vFakeableLeptons.at(1).id);
        if(sum_id == 26) {channel = 0;} //uu
        else if(sum_id == 24) {channel = 1;} //ue + eu
        else if(sum_id == 22) {channel = 2;} //ee -- not used in tHq2016
    }

    // ##########
    // # MET LD #
    // ##########
    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;
    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt, vSelectedJets.at(i).eta, vSelectedJets.at(i).phi, vSelectedJets.at(i).E);
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4.Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4.Py();
    }

    TLorentzVector taup4;
    for(int i=0; i<vSelectedTaus.size(); i++)
    {
        taup4.SetPtEtaPhiE(vSelectedTaus.at(i).pt, vSelectedTaus.at(i).eta, vSelectedTaus.at(i).phi, vSelectedTaus.at(i).E);
        tau_px = tau_px + taup4.Px();
        tau_py = tau_py + taup4.Py();
    }

    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt + 0.00265 * MHT;

    if(DEBUG) std::cout << " MHT =  " << MHT << "MET = " << vEvent->at(0).metpt << " met_ld = " << met_ld ;

    // #################################
    // # b-tagging nominal reweighting #
    // #################################
    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt                );
        jetEtas.push_back(    vSelectedJets.at(i).eta               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2             );
        // jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    // wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    // new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################
    // std::vector<double> weights_csv;
    // double wgt_csv_def_sys = 0;

    // for(int i=7; i<25; i++)
    // {
    //     wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
    //     weights_csv.push_back(wgt_csv_def_sys);
    // }

    // double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end());
    // double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end());
    // weight_csv_down = min_weight_csv;
    // weight_csv_up   = max_weight_csv;

// ##################################################################################################################################


    //--- Compute input variables (tHq 2016)
    //NB : "forward jet" = most forward non-CSV loose jet

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


    // double nJet25;
    // double maxEtaJet25;
    // double lepCharge;
    // double nJetEta1 ;
    // double dEtaFwdJetBJet;
    // double dEtaFwdJet2BJet;
    // double dEtaFwdJetClosestLep;
    // double dPhiHighestPtSSPair;
    // double minDRll;
    // double Lep3Pt;

	//--- Var1 : nof jets with pT>25 and |eta|<2.4
    nJet25 = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if(vSelectedJets.at(ijet).pt > 25 && fabs(vSelectedJets.at(ijet).eta ) < 2.4) {nJet25++;}
    }

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
        lepCharge+= vFakeableLeptons.at(ilep).charge;
    }

	//--- Var4 : nof 'non-csv-loose' jets with eta>1.0
    nJetEta1 = 0;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if( fabs(vLightJets.at(ijet).eta ) > 1.0) {nJetEta1++;}
    }

	//--- Var5 : dEta between forward light jet and hardest tagged jet
    dEtaFwdJetBJet = fabs( vLightJets.at(ijet_forward).eta - vLooseBTagJets.at(ijet_hardest_btag).eta );

	//--- Var6 : dEta between forward and 2nd hardest tagged jet
    if(ijet_2nd_hardest_btag < 0) {dEtaFwdJet2BJet = -1;}
    else {dEtaFwdJet2BJet = fabs( vLightJets.at(ijet_forward).eta - vLooseBTagJets.at(ijet_2nd_hardest_btag).eta );}

	//--- Var7 : dEta between forward light jet and closet lepton (angular dist.)
    dEtaFwdJetClosestLep = 0; tmp = 999;
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta - vFakeableLeptons.at(ilep).eta ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta - vFakeableLeptons.at(ilep).eta ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt, vFakeableLeptons.at(i).eta, vFakeableLeptons.at(i).phi, vFakeableLeptons.at(i).E );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt, vFakeableLeptons.at(j).eta, vFakeableLeptons.at(j).phi, vFakeableLeptons.at(j).E );
            if(vFakeableLeptons.at(i).charge==vFakeableLeptons.at(j).charge && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi - vFakeableLeptons.at(j).phi ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt, vFakeableLeptons.at(i).eta, vFakeableLeptons.at(i).phi, vFakeableLeptons.at(i).E );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt, vFakeableLeptons.at(j).eta, vFakeableLeptons.at(j).phi, vFakeableLeptons.at(j).E );
            if(lepi.DeltaR(lepj) < minDRll)
            {

                minDRll = lepi.DeltaR(lepj);
            }
        }
    }

    //--- Var10 : pT of 3rd hardest lepton
    if(region == "3l") {Lep3Pt = vFakeableLeptons.at(2).pt;}
    else {Lep3Pt = vFakeableLeptons.at(1).pt;}



    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt; lep2Pt = vFakeableLeptons.at(1).pt; lep3Pt = vFakeableLeptons.at(1).pt;
    if(region == "3l") {lep3Pt = vFakeableLeptons.at(2).pt;}
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt; hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta;
    fwdJetPt = vLightJets.at(ijet_forward).pt;

    if(vFakeableLeptons.size() > 2)
    {
        float mll_tmp = 0; inv_mll = 999;
        for(int i=0; i<vFakeableLeptons.size()-1; i++)
        {
            for(int j=i+1; j<vFakeableLeptons.size(); j++)
            {
                if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && vFakeableLeptons.at(i).charge == -vFakeableLeptons.at(j).charge )
                {
                    mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                    if(fabs(mll_tmp-91.188) < fabs(inv_mll-91.188)) {inv_mll = mll_tmp;}
                }
            }
        }
    }


    //--------------------------
    // PRINTOUT OF INFOS & INPUT VARS
    //--------------------------

    bool do_printout = false;

    if(do_printout)
    {
        cout<<endl<<endl<<BOLD(FBLU("------------ EVENT -----------"))<<endl;

        cout<<FYEL("--- Tagged Jets : ")<<endl;
        for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt<<" / eta = "<<vLooseBTagJets.at(ijet).eta<<" / phi = "<<vLooseBTagJets.at(ijet).phi<<" / CSV = "<<vLooseBTagJets.at(ijet).CSVv2<<endl;
        }
        cout<<"Hardest & 2nd hardest jets are "<<ijet_hardest_btag<<", "<<ijet_2nd_hardest_btag<<endl<<endl;

        cout<<FYEL("--- Forward Jets : ")<<endl;
        for(int ijet=0; ijet<vLightJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLightJets.at(ijet).pt<<" / eta = "<<vLightJets.at(ijet).eta<<" / phi = "<<vLightJets.at(ijet).phi<<" / CSV = "<<vLightJets.at(ijet).CSVv2<<endl;
        }
        cout<<"Forwardest jet is "<<ijet_forward<<endl<<endl;

        cout<<FYEL("--- Selected Leptons : ")<<endl;
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt<<" / eta = "<<vFakeableLeptons.at(ilep).eta<<" / phi = "<<vFakeableLeptons.at(ilep).phi<<" : Charge = "<<vFakeableLeptons.at(ilep).charge<<endl;
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


    //--- COMPUTE BDT OUTPUT from weight files
    if(region == "3l")
    {
        signal_3l_TT_MVA   = mva_3l_tt->EvaluateMVA("BDTG method");
        signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");
    }
    else
    {
        // signal_2lss_TT_MVA   = mva_2lss_tt->EvaluateMVA("BDTG method");
        // signal_2lss_TTV_MVA   = mva_2lss_ttV->EvaluateMVA("BDTG method");
    }
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

    // tOutput->Branch("is_debug",&is_debug,"is_debug/B");


	//-- Event main infos
	tOutput->Branch("channel",&channel,"channel/F");
    tOutput->Branch("weight",&weight,"weight/F");
	tOutput->Branch("weightfake",&weightfake,"weightfake/F");
	tOutput->Branch("weightflip",&weightflip,"weightflip/F");
    tOutput->Branch("event_id",&event_id,"event_id/F");
    tOutput->Branch("event_run",&event_run,"event_run/F");
	tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");
    // tOutput->Branch("is_trigger",&is_trigger,"is_trigger/O");


	//--- Categories & MVA
	tOutput->Branch("is_3l_THQ_SR",&is_3l_THQ_SR,"is_3l_THQ_SR/B");
    tOutput->Branch("is_2l_THQ_SR",&is_2l_THQ_SR,"is_2l_THQ_SR/B");
	tOutput->Branch("is_3l_THQ_Training",&is_3l_THQ_Training,"is_3l_THQ_Training/B");
    tOutput->Branch("is_2l_THQ_Training",&is_2l_THQ_Training,"is_2l_THQ_Training/B");
    tOutput->Branch("is_3l_Z_CR",&is_3l_Z_CR,"is_3l_Z_CR/B");
    tOutput->Branch("is_2l_EMU_CR",&is_2l_EMU_CR,"is_2l_EMU_CR/B");
    tOutput->Branch("is_3l_AppFakes_SR",&is_3l_AppFakes_SR,"is_3l_AppFakes_SR/B");
    tOutput->Branch("is_2l_AppFakes_SR",&is_2l_AppFakes_SR,"is_2l_AppFakes_SR/B");
    tOutput->Branch("is_2l_QFlip_SR",&is_2l_QFlip_SR,"is_2l_QFlip_SR/B");
    tOutput->Branch("is_3l_GammaConv_SR",&is_3l_GammaConv_SR,"is_3l_GammaConv_SR/B");
    tOutput->Branch("is_2l_GammaConv_SR",&is_2l_GammaConv_SR,"is_2l_GammaConv_SR/B");

	// tOutput->Branch("signal_3l_TT_MVA",&signal_3l_TT_MVA,"signal_3l_TT_MVA/F");
	// tOutput->Branch("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,"signal_3l_TTV_MVA/F");
    // tOutput->Branch("signal_2lss_TT_MVA",&signal_2lss_TT_MVA,"signal_2lss_TT_MVA/F");
    // tOutput->Branch("signal_2lss_TTV_MVA",&signal_2lss_TTV_MVA,"signal_2lss_TTV_MVA/F");

    //-- Input variables from tHq2016 analysis
    tOutput->Branch("nJet25",&nJet25,"nJet25/F");
    tOutput->Branch("maxEtaJet25",&maxEtaJet25,"maxEtaJet25/F");
    tOutput->Branch("lepCharge",&lepCharge,"lepCharge/F");
    tOutput->Branch("nJetEta1",&nJetEta1,"nJetEta1/F");
    tOutput->Branch("dEtaFwdJetBJet",&dEtaFwdJetBJet,"dEtaFwdJetBJet/F");
    tOutput->Branch("dEtaFwdJet2BJet",&dEtaFwdJet2BJet,"dEtaFwdJet2BJet/F");
    tOutput->Branch("dEtaFwdJetClosestLep",&dEtaFwdJetClosestLep,"dEtaFwdJetClosestLep/F");
    tOutput->Branch("dPhiHighestPtSSPair",&dPhiHighestPtSSPair,"dPhiHighestPtSSPair/F");
    tOutput->Branch("minDRll",&minDRll,"minDRll/F");
    tOutput->Branch("Lep3Pt",&Lep3Pt,"Lep3Pt/F");

	//-- More control vars
    tOutput->Branch("inv_mll",&inv_mll,"inv_mll/F");
    tOutput->Branch("hardestBjetPt",&hardestBjetPt,"hardestBjetPt/F");
    tOutput->Branch("hardestBjetEta",&hardestBjetEta,"hardestBjetEta/F");
    tOutput->Branch("fwdJetPt",&fwdJetPt,"fwdJetPt/F");
    tOutput->Branch("lep1Pt",&lep1Pt,"lep1Pt/F");
    tOutput->Branch("lep2Pt",&lep2Pt,"lep2Pt/F");
    tOutput->Branch("lep3Pt",&lep3Pt,"lep3Pt/F");
    tOutput->Branch("MET",&MET,"MET/F");

	//-- Cut variables
	tOutput->Branch("nLooseBJets",&nLooseBJets,"nLooseBJets/F");
	tOutput->Branch("nMediumBJets",&nMediumBJets,"nMediumBJets/F");
	tOutput->Branch("nLightJets",&nLightJets,"nLightJets/F");
	tOutput->Branch("nTightLep",&nTightLep,"nTightLep/F");
	tOutput->Branch("nFakeableLep",&nFakeableLep,"nFakeableLep/F");
	tOutput->Branch("nLightJets_Fwd40",&nLightJets_Fwd40,"nLightJets_Fwd40/F");


	//-- Other weights
    // tOutput->Branch("weight_scale_muF0p5",&weight_scale_muF0p5,"weight_scale_muF0p5/F");
    // tOutput->Branch("weight_scale_muF2",&weight_scale_muF2,"weight_scale_muF2/F");
    // tOutput->Branch("weight_scale_muR0p5",&weight_scale_muR0p5,"weight_scale_muR0p5/F");
    // tOutput->Branch("weight_scale_muR2",&weight_scale_muR2,"weight_scale_muR2/F");
    // tOutput->Branch("weight_csv_down",&weight_csv_down,"weight_csv_down/F");
    // tOutput->Branch("weight_csv_up",&weight_csv_up,"weight_csv_up/F");
    // tOutput->Branch("weights_pdf","std::vector<float>",&weights_pdf);
    // tOutput->Branch("ids_pdf","std::vector<std::string>",&ids_pdf);

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
 * @param vSelectedJets     Contains all selected jets (b & light )
 * @param ibsel1            Return bjet index 1
 * @param ibsel2            Return bjet index 2
 * @param doSelectOnlyBjets True if consider only bjets in vector
 */
void tHqMultileponAnalysis::SelectBjets(int &ibsel1, int &ibsel2, bool doSelectOnlyBjets=true)
{
    int ib1=-1, ib2=-1;

    double CSV_threshold = 0.5426; //Loose WP

    Float_t btag_max=-9999, btag_max2=-9999;
    for (int ib=0; ib<vSelectedJets.size(); ib++)
    {
        if (doSelectOnlyBjets && (vSelectedJets.at(ib).CSVv2 < CSV_threshold || fabs(vSelectedJets.at(ib).eta ) > 2.4) ) {continue;}

        if (vSelectedJets.at(ib).CSVv2>btag_max)
        {
          btag_max2 = btag_max;
          ib2 = ib1;
          btag_max = vSelectedJets.at(ib).CSVv2;
          ib1 = ib;
        }
        else if (vSelectedJets.at(ib).CSVv2<btag_max && vSelectedJets.at(ib).CSVv2>btag_max2)
        {
          btag_max2 = vSelectedJets.at(ib).CSVv2;
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
 * @param vSelectedJets [Contains all jets (btagged or not) passing the object selection]
 */
void tHqMultileponAnalysis::SelectOtherJets(const int ib1, const int ib2, int &ij1, int &ij2, int &ik1, int &ik2, int &ie1, int &ie2, int &il1, int &il2, int &im1, int &im2, int &io1, int &io2, int &ip1, int &ip2)
{
    TLorentzVector Pjet1, Pjet2;
    float pt_max=0, pt_max2=0;
    float diffmass_min = 10000, mass_min = 10000;
    float eta_max=0, eta_max2=0;

    for(int ij=0; ij<vSelectedJets.size(); ij++)
    {
        if (ij==ib1 || ij==ib2) {continue;} //Don't take bjets into account

        if(vSelectedJets.at(ij).pt > pt_max ) //Highest pT
        {
            pt_max2 = pt_max;
            ij2 = ij1;
            pt_max = vSelectedJets.at(ij).pt;
            ij1 = ij;
        }
        else if(vSelectedJets.at(ij).pt > pt_max2) //2nd Highest pT
        {
            pt_max2 = vSelectedJets.at(ij).pt;
            ij2 = ij;
        }

        if(fabs(vSelectedJets.at(ij).eta) > eta_max ) //Highest eta
        {
            eta_max2 = eta_max;
            ie2 = ie1;
            eta_max = fabs(vSelectedJets.at(ij).eta);
            ie1 = ij;
        }
        else if(fabs(vSelectedJets.at(ij).eta ) > eta_max2) //2nd Highest eta
        {
            eta_max2 = fabs(vSelectedJets.at(ij).eta);
            ie2 = ij;
        }

        for(int ik=ij+1; ik<vSelectedJets.size(); ik++) //Dijet w/ lowest (m - mW)
        {
            if (ik==ib1 || ik==ib2) {continue;} //Don't take bjets into account

            Pjet1.SetPtEtaPhiE(vSelectedJets.at(ij).pt, vSelectedJets.at(ij).eta, vSelectedJets.at(ij).phi, vSelectedJets.at(ij).E);
            Pjet2.SetPtEtaPhiE(vSelectedJets.at(ik).pt, vSelectedJets.at(ik).eta, vSelectedJets.at(ik).phi, vSelectedJets.at(ik).E);

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
    for(int ij=0; ij<vSelectedJets.size(); ij++)
    {
        if (ij==ib1 || ij==ib2 || ij==ik1 || ij==ik2) {continue;} //Don't take bjets and 1rst mW pair into account

        if(vSelectedJets.at(ij).pt > pt_max ) //Highest pT
        {
            pt_max2 = pt_max;
            im2 = im1;
            pt_max = vSelectedJets.at(ij).pt;
            im1 = ij;
        }
        else if(vSelectedJets.at(ij).pt > pt_max2) //2nd Highest pT
        {
            pt_max2 = vSelectedJets.at(ij).pt;
            im2 = ij;
        }

        for(int ik=ij+1; ik<vSelectedJets.size(); ik++) //Dijet w/ lowest (m - mW)
        {
            if (ik==ib1 || ik==ib2 || ik==ik1 || ik==ik2) {continue;} //Don't take bjets and 1rst pair into account

            Pjet1.SetPtEtaPhiE(vSelectedJets.at(ij).pt, vSelectedJets.at(ij).eta, vSelectedJets.at(ij).phi, vSelectedJets.at(ij).E);
            Pjet2.SetPtEtaPhiE(vSelectedJets.at(ik).pt, vSelectedJets.at(ik).eta, vSelectedJets.at(ik).phi, vSelectedJets.at(ik).E);

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
    if(_sampleName.Contains("ttZJets_13TeV_madgraphMLM") || _sampleName.Contains("ttWJets_13TeV_madgraphMLM") || _sampleName.Contains("THQ") || _sampleName.Contains("TTJets") )
    {
        return true;
    }

    return false;
}

//Decides if given sample is to be used for build gamma-conversion sample or not
bool tHqMultileponAnalysis::Sample_isUsed_forGammaConv()
{
    if(_sampleName.Contains("TTGJets") || _sampleName.Contains("WGToLNuG") || _sampleName.Contains("ZGTo2LG") || _sampleName.Contains("TGJets") )
    {
        return true;
    }

    return false;
}

/*
//Convert a double into a TString
// precision --> can choose if TString how many digits the TString should display
TString tHqMultileponAnalysis::Convert_Number_To_TString(double number, int precision=10)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}

//Convert a TString into a double
double tHqMultileponAnalysis::Convert_TString_To_Number(TString ts)
{
	double number = 0;
	string s = ts.Data();
	stringstream ss(s);
	ss >> number;
	return number;
}


//Use stat function (from library sys/stat) to check if a file exists
bool tHqMultileponAnalysis::Check_File_Existence(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
}
*/



























//--------------------------------------------
// ########  ######## ########  ##     ##  ######    ######   #### ##    ##  ######
// ##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ###   ## ##    ##
// ##     ## ##       ##     ## ##     ## ##        ##         ##  ####  ## ##
// ##     ## ######   ########  ##     ## ##   #### ##   ####  ##  ## ## ## ##   ####
// ##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ##  #### ##    ##
// ##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ##   ### ##    ##
// ########  ######## ########   #######   ######    ######   #### ##    ##  ######
//--------------------------------------------

/*

//DoubleEG 2016 event IDs
list_ids.push_back(49893899);
list_ids.push_back(533217326);
list_ids.push_back(468000459);
list_ids.push_back(249757834);
list_ids.push_back(1494760657);
list_ids.push_back(1661828372);
list_ids.push_back(539950871);
list_ids.push_back(115325936);
list_ids.push_back(2297159671);
list_ids.push_back(3042404601);
list_ids.push_back(427918880);
list_ids.push_back(2916149572);
list_ids.push_back(325972278);
list_ids.push_back(3995027221);
list_ids.push_back(416533726);
list_ids.push_back(344468491);
list_ids.push_back(1292476406);
list_ids.push_back(677435209);
list_ids.push_back(1773660800);
list_ids.push_back(1166598525);
list_ids.push_back(1779113700);
list_ids.push_back(360784971);
list_ids.push_back(1421624580);
list_ids.push_back(2453510338);
list_ids.push_back(1359708490);
list_ids.push_back(984984150);
list_ids.push_back(249008513);
list_ids.push_back(2338434137);
list_ids.push_back(2917515592);
list_ids.push_back(1153388564);
list_ids.push_back(2811039954);
list_ids.push_back(1181123914);
list_ids.push_back(67525210);
list_ids.push_back(11502597);
list_ids.push_back(73052667);
list_ids.push_back(67074121);
list_ids.push_back(292351561);
list_ids.push_back(2573612145);
list_ids.push_back(1424492475);
list_ids.push_back(987507944);
list_ids.push_back(1161863314);
list_ids.push_back(759720511);
list_ids.push_back(235298156);
list_ids.push_back(1033125580);
list_ids.push_back(2190755104);

//DoubleMuon 2016 event IDs
list_ids.push_back(693957090);
list_ids.push_back(605768208);
list_ids.push_back(76478780);
list_ids.push_back(353139236);
list_ids.push_back(613403042);
list_ids.push_back(820990721);
list_ids.push_back(3373504328);
list_ids.push_back(1160453054);
list_ids.push_back(1304884994);
list_ids.push_back(1859651621);
list_ids.push_back(612753841);
list_ids.push_back(856788860);
list_ids.push_back(795980629);
list_ids.push_back(1245388987);
list_ids.push_back(2080388751);
list_ids.push_back(1279943382);
list_ids.push_back(249576069);
list_ids.push_back(2580593845);
list_ids.push_back(782283930);
list_ids.push_back(1260747954);
list_ids.push_back(692735720);
list_ids.push_back(4795099248);
list_ids.push_back(1798328655);
list_ids.push_back(1621039363);
list_ids.push_back(77392441);
list_ids.push_back(1644425191);
list_ids.push_back(2011322588);
list_ids.push_back(1587940463);
list_ids.push_back(2918482755);
list_ids.push_back(372686400);
list_ids.push_back(1159352945);
list_ids.push_back(3209955538);
list_ids.push_back(816786352);
list_ids.push_back(213517976);
list_ids.push_back(918813656);
list_ids.push_back(356940857);
list_ids.push_back(232720721);
list_ids.push_back(2397102244);
list_ids.push_back(28566123);
list_ids.push_back(947375279);
list_ids.push_back(2120122022);
list_ids.push_back(2365406318);
list_ids.push_back(730968036);
list_ids.push_back(890523132);
list_ids.push_back(377641780);
list_ids.push_back(2067141857);
list_ids.push_back(3584452270);
list_ids.push_back(1295215576);
list_ids.push_back(24693238);
list_ids.push_back(152465675);
list_ids.push_back(260759050);
list_ids.push_back(324510783);
list_ids.push_back(1276360547);
list_ids.push_back(1629460524);
list_ids.push_back(1071622490);
list_ids.push_back(441283464);
list_ids.push_back(735122949);
list_ids.push_back(1066331089);
list_ids.push_back(328621184);
list_ids.push_back(2632762769);
list_ids.push_back(540024261);
list_ids.push_back(2136344002);
list_ids.push_back(1832390840);
list_ids.push_back(3587377843);
list_ids.push_back(266935169);
list_ids.push_back(3155971678);
list_ids.push_back(3014727543);
list_ids.push_back(1222324278);
list_ids.push_back(607443847);
list_ids.push_back(2217217561);
list_ids.push_back(714487283);
list_ids.push_back(557832244);
list_ids.push_back(508663945);
list_ids.push_back(2685100401);
list_ids.push_back(2461076082);
list_ids.push_back(99870196);



void tHqMultileponAnalysis::Debug_Selection(int evt)
{
    InitVariables();

    is_debug = true;

    //--- Var1 : nof jets with pT>25 and |eta|<2.4
    nJet25 = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if(vSelectedJets.at(ijet).pt > 25 && fabs(vSelectedJets.at(ijet).eta ) < 2.4) {nJet25++;}
    }


    //--------------------------
    // PRINTOUT OF INFOS & INPUT VARS
    //--------------------------

    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4 + vLooseLeptons.at(j).p4).M() ) < 12 ) {pass_cleanup = false;}
        }
    }

    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id) == fabs(vFakeableLeptons.at(j).id) && vFakeableLeptons.at(i).charge == -vFakeableLeptons.at(j).charge )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4 + vFakeableLeptons.at(j).p4 ).M();
                if( fabs(mll_tmp - 91.188) < 15) {pass_Zveto = false ;}
            }
        }
    }


    int nJet25_tmp = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if(vSelectedJets.at(ijet).pt > 25 && fabs(vSelectedJets.at(ijet).eta) < 2.4 ) nJet25_tmp++;
    }


    bool do_printout = true;
    if(do_printout)
    {
        cout<<endl<<endl<<endl<<endl;

        cout<<endl<<endl<<BOLD(FBLU("------------ EVENT "<<setprecision(15)<<event_id<<" (run "<<event_run<<") -----------"))<<endl;

        cout<<"Pass cleaning = "<<pass_cleanup<<" / Pass Zveto = "<<pass_Zveto<<endl<<endl;
        cout<<"isTrigger = "<<is_trigger<<" / weight = "<<weight<<endl;


		cout<<FYEL("--- Triggers : ")<<endl;
        cout<<"vEvent->at(0).trig_ee  "<<vEvent->at(0).trig_ee <<endl;
        cout<<"vEvent->at(0).trig_eee "<<vEvent->at(0).trig_eee<<endl;
        cout<<"vEvent->at(0).trig_mm  "<<vEvent->at(0).trig_mm <<endl;
        cout<<"vEvent->at(0).is_TRIGmmTk() "<<vEvent->at(0).is_TRIGmmTk()<<endl;
        cout<<"vEvent->at(0).trig_mmm() "<<vEvent->at(0).trig_mmm()<<endl<<endl;
		cout<<"vEvent->at(0).is_TRIGmTk() "<<vEvent->at(0).is_TRIGmTk()<<endl;
		cout<<"vEvent->at(0).trig_em  "<<vEvent->at(0).trig_em <<endl;
		cout<<"vEvent->at(0).is_TRIGem()  "<<vEvent->at(0).is_TRIGem() <<endl;
		cout<<"vEvent->at(0).trig_emm "<<vEvent->at(0).trig_emm<<endl;
		cout<<"vEvent->at(0).trig_eem "<<vEvent->at(0).trig_eem<<endl;
        cout<<"vEvent->at(0).trig_e   "<<vEvent->at(0).trig_e  <<endl;
        cout<<"vEvent->at(0).trig_m   "<<vEvent->at(0).trig_m  <<endl;


        cout<<FYEL("--- Loose Tagged Jets : ")<<endl;
        for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt<<" / eta = "<<vLooseBTagJets.at(ijet).eta<<" / phi = "<<vLooseBTagJets.at(ijet).phi<<" / CSV = "<<vLooseBTagJets.at(ijet).CSVv2<<endl;
        }

        cout<<FYEL("--- Light Jets : ")<<endl;
        for(int ijet=0; ijet<vLightJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLightJets.at(ijet).pt<<" / eta = "<<vLightJets.at(ijet).eta<<" / phi = "<<vLightJets.at(ijet).phi<<" / CSV = "<<vLightJets.at(ijet).CSVv2<<endl;
        }

        cout<<FYEL("--- FO Leptons : ")<<endl;
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            TString type = "";
            if(vFakeableLeptons[ilep].isElectron()) type = "Ele";
            else if(vFakeableLeptons[ilep].isMuon()) type = "Mu";

            if(ilep==0) lep1Pt = vFakeableLeptons.at(ilep).pt;
            else if(ilep==1) lep2Pt = vFakeableLeptons.at(ilep).pt;
            else if(ilep==2) lep3Pt = vFakeableLeptons.at(ilep).pt;

            cout<<ilep<<" "<<type<<": pT = "<<vFakeableLeptons.at(ilep).pt<<" / eta = "<<vFakeableLeptons.at(ilep).eta<<" / phi = "<<vFakeableLeptons.at(ilep).phi<<" : Charge = "<<vFakeableLeptons.at(ilep).charge<<" / ";
            cout<<"isTightTTH = "<<vFakeableLeptons.at(ilep).isTightTTH;
            // if(vFakeableLeptons[ilep].isElectron()) {cout<<" / CutEvtSel = "<<vFakeableLeptons.at(ilep).cutEventSel()<<" / noLostHit = "<<vFakeableLeptons.at(ilep).noLostHit();}
            cout<<endl<<endl;
        }

        cout<<endl<<"==> Njet25 = "<<nJet25_tmp<<" / nMediumBJets = "<<nMediumBJets<<endl;

        cout<<"------------------"<<endl<<endl;

        cout<<"Printing all leptons :"<<endl;
        // for(int ilep=0; ilep<vLooseLeptons.size(); ilep++)
        // {
        //     TString type = "";
        //     if(vLooseLeptons[ilep].isElectron()) type = "Ele";
        //     else if(vLooseLeptons[ilep].isMuon()) type = "Mu";
        //
        //     cout<<ilep<<" "<<type<<": pT = "<<vLooseLeptons.at(ilep).pt<<" / eta = "<<vLooseLeptons.at(ilep).eta<<" / phi = "<<vLooseLeptons.at(ilep).phi<<" : Charge = "<<vLooseLeptons.at(ilep).charge<<" / ";
        //     cout<<endl<<endl;
        // }

        cout<<FMAG("Muons : ")<<endl;
        for(int imu=0; imu<vMuon->size(); imu++)
        {
            cout<<"------------------"<<endl<<endl;
            cout<<"Muon "<<imu<<" : pT = "<<vMuon->at(imu).pt<<" / eta = "<<vMuon->at(imu).eta<<endl;
            cout<<FRED("(MVA = "<<vMuon->at(imu).lepMVA() <<") / MVA_ttH = "<<vMuon->at(imu).lepMVA_TTH<<" / lepMVA_jetCSV = "<<vMuon->at(imu).lepMVA_jetBTagCSV()<<" ")<<endl;
            cout<<"jetPtRatio = "<<vMuon->at(imu).lepMVA_jetPtRatio()<<" / isMedium = "<<vMuon->at(imu).isMedium()<<endl<<endl;
        }

        cout<<"------------------"<<endl<<endl;

        cout<<FMAG("Electrons : ")<<endl;
        for(int iele=0; iele<vElectron->size(); iele++)
        {
            cout<<"------------------"<<endl;
            cout<<"Ele "<<iele<<" : pT = "<<vElectron->at(iele).pt<<" / eta = "<<vElectron->at(iele).eta<<endl;
            cout<<FRED("(MVA = "<<vElectron->at(iele).lepMVA() <<") / MVA_ttH = "<<vElectron->at(iele).lepMVA_TTH<<" / lepMVA_jetCSV = "<<vElectron->at(iele).lepMVA_jetBTagCSV()<<" ")<<endl;

            cout<<"isLooseTTH : "<<vElectron->at(iele).isLooseTTH<<" / _mvaNonTrigV0 = "<<vElectron->at(iele).mvaNonTrigV0()<<" / nofLostHit "<<vElectron->at(iele).nlosthits()<<" / tightQ = "<<vElectron->at(iele).tightCharge<<" / CV = "<<vElectron->at(iele).cutEventSel()<<" / lepMVA_jetPtRatio = "<<vElectron->at(iele).lepMVA_jetPtRatio()<<endl;
            cout<<"Pass conditions = "<<vElectron->at(iele).passConditions()<<" / Pass MuOverlap = "<<vElectron->at(iele).passMuOverlap()<<endl;
            cout<<" / sigmaIetaIeta = "<<vElectron->at(iele).sigmaIetaIeta<<endl;
            cout<<" / hadronicOverEm = "<<vElectron->at(iele).hadronicOverEm()<<endl;
            cout<<" / deltaEtaSuperClusterTrackAtVtx = "<<vElectron->at(iele).deltaEtaSuperClusterTrackAtVtx()<<endl;
            cout<<" / deltaPhiSuperClusterTrackAtVtx = "<<vElectron->at(iele).deltaPhiSuperClusterTrackAtVtx()<<endl;
            cout<<" / superCluster_eta = "<<vElectron->at(iele).superCluster_eta<<endl;
            cout<<" / _eInvMinusPinv = "<<vElectron->at(iele).eInvMinusPinv()<<endl;
            cout<<" / _ecalEnergy = "<<vElectron->at(iele).ecalEnergy()<<endl;
            cout<<" / eSuperClusterOverP = "<<vElectron->at(iele).eSuperClusterOverP()<<endl;
        }

        cout<<"------------------"<<endl<<endl;
        cout<<"------------------"<<endl<<endl<<endl<<endl<<endl;
    }
    // ======================================================================================================


    fillOutputTree();

    return;
}
*/
