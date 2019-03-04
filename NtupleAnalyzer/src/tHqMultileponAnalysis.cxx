// ThreeLeptonSelection_THQ3l
// tHqMultileponAnalysis::fillOutputTree()


//----------  FIXMES  -------------------------
// - pT 25/15/15 in 3l ??
// - Pileup disactivated !!!
// - noisy jets ?
// - 4l SR, ZZ CR : differences with NTP code ?
//--------------------------------------------

#include "../include/tHqMultileponAnalysis.h"

#include "TSystem.h"
#include "SignalExtractionMVA.cxx"
#include "Helper.cxx"
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

using namespace std;

bool DEBUG = false;

// --- OPTIONS
//--------------------------------------------
bool do_tHq_analysis = false; // false <-> will save events passing >=1 ttH cat ; true <-> idem with tHq cat. Jet vectors will also be different, and thus other variables

bool add_JES_TTrees = true;
bool add_JER_TTrees = true; //adds 1 TTree with JER smearing applied (buggy?), and 2 JER variations TTrees

bool apply_JER_smearing = false; //NB : seems to have large effect, must xcheck ! true <-> smear pt of jets with JER correcting factor

bool make_ntuples_for_overlap_studies = true; // true <-> adds histograms containing overlaps between categs
bool add_orthogocal_cat = true; //true <-> adds new categs for overlap/ortho studies

bool write_branches_forMEM = true; //true <-> write inputs necessary for MEM code

bool write_allScale_Variations = true; //true <-> write all scale variations (+ sumWeights) separately, for studies

//WARNING : for tHq sample, storing LHE weights increased file space by 10 ! Must select weights to store ?
bool write_LHE_weights_allFiles = false; //true <-> store LHE weights for all samples (~1k weights per event for now...) ; else only THQ/THW

bool dump_synchro_info = false;
//--------------------------------------------

//--------------------------------------------
//--- DEBUGGING - global variables
// vector<double> v_cutflow(15); //Cutflow vector --> Count events for THQ_3l_SR selection, after each cut
// int counter = 0; //For debug checks






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
 * Destructor
 */
tHqMultileponAnalysis::~tHqMultileponAnalysis()
{
    cout<<"--- Destroy 'tHqMultileponAnalysis' object..."<<endl;

    delete vEvent;
    delete vElectronLoose; delete vElectronFakeable; delete vElectronTight;
    delete vMuonLoose; delete vMuonFakeable; delete vMuonTight;
    delete vTauLoose; delete vTauMedium; delete vTauFakeable;
    delete vJetLoose_original;
    delete vJetLoose_tmp;
    delete vTruth;

    // delete tOutput;
    for(int itree=0; itree<v_tOutput.size(); itree++)
    {
        delete v_tOutput[itree];
    }

    delete outputfile;

	delete f_QFwgt; delete f_FRwgt;

    delete h_PU_noCorr; delete h_PU_withCorr;
    delete h_nPV_noCorr; delete h_nPV_withCorr;
    delete h_overlap_ttH_tHq_cat; h_overlap_ttH_tHq_cat = NULL;
    delete h_totalYield_ttH_cat; h_totalYield_ttH_cat = NULL;
    delete h_totalYield_tHq_cat; h_totalYield_tHq_cat = NULL;

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
    _isdata = isdata;
    _xsec = xsec;
    _lumi = lumi;
    _nowe = nowe;
    _nmax = nmax;
    _outputFileName = outputFileName;
    _sampleName = sampleName;
    this->inputFileName = inputFileName;
    // _process = "toto";

    if(!isdata && xsec == 1) {cout<<BOLD(FRED("WARNING : xsec = 1 ! Most likely, the name of the process was not found/matched in the table.txt file ! The MC event weights are probably wrong !"))<<endl;}

    cout<<endl<<endl<<BOLD(UNDL(FBLU("* Sample name = "<<_sampleName<<"")))<<endl<<endl<<endl;

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

    bool debug_pileup = false; //true <-> print PU weights
    sf = new ScaleFactors(_sampleName, debug_pileup); //Opens all the SF files, fills SF histograms

    //Initialize once the sum of weights for all scale variations, before reading
    // sumWeights_nominal = -999;
    // sumWeights_scale_originalXWGTUP = -999;
    // sumWeights_scale_muF0p5 = -999;
    // sumWeights_scale_muF2 = -999;
    // sumWeights_scale_muR0p5 = -999;
    // sumWeights_scale_muR2 = -999;
    // sumWeights_scale_muR2muF2 = -999;
    // sumWeights_scale_muR0p5muF0p5 = -999;
    //Read/store the sum of weights for scale variations, if available (stored as histograms)
    //These correspond to *before preselection*
    // sf->Read_Scale_SumWeights(_sampleName, sumWeights_nominal, sumWeights_scale_originalXWGTUP, sumWeights_scale_muF0p5, sumWeights_scale_muF2, sumWeights_scale_muR0p5, sumWeights_scale_muR2, sumWeights_scale_muR2muF2, sumWeights_scale_muR0p5muF0p5);

    //PU reweighting
    h_PU_noCorr = new TH1F("", "", 100, 0, 100);
    h_PU_withCorr = new TH1F("", "", 100, 0, 100);
    h_nPV_noCorr = new TH1F("", "", 100, 0, 100);
    h_nPV_withCorr = new TH1F("", "", 100, 0, 100);

    //For overlap studies -- NB : events entering any of the 2 analyses are accounted for
    h_overlap_ttH_tHq_cat = new TH2F("", "", 25, 0, 25, 25, 0, 25); //Contains total yield of events overlapping between categories of both analyses
    h_totalYield_ttH_cat = new TH1F("", "", 25, 0, 25); //Contains total yield of events entering each ttH cat
    h_totalYield_tHq_cat = new TH1F("", "", 25, 0, 25); //Contains total yield of events entering each tHq cat

    InitFiles();

    Get_SumWeights();

    //Systematics TTrees
    v_systTree.push_back(""); //Default TTree
    if(!_isdata)
    {
        if(add_JES_TTrees) {v_systTree.push_back("JESUp"); v_systTree.push_back("JESDown");}
        if(add_JER_TTrees)
        {
            if(!apply_JER_smearing) {cout<<BOLD(FRED("--- You are not applying the JER smearing, but want to create a TTree corresponding to JER up/down variations. Are you sure ?"))<<endl;}

            v_systTree.push_back("JER"); //For now : don't apply JER in default tree, but only in a separate TTree for studies

            v_systTree.push_back("JERUp"); v_systTree.push_back("JERDown");
        }
    }

    //Create 2 output TTree per syst (+1 default)
    outputfile->cd();
    v_tOutput.resize(v_systTree.size());
    for(int itree=0; itree<v_tOutput.size(); itree++)
    {
        if(!itree) {v_tOutput[itree] = new TTree("Tree", "Tree");}
        else {v_tOutput[itree] = new TTree(v_systTree[itree], v_systTree[itree]);}
    }

    if(_sampleName.Contains("THQ_ctcvcp") || _sampleName.Contains("THW_ctcvcp")) {write_LHE_weights_allFiles = true;} //Always store LHE info for THQ/THW

    //FakeRates : define here all the FR variations that we want to compute/store
    //NB : each one calls a different FR histogram, cf FakeRate.cxx code
    //NB : nominal is stored in a dedicated variable, independently
    v_FR_type; //All types of shifted FR weight
    v_FR_type.push_back("FR_norm_elUp");
    v_FR_type.push_back("FR_norm_elDown");
    v_FR_type.push_back("FR_norm_muUp");
    v_FR_type.push_back("FR_norm_muDown");
    v_FR_type.push_back("FR_pt_elUp");
    v_FR_type.push_back("FR_pt_elDown");
    v_FR_type.push_back("FR_pt_muUp");
    v_FR_type.push_back("FR_pt_muDown");
    v_FR_type.push_back("FR_be_elUp");
    v_FR_type.push_back("FR_be_elDown");
    v_FR_type.push_back("FR_be_muUp");
    v_FR_type.push_back("FR_be_muDown");
    v_floats_FR_variations.resize(v_FR_type.size()); //reserve 1 float per FR variation
}


/**
 * Initialize vectors contained in input TTree
 */
void tHqMultileponAnalysis::InitTree()
{
    // Set branch addresses and branch pointers
    if (!fChain) return;

	//Object collections defined in NtupleProducer
    vEvent             = new std::vector<Event>();
    vElectronLoose     = new std::vector<Electron>();
    vElectronFakeable  = new std::vector<Electron>();
    vElectronTight     = new std::vector<Electron>();
    vMuonLoose         = new std::vector<Muon>();
    vMuonFakeable      = new std::vector<Muon>();
    vMuonTight         = new std::vector<Muon>();
    vTauLoose          = new std::vector<Tau>();
    vTauFakeable       = new std::vector<Tau>();
    vTauMedium         = new std::vector<Tau>();
    vJetLoose_original = new std::vector<Jet>();
    vJetLoose_tmp      = new std::vector<Jet>();
    vTruth             = new std::vector<Truth>();


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
    fChain->SetBranchAddress("JetLoose", &vJetLoose_original);
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
    cout<<"... Restored."<<endl;

    //Charge Flip Rate
    TString inputFileQF = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/ElectronChargeMisIdRates_2017_2018Jun22.root"; //HARD-CODED
    f_QFwgt    = new TFile (inputFileQF);
    fillQFhistos(f_QFwgt);

    //Fake Rate
    TString inputFileFR = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/FR_data_ttH_mva.root"; //HARD-CODED
    f_FRwgt    = new TFile (inputFileFR);
    Fill_FR_Histograms(f_FRwgt);

    //Synchro //By default, existing files are
    TString outdir = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/synchro/";
    if(dump_synchro_info) {mkdir(outdir.Data(), 0777);}

    //Get string in samplename
    TString inputFileName_tmp = inputFileName;
    int i = inputFileName_tmp.Index(".txt");
    inputFileName_tmp.Remove(i);
    i = inputFileName_tmp.Last('/');
    if(i>0) {inputFileName_tmp.Remove(0, i+1);}
    // cout<<"inputFileName = "<<inputFileName_tmp<<endl;

    if(dump_synchro_info)
    {
        outfile_2lSS_SR.open( (outdir+inputFileName_tmp+"_2lSS_SR.txt").Data() );
        outfile_3l_SR.open( (outdir+inputFileName_tmp+"_3l_SR.txt").Data() );
        outfile_ttW_CR.open( (outdir+inputFileName_tmp+"_ttW_CR.txt").Data() );
        outfile_ttZ_CR.open( (outdir+inputFileName_tmp+"_ttZ_CR.txt").Data() );
        outfile_WZ_CR.open( (outdir+inputFileName_tmp+"_WZ_CR.txt").Data() );
    }

    cout<<endl<<"-> FR/QFlip/etc. files opened"<<endl<<endl<<endl;

    return;
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

    vJetLoose.clear();
    vJetLoose_tHq.clear();
    vJetLoose_ttH.clear();
    vLightJets.clear();
    vLightJets_tHq.clear();
    vLightJets_ttH.clear();
	vLooseBTagJets.clear();

    LHEweights.clear();
    LHEweights_Ids.clear();

	return;
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
    is_tHq_3l_GammaConv       = 0;
    is_tHq_ttWctrl            = 0;
    is_tHq_ttWctrl_SR         = 0;
    is_tHq_ttWctrl_Fake       = 0;
    is_tHq_ttWctrl_Flip       = 0;
    is_tHq_ttWctrl_GammaConv  = 0;
    is_tHq_ttZctrl            = 0;
    is_tHq_ttZctrl_SR         = 0;
    is_tHq_ttZctrl_Fake       = 0;
    is_tHq_ttZctrl_GammaConv  = 0;
    is_tHq_WZctrl             = 0;
    is_tHq_WZctrl_SR          = 0;
    is_tHq_WZctrl_Fake        = 0;
    is_tHq_WZctrl_GammaConv   = 0;
    is_tHq_4l_SR              = 0;
    is_tHq_ZZctrl_SR          = 0;
    is_tHqFCNC_2lSS_SR        = 0;
    is_tHqFCNC_2lSS_Fake      = 0;
    is_tHqFCNC_2lSS_Flip      = 0;
    is_tHqFCNC_2lSS_GammaConv = 0;
    is_tHqFCNC_3l_SR          = 0;
    is_tHqFCNC_3l_Fake        = 0;
    is_tHqFCNC_3l_GammaConv   = 0;

    //ttH2017 categories
    is_ttH_2lSS              = 0;
    is_ttH_2lSS_Training     = 0;
    is_ttH_2lSS_SR           = 0;
    is_ttH_2lSS_Fake         = 0;
    is_ttH_2lSS_Flip         = 0;
    is_ttH_2lSS_GammaConv    = 0;
    is_ttH_3l                = 0;
    is_ttH_3l_Training       = 0;
    is_ttH_3l_SR             = 0;
    is_ttH_3l_Fake           = 0;
    is_ttH_3l_GammaConv      = 0;
    is_ttH_ttWctrl           = 0;
    is_ttH_ttWctrl_SR        = 0;
    is_ttH_ttWctrl_Fake      = 0;
    is_ttH_ttWctrl_Flip      = 0;
    is_ttH_ttWctrl_GammaConv = 0;
    is_ttH_ttZctrl           = 0;
    is_ttH_ttZctrl_SR        = 0;
    is_ttH_ttZctrl_Fake      = 0;
    is_ttH_ttZctrl_GammaConv = 0;
    is_ttH_WZctrl            = 0;
    is_ttH_WZctrl_SR         = 0;
    is_ttH_WZctrl_Fake       = 0;
    is_ttH_WZctrl_GammaConv  = 0;
    is_ttH_4l_SR             = 0;
    is_ttH_ZZctrl_SR         = 0;

    is_ttH_2lSS_SR_fwd    = 0;
    is_ttH_3l_SR_fwd      = 0;
    is_tHq_2lSS_SR_fwd    = 0;
    is_tHq_3l_SR_fwd      = 0;

    is_ttH_2lSS_SR_btag      = 0;
    is_ttH_3l_SR_btag        = 0;
    is_tHq_2lSS_SR_btag      = 0;
    is_tHq_3l_SR_btag        = 0;

    is_ttH_2lSS_SR_njet = 0;
    is_ttH_3l_SR_njet   = 0;
    is_tHq_2lSS_SR_njet = 0;
    is_tHq_3l_SR_njet   = 0;
    is_ttH_ttWctrl_SR_njet = 0;

    is_tHq_2lSS_SR_fwd2     = 0;
    is_ttH_ttWctrl_SR_fwd2  = 0;

    is_tHq_3l_SR_njet2 = 0;
    is_ttH_3l_SR_njet2 = 0;

    is_tHq_3l_SR_njet3 = 0;
    is_ttH_3l_SR_njet3 = 0;

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

    return;
}


/**
 * Initialize all variables and category : to be called between before each selection function (else event can activate 2 categories at once)
 */
void tHqMultileponAnalysis::InitVariables()
{
	// signal_3l_TT_MVA = -9; signal_3l_TTV_MVA = -9;
	// signal_2lss_TT_MVA = -9; signal_2lss_TTV_MVA = -9;

	channel = -1;

    //tHq2016 input vars
    nJet25               = -999;
    nJetLoose            = -999;
    maxEtaJet25          = -999;
    lepCharge            = -999;
    nJetEta1             = -999;
    dEtaFwdJetBJet       = -999;
    dEtaFwdJet2BJet      = -999;
    dEtaFwdJetClosestLep = -999;
    dPhiHighestPtSSPair  = -999;
    minDRll              = -999;
    Lep3Pt               = -999;

    //ttH2017 input vars
    lep1_conePt    = -999;
    lep2_conePt    = -999;
    lep3_conePt    = -999;
    mindr_lep1_jet = -999;
    mindr_lep2_jet = -999;
    mT_lep1        = -999;
    mT_lep2        = -999;
    max_lep_eta    = -999;

    //More input vars, to be tested
    minv_FwdJetBJet  = -999;
    FwdJetEta        = -999;
    FwdJetPt         = -999;
    LeadJetEta       = -999;
    LeadJetPt        = -999;
    dRjj_max         = -999;
    deepCSV_max      = -999;
    deepCSV_2nd      = -999;
    Mjj_max          = -999;
    dPhiLepBJet_max  = -999;
    dPhijj_max       = -999;
    m3l              = -999;
    dPhiLepLep_max   = -999;
    top_mass         = -999;
    mTW              = -999;
    lW_asym_mtop     = -999;
    dRBjetRecoilJet  = -999;
    dRLepWRecoilJet  = -999;
    RecoilJetPt      = -999;
    RecoilJetEta     = -999;
    LepWPt           = -999;
    LepWEta          = -999;
    top_Pt           = -999;
    HjTag_max        = -999;
    HjTag_mean       = -999;
    mass_LepBJet_min = -999;
    sum_jetPt        = -999;

    //Additional vars
    hardestBjetPt = -999; hardestBjetEta = -999;
    lep1Eta = -999; lep2Eta = -999; lep3Eta = -999; lep1Phi = -999; lep2Phi = -999; lep3Phi = -999;

    //Additional vars for FCNC analysis
    // nSoftJets                 = -999;
    min_dr_lep_bjet           = -999;
    min_dr_lep_lightjet       = -999;
    lW_asym                   = -999;
    ratio_lep3pt_closestJetPt = -999;
    dPhiLepLep_hardestOS  = -999;

	return;
}


void tHqMultileponAnalysis::Get_SumWeights()
{
    // cout<<"Entering Get_SumWeights()"<<endl;

    // sumWeights_nominal = 1;
    // sumWeights_scale_originalXWGTUP = 1;
    // sumWeights_scale_muF0p5 = 1;
    // sumWeights_scale_muF2 = 1;
    // sumWeights_scale_muR0p5 = 1;
    // sumWeights_scale_muR2 = 1;
    // sumWeights_scale_muR2muF2 = 1;
    // sumWeights_scale_muR0p5muF0p5 = 1;

    sumWeights_SMcoupling = 1;

    TString samplename_tmp = _sampleName;

    //HARD-CODED -- difference in automatic naming conventions for sample extension from FlatTrees to IPHCNtuple code...
    if(_sampleName == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8") {samplename_tmp = "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM";}
    else if(_sampleName == "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8") {samplename_tmp = "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_v1_MINIAODSIM";}

    TString path_sumweight_file = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/merged_histograms/"+samplename_tmp+".root";

    if(!Check_File_Existence(path_sumweight_file)) {cout<<BOLD(FRED("File "<<path_sumweight_file<<" not found ! Can not retrieve sum of weights to properly rescale events !"))<<endl; return;}
    TFile* f = TFile::Open(path_sumweight_file);
    TH1F* h_tmp = 0;

    //Read sums of weights
    hSumWeights = 0;
    if(!_isdata)
    {
        if(!f->GetListOfKeys()->Contains("hSumWeights")) {cout<<BOLD(FRED("hSumWeights histo not found !"))<<endl;}
        else
        {
            cout<<"-> Histogram containing sums of weights (nominal, scale, etc.) opened"<<endl<<endl<<endl;
            hSumWeights = (TH1F*) f->Get("hSumWeights")->Clone();
        }
    }

    //Read sums of weights for kT/kV weights (only SM for now ?)
    if(!_isdata)
    {
        h_tmp = (TH1F*) f->Get("hLHE");
        sumWeights_SMcoupling = h_tmp->GetBinContent(12);
    }

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

    // for(int i=0; i<15; i++)
    // {
    //     v_cutflow[i] = 0;
    // }

    Long64_t nentries = fChain->GetEntries();
    int nentries_max = nentries;
    if ( _nmax != -1 && _nmax < nentries ) nentries_max = _nmax;

    cout<<endl<<FYEL("--------------------------------------------")<<endl;
    cout<<BOLD(FYEL("Will process "<<nentries_max<<"/"<<nentries<<" events"))<<endl;
    cout<<BOLD(FYEL("\t * "<<v_systTree.size()<<" TTrees (default, JES, JER, ...)"))<<endl;
    cout<<FYEL("--------------------------------------------")<<endl<<endl;

    for(int itree=0; itree<v_systTree.size(); itree++)
    {
        cout<<endl<<FYEL("--------------------------------------------")<<endl;
        cout<<BOLD(FYEL("== TTREE NAME : '"<<v_systTree[itree]<<"'"))<<endl;
        cout<<FYEL("--------------------------------------------")<<endl<<endl;

        initializeOutputTree(itree);

        Long64_t nbytes = 0, nb = 0;
        for (Long64_t jentry=0; jentry<nentries_max;jentry++)
        {
            // cout<<endl<<jentry<< " / "<<nentries_max<<endl;

            Long64_t ientry = fChain->LoadTree(jentry);
            if (ientry < 0) break;

            if(jentry%200000 == 0) std::cout << "--- "<<jentry<< " / "<<nentries_max<<endl;

            // if(jentry > 100) {break;}

            nb = fChain->GetEntry(jentry);   nbytes += nb;

            InitCollections(); //Re-init object collections for each event
            InitCustomCategories(); //Re-init custom categories

            metpt = 0; metphi = 0; mHT = 0; metLD = 0;

            //Decay modes and jets flavours
            higgs_daughter_id = -1;
            wz_jetFlav_b   = 0;
            wz_jetFlav_c   = 0;
            wz_jetFlav_l   = 0;

            event_id  = vEvent->at(0).id;
            metpt       = vEvent->at(0).metpt;
            metphi       = vEvent->at(0).metphi;
            // metLD  = Compute_metLD_Alternative(metpt); //Helper func
            metLD     = vEvent->at(0).metLD;
            event_run = vEvent->at(0).run;
            event_lumi = vEvent->at(0).lumi;
            nPU       = vEvent->at(0).mc_pu_trueNumInt;
            nPV       = vEvent->at(0).pv_n;

            //LHE scale variation weights
            weight_originalXWGTUP = vEvent->at(0).weight_originalXWGTUP;
            weight_scale_muF0p5 = vEvent->at(0).weight_scale_muF0p5;
            weight_scale_muF2 = vEvent->at(0).weight_scale_muF2;
            weight_scale_muR0p5 = vEvent->at(0).weight_scale_muR0p5;
            weight_scale_muR2 = vEvent->at(0).weight_scale_muR2;
            weight_scale_muR2muF2 = vEvent->at(0).weight_scale_muR2muF2;
            weight_scale_muR0p5muF0p5 = vEvent->at(0).weight_scale_muR0p5muF0p5;

            //LHE kT/kV weights
            if(write_LHE_weights_allFiles && !itree) //only for nominal tree, and if user asks for it
            {
                LHEweights = vEvent->at(0).pdf_weights;
                LHEweights_Ids = vEvent->at(0).pdf_ids;

                //FIXME - sumWeights_SMcoupling
            }


            if(!_isdata )
            {
                mc_weight = vEvent->at(0).mc_weight;
                mc_weight_originalValue = vEvent->at(0).mc_weight_originalValue;

                //NB : old, obsolete weight (used to force gen weight to +-1, but can be different)
                //NB : shouldn't be used anymore ! Reads nof entries from 'table_xxx.txt' file, but not up-to-date anymore !!
                weight_old = mc_weight*_lumi*_xsec/_nowe;

                //NEW -- set all 'nowe' to 1 in table.txt, and instead read it from a root file obtained separately (from code NTProducer/src/Get_Merged_Histograms_From_FlatTrees.exe)
                if(!hSumWeights) {cout<<BOLD(FRED("Error : File containing sums of weights was not found, can not scale MC events. Abort ! (you need to produce it first)"))<<endl; return;}

                weight = mc_weight_originalValue*_lumi*_xsec/hSumWeights->GetBinContent(2);

                if(!hSumWeights->GetBinContent(2)) {weight = mc_weight_originalValue*_lumi*_xsec/hSumWeights->GetBinContent(1);} //Tmp fix : in WW_DPS sample, only the first sum of weights is filled for now...

                // cout<<"//--------------------------------------------"<<endl;
                // cout<<"_lumi = "<<_lumi<<endl;
                // cout<<"_xsec = "<<_xsec<<endl;
                // cout<<"_nowe = "<<_nowe<<endl;
                // cout<<"mc_weight = "<<mc_weight<<endl;
                // cout<<"mc_weight_originalValue = "<<mc_weight_originalValue<<endl;
                // cout<<"hSumWeights->GetBinContent(2) = "<<hSumWeights->GetBinContent(2)<<endl;
                // cout<<"=> weight_old = "<<weight_old<<endl;
                // cout<<"=> weight = "<<weight<<endl;
                // cout<<"//--------------------------------------------"<<endl;

                higgs_daughter_id = vTruth->at(0).higgs_daughter_id;
                if(_sampleName.Contains("WZTo")) {Get_HadronFlavour_WZsample();}
            }
            else
            {
                weight     = 1.;
                mc_weight  = 1.;
                lepton_SF  = 1.;
                trigger_SF = 1.;
                btag_SF    = 1.;
                PU_SF      = 1.;
                total_SF   = 1.;
            }

            weightflip = 1;
            weightfake = 1;

            //Fill pileup/nPV histos, before correction --only for nominal
            if(itree == 0)
            {
                h_PU_noCorr->Fill(nPU, weight);
                h_nPV_noCorr->Fill(nPV, weight);
            }


            //--------------------------------------------
            // #    # #    #  ####  #    #  ####
            // ##  ## #    # #    # ##   # #
            // # ## # #    # #    # # #  #  ####
            // #    # #    # #    # #  # #      #
            // #    # #    # #    # #   ## #    #
            // #    #  ####   ####  #    #  ####
            //--------------------------------------------

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



            //-- Or can directly use definitions from ttH2017 (coded in NtupleProducer)
            //NB : taus not included anywhere in my code !
            for(unsigned int itau=0; itau < vTauLoose->size() ; itau++)
            {
                Lepton l; l.setLepton(&vTauLoose->at(itau),itau,0,0);
            }
            for(unsigned int itau=0; itau < vTauFakeable->size() ; itau++)
            {
                Lepton l; l.setLepton(&vTauFakeable->at(itau),itau,0,0);
            }




            //--------------------------------------------
            //      # ###### #####  ####
            //      # #        #   #
            //      # #####    #    ####
            //      # #        #        #
            // #    # #        #   #    #
            //  ####  ######   #    ####
            //--------------------------------------------

            nLooseBJets  = 0;
            nMediumBJets = 0;
            nLightJets_eta_Gt2 = 0; //-- try CR enriched in ttbar (tHq2016)
            nLightJets_Fwd = 0;

            is_hasJetNoisyHCAL     = 0;
            is_hasManyJetNoisyHCAL = 0;
            JetNoisyHCALPt         = 0;
            JetNoisyHCALEta        = 0;

            //-- ttH btag definitions
            //NB : btag def must be consistent w/ SelectBjets()
            // isLooseBTag  =  (deepCSVb+deepCSVbb > 0.1522);
            // isMediumBTag =  (deepCSVb+deepCSVbb > 0.4941);
            // isTightBTag = (deepCSVb+deepCSVbb > 0.8001);


            //'vJetLoose_tmp' is a copy of 'vJetLoose_original', but here we can modify the pT of the jets (JES/JER variations)
            //Can also decide to apply the jet pt JER-smearing here
            //After this, the "original" jet vector is not used anymore
            vJetLoose_tmp->resize(0);
            for(unsigned int ijet=0; ijet < vJetLoose_original->size() ; ijet++)
            {
                vJetLoose_tmp->push_back(vJetLoose_original->at(ijet));

                if(apply_JER_smearing || v_systTree[itree] == "JER")
                {
                    vJetLoose_tmp->at(ijet).pt*= vJetLoose_tmp->at(ijet).JER_corr;
                }

                if(v_systTree[itree] == "JESUp")
                {
                    vJetLoose_tmp->at(ijet).pt = vJetLoose_original->at(ijet).pt_JES_up;
                    vJetLoose_tmp->at(ijet).E = vJetLoose_original->at(ijet).E_JES_up;
                }
                else if(v_systTree[itree] == "JESDown")
                {
                    vJetLoose_tmp->at(ijet).pt = vJetLoose_original->at(ijet).pt_JES_down;
                    vJetLoose_tmp->at(ijet).E = vJetLoose_original->at(ijet).E_JES_down;
                }
                else if(v_systTree[itree] == "JERUp")
                {
                    vJetLoose_tmp->at(ijet).pt = vJetLoose_original->at(ijet).pt_JER_up;
                    vJetLoose_tmp->at(ijet).E = vJetLoose_original->at(ijet).E_JER_up;
                }
                else if(v_systTree[itree] == "JERDown")
                {
                    vJetLoose_tmp->at(ijet).pt = vJetLoose_original->at(ijet).pt_JER_down;
                    vJetLoose_tmp->at(ijet).E = vJetLoose_original->at(ijet).E_JER_down;
                }
                else if(v_systTree[itree] != "" && v_systTree[itree] != "JER") {cout<<BOLD(FRED("ERROR ! systName '"<<v_systTree[itree]<<"' not known ! Abort"))<<endl; return;}

                if(itree>0 && vJetLoose_tmp->at(ijet).pt == 0)
                {
                    cout<<BOLD(FRED("ERROR ! JER or JES values not found in input TTree ! Abort"))<<endl; return;
                }
            }


            for(unsigned int ijet=0; ijet < vJetLoose_tmp->size() ; ijet++)
            {
                // cout<<"vJetLoose_tmp->at("<<ijet<<").pt "<<vJetLoose_tmp->at(ijet).pt<<" / eta = "<<vJetLoose_tmp->at(ijet).eta<<" / pT_JESup = "<<vJetLoose_tmp->at(ijet).pt_JES_up<<endl;

                //pT>25 cut. NB : for JES/JER variations, the pt variable was modified accordingly above
                if(vJetLoose_tmp->at(ijet).pt < 25) {continue;}
                if(fabs(vJetLoose_tmp->at(ijet).eta) > 5.0) {continue;}

                //Checks -- Noisy HCAL
                if(fabs(vJetLoose_tmp->at(ijet).eta) > 2.7 && fabs(vJetLoose_tmp->at(ijet).eta) < 3.0) //Noisy HCAL region -- for xchecks, cuts, ...
                {
                    if(is_hasJetNoisyHCAL==1) {is_hasManyJetNoisyHCAL = 1;}
                    else
                    {
                        is_hasJetNoisyHCAL = 1;
                        JetNoisyHCALPt = vJetLoose_tmp->at(ijet).pt;
                        JetNoisyHCALEta = vJetLoose_tmp->at(ijet).eta;
                    }
                }

                //--- Store B-tagged jets //Same for both analyses
                if(vJetLoose_tmp->at(ijet).isLooseBTag && fabs(vJetLoose_tmp->at(ijet).eta) < 2.4)
                {
                    vLooseBTagJets.push_back(vJetLoose_tmp->at(ijet)); nLooseBJets++;

                    if(vJetLoose_tmp->at(ijet).isMediumBTag) {nMediumBJets++;}

                    vJetLoose_tHq.push_back(vJetLoose_tmp->at(ijet));
                    vJetLoose_ttH.push_back(vJetLoose_tmp->at(ijet));
                }
                else //Store light jets //Differences between analyses
                {
                    //Central jets, used in both analyses
                    if(fabs(vJetLoose_tmp->at(ijet).eta) < 2.4)
                    {
                        vLightJets_tHq.push_back(vJetLoose_tmp->at(ijet));
                        vLightJets_ttH.push_back(vJetLoose_tmp->at(ijet));

                        vJetLoose_tHq.push_back(vJetLoose_tmp->at(ijet));
                        vJetLoose_ttH.push_back(vJetLoose_tmp->at(ijet));
                    }
                    //Forward hard jets, used in tHq analysis
                    // else if(vJetLoose_tmp->at(ijet).pt > 40 && fabs(vJetLoose_tmp->at(ijet).eta) > 2.4)
                    else if(fabs(vJetLoose_tmp->at(ijet).eta) < 2.7 || fabs(vJetLoose_tmp->at(ijet).eta) > 3.0)
                    {
                        vLightJets_tHq.push_back(vJetLoose_tmp->at(ijet));
                        vJetLoose_tHq.push_back(vJetLoose_tmp->at(ijet));
                        nLightJets_Fwd++;
                    }
                    //Forward soft jets -- removed completely from tHq2016 analysis. Keep ones outside noisy region for now
                    else //Problematic events : fwd jets with pT<40 -- only remove jets in transition region ? or all fwd jets with pt<40 ?
                    {
                        // if(fabs(vJetLoose_tmp->at(ijet).eta) > 2.7 && fabs(vJetLoose_tmp->at(ijet).eta) < 3.0) {continue;}
                        if(vJetLoose_tmp->at(ijet).pt < 60) {continue;} //Remove noise jets with pt<60 for now, like tZq 2017

                        vLightJets_tHq.push_back(vJetLoose_tmp->at(ijet));
                        vJetLoose_tHq.push_back(vJetLoose_tmp->at(ijet));
                        nLightJets_Fwd++;
                    }

                    if(fabs(vJetLoose_tmp->at(ijet).eta) > 2.0) {nLightJets_eta_Gt2++;}
                }
            } //jet loop

            nJets_tHq = vJetLoose_tHq.size();
            nLightJets_tHq = vLightJets_tHq.size();

            nJets_ttH = vJetLoose_ttH.size();
            nLightJets_ttH = vLightJets_ttH.size();

            //-- Select jet vector to be used, depending on analysis (NB : changes input vars, etc.)
            if(do_tHq_analysis) {vJetLoose = vJetLoose_tHq; vLightJets = vLightJets_tHq;}
            else {vJetLoose = vJetLoose_ttH; vLightJets = vLightJets_ttH;}

            nJets = vJetLoose.size();
            nLightJets = vLightJets.size();

            // cout<<"==> nJets_ttH = "<<nJets_ttH<<endl;


            //--------------------------------------------
            //  ####  #####  #####  ###### #####  # #    #  ####
            // #    # #    # #    # #      #    # # ##   # #    #
            // #    # #    # #    # #####  #    # # # #  # #
            // #    # #####  #    # #      #####  # #  # # #  ###
            // #    # #   #  #    # #      #   #  # #   ## #    #
            //  ####  #    # #####  ###### #    # # #    #  ####
            //--------------------------------------------

            //Sort object vectors by ConePt
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

            //Read trigger bits (15 trigger paths used in total, see NTP/src/EventExt.cxx)
            bool TRIGm   = vEvent->at(0).trig_m  ;
            bool TRIGe   = vEvent->at(0).trig_e  ;
            bool TRIGee  = vEvent->at(0).trig_ee ;
            bool TRIGmm  = vEvent->at(0).trig_mm ;
            bool TRIGem  = vEvent->at(0).trig_em ;
            bool TRIGeee = vEvent->at(0).trig_eee;
            bool TRIGmme = vEvent->at(0).trig_emm;
            bool TRIGeem = vEvent->at(0).trig_eem;
            bool TRIGmmm = vEvent->at(0).trig_mmm;

            //Updated
            bool E = false, M = false, EE = false, MM = false, EM = false;
            if(TRIGmmm || TRIGmm)               {MM = true;}
            if(TRIGeee || TRIGee)               {EE = true;}
            if(TRIGmme || TRIGeem || TRIGem)    {EM = true;}
            if(TRIGm)                           {M = true;}
            if(TRIGe)                           {E = true;}

            bool mmdataset = _sampleName.Contains("DoubleM");
            bool eedataset = _sampleName.Contains("DoubleE");
            bool emdataset = _sampleName.Contains("MuonE");
            bool mdataset  = _sampleName.Contains("SingleM");
            bool edataset  = _sampleName.Contains("SingleE");

            is_trigger_1lep = false;
            is_trigger_2lep = false;
            is_trigger_3lep = false;
            is_trigger     = false;
            is_trigger_2lss   = false;
            is_trigger_3l     = false;

            if(TRIGe || TRIGm) {is_trigger_1lep = true;}
            if(TRIGee || TRIGmm || TRIGem) {is_trigger_2lep = true;}
            if(TRIGeee || TRIGmme || TRIGeem || TRIGmmm) {is_trigger_3lep = true;}

            //From ttH2017 -- https://github.com/peruzzim/cmgtools-lite/blob/94X_dev_ttH/TTHAnalysis/python/tools/functionsTTH.py#L414
            //Require that events in 2lss categories pass 'is_trigger_2lss' (events in 3l cat. can pass any trigger)
            //Used for event selection
            is_trigger_2lss = (is_trigger_1lep || is_trigger_2lep);
            is_trigger_3l = (is_trigger_2lss || is_trigger_3lep);

            //Hypothesis : events satisfying EE will end up in eedataset, etc.
            if(_isdata) //For data
            {
                if      ( MM  &&                               (mmdataset) ) {is_trigger = true;}
                else if ( !MM && EE  &&                        (eedataset) ) {is_trigger = true;}
                else if ( !MM && !EE && EM  &&                 (emdataset) ) {is_trigger = true;}
                else if ( !MM && !EE && !EM && E  &&           (edataset ) ) {is_trigger = true;}
                else if ( !MM && !EE && !EM && !E && M &&      (mdataset ) ) {is_trigger = true;}
            }
            else //For MC
            {
                if(EM || MM || EE || M || E) {is_trigger = true;}
            }

            if(!is_trigger)
            {
                weight = 0;
                continue; //Skip data events failing trigger
            }



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
            pass_Zveto_ee = true;
            pass_cleanup = false;

            mllll = -1;
            inv_mll = 0;

            float mll_min = 10E+10;
            float mll_z_min = 10E+10;
            for(int il=0; il<vLeptonLoose.size(); il++)
            {
                TLorentzVector l1;
                TLorentzVector l1_cleanup;
                l1.SetPtEtaPhiE(vLeptonLoose.at(il).pt, vLeptonLoose.at(il).eta, vLeptonLoose.at(il).phi, vLeptonLoose.at(il).E);
                l1_cleanup.SetPtEtaPhiE(vLeptonLoose.at(il).pt, vLeptonLoose.at(il).eta, vLeptonLoose.at(il).phi, vLeptonLoose.at(il).E); //Using real pT for dilep mass cleaning

                for(int ill=il+1; ill<vLeptonLoose.size(); ill++)
                {
                    TLorentzVector l2;
                    TLorentzVector l2_cleanup; //Using real pT for dilep mass cleaning
                    l2.SetPtEtaPhiE(vLeptonLoose.at(ill).pt, vLeptonLoose.at(ill).eta, vLeptonLoose.at(ill).phi, vLeptonLoose.at(ill).E);
                    l2_cleanup.SetPtEtaPhiE(vLeptonLoose.at(ill).pt, vLeptonLoose.at(ill).eta, vLeptonLoose.at(ill).phi, vLeptonLoose.at(ill).E);

                    float mll_cleanup = (l1_cleanup+l2_cleanup).M();
                    if(mll_cleanup < mll_min) {mll_min = mll_cleanup;}

                    bool SFOS = (vLeptonLoose.at(il).id == -vLeptonLoose.at(ill).id);

                    if(!SFOS) {continue;}

                    float mll_z = fabs((l1+l2).M()-91.2);
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

                            float tmp = (l1+l2+l3+l4).M();

                            if(tmp > 0 && tmp < mllll) {mllll = tmp;}
                        }
                    }
                }
            }
            if(mll_min > 12) {pass_cleanup = true;}
            if(mll_z_min > 10) {pass_Zveto = true;}

            nSFOS = 0; //Among FO leptons
            for(int il=0; il<vLeptonFakeable.size(); il++)
            {
                for(int ill=il+1; ill<vLeptonFakeable.size(); ill++)
                {
                    bool SFOS = (vLeptonFakeable.at(il).id == -vLeptonFakeable.at(ill).id);
                    if(SFOS) {nSFOS++;}
                }
            }


            //--- Changed : will veto 2l events if two hardest FO are electrons within Z peak
            if(vLeptonFakeable.size() >= 2 && abs(vLeptonFakeable[0].id) == abs(vLeptonFakeable[1].id) && abs(vLeptonFakeable[0].id) == 11)
            {
                TLorentzVector e1;
                e1.SetPtEtaPhiE(vLeptonFakeable[0].pt, vLeptonFakeable[0].eta, vLeptonFakeable[0].phi, vLeptonFakeable[0].E);
                TLorentzVector e2;
                e2.SetPtEtaPhiE(vLeptonFakeable[1].pt, vLeptonFakeable[1].eta, vLeptonFakeable[1].phi, vLeptonFakeable[1].E);
                if(fabs(91.2 - (e1+e2).M() ) < 10) {pass_Zveto_ee = false;}
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
            ThreeLeptonSelection_THQ3l(jentry);
            TwoLeptonSelection_THQ2lSSl(jentry);
            FourLeptonSelection_THQ4l(jentry);

            if(Sample_isUsed_forTraining() ) //Looser training selection (training samples only)
            {
                ThreeLeptonSelection_THQ3l_TrainingSelection(jentry);
                TwoLeptonSelection_THQ2lSS_TrainingSelection(jentry);
            }

            if(add_orthogocal_cat) {Define_New_Categorization();} //Fill additional categ. booleans, for overlap studies
    //--------------------------------------------
            //Debug : make sure categories implemented in NTProd/NTAnalyzer codes are the same
            //Read NTProd categ., store them in the "is_tHq_xxx" booleans
            //Keep the "is_ttH_xxx" booleans as computed by this code ==> Compare both !
            //NB : only makes sense to do the comparison if the NTProd files which are used were produced with jet_eta<2.4 !
            //NB : else, the categories defined in NTProd are wrong (does not account properly for forward jets)
            // InitPredefinedCategories();
            // Read_Predefined_Categories();
            // Check_Disagreement_Category_NTP_NTA_Codes();

    //--------------------------------------------
            //Check if event should be saved, and how many leptons to consider for computations
            bool event_belongs_toAnalysis = false;

            int save_event_nlep = Check_If_Save_Event(do_tHq_analysis);
            if(save_event_nlep > 0) {event_belongs_toAnalysis = true;} //if part of considered analysis, will save full event
            else {save_event_nlep = Check_If_Save_Event(!do_tHq_analysis);} //else, if at least part of other analysis, will account for overlap study only

            //If event is (at least) to be included for overlap, must apply weights
            if(save_event_nlep > 0)
            {
                Apply_FakeRate_Weight();

                if(save_event_nlep == 2) {Compute_Variables("2l"); Apply_ScaleFactors(2, itree);}
                else if(save_event_nlep == 3 || save_event_nlep == 4) {Compute_Variables("3l"); Apply_ScaleFactors(3, itree);}
                else {cout<<"ERROR !"<<endl; continue;}

                //Account all events from both analyses for overlap study
                Fill_Overlap_Histo();

                //If part of considered analysis - or if making ntuples for overlap studies -, save full event
                if(event_belongs_toAnalysis || make_ntuples_for_overlap_studies) {fillOutputTree(itree);}
            }

            //Fill pileup/nPV histos, after correction -- only for nominal
            if(itree == 0)
            {
                h_PU_withCorr->Fill(nPU, weight*PU_SF);
                h_nPV_withCorr->Fill(nPV, weight*PU_SF);
            }



            // --------------------------------------------
            //  ####  #   # #    #  ####  #    # #####   ####
            // #       # #  ##   # #    # #    # #    # #    #
            //  ####    #   # #  # #      ###### #    # #    #
            //      #   #   #  # # #      #    # #####  #    #
            // #    #   #   #   ## #    # #    # #   #  #    #
            //  ####    #   #    #  ####  #    # #    #  ####
            //--------------------------------------------

            //Dump event infos into synchro text file
            if(dump_synchro_info && !do_tHq_analysis)
            {
                if(is_ttH_2lSS_SR) {Dump_EventInfo_Synchro(outfile_2lSS_SR);}
                if(is_ttH_3l_SR) {Dump_EventInfo_Synchro(outfile_3l_SR);}
                if(is_ttH_ttWctrl_SR) {Dump_EventInfo_Synchro(outfile_ttW_CR);}
                if(is_ttH_ttZctrl_SR) {Dump_EventInfo_Synchro(outfile_ttZ_CR);}
                if(is_ttH_WZctrl_SR) {Dump_EventInfo_Synchro(outfile_WZ_CR);}
            }



            //--------------------------------------------
            // #####  ###### #####  #    #  ####     #####  #####  # #    # #####  ####  #    # #####
            // #    # #      #    # #    # #    #    #    # #    # # ##   #   #   #    # #    #   #
            // #    # #####  #####  #    # #         #    # #    # # # #  #   #   #    # #    #   #
            // #    # #      #    # #    # #  ###    #####  #####  # #  # #   #   #    # #    #   #
            // #    # #      #    # #    # #    #    #      #   #  # #   ##   #   #    # #    #   #
            // #####  ###### #####   ####   ####     #      #    # # #    #   #    ####   ####    #
            //--------------------------------------------

            //-- events in CERN list but not mine
            //2lSS SR
            // if(event_id == 13583598 || event_id == 14142098 || event_id == 18797791 || event_id == 18166527 || event_id == 18169210 || event_id == 18170519 || event_id == 18836791 || event_id == 16730021 || event_id == 16798450 || event_id == 14862385 || event_id == 16796194 || event_id == 17844055 || event_id == 18840851)

            // if(is_ttH_WZctrl_SR) {DEBUG = true;}
            // if(event_id==3726793 || event_id==15991527 || event_id==14506738 || event_id==14507281 || event_id==14507493) {DEBUG = true;}


            if(DEBUG)
            {
                cout<<endl<<endl<<BOLD(FBLU("------------ EVENT "<<std::setprecision(15)<<event_id<<" (LS="<<event_lumi<<") -----------"))<<endl<<endl;

                cout<<BOLD(FBLU("### Loose Btag Jets : "))<<endl;
                for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
                {
                    cout<<FBLU(""<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt<<" / eta = "<<vLooseBTagJets.at(ijet).eta<<" / phi = "<<vLooseBTagJets.at(ijet).phi<<" / CSV = "<<vLooseBTagJets.at(ijet).DeepCSVbtag<<" / isMedium ? = "<<vLooseBTagJets.at(ijet).isMediumBTag<<"")<<endl;
                }
                cout<<endl;
                cout<<BOLD(FCYN("### Light Jets : "))<<endl;
                for(int ijet=0; ijet<vLightJets.size(); ijet++)
                {
                    cout<<FCYN("* "<<ijet<<" pT = "<<vLightJets.at(ijet).pt<<" / eta = "<<vLightJets.at(ijet).eta<<" / phi = "<<vLightJets.at(ijet).phi<<" / CSV = "<<vLightJets.at(ijet).DeepCSVbtag<<""<<" / pt_JES_up = "<<vLightJets.at(ijet).pt_JES_up<<"")<<endl;
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

                    cout<<FYEL(<<" "<<ilep<<" "<<type<<": FO ? "<<vLeptonFakeable[ilep].isFakeableTTH<<", Tight ? "<<vLeptonFakeable[ilep].isTightTTH<<" / conepT = "<<vLeptonFakeable.at(ilep).conept<<" / eta = "<<vLeptonFakeable.at(ilep).eta<<" / phi = "<<vLeptonFakeable.at(ilep).phi<<"")<<endl;
                }

                cout<<endl;
                cout<<BOLD(FYEL("### Loose taus : "))<<endl<<endl;
                for(int ilep=0; ilep<vTauLoose->size(); ilep++)
                {
                    cout<<"-- vTauLoose["<<ilep<<"].id = "<<vTauLoose->at(ilep).id<<endl;
                }

                cout<<BOLD(FBLU("------------------------------------"))<<endl<<endl;

            } //Debug condition

        } //end event loop



    	// ofstream file_out("cutflow.txt");
        // for(int i=0; i<v_cutflow.size(); i++)
        // {
    	// 	cout<<"v_cutflow "<<i<<" = "<<v_cutflow[i]<<endl;
        // }

        // cout<<"counter = "<<counter<<endl;

        cout<<"--- Writing output TTree"<<endl;

    	outputfile->cd();

        v_tOutput[itree]->Write();
        // v_tOutput[itree]->CloneTree()->Write(); //segfault ?

        // tOutput->Write();
    } //itree loop

    //Write histograms to output file
    h_PU_noCorr->Write("h_PU_noCorr");
    h_PU_withCorr->Write("h_PU_withCorr");
    h_nPV_noCorr->Write("h_nPV_noCorr");
    h_nPV_withCorr->Write("h_nPV_withCorr");
    h_overlap_ttH_tHq_cat->Write("h_overlap_ttH_tHq_cat");
    h_totalYield_ttH_cat->Write("h_totalYield_ttH_cat");
    h_totalYield_tHq_cat->Write("h_totalYield_tHq_cat");
    if(!_isdata) {hSumWeights->Write("hSumWeights");}

    return;
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
 * 3l event selection -- Common base selection
 */
bool tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l(int evt)
{
	if(DEBUG) {cout<<FYEL("-- ThreeLeptonSelection_THQ3l")<<endl;}

    InitVariables();
	// v_cutflow[0]++;

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return 0;}

	//At least 3 FO leptons
    if(vLeptonFakeable.size() < 3) {return 0;}
    // v_cutflow[1]++;
    if(DEBUG) {cout<<"- Passed nLep cut"<<endl;}

    //pT cuts
    if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15 || vLeptonFakeable.at(2).conept < 10) {return 0;}
    // if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15 || vLeptonFakeable.at(2).conept < 15) {return 0;}
    // v_cutflow[2]++;
    if(DEBUG) {cout<<"- Passed pT cut"<<endl;}

    //No pair of loose leptons with mll < 12
    if(!pass_cleanup) {return 0;}
    // v_cutflow[3]++;
    if(DEBUG) {cout<<"- Passed cleanup cut"<<endl;}

    if( fabs(vLeptonFakeable.at(0).charge + vLeptonFakeable.at(1).charge + vLeptonFakeable.at(2).charge) != 1) {return 0;}
    // v_cutflow[4]++;
    if(DEBUG) {cout<<"- Passed charge cut"<<endl;}

    if(nJets_tHq < 2) {return 0;} //All cat. contain at least 2 jets (NB : tHq jets looser than ttH)
    // v_cutflow[5]++;
    if(DEBUG) {cout<<"- Passed first nJet cut"<<endl;}

//--------------------------------------------

    if(!ThreeLeptonSelection_THQ3l_Regions(evt) ) {return 0;}

    return 1;
}



/**
 * 3l categories definitions //Add here all requirements which are not shared by all categories
 */
bool tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l_Regions(int evt)
{
    bool pass_cat = false;

    //Not shared among all categories
    bool pass_isTight       = false;
    bool pass_MCMatching    = false;
    bool is_GammaConv       = false;
    bool pass_m4l           = false;
    bool has_lessThan4Tight = false;
    bool has_noIsoTau       = false;

    if(vLeptonFakeable.at(0).isTightTTH && vLeptonFakeable.at(1).isTightTTH && vLeptonFakeable.at(2).isTightTTH) {pass_isTight = true;}
    if(_isdata || (vLeptonFakeable.at(0).hasMCMatch && vLeptonFakeable.at(1).hasMCMatch && vLeptonFakeable.at(2).hasMCMatch) ) {pass_MCMatching = true;}
    if(!_isdata && (vLeptonFakeable.at(0).hasPhotonMCMatch || vLeptonFakeable.at(1).hasPhotonMCMatch || vLeptonFakeable.at(2).hasPhotonMCMatch) ) {is_GammaConv = true;}
    if(mllll > 140 || mllll < 0) {pass_m4l = true;}
    if(vLeptonTight.size() < 4) {has_lessThan4Tight = true;}
    if(vTauLoose->size() == 0) {has_noIsoTau = true;} //NB : tight or loose taus ?
    bool pass_basic_SR_cuts(pass_m4l && has_lessThan4Tight && has_noIsoTau);

    //Depends on analysis type (jets)
    bool pass_njet_tth  = false;
    bool pass_nJets_tHq = false;
    bool pass_metLD     = true;

    if(nJets_ttH >= 2 && (nLooseBJets >= 2 || nMediumBJets >= 1) ) {pass_njet_tth = true;} //at least 2 b-loose or 1 b-medium
    if(nMediumBJets > 0 && nLightJets_tHq > 0) {pass_nJets_tHq = true;} //at least 1 b-medium & 1 fwd (light) jet

    if(nJets_ttH < 4)
    {
        if(nSFOS > 0 && metLD < 45) {pass_metLD = false;}
        else if(!nSFOS && metLD < 30) {pass_metLD = false;}
    }

    //Determine categories

    //-- 3l SR (ttH)
    if(pass_Zveto && pass_basic_SR_cuts && pass_metLD && pass_njet_tth)
    {
        is_ttH_3l = 1;
        if(pass_MCMatching && pass_isTight) {is_ttH_3l_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && !pass_isTight) {is_ttH_3l_Fake = 1; pass_cat = 1;}
        else if(is_GammaConv && pass_isTight) {is_ttH_3l_GammaConv = 1; pass_cat = 1;}
    }

    //-- 3l SR (tHq)
    if(pass_Zveto && pass_basic_SR_cuts && pass_nJets_tHq)
    {
        is_tHq_3l = 1;
        if(pass_MCMatching && pass_isTight) {is_tHq_3l_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && !pass_isTight) {is_tHq_3l_Fake = 1; pass_cat = 1;}
        else if(is_GammaConv && pass_isTight) {is_tHq_3l_GammaConv = 1; pass_cat = 1;}
    }

    //-- 3l SR (FCNC)
    if(pass_Zveto && pass_basic_SR_cuts && nJets_ttH >= 1 && nMediumBJets == 1)
    {
        if(pass_MCMatching && pass_isTight) {is_tHqFCNC_3l_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && !pass_isTight) {is_tHqFCNC_3l_Fake = 1; pass_cat = 1;}
        else if(is_GammaConv && pass_isTight) {is_tHqFCNC_3l_GammaConv = 1; pass_cat = 1;}
    }

    //ttZ CR (same for ttH/tHq)
    if(!pass_Zveto && pass_metLD && pass_njet_tth && pass_basic_SR_cuts)
    {
        is_ttH_ttZctrl = 1;
        if(pass_MCMatching && pass_isTight) {is_ttH_ttZctrl_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && !pass_isTight) {is_ttH_ttZctrl_Fake = 1; pass_cat = 1;}
        else if(is_GammaConv && pass_isTight) {is_ttH_ttZctrl_GammaConv = 1; pass_cat = 1;}

        is_tHq_ttZctrl = is_ttH_ttZctrl; is_tHq_ttZctrl_SR = is_ttH_ttZctrl_SR; is_tHq_ttZctrl_Fake = is_ttH_ttZctrl_Fake; is_tHq_ttZctrl_GammaConv = is_ttH_ttZctrl_GammaConv;
    }

    //WZ CR (same for ttH/tHq)
    if(!pass_Zveto && pass_basic_SR_cuts && pass_metLD && nJets_ttH >= 2 && !pass_njet_tth)
    {
        is_ttH_WZctrl = 1;
        if(pass_MCMatching && pass_isTight) {is_ttH_WZctrl_SR = 1; pass_cat = 1;} //removed basic cuts for now
        else if(pass_MCMatching && !pass_isTight) {is_ttH_WZctrl_Fake = 1; pass_cat = 1;}
        else if(is_GammaConv && pass_isTight) {is_ttH_WZctrl_GammaConv = 1; pass_cat = 1;}

        is_tHq_WZctrl = is_ttH_WZctrl; is_tHq_WZctrl_SR = is_ttH_WZctrl_SR; is_tHq_WZctrl_Fake = is_ttH_WZctrl_Fake; is_tHq_WZctrl_GammaConv = is_ttH_WZctrl_GammaConv;
    }


//--------------------------------------------

    return pass_cat;
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
	if(DEBUG) {cout<<FYEL("-- TwoLeptonSelection_THQ2lSS")<<endl;}

    InitVariables();

    //Must pass single or double lepton trigger
    if(weight==0 || !is_trigger_2lss) {return 0;}

    //At least 2 leptons
    if(vLeptonFakeable.size() < 2) {return 0;}
    if(vLeptonTight.size() > 2) {return 0;}
    if(DEBUG) {cout<<"- Passed nLep cut"<<endl;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    if(!vLeptonFakeable.at(0).passTightCharge || !vLeptonFakeable.at(1).passTightCharge ) {return 0;}
    if(DEBUG) {cout<<"- Passed TightQ cut"<<endl;}

    //pT cuts
    if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15) {return 0;}
    if(DEBUG) {cout<<"- Passed pT cut"<<endl;}

    //No pair of loose leptons with mll < 12
    if(!pass_cleanup) {return 0;}
    if(DEBUG) {cout<<"- Passed cleanup cut"<<endl;}

    //No loose tau
    if(vTauLoose->size() > 0) {return 0;}
    if(DEBUG) {cout<<"- Passed tau cut"<<endl;}

    //No ee pair with mll-mZ < 10
    if(!pass_Zveto_ee) {return 0;}
    if(DEBUG) {cout<<"- Passed Zee cut"<<endl;}

    if(fabs(vLeptonFakeable.at(0).id) == 11 && fabs(vLeptonFakeable.at(1).id) == 11 && metLD < 30) {return 0;}
    if(DEBUG) {cout<<"- Passed SS cut"<<endl;}

    if(nJets_tHq < 2) {return 0;} //All cat. contain at least 2 jets (NB : tHq jets looser than ttH)
    if(DEBUG) {cout<<"- Passed first jet cut"<<endl;}
//--------------------------------------------

    if(!TwoLeptonSelection_THQ2lSS_Regions(evt) ) {return 0;} //Check categories

    return 1;
}





/**
 * 2l categories selection
 */
bool tHqMultileponAnalysis::TwoLeptonSelection_THQ2lSS_Regions(int evt)
{
    bool pass_cat = false;

    //Not shared among all categories
    bool pass_isTight    = false;
    bool pass_MCMatching = false;
    bool is_GammaConv    = false;
    bool is_SS_pair      = false;
    bool contains_ele    = false;

    if(vLeptonFakeable.at(0).isTightTTH && vLeptonFakeable.at(1).isTightTTH) {pass_isTight = true;}
    if(_isdata || (vLeptonFakeable.at(0).hasChargeMCMatch && vLeptonFakeable.at(1).hasChargeMCMatch) ) {pass_MCMatching = true;} //In 2lSS, also require charge matching
    if(!_isdata && (vLeptonFakeable.at(0).hasPhotonMCMatch || vLeptonFakeable.at(1).hasPhotonMCMatch) ) {is_GammaConv = true;}
    if(vLeptonFakeable.at(0).charge * vLeptonFakeable.at(1).charge > 0) {is_SS_pair = true;}
    if( fabs(vLeptonFakeable.at(0).id) == 11 || fabs(vLeptonFakeable.at(1).id) == 11 ) {contains_ele = true;}

    //Depends on analysis type (jets)
    bool pass_njet_tth = false;
    bool pass_nJets_tHq = false;
    bool pass_njet_ttbar = false; //Inspired from tHq2016 CR ; store it as "ttW CR"

    if(nJets_ttH >= 3 && (nLooseBJets >= 2 || nMediumBJets >= 1) ) {pass_njet_tth = true;}
    if(nJets_tHq >= 2 && nMediumBJets > 0 && nLightJets_tHq > 0)
    {
        pass_nJets_tHq = true;

        if(nJets_tHq < 4 && !nLightJets_eta_Gt2)  {pass_njet_ttbar = true;}
    }

    //Determine categories

    //-- 2lSS SR (ttH)
    if(nJets_ttH >= 4 && pass_njet_tth)
    {
        is_ttH_2lSS = 1;
        if(pass_MCMatching && is_SS_pair && pass_isTight) {is_ttH_2lSS_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && is_SS_pair && !pass_isTight) {is_ttH_2lSS_Fake = 1; pass_cat = 1;}
        else if(pass_MCMatching && !is_SS_pair && contains_ele && pass_isTight) {is_ttH_2lSS_Flip = 1; pass_cat = 1; Apply_FlipRate_Weight();}
        else if(is_GammaConv && is_SS_pair && pass_isTight) {is_ttH_2lSS_GammaConv = 1; pass_cat = 1;}
    }
    //-- 2lSS SR (tHq)
    if(pass_nJets_tHq)
    {
        is_tHq_2lSS = 1;
        if(pass_MCMatching && is_SS_pair && pass_isTight) {is_tHq_2lSS_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && is_SS_pair && !pass_isTight) {is_tHq_2lSS_Fake = 1; pass_cat = 1;}
        else if(pass_MCMatching && !is_SS_pair && contains_ele && pass_isTight) {is_tHq_2lSS_Flip = 1; pass_cat = 1; Apply_FlipRate_Weight();}
        else if(is_GammaConv && is_SS_pair && pass_isTight) {is_tHq_2lSS_GammaConv = 1; pass_cat = 1;}
    }
    //-- 2lSS SR (FCNC)
    if(nJets_ttH >= 1 && nMediumBJets == 1)
    {
        if(pass_MCMatching && is_SS_pair && pass_isTight) {is_tHqFCNC_2lSS_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && is_SS_pair && !pass_isTight) {is_tHqFCNC_2lSS_Fake = 1; pass_cat = 1;}
        else if(pass_MCMatching && !is_SS_pair && contains_ele && pass_isTight) {is_tHqFCNC_2lSS_Flip = 1; pass_cat = 1; Apply_FlipRate_Weight();}
        else if(is_GammaConv && is_SS_pair && pass_isTight) {is_tHqFCNC_2lSS_GammaConv = 1; pass_cat = 1;}
    }

    //-- ttW CR (ttH)
    if(nJets_ttH == 3 && pass_njet_tth)
    {
        is_ttH_ttWctrl = 1;
        if(pass_MCMatching && is_SS_pair && pass_isTight) {is_ttH_ttWctrl_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && is_SS_pair && !pass_isTight) {is_ttH_ttWctrl_Fake = 1; pass_cat = 1;}
        else if(pass_MCMatching && !is_SS_pair && contains_ele && pass_isTight) {is_ttH_ttWctrl_Flip = 1; pass_cat = 1; Apply_FlipRate_Weight();}
        else if(is_GammaConv && is_SS_pair && pass_isTight) {is_ttH_ttWctrl_GammaConv = 1; pass_cat = 1;}
    }

    //ttW-ttbar CR (as in tHq2016 analysis)
    if(pass_njet_ttbar)
    {
        is_tHq_ttWctrl = 1;
        if(pass_MCMatching && is_SS_pair && pass_isTight) {is_tHq_ttWctrl_SR = 1; pass_cat = 1;}
        else if(pass_MCMatching && is_SS_pair && !pass_isTight) {is_tHq_ttWctrl_Fake = 1; pass_cat = 1;}
        else if(pass_MCMatching && !is_SS_pair && contains_ele && pass_isTight) {is_tHq_ttWctrl_Flip = 1; pass_cat = 1; Apply_FlipRate_Weight();}
        else if(is_GammaConv && is_SS_pair && pass_isTight) {is_tHq_ttWctrl_GammaConv = 1; pass_cat = 1;}
    }

//--------------------------------------------
	return pass_cat;
}












//--------------------------------------------
 // #
 // #    #        #
 // #    #        #
 // #    #        #
 // #######       #
 //      #        #
 //      #        ######
//--------------------------------------------

/**
 * 4l selection
 */
bool tHqMultileponAnalysis::FourLeptonSelection_THQ4l(int evt)
{
	if(DEBUG) {cout<<FYEL("-- FourLeptonSelection_THQ4l")<<endl;}

    InitVariables();

    //Must pass single/double/tri-lepton trigger
    if(weight==0 || !is_trigger_3l) {return 0;}

    //At least 4 FO leptons
    if(vLeptonFakeable.size() < 4) {return 0;}
    if(DEBUG) {cout<<"- Passed nLep cut"<<endl;}

    //pT cuts
    if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15 || vLeptonFakeable.at(2).conept < 15 || vLeptonFakeable.at(3).conept < 10) {return 0;}
    if(DEBUG) {cout<<"- Passed pT cut"<<endl;}

    //No pair of loose leptons with mll < 12
    if(!pass_cleanup) {return 0;}
    if(DEBUG) {cout<<"- Passed cleanup cut"<<endl;}

    if(nJets_ttH < 4 && ((nSFOS > 0 && metLD < 45) || (!nSFOS && metLD < 30)) ) {return 0;}
    if(DEBUG) {cout<<"- Passed metLD cut"<<endl;}

    if(vLeptonFakeable.at(0).charge + vLeptonFakeable.at(1).charge + vLeptonFakeable.at(2).charge + vLeptonFakeable.at(3).charge != 0) {return 0;}
    if(DEBUG) {cout<<"- Passed charge cut"<<endl;}

    if(mllll > 0 && mllll < 140) {return 0;}
    if(DEBUG) {cout<<"- Passed mllll cut"<<endl;}

    if(nJets_ttH < 2) {return 0;}
    if(DEBUG) {cout<<"- Passed first jet cut"<<endl;}
//--------------------------------------------

    if(!FourLeptonSelection_THQ4l_Regions(evt) ) {return 0;} //Check categories

    return 1;
}





/**
 * 4l categories selection
 */
bool tHqMultileponAnalysis::FourLeptonSelection_THQ4l_Regions(int evt)
{
    bool pass_cat = false;

    //Not shared among all categories
    bool pass_isTight = false;
    bool pass_MCMatching = false;
    bool pass_njet = false;
    if(vLeptonFakeable.at(0).isTightTTH && vLeptonFakeable.at(1).isTightTTH && vLeptonFakeable.at(2).isTightTTH && vLeptonFakeable.at(3).isTightTTH) {pass_isTight = true;}
    if(_isdata || (vLeptonFakeable.at(0).hasMCMatch && vLeptonFakeable.at(1).hasMCMatch && vLeptonFakeable.at(2).hasMCMatch && vLeptonFakeable.at(3).hasMCMatch) ) {pass_MCMatching = true;} //In 2lSS, also require charge matching
    if(nLooseBJets >= 2 || nMediumBJets >= 1) {pass_njet = true;}

    if(pass_njet && pass_Zveto && pass_isTight && pass_MCMatching) {is_tHq_4l_SR = true; pass_cat = 1;}

    if(!pass_njet && !pass_Zveto && pass_isTight && pass_MCMatching) {is_tHq_ZZctrl_SR = true; pass_cat = 1;}

//--------------------------------------------

    return pass_cat;
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
 * From : https://github.com/peruzzim/cmgtools-lite/blob/94X_dev_ttH/TTHAnalysis/macros/finalMVA/trainKinMVA.py#L118
 */
bool tHqMultileponAnalysis::ThreeLeptonSelection_THQ3l_TrainingSelection(int evt)
{
    if(DEBUG) {cout<<FYEL("-- ThreeLeptonSelection_THQ3l_TrainingSelection")<<endl;}

    InitVariables();
    // v_cutflow[10]++;

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return 0;}
    // v_cutflow[0]++;

	//At least 3 FO leptons
    if(vLeptonFakeable.size() < 3) {return 0;}
    // v_cutflow[1]++;

    //MC matching //REMOVED -- kills too much ttbar ?
    // if(!_isdata && (!vLeptonFakeable.at(0).hasMCMatch || !vLeptonFakeable.at(1).hasMCMatch || !vLeptonFakeable.at(2).hasMCMatch) ) {return 0;}

    if(!pass_Zveto) {return 0;}
    // v_cutflow[2]++;

    //metLD cut
    if(nJets_ttH < 4)
    {
        if(nSFOS > 0 && metLD < 45) {return 0;}
        else if(!nSFOS && metLD < 30) {return 0;}
    }

    //tHq training sel
    {
        is_tHq_3l_Training = 1;

        if(vLeptonFakeable.at(0).conept < 20 || vLeptonFakeable.at(1).conept < 10 || vLeptonFakeable.at(2).conept < 10) {is_tHq_3l_Training = 0;}
        // v_cutflow[3]++;

        // if(nLooseBJets < 2) {is_tHq_3l_Training = 0;}
        // if(!nLooseBJets || !nLightJets_tHq) {is_tHq_3l_Training = 0;}
        if(nLooseBJets < 2 && !nMediumBJets) {is_tHq_3l_Training = 0;}
    }

    // v_cutflow[4]++;

    //ttH training sel
    {
        is_ttH_3l_Training = 1;

        if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15 || vLeptonFakeable.at(2).conept < 10) {is_ttH_3l_Training = 0;} //NB -- says different in training script & preApp talk !
        if(nJets_ttH < 4)
        {
            if(nSFOS > 0 && metLD < 45) {is_ttH_3l_Training = 0;}
            else if(!nSFOS && metLD < 30) {is_ttH_3l_Training = 0;}
        }

        if(nLooseBJets < 2) {is_ttH_3l_Training = 0;}
    }

    //--------------------------------------------
    return 1;
}


/**
 *  * Selection function implementing the ttH2017 training selection, for 2lSS cat.
 *  From : https://github.com/peruzzim/cmgtools-lite/blob/94X_dev_ttH/TTHAnalysis/macros/finalMVA/trainKinMVA.py#L118
 */
bool tHqMultileponAnalysis::TwoLeptonSelection_THQ2lSS_TrainingSelection(int evt)
{
    if(DEBUG) {cout<<FYEL("-- TwoLeptonSelection_THQ2lSS_TrainingSelection")<<endl;}

    InitVariables();

    if(weight==0) {return 0;}

    //At least 2 FO leptons
    if(vLeptonFakeable.size() < 2) {return 0;}

    //MC matching //REMOVED -- kills too much ttbar ?
    // if(!_isdata && (!vLeptonFakeable.at(0).hasMCMatch || !vLeptonFakeable.at(1).hasMCMatch) ) {return 0;}

    //SS pair
    if(vLeptonFakeable.at(0).charge != vLeptonFakeable.at(1).charge) {return 0;}

    //No ee pair with mll-mZ < 10 //CHANGED -- removed ?
    // if(!pass_Zveto_ee) {return 0;}

    //tHq training sel -- Try looser cuts
    {
        is_tHq_2lSS_Training = 1;
        if(vLeptonFakeable.at(0).conept < 20 || vLeptonFakeable.at(1).conept < 10) {is_tHq_2lSS_Training = 0;}
        // if(nLooseBJets < 2 && nMediumBJets == 0) {is_tHq_2lSS_Training = 0;}
        // if(!nLooseBJets || !nLightJets_tHq) {is_tHq_2lSS_Training = 0;}
        if(nLooseBJets < 2 && !nMediumBJets) {is_tHq_2lSS_Training = 0;}
    }

    //ttH training sel
    {
        is_ttH_2lSS_Training = 1;
        if(vLeptonFakeable.at(0).conept < 25 || vLeptonFakeable.at(1).conept < 15) {is_ttH_2lSS_Training = 0;}
        if(nJets_ttH < 4) {is_ttH_2lSS_Training = 0;}
        if(nLooseBJets < 2 && nMediumBJets == 0) {is_ttH_2lSS_Training = 0;}
    }

    //--------------------------------------------
    return 1;
}





























//--------------------------------------------
// ########    ###    ##    ## ######## ########     ###    ######## ########
// ##         ## ##   ##   ##  ##       ##     ##   ## ##      ##    ##
// ##        ##   ##  ##  ##   ##       ##     ##  ##   ##     ##    ##
// ######   ##     ## #####    ######   ########  ##     ##    ##    ######
// ##       ######### ##  ##   ##       ##   ##   #########    ##    ##
// ##       ##     ## ##   ##  ##       ##    ##  ##     ##    ##    ##
// ##       ##     ## ##    ## ######## ##     ## ##     ##    ##    ########
//--------------------------------------------

/**
 * Set the value of "weightfake" variable. Weight read from ttH FR file
 * NB : nominal and variations stored separately
 */
void tHqMultileponAnalysis::Apply_FakeRate_Weight()
{
    weightfake = 1.; //nominal FR weight
    std::fill(v_floats_FR_variations.begin(), v_floats_FR_variations.end(), 1); //reset all FR variations to 1

//--------------------------------------------

    int nlep = 0;
    if(is_tHq_2lSS_Fake || is_tHq_ttWctrl_Fake || is_ttH_2lSS_Fake || is_ttH_ttWctrl_Fake) {nlep = 2;}
    else if(is_tHq_3l_Fake || is_tHq_ttZctrl_Fake || is_tHq_WZctrl_Fake || is_ttH_3l_Fake || is_ttH_ttZctrl_Fake || is_ttH_WZctrl_Fake) {nlep = 3;}
    if(!nlep) {return;}

    //If channel is ele-only or muon-only, can set some variations to nominal directly
    bool isChannel_withoutEl = false;
    bool isChannel_withoutMu = false;
    if(channel == 0) {isChannel_withoutEl = true;}
    if((nlep == 2 && channel == 2) || (nlep == 3 && channel == 3) ) {isChannel_withoutMu = true;}

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;
    //Consider only 2 or 3 hardest FO leptons
    for(int i=0; i<nlep; i++)
    {
        if( !vLeptonFakeable.at(i).isTightTTH ) //For each lepton failing the tight requirements, multiply event by a weight
        {
            leptonsPts.push_back(vLeptonFakeable.at(i).conept);
            leptonsEtas.push_back(vLeptonFakeable.at(i).eta);
            leptonsIds.push_back(vLeptonFakeable.at(i).id);
        }
    }

    //Store nominal FR weight in dedicated variable
    weightfake = Get_FR_Weight(leptonsPts, leptonsEtas, leptonsIds, "");

    for(int ivar=0; ivar<v_FR_type.size(); ivar++)
    {
        v_floats_FR_variations[ivar] = Get_FR_Weight(leptonsPts, leptonsEtas, leptonsIds, v_FR_type[ivar], isChannel_withoutEl, isChannel_withoutMu);
        if(!v_floats_FR_variations[ivar]) {v_floats_FR_variations[ivar] = weightfake;} //default case  set to nominal value
    }

    // cout<<"weightfake = "<<weightfake<<endl;

    return;
}




//--------------------------------------------
// ######## ##       #### ########  ########     ###    ######## ########
// ##       ##        ##  ##     ## ##     ##   ## ##      ##    ##
// ##       ##        ##  ##     ## ##     ##  ##   ##     ##    ##
// ######   ##        ##  ########  ########  ##     ##    ##    ######
// ##       ##        ##  ##        ##   ##   #########    ##    ##
// ##       ##        ##  ##        ##    ##  ##     ##    ##    ##
// ##       ######## #### ##        ##     ## ##     ##    ##    ########
//--------------------------------------------

/**
 * Set the value of "weightflip" variable. Weight read from ttH FR file
 */
void tHqMultileponAnalysis::Apply_FlipRate_Weight()
{
    weightflip = 1;

    weightflip = get_QF_weight(vLeptonFakeable.at(0).conept, fabs(vLeptonFakeable.at(0).eta), vLeptonFakeable.at(0).id, vLeptonFakeable.at(1).conept, fabs(vLeptonFakeable.at(1).eta), vLeptonFakeable.at(1).id);

    // cout<<"weightflip = "<<weightflip<<endl;

    return;
}







//--------------------------------------------
//  ######   ######     ###    ##       ########    ########    ###     ######  ########  #######  ########   ######
// ##    ## ##    ##   ## ##   ##       ##          ##         ## ##   ##    ##    ##    ##     ## ##     ## ##    ##
// ##       ##        ##   ##  ##       ##          ##        ##   ##  ##          ##    ##     ## ##     ## ##
//  ######  ##       ##     ## ##       ######      ######   ##     ## ##          ##    ##     ## ########   ######
//       ## ##       ######### ##       ##          ##       ######### ##          ##    ##     ## ##   ##         ##
// ##    ## ##    ## ##     ## ##       ##          ##       ##     ## ##    ##    ##    ##     ## ##    ##  ##    ##
//  ######   ######  ##     ## ######## ########    ##       ##     ##  ######     ##     #######  ##     ##  ######
//--------------------------------------------




/**
 * Apply scale factors (lepton, trigger, ...) to weight, using the ScaleFactors object
 * @param nlep : read different files depending if 2lss or 3l
 */
void tHqMultileponAnalysis::Apply_ScaleFactors(int nlep, int itree)
{
    if(_isdata) {return;}

    if(nlep != 2 && nlep != 3) {cout<<"Error : wrong nlep value !"<<endl; return;}

    //Init SFs
    lepton_SF  = 1;
    trigger_SF = 1;
    btag_SF    = 1;
    PU_SF      = 1;
    total_SF   = 1;

    //Init systematics weights with event MC weight
    TrigEffUp           = weight; TrigEffDown        = weight;
    LepEff_elLooseUp    = weight; LepEff_elLooseDown = weight;
    LepEff_muLooseUp    = weight; LepEff_muLooseDown = weight;
    LepEff_elTightUp    = weight; LepEff_elTightDown = weight;
    // LepEff_muTightUp = weight; LepEff_muTightDown = weight; //Use flat syst for muTight uncert.
    LFcontUp            = weight; LFcontDown         = weight;
    HFstats1Up          = weight; HFstats1Down       = weight;
    HFstats2Up          = weight; HFstats2Down       = weight;
    CFerr1Up            = weight; CFerr1Down         = weight;
    CFerr2Up            = weight; CFerr2Down         = weight;
    HFcontUp            = weight; HFcontDown         = weight;
    LFstats1Up          = weight; LFstats1Down       = weight;
    LFstats2Up          = weight; LFstats2Down       = weight;
    PUUp                = weight; PUDown             = weight;
    QCDscaleUp          = weight; QCDscaleDown       = weight;
    pdfUp               = weight; pdfDown            = weight;


    //Lepton SF (split)
    //--------------------------------------------
    for(int ilep=0; ilep<nlep; ilep++)
    {
        lepton_SF*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "");

        if(abs(vLeptonFakeable[ilep].id) == 11) //Ele
        {
            LepEff_elLooseUp*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "looseUp");
            LepEff_elLooseDown*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "looseDown");
            LepEff_elTightUp*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "tightUp");
            LepEff_elTightDown*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "tightDown");
        }
        else if(abs(vLeptonFakeable[ilep].id) == 13) //Mu
        {
            LepEff_muLooseUp*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "looseUp");
            LepEff_muLooseDown*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "looseDown");

            //Use flat syst for muTight uncert.
            // LepEff_muTightUp*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "tightUp");
            // LepEff_muTightDown*= sf->Get_Lepton_SF(nlep, vLeptonFakeable[ilep].id, vLeptonFakeable[ilep].conept, vLeptonFakeable[ilep].eta, "tightDown");
        }
    }

    // cout<<"lepton_SF = "<<lepton_SF<<endl;
    // cout<<"Weight : "<<weight<<" --> "<<weight * lepton_SF<<endl;

    // weight*= lepton_SF;
    //--------------------------------------------


    //Trigger SF
    //--------------------------------------------
    trigger_SF = sf->Get_Trigger_SF(nlep, vLeptonFakeable[0].id, vLeptonFakeable[0].conept, vLeptonFakeable[1].id, 0);

    TrigEffUp*= sf->Get_Trigger_SF(nlep, vLeptonFakeable[0].id, vLeptonFakeable[0].conept, vLeptonFakeable[1].id, +1);
    TrigEffDown*= sf->Get_Trigger_SF(nlep, vLeptonFakeable[0].id, vLeptonFakeable[0].conept, vLeptonFakeable[1].id, -1);

    // cout<<"trigger_SF = "<<trigger_SF<<endl;
    // cout<<"Weight : "<<weight<<" --> "<<weight * trigger_SF<<endl;

    // weight*= trigger_SF;
    //--------------------------------------------


   //Btag SF
   //--------------------------------------------
   //Using the BTagShapeCalibration method : https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
   //NB : assign flat scale factor of 1 (no correction) to charm jets
    for(int ijet=0; ijet<vJetLoose.size(); ijet++)
    {
        if(v_systTree[itree] == "JESUp")
        {
            btag_SF*= sf->Get_Btag_SF("up_jes", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
        }
        else if(v_systTree[itree] == "JESDown")
        {
            btag_SF*= sf->Get_Btag_SF("down_jes", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
        }
        else
        {
            btag_SF*= sf->Get_Btag_SF("central", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
        }

        if(v_systTree[itree] == "") //Only compute syst weights for nominal TTree
        {
            LFcontUp*= sf->Get_Btag_SF("up_lf", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            HFstats1Up*= sf->Get_Btag_SF("up_hfstats1", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            HFstats2Up*= sf->Get_Btag_SF("up_hfstats2", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            CFerr1Up*= sf->Get_Btag_SF("up_cferr1", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            CFerr2Up*= sf->Get_Btag_SF("up_cferr2", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            HFcontUp*= sf->Get_Btag_SF("up_hf", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            LFstats1Up*= sf->Get_Btag_SF("up_lfstats1", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            LFstats2Up*= sf->Get_Btag_SF("up_lfstats2", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);

            LFcontDown*= sf->Get_Btag_SF("down_lf", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            HFstats1Down*= sf->Get_Btag_SF("down_hfstats1", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            HFstats2Down*= sf->Get_Btag_SF("down_hfstats2", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            CFerr1Down*= sf->Get_Btag_SF("down_cferr1", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            CFerr2Down*= sf->Get_Btag_SF("down_cferr2", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            HFcontDown*= sf->Get_Btag_SF("down_hf", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            LFstats1Down*= sf->Get_Btag_SF("down_lfstats1", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
            LFstats2Down*= sf->Get_Btag_SF("down_lfstats2", abs(vJetLoose.at(ijet).jet_hadronFlavour), vJetLoose.at(ijet).pt, fabs(vJetLoose.at(ijet).eta), vJetLoose.at(ijet).DeepCSVbtag);
        }
    }

    // cout<<"btag_SF = "<<btag_SF<<endl;

    // weight*= btag_SF;
    //--------------------------------------------


    //PU SF
    //--------------------------------------------
    //The 'ScaleFactors' class possesses an object of the 'Pileup' class, used to get PU SF
    PU_SF = sf->Get_Pileup_SF(nPU, "nom");
    PUUp*= sf->Get_Pileup_SF(nPU, "up");
    PUDown*= sf->Get_Pileup_SF(nPU, "down");

    // cout<<"nPU = "<<nPU<<" => PU_SF = "<<PU_SF<<endl;
    // cout<<"PU_SF = "<<PU_SF<<" / PUUp = "<<PUUp<<" / PUDown = "<<PUDown<<endl;

    // weight*= PU_SF;
    //--------------------------------------------

    //Scale SF
    //--------------------------------------------
    Compute_Weight_ScaleSyst();
    //--------------------------------------------


    //Total SF
    //--------------------------------------------
    if(make_ntuples_for_overlap_studies) {total_SF = lepton_SF * trigger_SF * btag_SF;} //don't apply for now //FIXME
    // else {total_SF = lepton_SF * trigger_SF * btag_SF * PU_SF;}
    else {total_SF = lepton_SF * trigger_SF * btag_SF;} //FIXME

    weight*= total_SF; //Note that the "default" weight now contains the SF corrections
    weight_old*= total_SF;

    //Need to also multiply syst weights with total SF, so they are consistent with nominal weight
    //But e.g. need to divide triggerUp by trigger_SF, since we want to multiply it with upper variation, not nominal
    TrigEffUp*= total_SF/trigger_SF; TrigEffDown*= total_SF/trigger_SF;
    // EleEffUp*= total_SF/lepton_SF; EleEffDown*= total_SF/lepton_SF;
    // MuEffUp*= total_SF/lepton_SF; MuEffDown*= total_SF/lepton_SF;
    LepEff_elLooseUp*= total_SF/lepton_SF; LepEff_elLooseDown*= total_SF/lepton_SF;
    LepEff_muLooseUp*= total_SF/lepton_SF; LepEff_muLooseDown*= total_SF/lepton_SF;
    LepEff_elTightUp*= total_SF/lepton_SF; LepEff_elTightDown*= total_SF/lepton_SF;
    PUUp*= total_SF/PU_SF; PUDown*= total_SF/PU_SF;
    LFcontUp*= total_SF/btag_SF; LFcontDown*= total_SF/btag_SF;
    HFstats1Up*= total_SF/btag_SF; HFstats1Down*= total_SF/btag_SF;
    HFstats2Up*= total_SF/btag_SF; HFstats2Down*= total_SF/btag_SF;
    CFerr1Up*= total_SF/btag_SF; CFerr1Down*= total_SF/btag_SF;
    CFerr2Up*= total_SF/btag_SF; CFerr2Down*= total_SF/btag_SF;
    HFcontUp*= total_SF/btag_SF; HFcontDown*= total_SF/btag_SF;
    LFstats1Up*= total_SF/btag_SF; LFstats1Down*= total_SF/btag_SF;
    LFstats2Up*= total_SF/btag_SF; LFstats2Down*= total_SF/btag_SF;
    QCDscaleUp*= total_SF; QCDscaleDown*= total_SF;
    pdfUp*= total_SF; pdfDown*= total_SF;

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

//--------------------------------------------





/**
* Fill output tree. Compute lot of variables, needed either for MVA analysis of for MEM computation
*/
void tHqMultileponAnalysis::fillOutputTree(int itree)
{
    InitTLorentzVectors(); //Re-init the MEM inputs at each call of the function

    //--------------------------------------------
    //      # ###### #####  ####
    //      # #        #   #
    //      # #####    #    ####
    //      # #        #        #
    // #    # #        #   #    #
    //  ####  ######   #    ####
    //--------------------------------------------

    JetsPt = new vector<Float_t>;
    JetsEta = new vector<Float_t>;
    JetsPhi = new vector<Float_t>;
    JetsE = new vector<Float_t>;
    JetsCSV = new vector<Float_t>;

    //NEW -- store full jets collection, so MEM code can they select them automatically
    for(int ijet=0; ijet<nJets; ijet++)
    {
        JetsPt->push_back(vJetLoose.at(ijet).pt);
        JetsEta->push_back(vJetLoose.at(ijet).eta);
        JetsPhi->push_back(vJetLoose.at(ijet).phi);
        JetsE->push_back(vJetLoose.at(ijet).E);
        JetsCSV->push_back(vJetLoose.at(ijet).DeepCSVbtag);
    }


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
        // Bjet1.SetPtEtaPhiE(vJetLoose.at(ib1).pt, vJetLoose.at(ib1).eta, vJetLoose.at(ib1).phi, vJetLoose.at(ib1).E);
        // FillJetInfoOutputTree(&multilepton_Bjet1_Id, 5, &multilepton_Bjet1_P4, Bjet1, &multilepton_Bjet1_CSV, vJetLoose.at(ib1).CSVv2, &multilepton_Bjet1_JEC_Up, &multilepton_Bjet1_JEC_Down, vJetLoose.at(ib1).JES_uncert(), &multilepton_Bjet1_JER_Up, &multilepton_Bjet1_JER_Down, vJetLoose.at(ib1).pt_JER(), vJetLoose.at(ib1).pt_JER_up(), vJetLoose.at(ib1).pt_JER_down());
        multilepton_Bjet1_Id = 5;
        multilepton_Bjet1_CSV = vJetLoose.at(ib1).DeepCSVbtag;
        multilepton_Bjet1_P4.SetPtEtaPhiE(vJetLoose.at(ib1).pt, vJetLoose.at(ib1).eta, vJetLoose.at(ib1).phi, vJetLoose.at(ib1).E);

    }
    if (ib2!=-1)
    {
        // Bjet2.SetPtEtaPhiE(vJetLoose.at(ib2).pt, vJetLoose.at(ib2).eta, vJetLoose.at(ib2).phi, vJetLoose.at(ib2).E);
        // FillJetInfoOutputTree(&multilepton_Bjet2_Id, 5, &multilepton_Bjet2_P4, Bjet2, &multilepton_Bjet2_CSV, vJetLoose.at(ib2).CSVv2, &multilepton_Bjet2_JEC_Up, &multilepton_Bjet2_JEC_Down, vJetLoose.at(ib2).JES_uncert(), &multilepton_Bjet2_JER_Up, &multilepton_Bjet2_JER_Down, vJetLoose.at(ib2).pt_JER(), vJetLoose.at(ib2).pt_JER_up(), vJetLoose.at(ib2).pt_JER_down());
        multilepton_Bjet2_Id = 5;
        multilepton_Bjet2_CSV = vJetLoose.at(ib2).DeepCSVbtag;
        multilepton_Bjet2_P4.SetPtEtaPhiE(vJetLoose.at(ib2).pt, vJetLoose.at(ib2).eta, vJetLoose.at(ib2).phi, vJetLoose.at(ib2).E);
    }


    //Define jet category
    //2lss
    if (is_tHq_2lSS && ib1!=-1 && ib2!=-1 && vJetLoose.size()-2>=4) catJets = kCat_2lss_2b_4j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2==-1 && vJetLoose.size()-1>=4) catJets = kCat_2lss_1b_4j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2!=-1 && vJetLoose.size()-2==3) catJets = kCat_2lss_2b_3j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2==-1 && vJetLoose.size()-1==3) catJets = kCat_2lss_1b_3j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2!=-1 && vJetLoose.size()-2==2) catJets = kCat_2lss_2b_2j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2!=-1 && vJetLoose.size()-2==1) catJets = kCat_2lss_2b_1j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2==-1 && vJetLoose.size()-1==2) catJets = kCat_2lss_1b_2j;
    else if (is_tHq_2lSS && ib1!=-1 && ib2==-1 && vJetLoose.size()-1==1) catJets = kCat_2lss_1b_1j;
    //4l
    // else if (is4l && ib1!=-1 && ib2!=-1) catJets = kCat_4l_2b;
    // else if (is4l && ib1!=-1 && ib2==-1) catJets = kCat_4l_1b;
    //3l
    else if (is_tHq_3l && ib1!=-1 && ib2!=-1 && vJetLoose.size()-2>=2) catJets = kCat_3l_2b_2j;
    else if (is_tHq_3l && ib1!=-1 && ib2==-1 && vJetLoose.size()-1>=2) catJets = kCat_3l_1b_2j;
    else if (is_tHq_3l && ib1!=-1 && ib2!=-1 && vJetLoose.size()-2==1) catJets = kCat_3l_2b_1j;
    else if (is_tHq_3l && ib1!=-1 && ib2==-1 && vJetLoose.size()-1==1) catJets = kCat_3l_1b_1j;
    else if (is_tHq_3l && ib1!=-1 && ib2!=-1 && vJetLoose.size()-2==0) catJets = kCat_3l_2b_0j;
    else if (is_tHq_3l && ib1!=-1 && ib2==-1 && vJetLoose.size()-1==0) catJets = kCat_3l_1b_0j;
    else if (is_tHq_3l && ib1==-1 && ib2==-1 && vJetLoose.size()==1) catJets = kCat_3l_0b_1j;
    else if (is_tHq_3l && ib1==-1 && ib2==-1 && vJetLoose.size()==0) catJets = kCat_3l_0b_0j;

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
                    jet1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vJetLoose.at(ib1).eta, vJetLoose.at(ib1).phi );
                    if( jet1_dr_gen < jet1_dr_gen_min) {jet1_dr_gen_min = jet1_dr_gen;  jet1_matched = itruth;}
                }

                if(ib2!=-1)
                {
                    jet2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta.at(itruth),  vTruth->at(0).mc_truth_phi.at(itruth), vJetLoose.at(ib2).eta, vJetLoose.at(ib2).phi );
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
        // Jet1.SetPtEtaPhiE(vJetLoose.at(ij1).pt, vJetLoose.at(ij1).eta, vJetLoose.at(ij1).phi, vJetLoose.at(ij1).E);
        // FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vJetLoose.at(ij1).CSVv2, &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vJetLoose.at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vJetLoose.at(ij1).pt_JER(), vJetLoose.at(ij1).pt_JER_up(), vJetLoose.at(ij1).pt_JER_down());
        multilepton_JetHighestPt1_Id = 1;
        multilepton_JetHighestPt1_CSV = vJetLoose.at(ij1).DeepCSVbtag;
        multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vJetLoose.at(ij1).pt, vJetLoose.at(ij1).eta, vJetLoose.at(ij1).phi, vJetLoose.at(ij1).E);
    }
    if (ij2!=-1)
    {
        // Jet2.SetPtEtaPhiE(vJetLoose.at(ij2).pt, vJetLoose.at(ij2).eta, vJetLoose.at(ij2).phi, vJetLoose.at(ij2).E);
        // FillJetInfoOutputTree(&multilepton_JetHighestPt2_Id, 1, &multilepton_JetHighestPt2_P4, Jet2, &multilepton_JetHighestPt2_CSV, vJetLoose.at(ij2).CSVv2, &multilepton_JetHighestPt2_JEC_Up, &multilepton_JetHighestPt2_JEC_Down, vJetLoose.at(ij2).JES_uncert(), &multilepton_JetHighestPt2_JER_Up, &multilepton_JetHighestPt2_JER_Down, vJetLoose.at(ij2).pt_JER(), vJetLoose.at(ij2).pt_JER_up(), vJetLoose.at(ij2).pt_JER_down());
        multilepton_JetHighestPt2_Id = 1;
        multilepton_JetHighestPt2_CSV = vJetLoose.at(ij2).DeepCSVbtag;
        multilepton_JetHighestPt2_P4.SetPtEtaPhiE(vJetLoose.at(ij2).pt, vJetLoose.at(ij2).eta, vJetLoose.at(ij2).phi, vJetLoose.at(ij2).E);
    }

    if (ik1!=-1 && ik2!=-1){
        // Jet1.SetPtEtaPhiE(vJetLoose.at(ik1).pt, vJetLoose.at(ik1).eta, vJetLoose.at(ik1).phi, vJetLoose.at(ik1).E);
        // Jet2.SetPtEtaPhiE(vJetLoose.at(ik2).pt, vJetLoose.at(ik2).eta, vJetLoose.at(ik2).phi, vJetLoose.at(ik2).E);
        // FillJetInfoOutputTree(&multilepton_JetClosestMw1_Id, 2, &multilepton_JetClosestMw1_P4, Jet1, &multilepton_JetClosestMw1_CSV, vJetLoose.at(ik1).CSVv2, &multilepton_JetClosestMw1_JEC_Up, &multilepton_JetClosestMw1_JEC_Down, vJetLoose.at(ik1).JES_uncert(), &multilepton_JetClosestMw1_JER_Up, &multilepton_JetClosestMw1_JER_Down, vJetLoose.at(ik1).pt_JER(), vJetLoose.at(ik1).pt_JER_up(), vJetLoose.at(ik1).pt_JER_down());
        // FillJetInfoOutputTree(&multilepton_JetClosestMw2_Id, 2, &multilepton_JetClosestMw2_P4, Jet2, &multilepton_JetClosestMw2_CSV, vJetLoose.at(ik2).CSVv2, &multilepton_JetClosestMw2_JEC_Up, &multilepton_JetClosestMw2_JEC_Down, vJetLoose.at(ik2).JES_uncert(), &multilepton_JetClosestMw2_JER_Up, &multilepton_JetClosestMw2_JER_Down, vJetLoose.at(ik2).pt_JER(), vJetLoose.at(ik2).pt_JER_up(), vJetLoose.at(ik2).pt_JER_down());
        multilepton_JetClosestMw1_Id = 2;
        multilepton_JetClosestMw2_Id = 2;
        multilepton_JetClosestMw1_CSV = vJetLoose.at(ik1).DeepCSVbtag;
        multilepton_JetClosestMw2_CSV = vJetLoose.at(ik2).DeepCSVbtag;
        multilepton_JetClosestMw1_P4.SetPtEtaPhiE(vJetLoose.at(ik1).pt, vJetLoose.at(ik1).eta, vJetLoose.at(ik1).phi, vJetLoose.at(ik1).E);
        multilepton_JetClosestMw2_P4.SetPtEtaPhiE(vJetLoose.at(ik2).pt, vJetLoose.at(ik2).eta, vJetLoose.at(ik2).phi, vJetLoose.at(ik2).E);
    }
    if (il1!=-1 && il2!=-1){
        // Jet1.SetPtEtaPhiE(vJetLoose.at(il1).pt, vJetLoose.at(il1).eta, vJetLoose.at(il1).phi, vJetLoose.at(il1).E);
        // Jet2.SetPtEtaPhiE(vJetLoose.at(il2).pt, vJetLoose.at(il2).eta, vJetLoose.at(il2).phi, vJetLoose.at(il2).E);
        // FillJetInfoOutputTree(&multilepton_JetLowestMjj1_Id, 3, &multilepton_JetLowestMjj1_P4, Jet1, &multilepton_JetLowestMjj1_CSV, vJetLoose.at(il1).CSVv2, &multilepton_JetLowestMjj1_JEC_Up, &multilepton_JetLowestMjj1_JEC_Down, vJetLoose.at(il1).JES_uncert(), &multilepton_JetLowestMjj1_JER_Up, &multilepton_JetLowestMjj1_JER_Down, vJetLoose.at(il1).pt_JER(), vJetLoose.at(il1).pt_JER_up(), vJetLoose.at(il1).pt_JER_down());
        // FillJetInfoOutputTree(&multilepton_JetLowestMjj2_Id, 3, &multilepton_JetLowestMjj2_P4, Jet2, &multilepton_JetLowestMjj2_CSV, vJetLoose.at(il2).CSVv2, &multilepton_JetLowestMjj2_JEC_Up, &multilepton_JetLowestMjj2_JEC_Down, vJetLoose.at(il2).JES_uncert(), &multilepton_JetLowestMjj2_JER_Up, &multilepton_JetLowestMjj2_JER_Down, vJetLoose.at(il2).pt_JER(), vJetLoose.at(il2).pt_JER_up(), vJetLoose.at(il2).pt_JER_down());
        multilepton_JetLowestMjj1_Id = 3;
        multilepton_JetLowestMjj2_Id = 3;
        multilepton_JetLowestMjj1_CSV = vJetLoose.at(il1).DeepCSVbtag;
        multilepton_JetLowestMjj2_CSV = vJetLoose.at(il2).DeepCSVbtag;
        multilepton_JetLowestMjj1_P4.SetPtEtaPhiE(vJetLoose.at(il1).pt, vJetLoose.at(il1).eta, vJetLoose.at(il1).phi, vJetLoose.at(il1).E);
        multilepton_JetLowestMjj2_P4.SetPtEtaPhiE(vJetLoose.at(il2).pt, vJetLoose.at(il2).eta, vJetLoose.at(il2).phi, vJetLoose.at(il2).E);
    }
    if(ie1!=-1)
    {
        multilepton_JetHighestEta1_Id = 4;
        multilepton_JetHighestEta1_CSV = vJetLoose.at(ie1).DeepCSVbtag;
        multilepton_JetHighestEta1_P4.SetPtEtaPhiE(vJetLoose.at(ie1).pt, vJetLoose.at(ie1).eta, vJetLoose.at(ie1).phi, vJetLoose.at(ie1).E );
    }
    if(ie2!=-1) //2jets
    {
        multilepton_JetHighestEta2_Id = 4;
        multilepton_JetHighestEta2_CSV = vJetLoose.at(ie2).DeepCSVbtag;
        multilepton_JetHighestEta2_P4.SetPtEtaPhiE(vJetLoose.at(ie2).pt, vJetLoose.at(ie2).eta, vJetLoose.at(ie2).phi, vJetLoose.at(ie2).E );
    }


    //--- Fill 2nd pairs (first one is closest to mW) -- needed for 2l only (more jets)
    if(is_tHq_2lSS && ij1!=-1 && ij2!=-1)
    {
        if (im1!=-1)
        {
            // Jet1.SetPtEtaPhiE(vJetLoose.at(im1).pt, vJetLoose.at(im1).eta, vJetLoose.at(im1).phi, vJetLoose.at(im1).E);
            // FillJetInfoOutputTree(&multilepton_JetHighestPt1_2ndPair_Id, 1, &multilepton_JetHighestPt1_2ndPair_P4, Jet1, &multilepton_JetHighestPt1_2ndPair_CSV, vJetLoose.at(im1).CSVv2, &multilepton_JetHighestPt1_2ndPair_JEC_Up, &multilepton_JetHighestPt1_2ndPair_JEC_Down, vJetLoose.at(im1).JES_uncert(), &multilepton_JetHighestPt1_2ndPair_JER_Up, &multilepton_JetHighestPt1_2ndPair_JER_Down, vJetLoose.at(im1).pt_JER(), vJetLoose.at(im1).pt_JER_up(), vJetLoose.at(im1).pt_JER_down());
            multilepton_JetHighestPt1_2ndPair_Id = 1;
            multilepton_JetHighestPt1_2ndPair_CSV = vJetLoose.at(im1).DeepCSVbtag;
            multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vJetLoose.at(im1).pt, vJetLoose.at(im1).eta, vJetLoose.at(im1).phi, vJetLoose.at(im1).E);
        }
        if(im2!=-1){
            // Jet2.SetPtEtaPhiE(vJetLoose.at(im2).pt, vJetLoose.at(im2).eta, vJetLoose.at(im2).phi, vJetLoose.at(im2).E);
            // FillJetInfoOutputTree(&multilepton_JetHighestPt2_2ndPair_Id, 1, &multilepton_JetHighestPt2_2ndPair_P4, Jet2, &multilepton_JetHighestPt2_2ndPair_CSV, vJetLoose.at(im2).CSVv2, &multilepton_JetHighestPt2_2ndPair_JEC_Up, &multilepton_JetHighestPt2_2ndPair_JEC_Down, vJetLoose.at(im2).JES_uncert(), &multilepton_JetHighestPt2_2ndPair_JER_Up, &multilepton_JetHighestPt2_2ndPair_JER_Down, vJetLoose.at(im2).pt_JER(), vJetLoose.at(im2).pt_JER_up(), vJetLoose.at(im2).pt_JER_down());
            multilepton_JetHighestPt2_2ndPair_Id = 1;
            multilepton_JetHighestPt1_2ndPair_CSV = vJetLoose.at(im2).DeepCSVbtag;
            multilepton_JetHighestPt2_2ndPair_P4.SetPtEtaPhiE(vJetLoose.at(im2).pt, vJetLoose.at(im2).eta, vJetLoose.at(im2).phi, vJetLoose.at(im2).E);
        }
        if (io1!=-1 && io2!=-1){
            // Jet1.SetPtEtaPhiE(vJetLoose.at(ip1).pt, vJetLoose.at(ip1).eta, vJetLoose.at(ip1).phi, vJetLoose.at(ip1).E);
            // Jet2.SetPtEtaPhiE(vJetLoose.at(io2).pt, vJetLoose.at(io2).eta, vJetLoose.at(io2).phi, vJetLoose.at(io2).E);
            // FillJetInfoOutputTree(&multilepton_JetClosestMw1_2ndPair_Id, 2, &multilepton_JetClosestMw1_2ndPair_P4, Jet1, &multilepton_JetClosestMw1_2ndPair_CSV, vJetLoose.at(io1).CSVv2, &multilepton_JetClosestMw1_2ndPair_JEC_Up, &multilepton_JetClosestMw1_2ndPair_JEC_Down, vJetLoose.at(io1).JES_uncert(), &multilepton_JetClosestMw1_2ndPair_JER_Up, &multilepton_JetClosestMw1_2ndPair_JER_Down, vJetLoose.at(io1).pt_JER(), vJetLoose.at(io1).pt_JER_up(), vJetLoose.at(io1).pt_JER_down());
            // FillJetInfoOutputTree(&multilepton_JetClosestMw2_2ndPair_Id, 2, &multilepton_JetClosestMw2_2ndPair_P4, Jet2, &multilepton_JetClosestMw2_2ndPair_CSV, vJetLoose.at(io2).CSVv2, &multilepton_JetClosestMw2_2ndPair_JEC_Up, &multilepton_JetClosestMw2_2ndPair_JEC_Down, vJetLoose.at(io2).JES_uncert(), &multilepton_JetClosestMw2_2ndPair_JER_Up, &multilepton_JetClosestMw2_2ndPair_JER_Down, vJetLoose.at(io2).pt_JER(), vJetLoose.at(io2).pt_JER_up(), vJetLoose.at(io2).pt_JER_down());
            multilepton_JetClosestMw1_2ndPair_Id = 2;
            multilepton_JetClosestMw2_2ndPair_Id = 2;
            multilepton_JetClosestMw1_2ndPair_CSV = vJetLoose.at(io1).DeepCSVbtag;
            multilepton_JetClosestMw2_2ndPair_CSV = vJetLoose.at(io2).DeepCSVbtag;
            multilepton_JetClosestMw1_2ndPair_P4.SetPtEtaPhiE(vJetLoose.at(io1).pt, vJetLoose.at(io1).eta, vJetLoose.at(io1).phi, vJetLoose.at(io1).E);
            multilepton_JetClosestMw2_2ndPair_P4.SetPtEtaPhiE(vJetLoose.at(io2).pt, vJetLoose.at(io2).eta, vJetLoose.at(io2).phi, vJetLoose.at(io2).E);
        }
        if (ip1!=-1 && ip2!=-1){
            // Jet1.SetPtEtaPhiE(vJetLoose.at(ip1).pt, vJetLoose.at(ip1).eta, vJetLoose.at(ip1).phi, vJetLoose.at(ip1).E);
            // Jet2.SetPtEtaPhiE(vJetLoose.at(ip2).pt, vJetLoose.at(ip2).eta, vJetLoose.at(ip2).phi, vJetLoose.at(ip2).E);
            // FillJetInfoOutputTree(&multilepton_JetLowestMjj1_2ndPair_Id, 3, &multilepton_JetLowestMjj1_2ndPair_P4, Jet1, &multilepton_JetLowestMjj1_2ndPair_CSV, vJetLoose.at(ip1).CSVv2, &multilepton_JetLowestMjj1_2ndPair_JEC_Up, &multilepton_JetLowestMjj1_2ndPair_JEC_Down, vJetLoose.at(ip1).JES_uncert(), &multilepton_JetLowestMjj1_2ndPair_JER_Up, &multilepton_JetLowestMjj1_2ndPair_JER_Down, vJetLoose.at(ip1).pt_JER(), vJetLoose.at(ip1).pt_JER_up(), vJetLoose.at(ip1).pt_JER_down());
            // FillJetInfoOutputTree(&multilepton_JetLowestMjj2_2ndPair_Id, 3, &multilepton_JetLowestMjj2_2ndPair_P4, Jet2, &multilepton_JetLowestMjj2_2ndPair_CSV, vJetLoose.at(ip2).CSVv2, &multilepton_JetLowestMjj2_2ndPair_JEC_Up, &multilepton_JetLowestMjj2_2ndPair_JEC_Down, vJetLoose.at(ip2).JES_uncert(), &multilepton_JetLowestMjj2_2ndPair_JER_Up, &multilepton_JetLowestMjj2_2ndPair_JER_Down, vJetLoose.at(ip2).pt_JER(), vJetLoose.at(ip2).pt_JER_up(), vJetLoose.at(ip2).pt_JER_down());
            multilepton_JetLowestMjj1_2ndPair_Id = 3;
            multilepton_JetLowestMjj2_2ndPair_Id = 3;
            multilepton_JetLowestMjj1_2ndPair_CSV = vJetLoose.at(ip1).DeepCSVbtag;
            multilepton_JetLowestMjj2_2ndPair_CSV = vJetLoose.at(ip2).DeepCSVbtag;
            multilepton_JetLowestMjj1_2ndPair_P4.SetPtEtaPhiE(vJetLoose.at(ip1).pt, vJetLoose.at(ip1).eta, vJetLoose.at(ip1).phi, vJetLoose.at(ip1).E);
            multilepton_JetLowestMjj2_2ndPair_P4.SetPtEtaPhiE(vJetLoose.at(ip2).pt, vJetLoose.at(ip2).eta, vJetLoose.at(ip2).phi, vJetLoose.at(ip2).E);
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
    multilepton_sumET = vEvent->at(0).metsumet;

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

    v_tOutput[itree]->Fill();

    delete JetsPt;
    delete JetsEta;
    delete JetsPhi;
    delete JetsE;
    delete JetsCSV;

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

    TLorentzVector lep1, lep2, lep3;
    bool lep1_isSSpair = false;
    bool lep2_isSSpair = false;
    bool lep3_isSSpair = false;
    lep1.SetPtEtaPhiE(vLeptonFakeable.at(0).conept, vLeptonFakeable.at(0).eta, vLeptonFakeable.at(0).phi, vLeptonFakeable.at(0).E );
    lep2.SetPtEtaPhiE(vLeptonFakeable.at(1).conept, vLeptonFakeable.at(1).eta, vLeptonFakeable.at(1).phi, vLeptonFakeable.at(1).E );
    if(vLeptonFakeable.at(0).id == vLeptonFakeable.at(1).id) {lep1_isSSpair = true; lep2_isSSpair = true;}
    if(region==3l)
    {
        if(vLeptonFakeable.at(0).id == vLeptonFakeable.at(2).id) {lep1_isSSpair = true; lep3_isSSpair = true;}
        if(vLeptonFakeable.at(2).id == vLeptonFakeable.at(1).id) {lep3_isSSpair = true; lep2_isSSpair = true;}

        lep3.SetPtEtaPhiE(vLeptonFakeable.at(2).conept, vLeptonFakeable.at(2).eta, vLeptonFakeable.at(2).phi, vLeptonFakeable.at(2).E );
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
    double tmp = -999, tmp2=-999;

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
    for(int ijet=0; ijet<vJetLoose.size(); ijet++)
    {
        if(vJetLoose.at(ijet).pt > 25 && fabs(vJetLoose.at(ijet).eta) < 2.4) {nJet25++;}
    }

    //NB : nJet25 != vJetLoose.size(), because vJetLoose accounts for JES variations (pT could be slightly < 25 !)
    nJetLoose = vJetLoose.size();


    //--- Var2: max eta of any 'non-CSV-loose' jet
    tmp = -999;
    maxEtaJet25 = -999;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta) > tmp) {tmp = fabs(vLightJets.at(ijet).eta);}
    }
    maxEtaJet25 = tmp;

	//--- Var3 : sum of leptons charges
    lepCharge = 0;
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
    dPhiHighestPtSSPair = -999; tmp = -999; tmp2=-999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vLeptonFakeable.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vLeptonFakeable.at(i).conept, vLeptonFakeable.at(i).eta, vLeptonFakeable.at(i).phi, vLeptonFakeable.at(i).E );
        for(int j=i+1; j<vLeptonFakeable.size(); j++)
        {
            lepj.SetPtEtaPhiE(vLeptonFakeable.at(j).conept, vLeptonFakeable.at(j).eta, vLeptonFakeable.at(j).phi, vLeptonFakeable.at(j).E );

            if(vLeptonFakeable.at(i).charge==vLeptonFakeable.at(j).charge && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vLeptonFakeable.at(i).phi - vLeptonFakeable.at(j).phi ) );
            }

            if(fabs(Phi_MPi_Pi(vLeptonFakeable.at(i).phi - vLeptonFakeable.at(j).phi)) > tmp2)
            {
                dPhiLepLep_max = fabs(Phi_MPi_Pi(vLeptonFakeable.at(i).phi - vLeptonFakeable.at(j).phi));
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    lepi = TLorentzVector();
    lepj = TLorentzVector();
    for(int i=0; i<vLeptonFakeable.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vLeptonFakeable.at(i).conept, vLeptonFakeable.at(i).eta, vLeptonFakeable.at(i).phi, vLeptonFakeable.at(i).E );
        for(int j=i+1; j<vLeptonFakeable.size(); j++)
        {
            lepj.SetPtEtaPhiE(vLeptonFakeable.at(j).conept, vLeptonFakeable.at(j).eta, vLeptonFakeable.at(j).phi, vLeptonFakeable.at(j).E );
            if(lepi.DeltaR(lepj) < minDRll)
            {
                minDRll = lepi.DeltaR(lepj);
            }
        }
    }

    //--- Var10 : pT of 3rd hardest lepton
    if(region == "3l") {Lep3Pt = vLeptonFakeable.at(2).conept;}
    else {Lep3Pt = vLeptonFakeable.at(1).conept;}

    //--- Fill additionnal variables, used for control only
    lep1Eta = vLeptonFakeable.at(0).eta; lep2Eta = vLeptonFakeable.at(1).eta; lep3Eta = vLeptonFakeable.at(1).eta;
    lep1Phi = vLeptonFakeable.at(0).phi; lep2Phi = vLeptonFakeable.at(1).phi; lep3Phi = vLeptonFakeable.at(1).phi;
    if(region == "3l") {lep3Eta = vLeptonFakeable.at(2).eta; lep3Phi = vLeptonFakeable.at(2).phi;}
    if(ijet_hardest_btag >= 0)
    {
        hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt;
        hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta;
    }


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
    lep2_conePt = vLeptonFakeable.at(1).conept;
    if(nleptons == 3) {lep3_conePt = vLeptonFakeable.at(2).conept;} //'Trailing' lepton
    else {lep3_conePt = lep2_conePt;}

    mindr_lep1_jet = 999;
    mindr_lep2_jet = 999;
    for(int i=0; i<vJetLoose.size(); i++)
    {
        float dr1 = GetDeltaR(vLeptonFakeable.at(0).eta,vLeptonFakeable.at(0).phi,vJetLoose.at(i).eta,vJetLoose.at(i).phi);
        if(dr1 < mindr_lep1_jet) {mindr_lep1_jet = dr1;}

        float dr2 = GetDeltaR(vLeptonFakeable.at(1).eta,vLeptonFakeable.at(1).phi,vJetLoose.at(i).eta,vJetLoose.at(i).phi);
        if(dr2 < mindr_lep2_jet) {mindr_lep2_jet = dr2;}
    }

    // mT_lep1 = comp_MT_met_lep(lep_tmp, vEvent->at(0).metpt, vEvent->at(0).metphi);
    mT_lep1 = sqrt( 2*vLeptonFakeable.at(0).conept * metpt * (1 - cos(vLeptonFakeable.at(0).phi - vEvent->at(0).metphi) ) );
    // mT_lep2 = comp_MT_met_lep(lep_tmp, vEvent->at(0).metpt, vEvent->at(0).metphi);
    mT_lep2 = sqrt( 2*vLeptonFakeable.at(1).conept * metpt * (1 - cos(vLeptonFakeable.at(1).phi - vEvent->at(0).metphi) ) );

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

    minv_FwdJetBJet = -999;
    if(ijet_forward >= 0 && ijet_hardest_btag >= 0) {minv_FwdJetBJet = (fwdJet + BJet).M();}

    FwdJetEta = -999;
    if(ijet_forward >= 0) {FwdJetEta = fwdJet.Eta();}

    FwdJetPt = -999;
    if(ijet_forward >= 0) {FwdJetPt = fwdJet.Pt();}

    LeadJetEta = -999;
    LeadJetPt = -999;
    sum_jetPt = 0;
    // nSoftJets = 0;
    for(int j=0; j<vJetLoose.size(); j++)
    {
        // if(vJetLoose.at(j).pt > 15 && vJetLoose.at(j).pt < 25) {nSoftJets++;}

        sum_jetPt+= vJetLoose.at(j).pt;

        if(vJetLoose.at(j).pt > LeadJetPt)
        {
            LeadJetPt = vJetLoose.at(j).pt;
            LeadJetEta = vJetLoose.at(j).eta;
        }
    }

    dRjj_max = 0;
    dPhijj_max = 0;
    Mjj_max = 0;
    TLorentzVector jet_tmp;
    for(int i=0; i<vJetLoose.size(); i++)
    {
        jet_tmp.SetPtEtaPhiE(vJetLoose.at(i).pt, vJetLoose.at(i).eta, vJetLoose.at(i).phi, vJetLoose.at(i).E );
        for(int j=i+1; j<vJetLoose.size(); j++)
        {
            TLorentzVector jet_tmp_2;
            jet_tmp_2.SetPtEtaPhiE(vJetLoose.at(j).pt, vJetLoose.at(j).eta, vJetLoose.at(j).phi, vJetLoose.at(j).E );

            if( fabs(GetDeltaR(vJetLoose.at(i).eta, vJetLoose.at(i).phi, vJetLoose.at(j).eta, vJetLoose.at(j).phi) ) > fabs(dRjj_max))
            {
                dRjj_max = GetDeltaR(vJetLoose.at(i).eta, vJetLoose.at(i).phi, vJetLoose.at(j).eta, vJetLoose.at(j).phi);
            }
            if( (jet_tmp + jet_tmp_2).M() >  Mjj_max)
            {
                Mjj_max = (jet_tmp + jet_tmp_2).M();
            }
            if(Phi_MPi_Pi(vJetLoose.at(i).phi - vJetLoose.at(j).phi ) > fabs(dPhijj_max) )
            {
                dPhijj_max = Phi_MPi_Pi(vJetLoose.at(i).phi - vJetLoose.at(j).phi );
            }
        }
    }

    deepCSV_max = 0;
    deepCSV_2nd = 0;
    for(int i=0; i<vLooseBTagJets.size(); i++)
    {
        if(vLooseBTagJets.at(i).DeepCSVbtag > deepCSV_max) {deepCSV_2nd = deepCSV_max; deepCSV_max = vLooseBTagJets.at(i).DeepCSVbtag;}
        else if(vLooseBTagJets.at(i).DeepCSVbtag > deepCSV_2nd) {deepCSV_2nd = vLooseBTagJets.at(i).DeepCSVbtag;}
    }

    dPhiLepBJet_max = 0;
    if(ijet_hardest_btag >= 0)
    {
        for(int i=0; i<nleptons; i++)
        {
            if( fabs(Phi_MPi_Pi(vLooseBTagJets.at(ijet_hardest_btag).phi - vLeptonFakeable.at(i).phi)) > fabs(dPhiLepBJet_max) )
            {
                dPhiLepBJet_max = fabs(Phi_MPi_Pi(vLooseBTagJets.at(ijet_hardest_btag).phi - vLeptonFakeable.at(i).phi) );
            }
        }
    }

    m3l = -999;
    if(region == "3l") {m3l = (lep1+lep2+lep3).M();}
    else {m3l = (lep1+lep2).M();}

    mHT = Compute_mHT();

    Compute_Top_W_variables(); //TESTING -- Fill high-level reco variables

    //3l - compute charge asymmetry for lepton (not part of SS pair) with min dR with bjet
    lW_asym = -999;
    float dr_min = 999;
    tmp = 0;
    if(ijet_hardest_btag >= 0)
    {
        if(!lep1_isSSpair)
        {
            tmp = GetDeltaR(lep1.Eta(), lep1.Phi(), vLooseBTagJets.at(ijet_hardest_btag).eta, vLooseBTagJets.at(ijet_hardest_btag).phi);
            if(fabs(tmp) < fabs(dr_min) ) {lW_asym = vLeptonFakeable.at(0).charge * fabs(vLeptonFakeable.at(0).eta);}
        }
        if(!lep2_isSSpair)
        {
            tmp = GetDeltaR(lep2.Eta(), lep2.Phi(), vLooseBTagJets.at(ijet_hardest_btag).eta, vLooseBTagJets.at(ijet_hardest_btag).phi);
            if(fabs(tmp) < fabs(dr_min) ) {lW_asym = vLeptonFakeable.at(1).charge * fabs(vLeptonFakeable.at(1).eta);}
        }
        if(region == "3l" && !lep3_isSSpair)
        {
            tmp = GetDeltaR(lep3.Eta(), lep3.Phi(), vLooseBTagJets.at(ijet_hardest_btag).eta, vLooseBTagJets.at(ijet_hardest_btag).phi);
            if(fabs(tmp) < fabs(dr_min) ) {lW_asym = vLeptonFakeable.at(2).charge * fabs(vLeptonFakeable.at(2).eta);}
        }    
    }

    min_dr_lep_bjet = 999;
    min_dr_lep_lightjet = 999;
    for(int i=0; i<vLeptonFakeable.size(); i++)
    {
        lepi.SetPtEtaPhiE(vLeptonFakeable.at(i).conept, vLeptonFakeable.at(i).eta, vLeptonFakeable.at(i).phi, vLeptonFakeable.at(i).E );

        for(int j=0; j<vJetLoose.size(); j++)
        {
            jet_tmp.SetPtEtaPhiE(vJetLoose.at(j).pt, vJetLoose.at(j).eta, vJetLoose.at(j).phi, vJetLoose.at(j).E );

            tmp = GetDeltaR(lepi.Eta(), lepi.Phi(), vJetLoose.at(j).eta, vJetLoose.at(j).phi);

            if(fabs(tmp) < fabs(min_dr_lep_lightjet) && !vJetLoose.at(j).isLooseBTag) {min_dr_lep_lightjet = tmp;}
            else if(fabs(tmp) < fabs(min_dr_lep_lightjet) && vJetLoose.at(j).isLooseBTag) {min_dr_lep_bjet = tmp;}
        }
    }

    tmp = 999;
    ratio_lep3pt_closestJetPt = -999;
    if(region == "3l") {lepi.SetPtEtaPhiE(vLeptonFakeable.at(2).conept, vLeptonFakeable.at(2).eta, vLeptonFakeable.at(2).phi, vLeptonFakeable.at(2).E );}
    else {lepi.SetPtEtaPhiE(vLeptonFakeable.at(1).conept, vLeptonFakeable.at(1).eta, vLeptonFakeable.at(1).phi, vLeptonFakeable.at(1).E );}
    for(int j=0; j<vJetLoose.size(); j++)
    {
        jet_tmp.SetPtEtaPhiE(vJetLoose.at(j).pt, vJetLoose.at(j).eta, vJetLoose.at(j).phi, vJetLoose.at(j).E );

        if(fabs(GetDeltaR(lepi.Eta(), lepi.Phi(), vJetLoose.at(j).eta, vJetLoose.at(j).phi) ) < fabs(tmp) ) {ratio_lep3pt_closestJetPt = lepi.Pt() / vJetLoose.at(j).pt;}
    }

    dPhiLepLep_hardestOS = -999;
    if(region == "2l") {dPhiLepLep_hardestOS = fabs(Phi_MPi_Pi(vLeptonFakeable.at(1).phi - vLeptonFakeable.at(0).phi) );}
    else
    {
        if(vLeptonFakeable.at(0).id == - vLeptonFakeable.at(1).id) {dPhiLepLep_hardestOS =  fabs(Phi_MPi_Pi(vLeptonFakeable.at(1).phi - vLeptonFakeable.at(0).phi) );}
        else if(vLeptonFakeable.at(0).id == - vLeptonFakeable.at(2).id) {dPhiLepLep_hardestOS =  fabs(Phi_MPi_Pi(vLeptonFakeable.at(2).phi - vLeptonFakeable.at(0).phi) );}
        else if(vLeptonFakeable.at(1).id == - vLeptonFakeable.at(2).id) {dPhiLepLep_hardestOS =  fabs(Phi_MPi_Pi(vLeptonFakeable.at(1).phi - vLeptonFakeable.at(2).phi) );}
    }

    //--------------------------------------------
    // #    #      #    #####   ##    ####   ####  ###### #####
    // #    #      #      #    #  #  #    # #    # #      #    #
    // ######      #      #   #    # #      #      #####  #    #
    // #    #      #      #   ###### #  ### #  ### #      #####
    // #    # #    #      #   #    # #    # #    # #      #   #
    // #    #  ####       #   #    #  ####   ####  ###### #    #
    //--------------------------------------------

    //Input variables needed for Hj tagger
    Jet25_lepdrmin;
    Jet25_lepdrmax;
    Jet25_qg;
    Jet25_bDiscriminator;
    Jet25_pt;

    tmp = -999; tmp2 = 999;

    for(int ijet=0; ijet<vJetLoose.size(); ijet++)
    {
        for(int ilep=0; ilep<vLeptonFakeable.size(); ilep++)
        {
            float dr_lep_jet = fabs(GetDeltaR(vJetLoose.at(ijet).eta, vJetLoose.at(ijet).phi, vLeptonFakeable.at(ilep).eta, vLeptonFakeable.at(ilep).phi));
            if(dr_lep_jet > tmp) {tmp = dr_lep_jet;}
            if(dr_lep_jet < tmp2) {tmp2 = dr_lep_jet;}
        }

        Jet25_qg = vJetLoose.at(ijet).qgtag;  Jet25_qg = (Jet25_qg<0) ? 0 : Jet25_qg;
        Jet25_lepdrmin = tmp2;
        Jet25_lepdrmax = tmp;
        Jet25_bDiscriminator = vJetLoose.at(ijet).DeepCSVbtag; Jet25_bDiscriminator = (Jet25_bDiscriminator<0) ? 0 : Jet25_bDiscriminator;
        Jet25_pt = vJetLoose.at(ijet).pt;

        double HjTag_jet = mva_HjTagger->EvaluateMVA("BDTG method");

        if(ijet==0) {HjTag_max = HjTag_jet; HjTag_mean = HjTag_jet;}
        else
        {
            HjTag_mean+= HjTag_jet;

            if(HjTag_jet > HjTag_max) {HjTag_max = HjTag_jet;}
        }

        // cout<<"HjTag_jet = "<<HjTag_jet<<endl;
    }
    HjTag_mean/= vJetLoose.size();



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
            cout<<ilep<<" pT = "<<vLeptonFakeable.at(ilep).conept<<" / eta = "<<vLeptonFakeable.at(ilep).eta<<" / phi = "<<vLeptonFakeable.at(ilep).phi<<" : Charge = "<<vLeptonFakeable.at(ilep).charge<<endl;
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
        // signal_3l_TT_MVA   = mva_3l_tt->EvaluateMVA("BDTG method");
        // signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");//--------------------------------------------
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
 * systname : can be "" (default TTree), "JERUp", etc.
 */
void tHqMultileponAnalysis::initializeOutputTree(int itree)
{
    outputfile->cd();

    // tOutput = new TTree("Tree", "Tree");

	//-- Event main infos
	v_tOutput[itree]->Branch("channel",&channel,"channel/F");
    v_tOutput[itree]->Branch("weight",&weight,"weight/F");
    v_tOutput[itree]->Branch("weight_old",&weight_old,"weight_old/F");
	v_tOutput[itree]->Branch("weightfake",&weightfake,"weightfake/F");
    v_tOutput[itree]->Branch("weightflip",&weightflip,"weightflip/F");
    v_tOutput[itree]->Branch("event_id",&event_id,"event_id/I");
    v_tOutput[itree]->Branch("is_trigger_1lep",&is_trigger_1lep,"is_trigger_1lep/O");
    v_tOutput[itree]->Branch("is_trigger_2lep",&is_trigger_2lep,"is_trigger_2lep/O");
    v_tOutput[itree]->Branch("is_trigger_3lep",&is_trigger_3lep,"is_trigger_3lep/O");
    v_tOutput[itree]->Branch("is_trigger_2lss",&is_trigger_2lss,"is_trigger_2lss/O"); //not defined globally
    v_tOutput[itree]->Branch("is_trigger_3l",&is_trigger_3l,"is_trigger_3l/O");
    // v_tOutput[itree]->Branch("event_run",&event_run,"event_run/F");
    v_tOutput[itree]->Branch("mc_weight",&mc_weight,"mc_weight/F");

    v_tOutput[itree]->Branch("total_SF",&total_SF,"total_SF/F");
    v_tOutput[itree]->Branch("lepton_SF",&lepton_SF,"lepton_SF/F");
    v_tOutput[itree]->Branch("trigger_SF",&trigger_SF,"trigger_SF/F");
    v_tOutput[itree]->Branch("btag_SF",&btag_SF,"btag_SF/F");
    v_tOutput[itree]->Branch("PU_SF",&PU_SF,"PU_SF/F");
    v_tOutput[itree]->Branch("nPU",&nPU,"nPU/I");

    if(write_LHE_weights_allFiles && !itree) //only for nominal tree, and if user asks for it
    {
        v_tOutput[itree]->Branch("LHEweights", &LHEweights); //LHE weights
        v_tOutput[itree]->Branch("LHEweights_Ids", &LHEweights_Ids); //LHE weights Ids
        v_tOutput[itree]->Branch("sumWeights_SMcoupling",&sumWeights_SMcoupling,"sumWeights_SMcoupling/F");
    }


	//--- Categories & MVA

    //tHq 2017 categories (implemented in NTA)
    v_tOutput[itree]->Branch("is_tHq_2lSS",&is_tHq_2lSS,"is_tHq_2lSS/B");
	v_tOutput[itree]->Branch("is_tHq_2lSS_SR",&is_tHq_2lSS_SR,"is_tHq_2lSS_SR/B");
	v_tOutput[itree]->Branch("is_tHq_2lSS_Training",&is_tHq_2lSS_Training,"is_tHq_2lSS_Training/B");
    v_tOutput[itree]->Branch("is_tHq_2lSS_Fake",&is_tHq_2lSS_Fake,"is_tHq_2lSS_Fake/B");
    v_tOutput[itree]->Branch("is_tHq_2lSS_Flip",&is_tHq_2lSS_Flip,"is_tHq_2lSS_Flip/B");
    v_tOutput[itree]->Branch("is_tHq_2lSS_GammaConv",&is_tHq_2lSS_GammaConv,"is_tHq_2lSS_GammaConv/B");
    v_tOutput[itree]->Branch("is_tHq_3l",&is_tHq_3l,"is_tHq_3l/B");
	v_tOutput[itree]->Branch("is_tHq_3l_SR",&is_tHq_3l_SR,"is_tHq_3l_SR/B");
	v_tOutput[itree]->Branch("is_tHq_3l_Training",&is_tHq_3l_Training,"is_tHq_3l_Training/B");
    v_tOutput[itree]->Branch("is_tHq_3l_Fake",&is_tHq_3l_Fake,"is_tHq_3l_Fake/B");
    v_tOutput[itree]->Branch("is_tHq_3l_GammaConv",&is_tHq_3l_GammaConv,"is_tHq_3l_GammaConv/B");
    v_tOutput[itree]->Branch("is_tHq_ttWctrl",&is_tHq_ttWctrl,"is_tHq_ttWctrl/B");
    v_tOutput[itree]->Branch("is_tHq_ttWctrl_SR",&is_tHq_ttWctrl_SR,"is_tHq_ttWctrl_SR/B");
	v_tOutput[itree]->Branch("is_tHq_ttWctrl_Fake",&is_tHq_ttWctrl_Fake,"is_tHq_ttWctrl_Fake/B");
    v_tOutput[itree]->Branch("is_tHq_ttWctrl_Flip",&is_tHq_ttWctrl_Flip,"is_tHq_ttWctrl_Flip/B");
    v_tOutput[itree]->Branch("is_tHq_ttWctrl_GammaConv",&is_tHq_ttWctrl_GammaConv,"is_tHq_ttWctrl_GammaConv/B");
    v_tOutput[itree]->Branch("is_tHq_ttZctrl",&is_tHq_ttZctrl,"is_tHq_ttZctrl/B");
    v_tOutput[itree]->Branch("is_tHq_ttZctrl_SR",&is_tHq_ttZctrl_SR,"is_tHq_ttZctrl_SR/B");
    v_tOutput[itree]->Branch("is_tHq_ttZctrl_Fake",&is_tHq_ttZctrl_Fake,"is_tHq_ttZctrl_Fake/B");
    v_tOutput[itree]->Branch("is_tHq_ttZctrl_GammaConv",&is_tHq_ttZctrl_GammaConv,"is_tHq_ttZctrl_GammaConv/B");
    v_tOutput[itree]->Branch("is_tHq_WZctrl",&is_tHq_WZctrl,"is_tHq_WZctrl/B");
    v_tOutput[itree]->Branch("is_tHq_WZctrl_SR",&is_tHq_WZctrl_SR,"is_tHq_WZctrl_SR/B");
    v_tOutput[itree]->Branch("is_tHq_WZctrl_Fake",&is_tHq_WZctrl_Fake,"is_tHq_WZctrl_Fake/B");
    v_tOutput[itree]->Branch("is_tHq_WZctrl_GammaConv",&is_tHq_WZctrl_GammaConv,"is_tHq_WZctrl_GammaConv/B");
    v_tOutput[itree]->Branch("is_tHq_4l_SR",&is_tHq_4l_SR,"is_tHq_4l_SR/B");
    v_tOutput[itree]->Branch("is_tHq_ZZctrl_SR",&is_tHq_ZZctrl_SR,"is_tHq_ZZctrl_SR/B");
    v_tOutput[itree]->Branch("is_tHqFCNC_2lSS_SR",&is_tHqFCNC_2lSS_SR,"is_tHqFCNC_2lSS_SR/B");
    v_tOutput[itree]->Branch("is_tHqFCNC_2lSS_Fake",&is_tHqFCNC_2lSS_Fake,"is_tHqFCNC_2lSS_Fake/B");
    v_tOutput[itree]->Branch("is_tHqFCNC_2lSS_Flip",&is_tHqFCNC_2lSS_Flip,"is_tHqFCNC_2lSS_Flip/B");
    v_tOutput[itree]->Branch("is_tHqFCNC_2lSS_GammaConv",&is_tHqFCNC_2lSS_GammaConv,"is_tHqFCNC_2lSS_GammaConv/B");
    v_tOutput[itree]->Branch("is_tHqFCNC_3l_SR",&is_tHqFCNC_3l_SR,"is_tHqFCNC_3l_SR/B");
    v_tOutput[itree]->Branch("is_tHqFCNC_3l_Fake",&is_tHqFCNC_3l_Fake,"is_tHqFCNC_3l_Fake/B");
    v_tOutput[itree]->Branch("is_tHqFCNC_3l_GammaConv",&is_tHqFCNC_3l_GammaConv,"is_tHqFCNC_3l_GammaConv/B");


    //ttH2017 predefined categories (implemented in NTP)
    v_tOutput[itree]->Branch("is_ttH_2lSS",&is_ttH_2lSS,"is_ttH_2lSS/B");
    v_tOutput[itree]->Branch("is_ttH_2lSS_Training",&is_ttH_2lSS_Training,"is_ttH_2lSS_Training/B");
    v_tOutput[itree]->Branch("is_ttH_2lSS_SR",&is_ttH_2lSS_SR,"is_ttH_2lSS_SR/B");
    v_tOutput[itree]->Branch("is_ttH_2lSS_Fake",&is_ttH_2lSS_Fake,"is_ttH_2lSS_Fake/B");
    v_tOutput[itree]->Branch("is_ttH_2lSS_Flip",&is_ttH_2lSS_Flip,"is_ttH_2lSS_Flip/B");
    v_tOutput[itree]->Branch("is_ttH_2lSS_GammaConv",&is_ttH_2lSS_GammaConv,"is_ttH_2lSS_GammaConv/B");
    v_tOutput[itree]->Branch("is_ttH_3l",&is_ttH_3l,"is_ttH_3l/B");
    v_tOutput[itree]->Branch("is_ttH_3l_Training",&is_ttH_3l_Training,"is_ttH_3l_Training/B");
    v_tOutput[itree]->Branch("is_ttH_3l_SR",&is_ttH_3l_SR,"is_ttH_3l_SR/B");
    v_tOutput[itree]->Branch("is_ttH_3l_Fake",&is_ttH_3l_Fake,"is_ttH_3l_Fake/B");
    v_tOutput[itree]->Branch("is_ttH_3l_GammaConv",&is_ttH_3l_GammaConv,"is_ttH_3l_GammaConv/B");
    v_tOutput[itree]->Branch("is_ttH_ttWctrl",&is_ttH_ttWctrl,"is_ttH_ttWctrl/B");
    v_tOutput[itree]->Branch("is_ttH_ttWctrl_SR",&is_ttH_ttWctrl_SR,"is_ttH_ttWctrl_SR/B");
    v_tOutput[itree]->Branch("is_ttH_ttWctrl_Fake",&is_ttH_ttWctrl_Fake,"is_ttH_ttWctrl_Fake/B");
    v_tOutput[itree]->Branch("is_ttH_ttWctrl_Flip",&is_ttH_ttWctrl_Flip,"is_ttH_ttWctrl_Flip/B");
    v_tOutput[itree]->Branch("is_ttH_ttWctrl_GammaConv",&is_ttH_ttWctrl_GammaConv,"is_ttH_ttWctrl_GammaConv/B");
    v_tOutput[itree]->Branch("is_ttH_ttZctrl",&is_ttH_ttZctrl,"is_ttH_ttZctrl/B");
    v_tOutput[itree]->Branch("is_ttH_ttZctrl_SR",&is_ttH_ttZctrl_SR,"is_ttH_ttZctrl_SR/B");
    v_tOutput[itree]->Branch("is_ttH_ttZctrl_Fake",&is_ttH_ttZctrl_Fake,"is_ttH_ttZctrl_Fake/B");
    v_tOutput[itree]->Branch("is_ttH_ttZctrl_GammaConv",&is_ttH_ttZctrl_GammaConv,"is_ttH_ttZctrl_GammaConv/B");
    v_tOutput[itree]->Branch("is_ttH_WZctrl",&is_ttH_WZctrl,"is_ttH_WZctrl/B");
    v_tOutput[itree]->Branch("is_ttH_WZctrl_SR",&is_ttH_WZctrl_SR,"is_ttH_WZctrl_SR/B");
    v_tOutput[itree]->Branch("is_ttH_WZctrl_Fake",&is_ttH_WZctrl_Fake,"is_ttH_WZctrl_Fake/B");
    v_tOutput[itree]->Branch("is_ttH_WZctrl_GammaConv",&is_ttH_WZctrl_GammaConv,"is_ttH_WZctrl_GammaConv/B");
    v_tOutput[itree]->Branch("is_ttH_4l_SR",&is_ttH_4l_SR,"is_ttH_4l_SR/B");
    v_tOutput[itree]->Branch("is_ttH_ZZctrl_SR",&is_ttH_ZZctrl_SR,"is_ttH_ZZctrl_SR/B");



    if(add_orthogocal_cat)
    {
        v_tOutput[itree]->Branch("is_ttH_2lSS_SR_fwd",&is_ttH_2lSS_SR_fwd,"is_ttH_2lSS_SR_fwd/B");
        v_tOutput[itree]->Branch("is_ttH_3l_SR_fwd",&is_ttH_3l_SR_fwd,"is_ttH_3l_SR_fwd/B");
        v_tOutput[itree]->Branch("is_tHq_2lSS_SR_fwd",&is_tHq_2lSS_SR_fwd,"is_tHq_2lSS_SR_fwd/B");
        v_tOutput[itree]->Branch("is_tHq_3l_SR_fwd",&is_tHq_3l_SR_fwd,"is_tHq_3l_SR_fwd/B");

        v_tOutput[itree]->Branch("is_ttH_2lSS_SR_btag",&is_ttH_2lSS_SR_btag,"is_ttH_2lSS_SR_btag/B");
        v_tOutput[itree]->Branch("is_ttH_3l_SR_btag",&is_ttH_3l_SR_btag,"is_ttH_3l_SR_btag/B");
        v_tOutput[itree]->Branch("is_tHq_2lSS_SR_btag",&is_tHq_2lSS_SR_btag,"is_tHq_2lSS_SR_btag/B");
        v_tOutput[itree]->Branch("is_tHq_3l_SR_btag",&is_tHq_3l_SR_btag,"is_tHq_3l_SR_btag/B");

        v_tOutput[itree]->Branch("is_ttH_2lSS_SR_njet",&is_ttH_2lSS_SR_njet,"is_ttH_2lSS_SR_njet/B");
        v_tOutput[itree]->Branch("is_ttH_3l_SR_njet",&is_ttH_3l_SR_njet,"is_ttH_3l_SR_njet/B");
        v_tOutput[itree]->Branch("is_tHq_2lSS_SR_njet",&is_tHq_2lSS_SR_njet,"is_tHq_2lSS_SR_njet/B");
        v_tOutput[itree]->Branch("is_tHq_3l_SR_njet",&is_tHq_3l_SR_njet,"is_tHq_3l_SR_njet/B");
        v_tOutput[itree]->Branch("is_ttH_ttWctrl_SR_njet",&is_ttH_ttWctrl_SR_njet,"is_ttH_ttWctrl_SR_njet/B");

        v_tOutput[itree]->Branch("is_tHq_2lSS_SR_fwd2",&is_tHq_2lSS_SR_fwd2,"is_tHq_2lSS_SR_fwd2/B");
        v_tOutput[itree]->Branch("is_ttH_ttWctrl_SR_fwd2",&is_ttH_ttWctrl_SR_fwd2,"is_ttH_ttWctrl_SR_fwd2/B");

        v_tOutput[itree]->Branch("is_tHq_3l_SR_njet2",&is_tHq_3l_SR_njet2,"is_tHq_3l_SR_njet2/B");
        v_tOutput[itree]->Branch("is_ttH_3l_SR_njet2",&is_ttH_3l_SR_njet2,"is_ttH_3l_SR_njet2/B");

        v_tOutput[itree]->Branch("is_tHq_3l_SR_njet3",&is_tHq_3l_SR_njet3,"is_tHq_3l_SR_njet3/B");
        v_tOutput[itree]->Branch("is_ttH_3l_SR_njet3",&is_ttH_3l_SR_njet3,"is_ttH_3l_SR_njet3/B");
    }



	//-- Nof objects variables
	v_tOutput[itree]->Branch("nLooseBJets",&nLooseBJets,"nLooseBJets/F");
	v_tOutput[itree]->Branch("nMediumBJets",&nMediumBJets,"nMediumBJets/F");
	v_tOutput[itree]->Branch("nTightLep",&nTightLep,"nTightLep/F");
	v_tOutput[itree]->Branch("nFakeableLep",&nFakeableLep,"nFakeableLep/F");

    v_tOutput[itree]->Branch("nJets",&nJets,"nJets/I");
    v_tOutput[itree]->Branch("nLightJets",&nLightJets,"nLightJets/F");
    v_tOutput[itree]->Branch("nLightJets_Fwd",&nLightJets_Fwd,"nLightJets_Fwd/F");

    v_tOutput[itree]->Branch("is_hasJetNoisyHCAL",&is_hasJetNoisyHCAL,"is_hasJetNoisyHCAL/O");
    v_tOutput[itree]->Branch("is_hasManyJetNoisyHCAL",&is_hasManyJetNoisyHCAL,"is_hasManyJetNoisyHCAL/O");
    v_tOutput[itree]->Branch("JetNoisyHCALPt",&JetNoisyHCALPt,"JetNoisyHCALPt/F");
    v_tOutput[itree]->Branch("JetNoisyHCALEta",&JetNoisyHCALEta,"JetNoisyHCALEta/F");


	// v_tOutput[itree]->Branch("signal_3l_TT_MVA",&signal_3l_TT_MVA,"signal_3l_TT_MVA/F");
	// v_tOutput[itree]->Branch("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,"signal_3l_TTV_MVA/F");
    // v_tOutput[itree]->Branch("signal_2lss_TT_MVA",&signal_2lss_TT_MVA,"signal_2lss_TT_MVA/F");
    // v_tOutput[itree]->Branch("signal_2lss_TTV_MVA",&signal_2lss_TTV_MVA,"signal_2lss_TTV_MVA/F");

    //Default variables "nJets" & "nLightJets" correspond to the chosen analysis (tHq/ttH) and are always stored
    //Therefore, we will also store these variables for the other analysis
    if(do_tHq_analysis) //if considering tHq analysis, also store ttH variables
    {
        v_tOutput[itree]->Branch("nJets_ttH",&nJets_ttH,"nJets_ttH/I");
        v_tOutput[itree]->Branch("nLightJets_ttH",&nLightJets_ttH,"nLightJets_ttH/F");
    }
    else //if considering ttH analysis, also store tHq variables
    {
        v_tOutput[itree]->Branch("nJets_tHq",&nJets_tHq,"nJets_tHq/I");
        v_tOutput[itree]->Branch("nLightJets_tHq",&nLightJets_tHq,"nLightJets_tHq/F");
    }

    // if(!make_ntuples_for_overlap_studies)
    {
        //-- Input variables from tHq2016 analysis
        v_tOutput[itree]->Branch("nJet25",&nJet25,"nJet25/F");
        v_tOutput[itree]->Branch("nJetLoose",&nJetLoose,"nJetLoose/F");
        v_tOutput[itree]->Branch("maxEtaJet25",&maxEtaJet25,"maxEtaJet25/F");
        v_tOutput[itree]->Branch("lepCharge",&lepCharge,"lepCharge/F");
        v_tOutput[itree]->Branch("nJetEta1",&nJetEta1,"nJetEta1/F");
        v_tOutput[itree]->Branch("dEtaFwdJetBJet",&dEtaFwdJetBJet,"dEtaFwdJetBJet/F");
        v_tOutput[itree]->Branch("dEtaFwdJet2BJet",&dEtaFwdJet2BJet,"dEtaFwdJet2BJet/F");
        v_tOutput[itree]->Branch("dEtaFwdJetClosestLep",&dEtaFwdJetClosestLep,"dEtaFwdJetClosestLep/F");
        v_tOutput[itree]->Branch("dPhiHighestPtSSPair",&dPhiHighestPtSSPair,"dPhiHighestPtSSPair/F");
        v_tOutput[itree]->Branch("minDRll",&minDRll,"minDRll/F");
        v_tOutput[itree]->Branch("Lep3Pt",&Lep3Pt,"Lep3Pt/F");

        //-- Input variables from ttH2017 analysis
        v_tOutput[itree]->Branch("lep1_conePt",&lep1_conePt,"lep1_conePt/F");
        v_tOutput[itree]->Branch("lep2_conePt",&lep2_conePt,"lep2_conePt/F");
        v_tOutput[itree]->Branch("lep3_conePt",&lep3_conePt,"lep3_conePt/F");
        v_tOutput[itree]->Branch("mindr_lep1_jet",&mindr_lep1_jet,"mindr_lep1_jet/F");
        v_tOutput[itree]->Branch("mindr_lep2_jet",&mindr_lep2_jet,"mindr_lep2_jet/F");
        v_tOutput[itree]->Branch("mT_lep1",&mT_lep1,"mT_lep1/F");
        v_tOutput[itree]->Branch("mT_lep2",&mT_lep2,"mT_lep2/F");
        v_tOutput[itree]->Branch("max_lep_eta",&max_lep_eta,"max_lep_eta/F");

        //-- new variables
        v_tOutput[itree]->Branch("minv_FwdJetBJet",&minv_FwdJetBJet,"minv_FwdJetBJet/F");
        v_tOutput[itree]->Branch("FwdJetEta",&FwdJetEta,"FwdJetEta/F");
        v_tOutput[itree]->Branch("FwdJetPt",&FwdJetPt,"FwdJetPt/F");
        v_tOutput[itree]->Branch("LeadJetEta",&LeadJetEta,"LeadJetEta/F");
        v_tOutput[itree]->Branch("LeadJetPt",&LeadJetPt,"LeadJetPt/F");
        v_tOutput[itree]->Branch("dRjj_max",&dRjj_max,"dRjj_max/F");
        v_tOutput[itree]->Branch("deepCSV_max",&deepCSV_max,"deepCSV_max/F");
        v_tOutput[itree]->Branch("deepCSV_2nd",&deepCSV_2nd,"deepCSV_2nd/F");
        v_tOutput[itree]->Branch("Mjj_max",&Mjj_max,"Mjj_max/F");
        v_tOutput[itree]->Branch("dPhiLepBJet_max",&dPhiLepBJet_max,"dPhiLepBJet_max/F");
        v_tOutput[itree]->Branch("dPhijj_max",&dPhijj_max,"dPhijj_max/F");
        v_tOutput[itree]->Branch("m3l",&m3l,"m3l/F");
        v_tOutput[itree]->Branch("dPhiLepLep_max",&dPhiLepLep_max,"dPhiLepLep_max/F");

        v_tOutput[itree]->Branch("top_mass",&top_mass,"top_mass/F");
        v_tOutput[itree]->Branch("mTW",&mTW,"mTW/F");
        v_tOutput[itree]->Branch("lW_asym_mtop",&lW_asym_mtop,"lW_asym_mtop/F");
        v_tOutput[itree]->Branch("dRBjetRecoilJet",&dRBjetRecoilJet,"dRBjetRecoilJet/F");
        v_tOutput[itree]->Branch("dRLepWRecoilJet",&dRLepWRecoilJet,"dRLepWRecoilJet/F");
        v_tOutput[itree]->Branch("RecoilJetPt",&RecoilJetPt,"RecoilJetPt/F");
        v_tOutput[itree]->Branch("RecoilJetEta",&RecoilJetEta,"RecoilJetEta/F");
        v_tOutput[itree]->Branch("LepWPt",&LepWPt,"LepWPt/F");
        v_tOutput[itree]->Branch("LepWEta",&LepWEta,"LepWEta/F");
        v_tOutput[itree]->Branch("top_Pt",&top_Pt,"top_Pt/F");
        v_tOutput[itree]->Branch("mass_LepBJet_min",&mass_LepBJet_min,"mass_LepBJet_min/F");
        v_tOutput[itree]->Branch("sum_jetPt",&sum_jetPt,"sum_jetPt/F");

        //not used for now (would need to rerun NTP with qgjet info, ...)
        v_tOutput[itree]->Branch("HjTag_max",&HjTag_max,"HjTag_max/F");
        v_tOutput[itree]->Branch("HjTag_mean",&HjTag_mean,"HjTag_mean/F");

    	//-- More control vars
        v_tOutput[itree]->Branch("inv_mll",&inv_mll,"inv_mll/F");
        v_tOutput[itree]->Branch("hardestBjetPt",&hardestBjetPt,"hardestBjetPt/F");
        v_tOutput[itree]->Branch("hardestBjetEta",&hardestBjetEta,"hardestBjetEta/F");
        v_tOutput[itree]->Branch("lep1Eta",&lep1Eta,"lep1Eta/F");
        v_tOutput[itree]->Branch("lep2Eta",&lep2Eta,"lep2Eta/F");
        v_tOutput[itree]->Branch("lep3Eta",&lep3Eta,"lep3Eta/F");
        v_tOutput[itree]->Branch("lep1Phi",&lep1Phi,"lep1Phi/F");
        v_tOutput[itree]->Branch("lep2Phi",&lep2Phi,"lep2Phi/F");
        v_tOutput[itree]->Branch("lep3Phi",&lep3Phi,"lep3Phi/F");
        v_tOutput[itree]->Branch("metpt",&metpt,"metpt/F");
        v_tOutput[itree]->Branch("metphi",&metphi,"metphi/F");
        v_tOutput[itree]->Branch("metLD",&metLD,"metLD/F");
        v_tOutput[itree]->Branch("mHT",&mHT,"mHT/F");

        //Additonal variables for FCNC analysis
        // v_tOutput[itree]->Branch("nSoftJets",&nSoftJets,"nSoftJets/F");
        v_tOutput[itree]->Branch("min_dr_lep_bjet",&min_dr_lep_bjet,"min_dr_lep_bjet/F");
        v_tOutput[itree]->Branch("min_dr_lep_lightjet",&min_dr_lep_lightjet,"min_dr_lep_lightjet/F");
        v_tOutput[itree]->Branch("lW_asym",&lW_asym,"lW_asym/F");
        v_tOutput[itree]->Branch("ratio_lep3pt_closestJetPt",&ratio_lep3pt_closestJetPt,"ratio_lep3pt_closestJetPt/F");
        v_tOutput[itree]->Branch("dPhiLepLep_hardestOS",&dPhiLepLep_hardestOS,"dPhiLepLep_hardestOS/F");

        //Decay modes
        v_tOutput[itree]->Branch("higgs_daughter_id",&higgs_daughter_id,"higgs_daughter_id/I");

        if(_sampleName.Contains("WZ"))
        {
            v_tOutput[itree]->Branch("wz_jetFlav_b",&wz_jetFlav_b,"wz_jetFlav_b/B");
            v_tOutput[itree]->Branch("wz_jetFlav_c",&wz_jetFlav_c,"wz_jetFlav_c/B");
            v_tOutput[itree]->Branch("wz_jetFlav_l",&wz_jetFlav_l,"wz_jetFlav_l/B");
        }

        if(!_isdata && v_systTree[itree] == "") //Save for default TTree only
        {
            //Systematics weights
            v_tOutput[itree]->Branch("TrigEffUp",&TrigEffUp,"TrigEffUp/F");
            v_tOutput[itree]->Branch("TrigEffDown",&TrigEffDown,"TrigEffDown/F");
            // v_tOutput[itree]->Branch("MuEffUp",&MuEffUp,"MuEffUp/F");
            // v_tOutput[itree]->Branch("MuEffDown",&MuEffDown,"MuEffDown/F");
            // v_tOutput[itree]->Branch("EleEffUp",&EleEffUp,"EleEffUp/F");
            // v_tOutput[itree]->Branch("EleEffDown",&EleEffDown,"EleEffDown/F");
            v_tOutput[itree]->Branch("LepEff_elLooseUp",&LepEff_elLooseUp,"LepEff_elLooseUp/F");
            v_tOutput[itree]->Branch("LepEff_elLooseDown",&LepEff_elLooseDown,"LepEff_elLooseDown/F");
            v_tOutput[itree]->Branch("LepEff_muLooseUp",&LepEff_muLooseUp,"LepEff_muLooseUp/F");
            v_tOutput[itree]->Branch("LepEff_muLooseDown",&LepEff_muLooseDown,"LepEff_muLooseDown/F");
            v_tOutput[itree]->Branch("LepEff_elTightUp",&LepEff_elTightUp,"LepEff_elTightUp/F");
            v_tOutput[itree]->Branch("LepEff_elTightDown",&LepEff_elTightDown,"LepEff_elTightDown/F");
            v_tOutput[itree]->Branch("LFcontUp",&LFcontUp,"LFcontUp/F");
            v_tOutput[itree]->Branch("LFcontDown",&LFcontDown,"LFcontDown/F");
            v_tOutput[itree]->Branch("HFstats1Up",&HFstats1Up,"HFstats1Up/F");
            v_tOutput[itree]->Branch("HFstats1Down",&HFstats1Down,"HFstats1Down/F");
            v_tOutput[itree]->Branch("HFstats2Up",&HFstats2Up,"HFstats2Up/F");
            v_tOutput[itree]->Branch("HFstats2Down",&HFstats2Down,"HFstats2Down/F");
            v_tOutput[itree]->Branch("CFerr1Up",&CFerr1Up,"CFerr1Up/F");
            v_tOutput[itree]->Branch("CFerr1Down",&CFerr1Down,"CFerr1Down/F");
            v_tOutput[itree]->Branch("CFerr2Up",&CFerr2Up,"CFerr2Up/F");
            v_tOutput[itree]->Branch("CFerr2Down",&CFerr2Down,"CFerr2Down/F");
            v_tOutput[itree]->Branch("HFcontUp",&HFcontUp,"HFcontUp/F");
            v_tOutput[itree]->Branch("HFcontDown",&HFcontDown,"HFcontDown/F");
            v_tOutput[itree]->Branch("LFstats1Up",&LFstats1Up,"LFstats1Up/F");
            v_tOutput[itree]->Branch("LFstats1Down",&LFstats1Down,"LFstats1Down/F");
            v_tOutput[itree]->Branch("LFstats2Up",&LFstats2Up,"LFstats2Up/F");
            v_tOutput[itree]->Branch("LFstats2Down",&LFstats2Down,"LFstats2Down/F");
            v_tOutput[itree]->Branch("PUUp",&PUUp,"PUUp/F");
            v_tOutput[itree]->Branch("PUDown",&PUDown,"PUDown/F");
            v_tOutput[itree]->Branch("QCDscaleUp",&QCDscaleUp,"QCDscaleUp/F");
            v_tOutput[itree]->Branch("QCDscaleDown",&QCDscaleDown,"QCDscaleDown/F");
            v_tOutput[itree]->Branch("pdfUp",&pdfUp,"pdfUp/F");
            v_tOutput[itree]->Branch("pdfDown",&pdfDown,"pdfDown/F");

            //FakeRate shape variations
            for(int ivar=0; ivar<v_FR_type.size(); ivar++)
            {
                v_tOutput[itree]->Branch(v_FR_type[ivar],&v_floats_FR_variations[ivar],v_FR_type[ivar]+"/F");
            }


            if(write_allScale_Variations)
            {
                // v_tOutput[itree]->Branch("sumWeights_nominal",&sumWeights_nominal,"sumWeights_nominal/F");
                // v_tOutput[itree]->Branch("sumWeights_scale_originalXWGTUP",&sumWeights_scale_originalXWGTUP,"sumWeights_scale_originalXWGTUP/F");
                // v_tOutput[itree]->Branch("sumWeights_scale_muF0p5",&sumWeights_scale_muF0p5,"sumWeights_scale_muF0p5/F");
                // v_tOutput[itree]->Branch("sumWeights_scale_muF2",&sumWeights_scale_muF2,"sumWeights_scale_muF2/F");
                // v_tOutput[itree]->Branch("sumWeights_scale_muR0p5",&sumWeights_scale_muR0p5,"sumWeights_scale_muR0p5/F");
                // v_tOutput[itree]->Branch("sumWeights_scale_muR2",&sumWeights_scale_muR2,"sumWeights_scale_muR2/F");
                // v_tOutput[itree]->Branch("sumWeights_scale_muR2muF2",&sumWeights_scale_muR2muF2,"sumWeights_scale_muR2muF2/F");
                // v_tOutput[itree]->Branch("sumWeights_scale_muR0p5muF0p5",&sumWeights_scale_muR0p5muF0p5,"sumWeights_scale_muR0p5muF0p5/F");

                v_tOutput[itree]->Branch("weight_originalXWGTUP",&weight_originalXWGTUP,"weight_originalXWGTUP/F");
                v_tOutput[itree]->Branch("weight_scale_muF0p5",&weight_scale_muF0p5,"weight_scale_muF0p5/F");
                v_tOutput[itree]->Branch("weight_scale_muF2",&weight_scale_muF2,"weight_scale_muF2/F");
                v_tOutput[itree]->Branch("weight_scale_muR0p5",&weight_scale_muR0p5,"weight_scale_muR0p5/F");
                v_tOutput[itree]->Branch("weight_scale_muR2",&weight_scale_muR2,"weight_scale_muR2/F");
                v_tOutput[itree]->Branch("weight_scale_muR2muF2",&weight_scale_muR2muF2,"weight_scale_muR2muF2/F");
                v_tOutput[itree]->Branch("weight_scale_muR0p5muF0p5",&weight_scale_muR0p5muF0p5,"weight_scale_muR0p5muF0p5/F");
            } //if save scale variations independently
        } //don't save if JES/JER ttree
    } //don't save if ntuple is for overlap studies only


	//-- Other, MEM necessary vars, ...
    if(write_branches_forMEM && !itree) //only for nominal tree, if user asks for it
    {
        //-- NEW : added vector of LooseJets, so that Jets are auto selected within MEM code
        // v_tOutput[itree]->Branch("nJets", &nJets, "nJets/I");
        v_tOutput[itree]->Branch("JetsPt", &JetsPt);
        v_tOutput[itree]->Branch("JetsEta", &JetsEta);
        v_tOutput[itree]->Branch("JetsPhi", &JetsPhi);
        v_tOutput[itree]->Branch("JetsE", &JetsE);
        v_tOutput[itree]->Branch("JetsCSV", &JetsCSV);

        v_tOutput[itree]->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/B");
        v_tOutput[itree]->Branch("catJets",&catJets,"catJets/I");

        v_tOutput[itree]->Branch("multilepton_Lepton1_Id",               &multilepton_Lepton1_Id,                "multilepton_Lepton1_Id/I");
        v_tOutput[itree]->Branch("multilepton_Lepton1_P4",               "TLorentzVector",                       &multilepton_Lepton1_P4);
        v_tOutput[itree]->Branch("multilepton_Lepton1_DeltaR_Matched",   &multilepton_Lepton1_DeltaR_Matched,    "multilepton_Lepton1_DeltaR_Matched/F");
        v_tOutput[itree]->Branch("multilepton_Lepton1_Label_Matched",    &multilepton_Lepton1_Label_Matched,     "multilepton_Lepton1_Label_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Lepton1_Id_Matched",       &multilepton_Lepton1_Id_Matched,        "multilepton_Lepton1_Id_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Lepton1_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton1_P4_Matched);
        v_tOutput[itree]->Branch("multilepton_Lepton2_Id",               &multilepton_Lepton2_Id,                "multilepton_Lepton2_Id/I");
        v_tOutput[itree]->Branch("multilepton_Lepton2_P4",               "TLorentzVector",                       &multilepton_Lepton2_P4);
        v_tOutput[itree]->Branch("multilepton_Lepton2_DeltaR_Matched",   &multilepton_Lepton2_DeltaR_Matched,    "multilepton_Lepton2_DeltaR_Matched/F");
        v_tOutput[itree]->Branch("multilepton_Lepton2_Label_Matched",    &multilepton_Lepton2_Label_Matched,     "multilepton_Lepton2_Label_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Lepton2_Id_Matched",       &multilepton_Lepton2_Id_Matched,        "multilepton_Lepton2_Id_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Lepton2_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton2_P4_Matched);
        v_tOutput[itree]->Branch("multilepton_Lepton3_Id",               &multilepton_Lepton3_Id,                "multilepton_Lepton3_Id/I");
        v_tOutput[itree]->Branch("multilepton_Lepton3_P4",               "TLorentzVector",                       &multilepton_Lepton3_P4);
        v_tOutput[itree]->Branch("multilepton_Lepton3_DeltaR_Matched",   &multilepton_Lepton3_DeltaR_Matched,    "multilepton_Lepton3_DeltaR_Matched/F");
        v_tOutput[itree]->Branch("multilepton_Lepton3_Label_Matched",    &multilepton_Lepton3_Label_Matched,     "multilepton_Lepton3_Label_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Lepton3_Id_Matched",       &multilepton_Lepton3_Id_Matched,        "multilepton_Lepton3_Id_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Lepton3_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton3_P4_Matched);
        v_tOutput[itree]->Branch("multilepton_Lepton4_Id",               &multilepton_Lepton4_Id,                "multilepton_Lepton4_Id/I");
        v_tOutput[itree]->Branch("multilepton_Lepton4_P4",               "TLorentzVector",                       &multilepton_Lepton4_P4);
        v_tOutput[itree]->Branch("multilepton_Lepton4_DeltaR_Matched",   &multilepton_Lepton4_DeltaR_Matched,    "multilepton_Lepton4_DeltaR_Matched/F");
        v_tOutput[itree]->Branch("multilepton_Lepton4_Label_Matched",    &multilepton_Lepton4_Label_Matched,     "multilepton_Lepton4_Label_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Lepton4_Id_Matched",       &multilepton_Lepton4_Id_Matched,        "multilepton_Lepton4_Id_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Lepton4_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton4_P4_Matched);

        v_tOutput[itree]->Branch("multilepton_Bjet1_Id",                 &multilepton_Bjet1_Id,                  "multilepton_Bjet1_Id/I");
        v_tOutput[itree]->Branch("multilepton_Bjet1_P4",                 "TLorentzVector",                       &multilepton_Bjet1_P4);
        v_tOutput[itree]->Branch("multilepton_Bjet1_CSV",                &multilepton_Bjet1_CSV,                 "multilepton_Bjet1_CSV/F");
        v_tOutput[itree]->Branch("multilepton_Bjet1_JEC_Up",             &multilepton_Bjet1_JEC_Up,              "multilepton_Bjet1_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_Bjet1_JEC_Down",           &multilepton_Bjet1_JEC_Down,            "multilepton_Bjet1_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_Bjet1_JER_Up",             &multilepton_Bjet1_JER_Up,              "multilepton_Bjet1_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_Bjet1_JER_Down",           &multilepton_Bjet1_JER_Down,            "multilepton_Bjet1_JER_Down/F");
        v_tOutput[itree]->Branch("multilepton_Bjet1_DeltaR_Matched",     &multilepton_Bjet1_DeltaR_Matched,      "multilepton_Bjet1_DeltaR_Matched/F");
        v_tOutput[itree]->Branch("multilepton_Bjet1_Label_Matched",      &multilepton_Bjet1_Label_Matched,       "multilepton_Bjet1_Label_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Bjet1_Id_Matched",         &multilepton_Bjet1_Id_Matched,          "multilepton_Bjet1_Id_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Bjet1_P4_Matched",         "TLorentzVector",                       &multilepton_Bjet1_P4_Matched);

        v_tOutput[itree]->Branch("multilepton_Bjet2_Id",                 &multilepton_Bjet2_Id,                  "multilepton_Bjet2_Id/I");
        v_tOutput[itree]->Branch("multilepton_Bjet2_P4",                 "TLorentzVector",                       &multilepton_Bjet2_P4);
        v_tOutput[itree]->Branch("multilepton_Bjet2_CSV",                &multilepton_Bjet2_CSV,                 "multilepton_Bjet2_CSV/F");
        v_tOutput[itree]->Branch("multilepton_Bjet2_JEC_Up",             &multilepton_Bjet2_JEC_Up,              "multilepton_Bjet2_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_Bjet2_JEC_Down",           &multilepton_Bjet2_JEC_Down,            "multilepton_Bjet2_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_Bjet2_JER_Up",             &multilepton_Bjet2_JER_Up,              "multilepton_Bjet2_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_Bjet2_JER_Down",           &multilepton_Bjet2_JER_Down,            "multilepton_Bjet2_JER_Down/F");
        v_tOutput[itree]->Branch("multilepton_Bjet2_DeltaR_Matched",     &multilepton_Bjet2_DeltaR_Matched,      "multilepton_Bjet2_DeltaR_Matched/F");
        v_tOutput[itree]->Branch("multilepton_Bjet2_Label_Matched",      &multilepton_Bjet2_Label_Matched,       "multilepton_Bjet2_Label_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Bjet2_Id_Matched",         &multilepton_Bjet2_Id_Matched,          "multilepton_Bjet2_Id_Matched/I");
        v_tOutput[itree]->Branch("multilepton_Bjet2_P4_Matched",         "TLorentzVector",                       &multilepton_Bjet2_P4_Matched);

        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_CSV",&multilepton_JetHighestPt1_CSV,"multilepton_JetHighestPt1_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_JEC_Up",&multilepton_JetHighestPt1_JEC_Up,"multilepton_JetHighestPt1_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_JEC_Down",&multilepton_JetHighestPt1_JEC_Down,"multilepton_JetHighestPt1_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_JER_Up",&multilepton_JetHighestPt1_JER_Up,"multilepton_JetHighestPt1_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_JER_Down",&multilepton_JetHighestPt1_JER_Down,"multilepton_JetHighestPt1_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_CSV",&multilepton_JetHighestPt2_CSV,"multilepton_JetHighestPt2_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_JEC_Up",&multilepton_JetHighestPt2_JEC_Up,"multilepton_JetHighestPt2_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_JEC_Down",&multilepton_JetHighestPt2_JEC_Down,"multilepton_JetHighestPt2_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_JER_Up",&multilepton_JetHighestPt2_JER_Up,"multilepton_JetHighestPt2_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_JER_Down",&multilepton_JetHighestPt2_JER_Down,"multilepton_JetHighestPt2_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_CSV",&multilepton_JetClosestMw1_CSV,"multilepton_JetClosestMw1_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_JEC_Up",&multilepton_JetClosestMw1_JEC_Up,"multilepton_JetClosestMw1_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_JEC_Down",&multilepton_JetClosestMw1_JEC_Down,"multilepton_JetClosestMw1_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_JER_Up",&multilepton_JetClosestMw1_JER_Up,"multilepton_JetClosestMw1_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_JER_Down",&multilepton_JetClosestMw1_JER_Down,"multilepton_JetClosestMw1_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_CSV",&multilepton_JetClosestMw2_CSV,"multilepton_JetClosestMw2_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_JEC_Up",&multilepton_JetClosestMw2_JEC_Up,"multilepton_JetClosestMw2_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_JEC_Down",&multilepton_JetClosestMw2_JEC_Down,"multilepton_JetClosestMw2_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_JER_Up",&multilepton_JetClosestMw2_JER_Up,"multilepton_JetClosestMw2_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_JER_Down",&multilepton_JetClosestMw2_JER_Down,"multilepton_JetClosestMw2_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_CSV",&multilepton_JetLowestMjj1_CSV,"multilepton_JetLowestMjj1_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_JEC_Up",&multilepton_JetLowestMjj1_JEC_Up,"multilepton_JetLowestMjj1_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_JEC_Down",&multilepton_JetLowestMjj1_JEC_Down,"multilepton_JetLowestMjj1_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_JER_Up",&multilepton_JetLowestMjj1_JER_Up,"multilepton_JetLowestMjj1_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_JER_Down",&multilepton_JetLowestMjj1_JER_Down,"multilepton_JetLowestMjj1_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_CSV",&multilepton_JetLowestMjj2_CSV,"multilepton_JetLowestMjj2_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_JEC_Up",&multilepton_JetLowestMjj2_JEC_Up,"multilepton_JetLowestMjj2_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_JEC_Down",&multilepton_JetLowestMjj2_JEC_Down,"multilepton_JetLowestMjj2_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_JER_Up",&multilepton_JetLowestMjj2_JER_Up,"multilepton_JetLowestMjj2_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_JER_Down",&multilepton_JetLowestMjj2_JER_Down,"multilepton_JetLowestMjj2_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetHighestEta1_Id",&multilepton_JetHighestEta1_Id,"multilepton_JetHighestEta1_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta1_P4","TLorentzVector",&multilepton_JetHighestEta1_P4);
        v_tOutput[itree]->Branch("multilepton_JetHighestEta1_CSV",&multilepton_JetHighestEta1_CSV,"multilepton_JetHighestEta1_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta2_Id",&multilepton_JetHighestEta2_Id,"multilepton_JetHighestEta2_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta2_P4","TLorentzVector",&multilepton_JetHighestEta2_P4);
        v_tOutput[itree]->Branch("multilepton_JetHighestEta2_CSV",&multilepton_JetHighestEta2_CSV,"multilepton_JetHighestEta2_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta1_JEC_Up",&multilepton_JetHighestEta1_JEC_Up,"multilepton_JetHighestEta1_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta1_JEC_Down",&multilepton_JetHighestEta1_JEC_Down,"multilepton_JetHighestEta1_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta1_JER_Up",&multilepton_JetHighestEta1_JER_Up,"multilepton_JetHighestEta1_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta1_JER_Down",&multilepton_JetHighestEta1_JER_Down,"multilepton_JetHighestEta1_JER_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta2_JEC_Up",&multilepton_JetHighestEta2_JEC_Up,"multilepton_JetHighestEta2_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta2_JEC_Down",&multilepton_JetHighestEta2_JEC_Down,"multilepton_JetHighestEta2_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta2_JER_Up",&multilepton_JetHighestEta2_JER_Up,"multilepton_JetHighestEta2_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestEta2_JER_Down",&multilepton_JetHighestEta2_JER_Down,"multilepton_JetHighestEta2_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_2ndPair_Id",&multilepton_JetHighestPt1_2ndPair_Id,"multilepton_JetHighestPt1_2ndPair_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt1_2ndPair_P4);
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_2ndPair_CSV",&multilepton_JetHighestPt1_2ndPair_CSV,"multilepton_JetHighestPt1_2ndPair_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Up",&multilepton_JetHighestPt1_2ndPair_JEC_Up,"multilepton_JetHighestPt1_2ndPair_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Down",&multilepton_JetHighestPt1_2ndPair_JEC_Down,"multilepton_JetHighestPt1_2ndPair_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_2ndPair_JER_Up",&multilepton_JetHighestPt1_2ndPair_JER_Up,"multilepton_JetHighestPt1_2ndPair_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt1_2ndPair_JER_Down",&multilepton_JetHighestPt1_2ndPair_JER_Down,"multilepton_JetHighestPt1_2ndPair_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_2ndPair_Id",&multilepton_JetHighestPt2_2ndPair_Id,"multilepton_JetHighestPt2_2ndPair_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt2_2ndPair_P4);
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_2ndPair_CSV",&multilepton_JetHighestPt2_2ndPair_CSV,"multilepton_JetHighestPt2_2ndPair_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Up",&multilepton_JetHighestPt2_2ndPair_JEC_Up,"multilepton_JetHighestPt2_2ndPair_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Down",&multilepton_JetHighestPt2_2ndPair_JEC_Down,"multilepton_JetHighestPt2_2ndPair_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_2ndPair_JER_Up",&multilepton_JetHighestPt2_2ndPair_JER_Up,"multilepton_JetHighestPt2_2ndPair_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetHighestPt2_2ndPair_JER_Down",&multilepton_JetHighestPt2_2ndPair_JER_Down,"multilepton_JetHighestPt2_2ndPair_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_2ndPair_Id",&multilepton_JetClosestMw1_2ndPair_Id,"multilepton_JetClosestMw1_2ndPair_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw1_2ndPair_P4);
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_2ndPair_CSV",&multilepton_JetClosestMw1_2ndPair_CSV,"multilepton_JetClosestMw1_2ndPair_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Up",&multilepton_JetClosestMw1_2ndPair_JEC_Up,"multilepton_JetClosestMw1_2ndPair_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Down",&multilepton_JetClosestMw1_2ndPair_JEC_Down,"multilepton_JetClosestMw1_2ndPair_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_2ndPair_JER_Up",&multilepton_JetClosestMw1_2ndPair_JER_Up,"multilepton_JetClosestMw1_2ndPair_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw1_2ndPair_JER_Down",&multilepton_JetClosestMw1_2ndPair_JER_Down,"multilepton_JetClosestMw1_2ndPair_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_2ndPair_Id",&multilepton_JetClosestMw2_2ndPair_Id,"multilepton_JetClosestMw2_2ndPair_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw2_2ndPair_P4);
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_2ndPair_CSV",&multilepton_JetClosestMw2_2ndPair_CSV,"multilepton_JetClosestMw2_2ndPair_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Up",&multilepton_JetClosestMw2_2ndPair_JEC_Up,"multilepton_JetClosestMw2_2ndPair_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Down",&multilepton_JetClosestMw2_2ndPair_JEC_Down,"multilepton_JetClosestMw2_2ndPair_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_2ndPair_JER_Up",&multilepton_JetClosestMw2_2ndPair_JER_Up,"multilepton_JetClosestMw2_2ndPair_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetClosestMw2_2ndPair_JER_Down",&multilepton_JetClosestMw2_2ndPair_JER_Down,"multilepton_JetClosestMw2_2ndPair_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_2ndPair_Id",&multilepton_JetLowestMjj1_2ndPair_Id,"multilepton_JetLowestMjj1_2ndPair_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj1_2ndPair_P4);
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_2ndPair_CSV",&multilepton_JetLowestMjj1_2ndPair_CSV,"multilepton_JetLowestMjj1_2ndPair_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Up",&multilepton_JetLowestMjj1_2ndPair_JEC_Up,"multilepton_JetLowestMjj1_2ndPair_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Down",&multilepton_JetLowestMjj1_2ndPair_JEC_Down,"multilepton_JetLowestMjj1_2ndPair_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Up",&multilepton_JetLowestMjj1_2ndPair_JER_Up,"multilepton_JetLowestMjj1_2ndPair_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Down",&multilepton_JetLowestMjj1_2ndPair_JER_Down,"multilepton_JetLowestMjj1_2ndPair_JER_Down/F");

        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_2ndPair_Id",&multilepton_JetLowestMjj2_2ndPair_Id,"multilepton_JetLowestMjj2_2ndPair_Id/I");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj2_2ndPair_P4);
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_2ndPair_CSV",&multilepton_JetLowestMjj2_2ndPair_CSV,"multilepton_JetLowestMjj2_2ndPair_CSV/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Up",&multilepton_JetLowestMjj2_2ndPair_JEC_Up,"multilepton_JetLowestMjj2_2ndPair_JEC_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Down",&multilepton_JetLowestMjj2_2ndPair_JEC_Down,"multilepton_JetLowestMjj2_2ndPair_JEC_Down/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Up",&multilepton_JetLowestMjj2_2ndPair_JER_Up,"multilepton_JetLowestMjj2_2ndPair_JER_Up/F");
        v_tOutput[itree]->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Down",&multilepton_JetLowestMjj2_2ndPair_JER_Down,"multilepton_JetLowestMjj2_2ndPair_JER_Down/F");

        // Test adding truth information

        v_tOutput[itree]->Branch("multilepton_h0_Id",                    &multilepton_h0_Id,                 "multilepton_h0_Id/I");
        v_tOutput[itree]->Branch("multilepton_h0_P4",                    "TLorentzVector",                   &multilepton_h0_P4);
        v_tOutput[itree]->Branch("multilepton_h0_Label",                    &multilepton_h0_Label,                 "multilepton_h0_Label/I");
        v_tOutput[itree]->Branch("multilepton_t1_Id",                    &multilepton_t1_Id,                 "multilepton_t1_Id/I");
        v_tOutput[itree]->Branch("multilepton_t1_P4",                    "TLorentzVector",                   &multilepton_t1_P4);
        v_tOutput[itree]->Branch("multilepton_t1_Label",                    &multilepton_t1_Label,                 "multilepton_h0_Label/I");
        v_tOutput[itree]->Branch("multilepton_t2_Id",                    &multilepton_t2_Id,                 "multilepton_t2_Id/I");
        v_tOutput[itree]->Branch("multilepton_t2_P4",                    "TLorentzVector",                   &multilepton_t2_P4);
        v_tOutput[itree]->Branch("multilepton_t2_Label",                    &multilepton_t2_Label,                 "multilepton_h0_Label/I");

        // End test adding truth information

        v_tOutput[itree]->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
        v_tOutput[itree]->Branch("multilepton_mETcov00",&multilepton_mETcov00,"multilepton_mETcov00/D");
        v_tOutput[itree]->Branch("multilepton_mETcov01",&multilepton_mETcov01,"multilepton_mETcov01/D");
        v_tOutput[itree]->Branch("multilepton_mETcov10",&multilepton_mETcov10,"multilepton_mETcov10/D");
        v_tOutput[itree]->Branch("multilepton_mETcov11",&multilepton_mETcov11,"multilepton_mETcov11/D");
        v_tOutput[itree]->Branch("multilepton_sumET",&multilepton_sumET,"multilepton_sumET/F");
        v_tOutput[itree]->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);
    }

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




void tHqMultileponAnalysis::FillJetInfoOutputTree(int* tree_Id, int Id, TLorentzVector* tree_P4, TLorentzVector P4, float* tree_CSV, float CSV, float* tree_JEC_Up, float* tree_JEC_Down, float JEC_value, float* tree_JER_Up, float* tree_JER_Down, float JER, float JERUp, float JERDown)
{
    *tree_Id = Id;
    *tree_P4 = P4;

    *tree_CSV = CSV;

    *tree_JEC_Up = P4.E()*(1.+JEC_value);
    *tree_JEC_Down = P4.E()*(1.-JEC_value);

    *tree_JER_Up = P4.E()*JERUp/JER;
    *tree_JER_Down = P4.E()*JERDown/JER;

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
    for (int ib=0; ib<vJetLoose.size(); ib++)

    {
        if(doSelectOnlyBjets && !vJetLoose.at(ib).isLooseBTag) {continue;}

        if(vJetLoose.at(ib).DeepCSVbtag>btag_max)
        {
          btag_max2 = btag_max;
          ib2 = ib1;
          btag_max = vJetLoose.at(ib).DeepCSVbtag;
          ib1 = ib;
        }
        else if (vJetLoose.at(ib).DeepCSVbtag<btag_max && vJetLoose.at(ib).DeepCSVbtag>btag_max2)
        {
          btag_max2 = vJetLoose.at(ib).DeepCSVbtag;
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

    for(int ij=0; ij<vJetLoose.size(); ij++)
    {
        if (ij==ib1 || ij==ib2) {continue;} //Don't take bjets into account

        if(vJetLoose.at(ij).pt > pt_max ) //Highest pT
        {
            pt_max2 = pt_max;
            ij2 = ij1;
            pt_max = vJetLoose.at(ij).pt;
            ij1 = ij;
        }
        else if(vJetLoose.at(ij).pt > pt_max2) //2nd Highest pT
        {
            pt_max2 = vJetLoose.at(ij).pt;
            ij2 = ij;
        }

        if(fabs(vJetLoose.at(ij).eta) > eta_max ) //Highest eta
        {
            eta_max2 = eta_max;
            ie2 = ie1;
            eta_max = fabs(vJetLoose.at(ij).eta);
            ie1 = ij;
        }
        else if(fabs(vJetLoose.at(ij).eta ) > eta_max2) //2nd Highest eta
        {
            eta_max2 = fabs(vJetLoose.at(ij).eta);
            ie2 = ij;
        }

        for(int ik=ij+1; ik<vJetLoose.size(); ik++) //Dijet w/ lowest (m - mW)
        {
            if (ik==ib1 || ik==ib2) {continue;} //Don't take bjets into account

            Pjet1.SetPtEtaPhiE(vJetLoose.at(ij).pt, vJetLoose.at(ij).eta, vJetLoose.at(ij).phi, vJetLoose.at(ij).E);
            Pjet2.SetPtEtaPhiE(vJetLoose.at(ik).pt, vJetLoose.at(ik).eta, vJetLoose.at(ik).phi, vJetLoose.at(ik).E);

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
    for(int ij=0; ij<vJetLoose.size(); ij++)
    {
        if (ij==ib1 || ij==ib2 || ij==ik1 || ij==ik2) {continue;} //Don't take bjets and 1rst mW pair into account

        if(vJetLoose.at(ij).pt > pt_max ) //Highest pT
        {
            pt_max2 = pt_max;
            im2 = im1;
            pt_max = vJetLoose.at(ij).pt;
            im1 = ij;
        }
        else if(vJetLoose.at(ij).pt > pt_max2) //2nd Highest pT
        {
            pt_max2 = vJetLoose.at(ij).pt;
            im2 = ij;
        }

        for(int ik=ij+1; ik<vJetLoose.size(); ik++) //Dijet w/ lowest (m - mW)
        {
            if (ik==ib1 || ik==ib2 || ik==ik1 || ik==ik2) {continue;} //Don't take bjets and 1rst pair into account

            Pjet1.SetPtEtaPhiE(vJetLoose.at(ij).pt, vJetLoose.at(ij).eta, vJetLoose.at(ij).phi, vJetLoose.at(ij).E);
            Pjet2.SetPtEtaPhiE(vJetLoose.at(ik).pt, vJetLoose.at(ik).eta, vJetLoose.at(ik).phi, vJetLoose.at(ik).E);

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
    if(_sampleName.Contains("ttZ", TString::kIgnoreCase) || (_sampleName.Contains("ttW", TString::kIgnoreCase) && !_sampleName.Contains("ttWH", TString::kIgnoreCase) && !_sampleName.Contains("ttWW", TString::kIgnoreCase)) || _sampleName.Contains("THQ", TString::kIgnoreCase) || _sampleName.Contains("THW", TString::kIgnoreCase) || _sampleName.Contains("ttH", TString::kIgnoreCase) || _sampleName.Contains("TTJet", TString::kIgnoreCase) || _sampleName.Contains("TTTo2", TString::kIgnoreCase) || _sampleName.Contains("TTToSemi", TString::kIgnoreCase) || _sampleName.Contains("TTbar", TString::kIgnoreCase) || _sampleName.Contains("FCNC", TString::kIgnoreCase) || _sampleName.Contains("Tprime", TString::kIgnoreCase) )
    {
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

/**
 * mHT: vector sum of fakeable leptons, taus and selected jets
 * Worse resolution but more robust against pile-up
 */
double tHqMultileponAnalysis::Compute_mHT()
{
    TLorentzVector obj;
    float obj_px = 0;
    float obj_py = 0;
    for(int ij=0;ij<vJetLoose.size();ij++)
    {
        obj.SetPtEtaPhiE(vJetLoose.at(ij).pt, vJetLoose.at(ij).eta, vJetLoose.at(ij).phi, vJetLoose.at(ij).E);
        obj_px += obj.Px();
        obj_py += obj.Py();
    }
    for(int ij=0;ij<vLeptonFakeable.size();ij++)
    {
        obj.SetPtEtaPhiE(vLeptonFakeable.at(ij).pt, vLeptonFakeable.at(ij).eta, vLeptonFakeable.at(ij).phi, vLeptonFakeable.at(ij).E);
        obj_px += obj.Px();
        obj_py += obj.Py();
    }
    for(int ij=0;ij<vTauLoose->size();ij++)
    {
        obj.SetPtEtaPhiE(vTauLoose->at(ij).pt, vTauLoose->at(ij).eta, vTauLoose->at(ij).phi, vTauLoose->at(ij).E);
        obj_px += obj.Px();
        obj_py += obj.Py();
    }

    return sqrt( (obj_px*obj_px) + (obj_py*obj_py) );
}





/**
 * Try to co"mpute some high-level reco variables related to W boson or top
 * //-- could elaborate algorithm by adding result from Hj tagger ?
 */
void tHqMultileponAnalysis::Compute_Top_W_variables()
{
    double mW = 80.4;
    double mTop = 173;

    int i_lep_fromTopW=-1;
    int i_bJet_fromTop=-1;
    int i_jet_recoiling=-1;

    TLorentzVector lep_fromTopW; int charge_lepFromTopW;
    TLorentzVector Bjet_fromTop;
    TLorentzVector RecoilJet;
    TLorentzVector Neutrino;

    double tmp = 999; double mass_LepBJet_min_tmp = 999;
    for(int ilep=0; ilep<vLeptonFakeable.size(); ilep++)
    {
        // cout<<endl<<"ilep = "<<ilep<<endl;

        TLorentzVector lep_tmp;
        lep_tmp.SetPtEtaPhiE(vLeptonFakeable.at(ilep).conept, vLeptonFakeable.at(ilep).eta, vLeptonFakeable.at(ilep).phi, vLeptonFakeable.at(ilep).E );

        double neutrino_Px, neutrino_Py, neutrino_Pz, neutrino_E, k, a, b, c;
        neutrino_Px = metpt*sin(metphi);
        neutrino_Py = metpt*cos(metphi);
        a = lep_tmp.Pt()*lep_tmp.Pt();

        vector<double> solutions_poly(2);
        double metpt_tmp = metpt;
        bool found_solution = false;
        while(!found_solution)
        {
            k = (mW*mW/2) + lep_tmp.Pt() * metpt_tmp;
            b = -2 * lep_tmp.Pz() * k;
            c = metpt_tmp*metpt_tmp * fabs(lep_tmp.E()*lep_tmp.E()) - k*k;

            found_solution = Get_Polynom2_Root(a,b,c, solutions_poly);

            if(!found_solution) // -- how to deal w/ imaginary sols ?
            {
                // cout<<"//--------------------------------------------"<<endl;
                // cout<<"vLeptonFakeable.at(ilep).conept = "<<vLeptonFakeable.at(ilep).conept<<endl;
                // cout<<"vLeptonFakeable.at(ilep).pt = "<<vLeptonFakeable.at(ilep).pt<<endl;

                // metpt_tmp-= 5;
                // cout<<"No solution found for neutrino reco ! Try : MET - 5GEV = "<<metpt_tmp<<endl;
            }
        }

        int n_sol = 1;
        if(solutions_poly[1] != 0) {n_sol = 2;}
        for(int isol=0; isol<n_sol; isol++)
        {
            neutrino_Pz = fabs(solutions_poly[isol]);
            // cout<<"neutrino_Pz = "<<neutrino_Pz<<endl;

            neutrino_E = sqrt(neutrino_Px*neutrino_Px+neutrino_Py*neutrino_Py+neutrino_Pz*neutrino_Pz);

            TLorentzVector neutrino;
            neutrino.SetPxPyPzE(neutrino_Px, neutrino_Py, neutrino_Pz, neutrino_E);

            for(int ijet=0; ijet<vJetLoose.size(); ijet++)
            {
                // cout<<endl<<"ijet = "<<ijet<<endl;

                TLorentzVector Bjet_tmp;

                // if(vLooseBTagJets.at(ijet).DeepCSVbtag > deepCSV_max) {i_bJet_fromTop = ijet;}
                if(!vJetLoose.at(ijet).isMediumBTag) {continue;} //Only consider b-tagged jets

                Bjet_tmp.SetPtEtaPhiE(vJetLoose.at(ijet).pt, vJetLoose.at(ijet).eta, vJetLoose.at(ijet).phi, vJetLoose.at(ijet).E);

                double mtop_tmp = (lep_tmp+neutrino+Bjet_tmp).M();
                // cout<<"mtop_tmp = "<<mtop_tmp<<endl;

                if(fabs(mtop_tmp - mTop) < tmp) //If config is closest to mtop, store values/indices
                {
                    tmp = fabs(mtop_tmp - mTop);
                    top_mass = mtop_tmp;
                    top_Pt = mtop_tmp;
                    i_lep_fromTopW = ilep; charge_lepFromTopW = vLeptonFakeable.at(ilep).charge;
                    i_bJet_fromTop = ijet;

                    Neutrino.SetPxPyPzE(neutrino_Px, neutrino_Py, neutrino_Pz, neutrino_E);
                }

                //Also compute minimal inv mass of (lep+bjet) system
                mass_LepBJet_min_tmp = (lep_tmp + Bjet_tmp).M();
                if(mass_LepBJet_min_tmp < mass_LepBJet_min)
                {
                    mass_LepBJet_min = mass_LepBJet_min_tmp;
                }
            } //loop on bjets
        } //loop on neutrino solutions
    } //loop on leptons

    tmp = -999;
    for(int ijet=0; ijet<vJetLoose.size(); ijet++)
    {
        if(ijet == i_bJet_fromTop) {continue;}

        if(vJetLoose.at(ijet).pt > tmp) {i_jet_recoiling = ijet;}
    }

    if(i_lep_fromTopW < 0 || i_bJet_fromTop < 0 || i_jet_recoiling < 0) {return;}
    lep_fromTopW.SetPtEtaPhiE(vLeptonFakeable.at(i_lep_fromTopW).conept, vLeptonFakeable.at(i_lep_fromTopW).eta, vLeptonFakeable.at(i_lep_fromTopW).phi, vLeptonFakeable.at(i_lep_fromTopW).E );
    Bjet_fromTop.SetPtEtaPhiE(vJetLoose.at(i_bJet_fromTop).pt, vJetLoose.at(i_bJet_fromTop).eta, vJetLoose.at(i_bJet_fromTop).phi, vJetLoose.at(i_bJet_fromTop).E);
    RecoilJet.SetPtEtaPhiE(vJetLoose.at(i_jet_recoiling).pt, vJetLoose.at(i_jet_recoiling).eta, vJetLoose.at(i_jet_recoiling).phi, vJetLoose.at(i_jet_recoiling).E);

    mTW = sqrt(2*lep_fromTopW.Pt()*metpt*(1-cos(Phi_MPi_Pi(lep_fromTopW.Phi()-metphi)) ) );
    // cout<<"mTW = "<<mTW<<endl;

    lW_asym_mtop = charge_lepFromTopW * fabs(lep_fromTopW.Eta() );

    dRBjetRecoilJet = GetDeltaR(Bjet_fromTop.Eta(), Bjet_fromTop.Phi(), RecoilJet.Eta(), RecoilJet.Phi());
    dRLepWRecoilJet = GetDeltaR(lep_fromTopW.Eta(), lep_fromTopW.Phi(), RecoilJet.Eta(), RecoilJet.Phi());

    RecoilJetPt = RecoilJet.Pt();
    RecoilJetEta = RecoilJet.Eta();
    LepWPt = lep_fromTopW.Pt();
    LepWEta = lep_fromTopW.Eta();

    return;
}





/**
 * Return smaller root of order 2 polynom, when solution exists
 * //NB - "In case of an imaginary solution, the p Tmiss value is varied until one real solutionis found." ==> correct ??
 */
bool tHqMultileponAnalysis::Get_Polynom2_Root(double a, double b, double c, vector<double>& result)
{
    result.resize(2);
    result[0] = 0;
    result[1] = 0;

    double delta = b*b - 4*a*c;

    if(delta > 0)
    {
        double res1 = (-b - sqrt(delta)) / (2*a);
        double res2 = (-b + sqrt(delta)) / (2*a);

        result[0] = res1;
        result[1] = res2;

        // return (fabs(res1)<fabs(res2)) ? res1 : res2;
        return true;
    }
    else if(delta == 0) {result[0] = -b / (2*a); return true;} //1 solution
    else if(delta < 0) //Imaginary solutions //NB : how to deal with it ?
    {
        //For now, drop the imaginary part (see https://indico.cern.ch/event/181837/contributions/1449128/attachments/245620/343594/Reconstruction_of_the_Single_Top_v2.pdf)
        // return false;

        result[0] = fabs(-b / 2*a);

        return true;
    }

    return false;
}


/**
 * Combine metpt and mHT into a metpt linear discriminator “metLD”
 * metpt and mHT are less correlated in events due to instrumental missing energy compared to the real one
 * From PreApp slides : metLD = 0.6 * metpt+ 0.4 * mHT
 * HERE : using alternative metLD def, used by ttH2017 for training ?? to be checked, don't use now
 */
double tHqMultileponAnalysis::Compute_metLD_Alternative(double metpt)
{
    double mHT = Compute_mHT();

    return metpt*0.00397 + mHT*0.00265; //Alternative def, used for training sel ?
    // return metpt*0.6 + mHT*0.4;
}






double tHqMultileponAnalysis::Get_Distance(double rec_pt, double rec_eta, double rec_phi, double gen_pt, double gen_eta, double gen_phi)
{
    double dr = GetDeltaR(rec_eta, rec_phi, gen_eta, gen_phi);
    double dptRel = fabs(rec_pt - gen_pt) / gen_pt;

    return dr + 0.2 * dptRel;
}




void tHqMultileponAnalysis::Dump_EventInfo_Synchro(std::ofstream& outfile)
{
    // cout<<setprecision(15)<<event_run<<":"<<event_lumi<<":"<<event_id<<endl;
    outfile<<setprecision(15)<<event_run<<":"<<event_lumi<<":"<<event_id<<endl;
}





/**
 * Read the ttH2017 predefined categories values (implemented at NTP level), and fill the "tHq" cat. with their values
 * This way we can compare the "ttH" cat. values between NTP/NTA codes
 */
void tHqMultileponAnalysis::Read_Predefined_Categories()
{
    //Can't convert directly boolean into Char_t
    //NB : using Char_t because booleans should not be used in vectors (later)
    if(vEvent->at(0).is_2lSS) {is_tHq_2lSS                           = 1;}
    if(vEvent->at(0).is_2lSS_Training) {is_tHq_2lSS_Training         = 1;}
    if(vEvent->at(0).is_2lSS_SR) {is_tHq_2lSS_SR                     = 1;}
    if(vEvent->at(0).is_2lSS_Fake) {is_tHq_2lSS_Fake                 = 1;}
    if(vEvent->at(0).is_2lSS_Flip) {is_tHq_2lSS_Flip                 = 1;}
    if(vEvent->at(0).is_3l) {is_tHq_3l                               = 1;}
    if(vEvent->at(0).is_3l_Training) {is_tHq_3l_Training             = 1;}
    if(vEvent->at(0).is_3l_SR) {is_tHq_3l_SR                         = 1;}
    if(vEvent->at(0).is_3l_Fake) {is_tHq_3l_Fake                     = 1;}
    if(vEvent->at(0).is_ttWctrl) {is_tHq_ttWctrl                     = 1;}
    if(vEvent->at(0).is_ttWctrl_SR) {is_tHq_ttWctrl_SR               = 1;}
    if(vEvent->at(0).is_ttWctrl_Fake) {is_tHq_ttWctrl_Fake           = 1;}
    if(vEvent->at(0).is_ttWctrl_Flip) {is_tHq_ttWctrl_Flip           = 1;}
    if(vEvent->at(0).is_ttZctrl) {is_tHq_ttZctrl                     = 1;}
    if(vEvent->at(0).is_ttZctrl_SR) {is_tHq_ttZctrl_SR               = 1;}
    if(vEvent->at(0).is_ttZctrl_Fake) {is_tHq_ttZctrl_Fake           = 1;}
    if(vEvent->at(0).is_WZctrl) {is_tHq_WZctrl                       = 1;}
    if(vEvent->at(0).is_WZctrl_SR) {is_tHq_WZctrl_SR                 = 1;}
    if(vEvent->at(0).is_WZctrl_Fake) {is_tHq_WZctrl_Fake             = 1;}
    if(vEvent->at(0).is_4l_SR) {is_tHq_4l_SR                         = 1;}
    if(vEvent->at(0).is_ZZctrl_SR) {is_tHq_ZZctrl_SR                 = 1;}

    return;
}





/**
 * The default "procedure" is to fill both the "tHq" cat. and "ttH" cat in this NTAnalysis code
 * However the "ttH" cat. are also filled at the NTProducer level (for sync, skim, etc)
 * This debugging function reads the "ttH" cat. from NTP files and fills the "tHq" booleans with their values
 * This allows to then compare the "ttH" cat. values from NTP and NTA codes
 */
void tHqMultileponAnalysis::InitPredefinedCategories()
{
    is_tHq_2lSS              = 0;
    is_tHq_2lSS_Training     = 0;
    is_tHq_2lSS_SR           = 0;
    is_tHq_2lSS_Fake         = 0;
    is_tHq_2lSS_Flip         = 0;
    is_tHq_3l                = 0;
    is_tHq_3l_Training       = 0;
    is_tHq_3l_SR             = 0;
    is_tHq_3l_Fake           = 0;
    is_tHq_ttWctrl           = 0;
    is_tHq_ttWctrl_SR        = 0;
    is_tHq_ttWctrl_Fake      = 0;
    is_tHq_ttWctrl_Flip      = 0;
    is_tHq_ttZctrl           = 0;
    is_tHq_ttZctrl_SR        = 0;
    is_tHq_ttZctrl_Fake      = 0;
    is_tHq_WZctrl            = 0;
    is_tHq_WZctrl_SR         = 0;
    is_tHq_WZctrl_Fake       = 0;
    is_tHq_4l_SR             = 0;
    is_tHq_ZZctrl_SR         = 0;

    return;
}




/**
 * Check if event should be saved or not, depending on the type of anlysis we're doing
 * If yes, returns the number of leptons to consider to compute the variables (as some 2l/3l categ may not be perfectly orthogonal, need to choose)
 */
int tHqMultileponAnalysis::Check_If_Save_Event(bool do_tHq_analysis)
{
    int nlep = 0;

    if(do_tHq_analysis)
    {
        if(is_tHq_2lSS_SR
        || is_tHq_2lSS_Training
        || is_tHq_2lSS_Fake
        || is_tHq_2lSS_Flip
        || is_tHq_2lSS_GammaConv
        || is_tHq_ttWctrl_SR
        || is_tHq_ttWctrl_Fake
        || is_tHq_ttWctrl_Flip
        || is_tHq_ttWctrl_GammaConv
        || is_tHqFCNC_2lSS_SR
        || is_tHqFCNC_2lSS_Fake
        || is_tHqFCNC_2lSS_Flip
        || is_tHqFCNC_2lSS_GammaConv) {nlep = 2;}
        else if(is_tHq_3l_SR
        || is_tHq_3l_Training
        || is_tHq_3l_Fake
        || is_tHq_3l_GammaConv
        || is_tHq_ttZctrl_SR
        || is_tHq_ttZctrl_Fake
        || is_tHq_ttZctrl_GammaConv
        || is_tHq_WZctrl_SR
        || is_tHq_WZctrl_Fake
        || is_tHq_WZctrl_GammaConv
        || is_tHqFCNC_3l_SR
        || is_tHqFCNC_3l_Fake
        || is_tHqFCNC_3l_GammaConv) {nlep = 3;}
        else if(is_tHq_4l_SR
        || is_tHq_ZZctrl_SR) {nlep = 4;}
    }
    else
    {
        if(is_ttH_2lSS_SR
        || is_ttH_2lSS_Training
        || is_ttH_2lSS_Fake
        || is_ttH_2lSS_Flip
        || is_ttH_2lSS_GammaConv
        || is_ttH_ttWctrl_SR
        || is_ttH_ttWctrl_Fake
        || is_ttH_ttWctrl_Flip
        || is_ttH_ttWctrl_GammaConv) {nlep = 2;}
        else if(is_ttH_3l_SR
        || is_ttH_3l_Training
        || is_ttH_3l_Fake
        || is_ttH_3l_GammaConv
        || is_ttH_ttZctrl_SR
        || is_ttH_ttZctrl_Fake
        || is_ttH_WZctrl_SR
        || is_ttH_WZctrl_Fake
        || is_ttH_WZctrl_GammaConv) {nlep = 3;}
        else if(is_ttH_4l_SR
        || is_ttH_ZZctrl_SR) {nlep = 4;}
    }

    if(make_ntuples_for_overlap_studies) //also check ortho categories (events should normally also belong to another categ, but...)
    {
        if(is_ttH_2lSS_SR_fwd
        || is_tHq_2lSS_SR_fwd

        || is_ttH_2lSS_SR_btag
        || is_tHq_2lSS_SR_btag

        || is_tHq_2lSS_SR_njet
        || is_tHq_2lSS_SR_njet
        || is_ttH_ttWctrl_SR_njet

        || is_tHq_2lSS_SR_fwd2
        || is_ttH_ttWctrl_SR_fwd2

        ) {nlep = 2;}

        else if(is_ttH_3l_SR_fwd
        || is_tHq_3l_SR_fwd

        || is_ttH_3l_SR_btag
        || is_tHq_3l_SR_btag

        || is_ttH_3l_SR_njet
        || is_tHq_3l_SR_njet

        || is_tHq_3l_SR_njet2
        || is_ttH_3l_SR_njet2

        || is_tHq_3l_SR_njet3
        || is_ttH_3l_SR_njet3

        ) {nlep = 3;}
    }

    return nlep;
}


/**
 * Fill overlap plot with event weight, depending on which tHq & ttH categories are activated
 */
void tHqMultileponAnalysis::Fill_Overlap_Histo()
{
    vector<Char_t> v_tHq_cat(0);
    // v_tHq_cat.push_back(is_tHq_2lSS);
    v_tHq_cat.push_back(is_tHq_2lSS_SR);
    // v_tHq_cat.push_back(is_tHq_2lSS_Training);
    // v_tHq_cat.push_back(is_tHq_2lSS_Fake);
    // v_tHq_cat.push_back(is_tHq_2lSS_Flip);
    // v_tHq_cat.push_back(is_tHq_2lSS_GammaConv);
    // v_tHq_cat.push_back(is_tHq_ttWctrl);
    v_tHq_cat.push_back(is_tHq_ttWctrl_SR);
    // v_tHq_cat.push_back(is_tHq_ttWctrl_Fake);
    // v_tHq_cat.push_back(is_tHq_ttWctrl_Flip);
    // v_tHq_cat.push_back(is_tHq_ttWctrl_GammaConv);
    // v_tHq_cat.push_back(is_ttH_3l);
    v_tHq_cat.push_back(is_tHq_3l_SR);
    // v_tHq_cat.push_back(is_tHq_3l_Training);
    // v_tHq_cat.push_back(is_tHq_3l_Fake);
    // v_tHq_cat.push_back(is_tHq_ttZctrl);
    v_tHq_cat.push_back(is_tHq_ttZctrl_SR);
    // v_tHq_cat.push_back(is_tHq_ttZctrl_Fake);
    // v_tHq_cat.push_back(is_tHq_WZctrl);
    v_tHq_cat.push_back(is_tHq_WZctrl_SR);
    // v_tHq_cat.push_back(is_tHq_WZctrl_Fake);

    vector<Char_t> v_ttH_cat(0);
    // v_ttH_cat.push_back(is_ttH_2lSS);
    v_ttH_cat.push_back(is_ttH_2lSS_SR);
    // v_ttH_cat.push_back(is_ttH_2lSS_Training);
    // v_ttH_cat.push_back(is_ttH_2lSS_Fake);
    // v_ttH_cat.push_back(is_ttH_2lSS_Flip);
    // v_ttH_cat.push_back(is_ttH_ttWctrl);
    v_ttH_cat.push_back(is_ttH_ttWctrl_SR);
    // v_ttH_cat.push_back(is_ttH_ttWctrl_Fake);
    // v_ttH_cat.push_back(is_ttH_ttWctrl_Flip);
    // v_ttH_cat.push_back(is_ttH_ttWctrl_GammaConv);
    // v_ttH_cat.push_back(is_ttH_3l);
    v_ttH_cat.push_back(is_ttH_3l_SR);
    // v_ttH_cat.push_back(is_ttH_3l_Training);
    // v_ttH_cat.push_back(is_ttH_3l_Fake);
    // v_ttH_cat.push_back(is_ttH_ttZctrl);
    v_ttH_cat.push_back(is_ttH_ttZctrl_SR);
    // v_ttH_cat.push_back(is_ttH_ttZctrl_Fake);
    // v_ttH_cat.push_back(is_ttH_WZctrl);
    v_ttH_cat.push_back(is_ttH_WZctrl_SR);
    // v_ttH_cat.push_back(is_ttH_WZctrl_Fake);

    //Store total yield of events entering each tHq cat
    for(int icat=0; icat<v_tHq_cat.size(); icat++)
    {
        if(v_tHq_cat[icat] == 1) //Add event yield to total
        {
            float tmp = h_totalYield_tHq_cat->GetBinContent(icat+1);
            h_totalYield_tHq_cat->SetBinContent(icat+1, tmp+weight);
        }
    }

    //Store total yield of events entering each ttH cat
    for(int icat=0; icat<v_ttH_cat.size(); icat++)
    {
        if(v_ttH_cat[icat] == 1) //Add event yield to total
        {
            float tmp = h_totalYield_ttH_cat->GetBinContent(icat+1);
            h_totalYield_ttH_cat->SetBinContent(icat+1, tmp+weight);
        }
    }

    //Save total yield of events overlapping between tHq/ttH categories
    for(int i=0; i<v_tHq_cat.size(); i++)
    {
        if(v_tHq_cat[i] == 1)
        {
            for(int j=0; j<v_ttH_cat.size(); j++)
            {
                if(v_ttH_cat[j] == 1) //Add event yield to total
                {
                    float tmp = h_overlap_ttH_tHq_cat->GetBinContent(i+1, j+1);
                    h_overlap_ttH_tHq_cat->SetBinContent(i+1, j+1, tmp+weight);
                }
            }
        }
    }

    return;
}



/**
 * Check if a "tHq" cat. is different than its "ttH" counterpart (defined at NTP code level)
 * If there is a difference, it means one code needs debugging (NTA or NTP)
 * /!\ NB : this is not useful (as is) because in NTP code the jet collection always includes forward jet, not this code...!
 */
void tHqMultileponAnalysis::Check_Disagreement_Category_NTP_NTA_Codes()
{
    TString catname = ""; Char_t NTA_cat_value = 0;

    if(is_tHq_2lSS_SR != is_ttH_2lSS_SR) {catname =  "is_tHq_2lSS_SR"; NTA_cat_value = is_tHq_2lSS_SR;}
    if(is_tHq_2lSS_Fake != is_ttH_2lSS_Fake) {catname =  "is_tHq_2lSS_Fake"; NTA_cat_value = is_tHq_2lSS_Fake;}
    if(is_tHq_2lSS_Flip != is_ttH_2lSS_Flip) {catname =  "is_tHq_2lSS_Flip"; NTA_cat_value = is_tHq_2lSS_Flip;}
    if(is_tHq_3l_SR != is_ttH_3l_SR) {catname =  "is_tHq_3l_SR"; NTA_cat_value = is_tHq_3l_SR;}
    if(is_tHq_3l_Fake != is_ttH_3l_Fake) {catname =  "is_tHq_3l_Fake"; NTA_cat_value = is_tHq_3l_Fake;}
    if(is_tHq_ttWctrl_SR != is_ttH_ttWctrl_SR) {catname =  "is_tHq_ttWctrl_SR"; NTA_cat_value = is_tHq_ttWctrl_SR;}
    if(is_tHq_ttWctrl_Fake != is_ttH_ttWctrl_Fake) {catname =  "is_tHq_ttWctrl_Fake"; NTA_cat_value = is_tHq_ttWctrl_Fake;}
    if(is_tHq_ttWctrl_Flip != is_ttH_ttWctrl_Flip) {catname =  "is_tHq_ttWctrl_Flip"; NTA_cat_value = is_tHq_ttWctrl_Flip;}
    if(is_tHq_ttZctrl_SR != is_ttH_ttZctrl_SR) {catname =  "is_tHq_ttZctrl_SR"; NTA_cat_value = is_tHq_ttZctrl_SR;}
    if(is_tHq_ttZctrl_Fake != is_ttH_ttZctrl_Fake) {catname =  "is_tHq_ttZctrl_Fake"; NTA_cat_value = is_tHq_ttZctrl_Fake;}
    if(is_tHq_WZctrl_SR != is_ttH_WZctrl_SR) {catname =  "is_tHq_WZctrl_SR"; NTA_cat_value = is_tHq_WZctrl_SR;}
    if(is_tHq_WZctrl_Fake != is_ttH_WZctrl_Fake) {catname =  "is_tHq_WZctrl_Fake"; NTA_cat_value = is_tHq_WZctrl_Fake;}
    if(is_tHq_4l_SR != is_ttH_4l_SR) {catname =  "is_tHq_4l_SR"; NTA_cat_value = is_tHq_4l_SR;}
    if(is_tHq_ZZctrl_SR != is_ttH_ZZctrl_SR) {catname =  "is_tHq_ZZctrl_SR"; NTA_cat_value = is_tHq_ZZctrl_SR;}

    if(catname != "")
    {
        cout<<"//--------------------------------------------"<<endl;
        cout<<"Event ID : "<<setprecision(15)<<event_id<<endl;
        cout<<"* ==> Cat. "<<catname<<" is defined differently in NTP & NTA code. CHECK !"<<endl;
        cout<<"(the categ. filled by Sync.cxx code has value "<<(bool) NTA_cat_value<<")"<<endl;

        cout<<"nJets = "<<nJets<<endl;
        cout<<"nLightJets = "<<nLightJets<<endl;
        cout<<"nLooseBJets = "<<nLooseBJets<<endl;
        cout<<"nMediumBJets = "<<nMediumBJets<<endl;
        cout<<"nFakeableLep = "<<nFakeableLep<<endl;

        DEBUG = true;
    }
    else {DEBUG = false;}

    return;
}

/**
 * Check if event enters 1 of the new categories -- preliminary studies
 * Logic FOR 3l CHANNEL : to enter a "new categ", event must 1) enter the old categ, and 2) either a) satisfy the additional requirement or b) not enter the category of the other analysis
 * 2l channel : more complicated, there is also overlap with ttW CR. Try 2 approaches (deal or do not deal with this overlap)
 */
void tHqMultileponAnalysis::Define_New_Categorization()
{
//--------------------------------------------
//--------------------------------------------
//DEFAULT -- STILL OVERLAP BETWEEN THQ 2LSS SR AND TTW CR
    if(is_ttH_3l_SR) //ttH 3l SR
    {
        if(!is_tHq_3l_SR)
        {
            is_ttH_3l_SR_btag = 1;
            is_ttH_3l_SR_fwd = 1;
        }
        else
        {
            if(nLooseBJets != 1 || nMediumBJets != 1) {is_ttH_3l_SR_btag = 1;}
            if(nLightJets_Fwd == 0) {is_ttH_3l_SR_fwd = 1;}
        }
    }

    if(is_tHq_3l_SR) //tHq 3l SR
    {
        if(!is_ttH_3l_SR)
        {
            is_tHq_3l_SR_btag = 1;
            is_tHq_3l_SR_fwd = 1;
        }
        else
        {
            if(nLooseBJets == 1 && nMediumBJets == 1) {is_tHq_3l_SR_btag = 1;}
            if(nLightJets_Fwd > 0) {is_tHq_3l_SR_fwd = 1;}
        }
    }

    if(is_ttH_2lSS_SR) //ttH 2lSS SR
    {
        if(!is_tHq_2lSS_SR)
        {
            is_ttH_2lSS_SR_btag = 1;
            is_ttH_2lSS_SR_fwd = 1;
        }
        else
        {
            if(nLooseBJets != 1 || nMediumBJets != 1) {is_ttH_2lSS_SR_btag = 1;}
            if(nLightJets_Fwd == 0) {is_ttH_2lSS_SR_fwd = 1;}
        }
    }

    if(is_tHq_2lSS_SR) //tHq 2lSS SR
    {
        if(!is_ttH_2lSS_SR)
        {
            is_tHq_2lSS_SR_btag = 1;
            is_tHq_2lSS_SR_fwd = 1;
        }
        else
        {
            if(nLooseBJets == 1 && nMediumBJets == 1) {is_tHq_2lSS_SR_btag = 1;}
            if(nLightJets_Fwd > 0) {is_tHq_2lSS_SR_fwd = 1;}
        }

        //Additional requirement to make tHq 2lSS SR ortho with ttW CR
        is_tHq_2lSS_SR_fwd2 = is_tHq_2lSS_SR_fwd;
        if(is_ttH_ttWctrl_SR && nLightJets_Fwd == 0) {is_tHq_2lSS_SR_fwd2 = 0;}
    }

    if(is_ttH_ttWctrl_SR)
    {
        if(!is_tHq_2lSS_SR) {is_ttH_ttWctrl_SR_fwd2 = 1;}
        else if(nLightJets_Fwd == 0) {is_ttH_ttWctrl_SR_fwd2 = 1;}
    }


//--------------------------------------------
//--------------------------------------------
//Also remove overlap with ttW CR

//--- 2l

    //Orthogonalize SRs  -> all overlap events go to ttH (most tHq signal is not in overlap anyway)
    //Orthogonalize ttW CR -> use nFwdJet
    if(is_ttH_2lSS_SR)
    {
        if(!is_tHq_2lSS_SR) {is_ttH_2lSS_SR_njet = 1;}
        else if(nJets_tHq >= 4) {is_ttH_2lSS_SR_njet = 1;} //needed for 3l SR, not 2lSS
    }

    if(is_tHq_2lSS_SR)
    {
        if(!is_ttH_2lSS_SR && !is_ttH_ttWctrl_SR) {is_tHq_2lSS_SR_njet = 1;}
        else if(is_ttH_2lSS_SR && nJets_tHq < 4) {is_tHq_2lSS_SR_njet = 1;}
        else if(is_ttH_ttWctrl_SR && nLightJets_Fwd > 0) {is_tHq_2lSS_SR_njet = 1;}
    }

    if(is_ttH_ttWctrl_SR)
    {
        if(!is_tHq_2lSS_SR) {is_ttH_ttWctrl_SR_njet = 1;}
        else if(nLightJets_Fwd == 0) {is_ttH_ttWctrl_SR_njet = 1;}
    }

    //--- 3l
    if(is_ttH_3l_SR)
    {
        is_ttH_3l_SR_njet2 = 1; //All overlap goes to ttH

        if(!is_tHq_3l_SR || nJets_tHq >= 4) {is_ttH_3l_SR_njet = 1;}
        if(!is_tHq_3l_SR || nJets_tHq >= 3) {is_ttH_3l_SR_njet3 = 1;}
    }

    if(is_tHq_3l_SR)
    {
        if(!is_ttH_3l_SR_njet2) {is_tHq_3l_SR_njet2 = 1;} //All overlap goes to ttH

        if(!is_ttH_3l_SR || nJets_tHq < 4) {is_tHq_3l_SR_njet = 1;}
        if(!is_ttH_3l_SR || nJets_tHq < 3) {is_tHq_3l_SR_njet3 = 1;}
    }

    return;
}


































//--------------------------------------------
//  #######  ########   ######   #######  ##       ######## ######## ########
// ##     ## ##     ## ##    ## ##     ## ##       ##          ##    ##
// ##     ## ##     ## ##       ##     ## ##       ##          ##    ##
// ##     ## ########   ######  ##     ## ##       ######      ##    ######
// ##     ## ##     ##       ## ##     ## ##       ##          ##    ##
// ##     ## ##     ## ##    ## ##     ## ##       ##          ##    ##
//  #######  ########   ######   #######  ######## ########    ##    ########
//--------------------------------------------

/*
bool tHqMultileponAnalysis::is_GammaConv_Event(int nlep)
{
    if(_isdata) {return false;}

    if(abs(vLeptonFakeable[0].id) == 11 && vLeptonFakeable[0].hasPhotonMCMatch) {return true;}
    if(abs(vLeptonFakeable[1].id) == 11 && vLeptonFakeable[1].hasPhotonMCMatch) {return true;}

    if(nlep==3)
    {
        if(abs(vLeptonFakeable[2].id) == 11 && vLeptonFakeable[2].hasPhotonMCMatch) {return true;}
    }

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
        double dlep = Get_Distance(vLeptonFakeable[lepRec_idx].conept, vLeptonFakeable[lepRec_idx].eta, vLeptonFakeable[lepRec_idx].phi, vTruth->at(0).gen_pt.at(genLep_idx), vTruth->at(0).gen_eta.at(genLep_idx), vTruth->at(0).gen_phi.at(genLep_idx) );
        double dpho = Get_Distance(vLeptonFakeable[lepRec_idx].conept, vLeptonFakeable[lepRec_idx].eta, vLeptonFakeable[lepRec_idx].phi, vTruth->at(0).gen_pt.at(genPho_idx), vTruth->at(0).gen_eta.at(genPho_idx), vTruth->at(0).gen_phi.at(genPho_idx) );

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
        // else if(vLeptonFakeable[lepRec_idx].conept < 10 && abs(vLeptonFakeable[lepRec_idx].id) == 13 && vLeptonFakeable[lepRec_idx].id != vTruth->at(0).gen_id.at(itruth) ) {return -1;}
        // else if(dR < 0.7) {return itruth;}
        // else if(std::min(vLeptonFakeable[lepRec_idx].conept, vTruth->at(0).gen_pt.at(itruth)) / std::max(vLeptonFakeable[lepRec_idx].conept, vTruth->at(0).gen_pt.at(itruth)) < 0.3) {return -1;}
        // else if(dR < 1.2) {return itruth;}
    }

    return idx_bestMatch;
}
*/























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
 * Compute the event weight associated with the scale uncertainty
 * See : https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#Factorization_and_renormalizatio
 * FIXME : still need to check the recipe : take enveloppe of all scale variations ? Normalization ?
 * NB : this is a try to get the QCDscaleUp/QCDscaleDown systematics : renorm. each variation, and take the largest of the 3 up/down variations
 * NB : for scale studies : also save each variation+sumWeights independently
 */
void tHqMultileponAnalysis::Compute_Weight_ScaleSyst()
{
    if(_isdata) {return;}

    double QCDscaleDown_tmp = 0, QCDscaleUp_tmp = 0;
    double tmp = 0;

//--------------------------------------------
//Scale DOWN variations => Take enveloppe
    //muR = 0.5 / muF = 1
    tmp = _lumi*_xsec*weight_scale_muR0p5/hSumWeights->GetBinContent(4);
    if(fabs(tmp) > fabs(QCDscaleDown_tmp)) {QCDscaleDown_tmp = tmp;}

    //muR = 1 / muF = 0.5
    tmp = _lumi*_xsec*weight_scale_muF0p5/hSumWeights->GetBinContent(5);
    if(fabs(tmp) > fabs(QCDscaleDown_tmp)) {QCDscaleDown_tmp = tmp;}

    //muR = 0.5 / muF = 0.5
    tmp = _lumi*_xsec*weight_scale_muR0p5muF0p5/hSumWeights->GetBinContent(6);
    if(fabs(tmp) > fabs(QCDscaleDown_tmp)) {QCDscaleDown_tmp = tmp;}

//--------------------------------------------
//Scale UP variations => Take enveloppe
    //muR = 2 / muF = 1
    tmp = _lumi*_xsec*weight_scale_muR2/hSumWeights->GetBinContent(7);
    if(fabs(tmp) > fabs(QCDscaleUp_tmp)) {QCDscaleUp_tmp = tmp;}

    //muR = 1 / muF = 2
    tmp = _lumi*_xsec*weight_scale_muF2/hSumWeights->GetBinContent(8);
    if(fabs(tmp) > fabs(QCDscaleUp_tmp)) {QCDscaleUp_tmp = tmp;}

    //muR = 2 / muF = 2
    tmp = _lumi*_xsec*weight_scale_muR2muF2/hSumWeights->GetBinContent(9);
    if(fabs(tmp) > fabs(QCDscaleUp_tmp)) {QCDscaleUp_tmp = tmp;}

    //Need protection for now, some weird values
    if(fabs(QCDscaleUp_tmp) > 2) {QCDscaleUp_tmp = 1;}
    if(fabs(QCDscaleDown_tmp) > 2) {QCDscaleDown_tmp = 1;}

    // -> weight * total_SF * QCDscale
    QCDscaleDown*= QCDscaleDown_tmp;
    QCDscaleUp*= QCDscaleUp_tmp;

    return;
}


/**
 * Check if WZ event contains at least 1bjet, or at least 1 cjet, else only light jets => split sample
 */
void tHqMultileponAnalysis::Get_HadronFlavour_WZsample()
{
    if(_isdata) {return;}

    for(int ijet=0; ijet<vJetLoose.size(); ijet++)
    {
        if(fabs(vJetLoose.at(ijet).jet_hadronFlavour) == 5) {wz_jetFlav_b = 1;}
        if(fabs(vJetLoose.at(ijet).jet_hadronFlavour) == 4) {wz_jetFlav_c = 1;}
    }

    //Orthogonal Char_t
    if(wz_jetFlav_b == 1) //first check bjets
    {
        wz_jetFlav_c = 0;
        wz_jetFlav_l = 0;
    }
    else if(wz_jetFlav_c == 1) //then check cjets
    {
        wz_jetFlav_l = 0;
    }
    else {wz_jetFlav_l = 1;} //else, only light jets

    return;
}
