#ifndef tHqMultileponAnalysis_H
#define tHqMultileponAnalysis_H

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

#include <TString.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TChain.h>
#include <TFile.h>
#include <TObject.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

//NTProducer
#include "Event.h"
#include "Muon.h"
#include "Tau.h"
#include "Electron.h"
#include "Jet.h"
#include "Truth.h"

//NTAnalyzer
#include "Lepton.h"
#include "ScaleFactors.h" //Class for Scale factors

//ttH Kinematic Fitter (needed to compute resHTT tagger)
#include "../test/HTT_kinfit/HadTopKinFit.h" // HadTopKinFit
#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

class tHqMultileponAnalysis
{

    public:

        // tHqMultileponAnalysis(); //Default constructor
        tHqMultileponAnalysis(TString inputFileName, TString sampleName, TString treeName, TString outputFileName, bool isdata, bool doSystCombine, float xsec, float lumi, int nowe, int nmax); //Constructor
        ~tHqMultileponAnalysis(); //Destructor

        //-- ttH basic selections
		bool ThreeLeptonSelection_THQ3l(int evt);
        bool TwoLeptonSelection_THQ2lSSl(int evt);
        bool FourLeptonSelection_THQ4l(int evt);

		//-- ttH regions
        bool ThreeLeptonSelection_THQ3l_Regions(int evt);
		bool TwoLeptonSelection_THQ2lSS_Regions(int evt);
        bool FourLeptonSelection_THQ4l_Regions(int evt);
		bool TwoLeptonSelection_THQ2lSS_TrainingSelection(int evt);
		bool ThreeLeptonSelection_THQ3l_TrainingSelection(int evt);

		void InitFiles();
		void InitCollections();
		void InitPredefinedCategories();
		void InitCustomCategories();
		void InitTLorentzVectors();
		void InitVariables();
		void InitTree();
		void initializeOutputTree(int);

		void Loop();

		void Compute_Variables(TString);
		void fillOutputTree(int);
		void FillJetInfoOutputTree(int*, int, TLorentzVector*, TLorentzVector, float*, float, float*, float*, float, float*, float*, float, float, float);

        void SelectBjets(int&, int&, bool);
        void SelectOtherJets(const int, const int, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&);
        float Phi_0_2Pi(float phi);
		float Phi_MPi_Pi(float phi);
		float GetDeltaR(float eta1,float phi1,float eta2,float phi2);
        bool Sample_isUsed_forTraining();
		bool Sample_isUsed_forGammaConv();
		double comp_MT_met_lep(TLorentzVector, double, double);
        double Compute_mHT();
        double Compute_metLD_Alternative(double);
		void Compute_Top_W_variables();
		bool Get_Polynom2_Root(double, double, double, vector<double>&);
		void Read_Predefined_Categories();
		void Apply_FakeRate_Weight();
		void Apply_FlipRate_Weight();
		double Get_Distance(double, double, double, double, double, double);
		bool is_GammaConv_Event(int);
		// bool is_Electron_Matched_To_Gamma(int);
		// int Ele_To_Lepton_Matching(int);
		// int Ele_To_Gamma_Matching(int);
		void Dump_EventInfo_Synchro(std::ofstream&);
		int Check_If_Save_Event(bool);
		void Apply_ScaleFactors(int, int);
		void Compute_Weight_ScaleSyst();
		void Fill_Overlap_Histo();
        void Check_Disagreement_Category_NTP_NTA_Codes();
        void Define_New_Categorization();
        void Get_HadronFlavour_WZsample();
		void Get_SumWeights();
		void Get_Mapping_LHE_PDF();
		void Read_LHE_ScaleVariations();
		void Get_PDFset_Weights();

		//TESTING
		void Modify_DefaultCategories_Orthogonal(TString);

		ScaleFactors* sf;

        TChain *fChain;   //!pointer to the analyzed TTree or TChain

		TH1F* h_PU_noCorr;
        TH1F* h_PU_withCorr;
		TH1F* h_nPV_noCorr;
        TH1F* h_nPV_withCorr;
		Int_t nPU;
		Int_t nPV;

        std::vector<Event>    *vEvent;
		std::vector<Electron> *vElectronLoose;
		std::vector<Electron> *vElectronFakeable;
		std::vector<Electron> *vElectronTight;
		std::vector<Muon>     *vMuonLoose;
		std::vector<Muon>     *vMuonFakeable;
		std::vector<Muon>     *vMuonTight;
		std::vector<Tau>      *vTauLoose;
		std::vector<Tau>      *vTauFakeable;
		std::vector<Tau>      *vTauMedium;
		std::vector<Jet>      *vJetLoose_original; //This is the 'original' vector of loose jets, directly taken from NtupleProducer output TTree
		std::vector<Jet>      *vJetSoft; //new vector, containing only soft jets (pt b/w 15 and 25)
		std::vector<Truth>    *vTruth;
		// std::vector<Truth>    *vGenJet;

		std::vector<Lepton>   vLeptonLoose;
		std::vector<Lepton>   vLeptonFakeable;
		std::vector<Lepton>   vLeptonTight;

        //New categories
		std::vector<Jet>      *vJetLoose_tmp; //Exact copy of 'vJetLoose_original', but we use this temporary collec to be able to modify pT by JES or JER variations

		//Then, fill several different of jets from there, that later serve the purposes of the analysis
		std::vector<Jet>	  vJetLoose; //Jets to be used
		std::vector<Jet>	  vJetLoose_tHq; //Jets included in tHq analysis
		std::vector<Jet>	  vJetLoose_ttH; //Jets included in ttH analysis
		std::vector<Jet>      vLightJets; //Non-loose CSV jets to be used
		std::vector<Jet>      vLightJets_tHq; //Non-loose CSV jets included in tHq analysis
		std::vector<Jet>      vLightJets_ttH; //Non-loose CSV jets included in ttH analysis
		// std::vector<Jet>      vLightJets_FwdPtCut; //Non-loose CSV jets (with pT>40 cut on forward jets, SR)
		std::vector<Jet>	  vLooseBTagJets; //Loose-CSV jets (+pT/eta cuts)

		//Nof objects
        float nLooseBJets;
        float nMediumBJets;
		float nLightJets, nLightJets_tHq, nLightJets_ttH, nLightJets_Fwd, nJetLoose;
		float nTightLep;
		float nFakeableLep;
		float nLightJets_eta_Gt2;

		//NB -- input vars (tHq2016, ttH2017) are declared in SignalExtraction.h !

		//Additional vars
		bool pass_Zveto;
		bool pass_Zveto_ee;
		bool pass_cleanup;
		int nSFOS;
		double mllll;

        //NB : all categories (used for cuts) should start with "is_" <-> automated in analysis code
		//tHq2017 categories : based on ttH2017, defined at NTA level
		Char_t is_tHq_2lSS;
		Char_t is_tHq_2lSS_SR;
		Char_t is_tHq_2lSS_Training;
		Char_t is_tHq_2lSS_Fake;
        Char_t is_tHq_2lSS_Flip;
        Char_t is_tHq_2lSS_GammaConv;
        Char_t is_tHq_3l;
		Char_t is_tHq_3l_SR;
		Char_t is_tHq_3l_Training;
		Char_t is_tHq_3l_Fake;
		Char_t is_tHq_3l_GammaConv;
		Char_t is_tHq_ttWctrl;
		Char_t is_tHq_ttWctrl_SR;
		Char_t is_tHq_ttWctrl_Fake;
        Char_t is_tHq_ttWctrl_Flip;
        Char_t is_tHq_ttWctrl_GammaConv;
		Char_t is_tHq_ttZctrl;
		Char_t is_tHq_ttZctrl_SR;
		Char_t is_tHq_ttZctrl_Fake;
		Char_t is_tHq_ttZctrl_GammaConv;
		Char_t is_tHq_WZctrl;
		Char_t is_tHq_WZctrl_SR;
		Char_t is_tHq_WZctrl_Fake;
		Char_t is_tHq_WZctrl_GammaConv;
		Char_t is_tHq_4l_SR;
        Char_t is_tHq_ZZctrl_SR;
		Char_t is_tHqFCNC_2lSS_SR;
		Char_t is_tHqFCNC_2lSS_Fake;
		Char_t is_tHqFCNC_2lSS_Flip;
		Char_t is_tHqFCNC_2lSS_GammaConv;
		Char_t is_tHqFCNC_3l_SR;
		Char_t is_tHqFCNC_3l_Fake;
		Char_t is_tHqFCNC_3l_GammaConv;

		//-- ttH 2017 categories
		Char_t is_ttH_2lSS;
		Char_t is_ttH_2lSS_Training;
        Char_t is_ttH_2lSS_SR;
		Char_t is_ttH_2lSS_Fake;
        Char_t is_ttH_2lSS_Flip;
		Char_t is_ttH_2lSS_GammaConv;
		Char_t is_ttH_3l;
		Char_t is_ttH_3l_Training;
		Char_t is_ttH_3l_SR;
		Char_t is_ttH_3l_GammaConv;
		Char_t is_ttH_3l_Fake;
		Char_t is_ttH_ttWctrl;
		Char_t is_ttH_ttWctrl_SR;
		Char_t is_ttH_ttWctrl_Fake;
		Char_t is_ttH_ttWctrl_Flip;
		Char_t is_ttH_ttWctrl_GammaConv;
		Char_t is_ttH_ttZctrl;
		Char_t is_ttH_ttZctrl_SR;
		Char_t is_ttH_ttZctrl_Fake;
		Char_t is_ttH_ttZctrl_GammaConv;
		Char_t is_ttH_WZctrl;
		Char_t is_ttH_WZctrl_SR;
		Char_t is_ttH_WZctrl_Fake;
		Char_t is_ttH_WZctrl_GammaConv;
		Char_t is_ttH_4l_SR;
		Char_t is_ttH_ZZctrl_SR;

        //-- Test orthogonal regions

        //Using nof btags (NB : ttW CR and tHq_2lSS_SR not orthogonalized)
        Char_t is_ttH_2lSS_SR_btag;
        Char_t is_ttH_3l_SR_btag;
        Char_t is_tHq_2lSS_SR_btag;
        Char_t is_tHq_3l_SR_btag;

        //Using presence of fwd jet
        Char_t is_ttH_2lSS_SR_fwd;
        Char_t is_ttH_3l_SR_fwd;
        Char_t is_tHq_2lSS_SR_fwd;
        Char_t is_tHq_3l_SR_fwd;
        //Modified to account for overlap with ttW CR (made orthog.)
        Char_t is_tHq_2lSS_SR_fwd2;
        Char_t is_ttH_ttWctrl_SR_fwd2;

        //Using total nof jets (including fwd)
        //If nJets < 4, goes to tHq SR
        //Then, if overlap with ttW CR : check presence of fwd jet
        //NB : in 2lSS <-> all overlap events go to ttH SR
        Char_t is_ttH_2lSS_SR_njet4;
        Char_t is_ttH_3l_SR_njet4;
        Char_t is_tHq_2lSS_SR_njet4;
        Char_t is_tHq_3l_SR_njet4;
		Char_t is_ttH_ttWctrl_SR_njet4;

        //Modified nJets<3 criterion
        Char_t is_tHq_3l_SR_njet3;
        Char_t is_ttH_3l_SR_njet3;

        //For 3l, also put all overlap in ttH SR
        Char_t is_tHq_3l_SR_ttHfirst;
		Char_t is_ttH_3l_SR_ttHfirst;



        Char_t is_trigger_1lep;
        Char_t is_trigger_2lep;
        Char_t is_trigger_3lep;
        Char_t is_trigger;
		Char_t is_trigger_2lss;
		Char_t is_trigger_3l;

        //Xcheck HCAL noisy region
		Char_t is_hasJetNoisyHCAL;
		Char_t is_hasManyJetNoisyHCAL;
		Float_t JetNoisyHCALPt;
		Float_t JetNoisyHCALEta;

		//3l : 0 uuu, 1 uue, 2 eeu, 3 eee
		//2l : 0 uu, 1 eu+ue, 2 ee
		float channel;

		// TTree* tOutput;
		vector<TTree*> v_tOutput;
        Int_t event_id;
		Float_t event_run;
		Int_t event_lumi;




 // # #    # #####  #    # #####    #    #   ##   #####  #   ##   #####  #      ######  ####
 // # ##   # #    # #    #   #      #    #  #  #  #    # #  #  #  #    # #      #      #
 // # # #  # #    # #    #   #      #    # #    # #    # # #    # #####  #      #####   ####
 // # #  # # #####  #    #   #      #    # ###### #####  # ###### #    # #      #           #
 // # #   ## #      #    #   #       #  #  #    # #   #  # #    # #    # #      #      #    #
 // # #    # #       ####    #        ##   #    # #    # # #    # #####  ###### ######  ####

		//NB -- all variables that may be needed to compute BDT output at NTA level (i.e. variables used in ttH2017 or tHq2016 BDTs) are defined in SignalExtraction.h

        //Input variables (some for testing only) === See definitions in function 'Compute_Variables()'
        //NB1 : "bjet" <-> loose WP, unless specified
        //NB2 : "fwd jet" <-> most forward light jet, unless specified
        //NB3 : jets are taken from 'vJetLoose'/'vLooseBTagJets' vectors, and leptons from 'vLeptonFakeable' (see .cxx code)
        //--------------------------------------------
        Float_t inv_mll; //Invariant mass of SFOS lepton pair, closest to mZ
		Float_t metpt, metphi; //MET (pt & phi)
        Float_t mHT; // mHT = sqrt(Sum over squared pT of all jets and leptons)
        Float_t metLD; //metLD = 0.6 * metpt + 0.4 * mHT (used by ttH, more robust against pileup)
		Float_t hardestBjetPt; //pT of hardest bjet
		Float_t hardestBjetEta; //eta of hardest bjet
		Float_t minv_FwdJetBJet; //Invariant mass of most forward light jet and hardest bjet
        // Float_t FwdJetPt; //pt of most forward light jet //declared in SignalExtraction.h, used in BDT at NTA level
        Float_t FwdJetEta; //eta of most forward light jet
        Float_t LeadJetPt; //pt of leading jet (any jet)
        Float_t LeadJetEta; //eta of leading jet (any jet)
        Float_t dRjj_max; //Max dR among any 2 jets
        Float_t deepCSV_max; //Max deepCSV score among jets
        Float_t deepCSV_2nd; //Second max deepCSV score among jets
        Float_t Mjj_max; //Max. invariant mass of any 2 jets
        Float_t dPhiLepBJet_max; //Max. dPhi between hardest btag and any lepton
        Float_t dPhijj_max; //Max dPhi between any 2 jets
        Float_t m3l; //Invariant mass of all fakeable leptons (2 or 3)
        Float_t dPhiLepLep_max; //Max. dPhi between any 2 leptons
        Float_t sum_jetPt; //Sum of the pt of all jets

        //These variables are filled in the testing function 'Compute_Top_W_variables()'
        // == !! May be bugged !! ==
        //Try to infer neutrino pT from MET and lepton infos
        //Lepton associated with W decay => reconstructed system has mass closest to 173 GeV
        //Recoiling jet defined has jet with highest pT, not associated to reconstructed top decay
        Float_t top_mass; //Reconstructed top mass
        Float_t mTW; //Reconstructed mTW
		Float_t lW_asym_mtop; //Charge asymmetry from the lepton assigned to top (<-> reconstructed system has mass closest to 173 GeV)
		Float_t dRBjetRecoilJet; //dR of hardest bjet and recoil jet
		Float_t dRLepWRecoilJet; //dR of lepton from W and recoil jet
		Float_t RecoilJetPt; //pt of recoil jet
		Float_t RecoilJetEta; //eta of recoil jet
		Float_t LepWPt; //pt of lepton from W
		Float_t LepWEta; //eta of lepton from W
        Float_t top_Pt; //pt of reconstructed top
        Float_t mass_LepBJet_min; //Min inv mass between any lepton and any *medium* bjet

        //HjTagger (testing)
		Float_t HjTagger; //Max value of the HjTag classifier among all jets

        //Low-level variables (for DNN tests, ...)
        // Float_t lep1_conePt, lep2_conePt, lep3_conePt; //Already declared in SignalExtraction.h, since may be used to book/evaluate MVA !
        Float_t lep1Eta, lep2Eta, lep3Eta, lep1Phi, lep2Phi, lep3Phi;

        //Additional variables for FCNC analysis (consider central jets only !)
		Float_t nSoftJets; //Nof 'ttH loose jets' with pT between 15 and 25
		Float_t min_dr_lep_bjet; //Min. dR between any lepton and any bjet
		Float_t min_dr_lep_lightjet; //Min. dR between any lepton and any light jet
		Float_t lW_asym; //Charge asymmetry of the lepton (not part of OS pair) having min dR with hardest bjet
		Float_t ratio_lep3pt_closestJetPt; //Ratio of pt of softest lepton and pt of closest jet (any jet)
		Float_t dPhiLepLep_hardestOS; //dPhi between the 2 leptons in 2lSS, and b/w the 2 leptons forming hardest OS pair in 3l
		Float_t jet1_pt, jet2_pt, jet3_pt, jet4_pt; //pt of the 4 leading jets -- INCLUDE SOFT JETS WITH PT > 15

/*
* le nombre de jets "soft", ayant un pT en entre 15 et 25 GeV,
* le deltaR de la pair de lepton-bjet les plus proches,
* le deltaR de la pair de lepton-ljet les plus proches,
* la somme de la charge de leptons selectionnés,
* dans le cas 3l, asymmetry de charge q*|eta_L| pour le candidat lepton venant du top  (celui ayant le deltaR avec le jet b le plus petit et qui ne donne pas de paire de leptons de même signe),
* le pT du lepton le plus soft,
* le pT du lepton le plus soft divisé par le pT du jet le plus proche,
* dans le 2LSS; Delta Phi entre les leptons, dans 3L, Delta Phi entre les deux leptons de chargé opposé avec le pT le plus élevé
 */


        //Event weights
        Float_t mc_weight; //set to +/-1 (old)
        Float_t weight_old; //=> weight_old = +-1 * xsec * lumi / swe_old

        Float_t mc_weight_originalValue; //'true' mc weight, can different from +/-1 (new)
        Float_t weight; //=> weight = real_gen_weight * xsec * lumi / swe_real

        //FakeRate weights
        Float_t weightfake; //store nominal FR weight separately
		vector<TString> v_FR_type; //name of FR variations
		vector<Float_t> v_floats_FR_variations; //store all FR variations
        // Float_t FR_norm_elUp;
        // Float_t FR_norm_elDown;
        // Float_t FR_norm_muUp;
        // Float_t FR_norm_muDown;
        // Float_t FR_pt_elUp;
        // Float_t FR_pt_elDown;
        // Float_t FR_pt_muUp;
        // Float_t FR_pt_muDown;
        // Float_t FR_be_elUp;
        // Float_t FR_be_elDown;
        // Float_t FR_be_muUp;
        // Float_t FR_be_muDown;


        //QFlip weight
        Float_t weightflip;

		//Scale factors
        Float_t lepton_SF;
		Float_t trigger_SF;
		Float_t btag_SF;
		Float_t PU_SF;
		Float_t total_SF;


        std::vector<Float_t> weights_pdf;
        std::vector<std::string> ids_pdf;

		//Systematic weights
		Float_t TrigEffUp, TrigEffDown;
		// Float_t MuEffUp, MuEffDown; //split lepton eff
		// Float_t EleEffUp, EleEffDown; //split lepton eff
        Float_t LepEff_elLooseUp, LepEff_elLooseDown;
		Float_t LepEff_muLooseUp, LepEff_muLooseDown;
        Float_t LepEff_elTightUp, LepEff_elTightDown;
		// Float_t LepEff_muTightUp, LepEff_muTightDown; //Use flat syst for muTight uncert.

		Float_t LFcontUp, LFcontDown;
		Float_t HFstats1Up, HFstats1Down;
		Float_t HFstats2Up, HFstats2Down;
		Float_t CFerr1Up, CFerr1Down;
		Float_t CFerr2Up, CFerr2Down;
		Float_t HFcontUp, HFcontDown;
		Float_t LFstats1Up, LFstats1Down;
		Float_t LFstats2Up, LFstats2Down;
		Float_t PUUp, PUDown;
        Float_t scaleShapeAccUp, scaleShapeAccDown;
		Float_t pdfUp, pdfDown;
		Float_t fwdJetUp, fwdJetDown;



		//--- Lots of variables, used for MEM computation
        Int_t mc_ttZhypAllowed;
        Int_t catJets;

        //NEW -- add collections of jets without selection -- store all infos in vectors of float
        Int_t nJets, nJets_tHq, nJets_ttH;
		vector<Float_t>* JetsPt;
        vector<Float_t>* JetsEta;
        vector<Float_t>* JetsPhi;
        vector<Float_t>* JetsE;
        vector<Float_t>* JetsCSV;

        Int_t           multilepton_Lepton1_Id,             multilepton_Lepton2_Id,                 multilepton_Lepton3_Id,             multilepton_Lepton4_Id;
        TLorentzVector  multilepton_Lepton1_P4,             multilepton_Lepton2_P4,                 multilepton_Lepton3_P4,             multilepton_Lepton4_P4;

        Float_t         multilepton_Lepton1_DeltaR_Matched,  multilepton_Lepton2_DeltaR_Matched,    multilepton_Lepton3_DeltaR_Matched,  multilepton_Lepton4_DeltaR_Matched;
        Int_t           multilepton_Lepton1_Label_Matched,   multilepton_Lepton2_Label_Matched,     multilepton_Lepton3_Label_Matched,   multilepton_Lepton4_Label_Matched;
        Int_t           multilepton_Lepton1_Id_Matched,      multilepton_Lepton2_Id_Matched,        multilepton_Lepton3_Id_Matched,      multilepton_Lepton4_Id_Matched;
        TLorentzVector  multilepton_Lepton1_P4_Matched,      multilepton_Lepton2_P4_Matched,        multilepton_Lepton3_P4_Matched,      multilepton_Lepton4_P4_Matched;

        Int_t multilepton_Bjet1_Id, multilepton_Bjet2_Id;
        TLorentzVector multilepton_Bjet1_P4, multilepton_Bjet2_P4;
        Float_t multilepton_Bjet1_CSV, multilepton_Bjet2_CSV;
        Float_t multilepton_Bjet1_JEC_Up, multilepton_Bjet1_JEC_Down, multilepton_Bjet2_JEC_Up, multilepton_Bjet2_JEC_Down;
        Float_t multilepton_Bjet1_JER_Up, multilepton_Bjet1_JER_Down, multilepton_Bjet2_JER_Up, multilepton_Bjet2_JER_Down;

        Float_t         multilepton_Bjet1_DeltaR_Matched,   multilepton_Bjet2_DeltaR_Matched;
        Int_t           multilepton_Bjet1_Label_Matched,    multilepton_Bjet2_Label_Matched;
        Int_t           multilepton_Bjet1_Id_Matched,       multilepton_Bjet2_Id_Matched;
        TLorentzVector  multilepton_Bjet1_P4_Matched,       multilepton_Bjet2_P4_Matched;

        Int_t multilepton_JetHighestPt1_Id, multilepton_JetHighestPt2_Id, multilepton_JetClosestMw1_Id, multilepton_JetClosestMw2_Id, multilepton_JetLowestMjj1_Id, multilepton_JetLowestMjj2_Id;
        TLorentzVector multilepton_JetHighestPt1_P4, multilepton_JetHighestPt2_P4, multilepton_JetClosestMw1_P4, multilepton_JetClosestMw2_P4, multilepton_JetLowestMjj1_P4, multilepton_JetLowestMjj2_P4;
        Float_t multilepton_JetHighestPt1_CSV, multilepton_JetHighestPt2_CSV, multilepton_JetClosestMw1_CSV, multilepton_JetClosestMw2_CSV, multilepton_JetLowestMjj1_CSV, multilepton_JetLowestMjj2_CSV;
        Float_t multilepton_JetHighestPt1_JEC_Up, multilepton_JetHighestPt2_JEC_Up, multilepton_JetClosestMw1_JEC_Up, multilepton_JetClosestMw2_JEC_Up, multilepton_JetLowestMjj1_JEC_Up, multilepton_JetLowestMjj2_JEC_Up;
        Float_t multilepton_JetHighestPt1_JEC_Down, multilepton_JetHighestPt2_JEC_Down, multilepton_JetClosestMw1_JEC_Down, multilepton_JetClosestMw2_JEC_Down, multilepton_JetLowestMjj1_JEC_Down, multilepton_JetLowestMjj2_JEC_Down;
        Float_t multilepton_JetHighestPt1_JER_Up, multilepton_JetHighestPt2_JER_Up, multilepton_JetClosestMw1_JER_Up, multilepton_JetClosestMw2_JER_Up, multilepton_JetLowestMjj1_JER_Up, multilepton_JetLowestMjj2_JER_Up;
        Float_t multilepton_JetHighestPt1_JER_Down, multilepton_JetHighestPt2_JER_Down, multilepton_JetClosestMw1_JER_Down, multilepton_JetClosestMw2_JER_Down, multilepton_JetLowestMjj1_JER_Down, multilepton_JetLowestMjj2_JER_Down;

        Int_t multilepton_JetHighestPt1_2ndPair_Id, multilepton_JetHighestPt2_2ndPair_Id, multilepton_JetClosestMw1_2ndPair_Id, multilepton_JetClosestMw2_2ndPair_Id, multilepton_JetLowestMjj1_2ndPair_Id, multilepton_JetLowestMjj2_2ndPair_Id;
        TLorentzVector  multilepton_JetHighestPt1_2ndPair_P4, multilepton_JetHighestPt2_2ndPair_P4, multilepton_JetClosestMw1_2ndPair_P4, multilepton_JetClosestMw2_2ndPair_P4, multilepton_JetLowestMjj1_2ndPair_P4, multilepton_JetLowestMjj2_2ndPair_P4;
        Float_t multilepton_JetHighestPt1_2ndPair_CSV, multilepton_JetHighestPt2_2ndPair_CSV, multilepton_JetClosestMw1_2ndPair_CSV, multilepton_JetClosestMw2_2ndPair_CSV, multilepton_JetLowestMjj1_2ndPair_CSV, multilepton_JetLowestMjj2_2ndPair_CSV;
        Float_t multilepton_JetHighestPt1_2ndPair_JEC_Up, multilepton_JetHighestPt2_2ndPair_JEC_Up, multilepton_JetClosestMw1_2ndPair_JEC_Up, multilepton_JetClosestMw2_2ndPair_JEC_Up, multilepton_JetLowestMjj1_2ndPair_JEC_Up, multilepton_JetLowestMjj2_2ndPair_JEC_Up;
        Float_t multilepton_JetHighestPt1_2ndPair_JEC_Down, multilepton_JetHighestPt2_2ndPair_JEC_Down, multilepton_JetClosestMw1_2ndPair_JEC_Down, multilepton_JetClosestMw2_2ndPair_JEC_Down, multilepton_JetLowestMjj1_2ndPair_JEC_Down, multilepton_JetLowestMjj2_2ndPair_JEC_Down;
        Float_t multilepton_JetHighestPt1_2ndPair_JER_Up, multilepton_JetHighestPt2_2ndPair_JER_Up, multilepton_JetClosestMw1_2ndPair_JER_Up, multilepton_JetClosestMw2_2ndPair_JER_Up, multilepton_JetLowestMjj1_2ndPair_JER_Up, multilepton_JetLowestMjj2_2ndPair_JER_Up;
        Float_t multilepton_JetHighestPt1_2ndPair_JER_Down, multilepton_JetHighestPt2_2ndPair_JER_Down, multilepton_JetClosestMw1_2ndPair_JER_Down, multilepton_JetClosestMw2_2ndPair_JER_Down, multilepton_JetLowestMjj1_2ndPair_JER_Down, multilepton_JetLowestMjj2_2ndPair_JER_Down;

        Int_t multilepton_JetHighestEta1_Id, multilepton_JetHighestEta2_Id;
        TLorentzVector multilepton_JetHighestEta1_P4, multilepton_JetHighestEta2_P4;
		Float_t multilepton_JetHighestEta1_CSV, multilepton_JetHighestEta2_CSV;
		Float_t multilepton_JetHighestEta1_JEC_Up, multilepton_JetHighestEta1_JEC_Down, multilepton_JetHighestEta2_JEC_Up, multilepton_JetHighestEta2_JEC_Down;
		Float_t multilepton_JetHighestEta1_JER_Up, multilepton_JetHighestEta1_JER_Down, multilepton_JetHighestEta2_JER_Up, multilepton_JetHighestEta2_JER_Down;

        Int_t           multilepton_h0_Label,   multilepton_t1_Label,   multilepton_t2_Label;
        Int_t           multilepton_h0_Id,      multilepton_t1_Id,      multilepton_t2_Id;
        TLorentzVector  multilepton_h0_P4,      multilepton_t1_P4,      multilepton_t2_P4;

        TLorentzVector multilepton_mET, multilepton_Ptot;
        Float_t multilepton_sumET;
        Double_t multilepton_mETcov00;
        Double_t multilepton_mETcov01;
        Double_t multilepton_mETcov10;
        Double_t multilepton_mETcov11;


    private:
        TFile * outputfile;

        // TFile *f_CSVwgt_HF, *f_CSVwgt_LF, *f_BTag_eff, *f_QFwgt, *f_FRwgt;
		TFile *f_QFwgt, *f_FRwgt;

        TString _sampleName;
        TString _process; // Needed for Combine
        TString _outputFileName;

        bool _isdata;
        float _xsec;
        float _lumi;
        int   _nowe; // number of weighted events
        int   _nmax; // max number of events to process

        //SYNCHRO
		std::ofstream outfile_2lSS_SR;
		std::ofstream outfile_3l_SR;
		std::ofstream outfile_ttW_CR;
		std::ofstream outfile_ttZ_CR;
		std::ofstream outfile_WZ_CR;

		TString inputFileName;

        TH2F* h_overlap_ttH_tHq_cat;
        TH1F* h_totalYield_tHq_cat;
		TH1F* h_totalYield_ttH_cat;

		vector<TString> v_systTree; //can choose to store JES/JER TTrees, in addition to default

        //Higgs decay modes
        int higgs_daughter_id;

		//Flavour of jet, for splitting of WZ sample only
		Char_t wz_jetFlav_b;
		Char_t wz_jetFlav_c;
		Char_t wz_jetFlav_l;

		//Scale variation weights, stored as LHE weights
		float weight_originalXWGTUP; //can differ from mc_weight_originalValue
		float weight_scale_muR0p5;
		float weight_scale_muF0p5;
		float weight_scale_muR0p5muF0p5;
		float weight_scale_muR2;
		float weight_scale_muF2;
		float weight_scale_muR2muF2;

		//Sums of weights (for nominal and scales)
		float sumWeights_mc_weight_originalValue;
		float sumWeights_originalXWGTUP;
		float sumWeights_scale_muR0p5;
		float sumWeights_scale_muF0p5;
		float sumWeights_scale_muR0p5muF0p5;
		float sumWeights_scale_muR2;
		float sumWeights_scale_muF2;
		float sumWeights_scale_muR2muF2;

		TH1D* hSumWeights; //Histo containing sums of weights for nominal + all scale variations //NB : for scale variations, now reading from hLHE ; but still use this nominal
		TH1D* hLHE; //Histo containing sums of weights for all LHE weights
		// TH1D* hktkv; //Histo containing sums of weights for kT/kV reweights //not present in FlatTrees yet, called "hLHE" for now

		vector<Float_t> v_sums_LHEweights; //Sum of weights for the LHE kT/kV reweights
		vector<Float_t> v_couplings_SF; //Scale factors for the LHE kT/kV reweights (kinematics/acceptance only ! Not xsec)
		Float_t SMcoupling_SF; //SF for SM coupling (kin + acceptance + xsec !)
		Float_t xsec_SM; //SM xsec (hardcoded, for tHq and tHW)

        //--- LHE weights
		//Full vectors
		std::vector<Float_t> LHEweights; //LHE weights
		std::vector<std::string> LHEweights_Ids; //LHE weights ids

		//--- Mapping LHE <-> PDF
		vector<int> v_LHE_ids;
		vector<int> v_PDF_ids;
		vector<TString> v_scale_ids;

		//Subvector, only to store given PDF sets
        std::vector<Float_t> v_PDF_weights_1;
        std::vector<Float_t> v_PDF_SumWeights_1;
		std::vector<Float_t> v_PDF_weights_2;
		std::vector<Float_t> v_PDF_SumWeights_2;
		std::vector<Float_t> v_PDF_weights_3;
		std::vector<Float_t> v_PDF_SumWeights_3;
		std::vector<Float_t> v_PDF_weights_4;
		std::vector<Float_t> v_PDF_SumWeights_4;
		Int_t minIdx_PDFset1, maxIdx_PDFset1; //Min/max indices (in LHE vector) of PDF set 1
		Int_t minIdx_PDFset2, maxIdx_PDFset2; //Min/max indices (in LHE vector) of PDF set 2
		Int_t minIdx_PDFset3, maxIdx_PDFset3; //Min/max indices (in LHE vector) of PDF set 3
		Int_t minIdx_PDFset4, maxIdx_PDFset4; //Min/max indices (in LHE vector) of PDF set 4

		//NEW -- Kin fit for resHTT
		// class HadTopKinFit; // forward declaration
		HadTopKinFit * kinFit_;
		Float_t resHTT;
};

#endif
