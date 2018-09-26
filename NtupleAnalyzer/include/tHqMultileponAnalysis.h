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
#include <fstream>

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

// #include "BTagCalibrationXStandalone.h"

class tHqMultileponAnalysis
{

    public:

        // tHqMultileponAnalysis(); //Default constructor
        tHqMultileponAnalysis(TString inputFileName, TString sampleName, TString treeName, TString outputFileName, bool isdata, bool doSystCombine, float xsec, float lumi, int nowe, int nmax); //Constructor
        ~tHqMultileponAnalysis(); //Destructor

        // void Debug_Selection(int evt);

        //-- ttH basic selections
		bool ThreeLeptonSelection_THQ3l(int evt);
		bool TwoLeptonSelection_THQ2lSSl(int evt);

		//-- ttH regions
		bool ThreeLeptonSelection_THQ3l_Regions(int evt);
		bool TwoLeptonSelection_THQ2lSS_Regions(int evt);

		void InitFiles();
		void InitCollections();
		void InitPredefinedCategories();
		void InitCustomCategories();
		void InitTLorentzVectors();
		void InitVariables();
		void InitTree();
		void initializeOutputTree();

		void Loop();

		void Compute_Variables(TString);
		void fillOutputTree();
		void FillJetInfoOutputTree(int*, int, TLorentzVector*, TLorentzVector, float*, float, float*, float*, float, float*, float*, float, float, float);

        void SelectBjets(int&, int&, bool);
        void SelectOtherJets(const int, const int, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&);
        float Phi_0_2Pi(float phi);
		float Phi_MPi_Pi(float phi);
		float GetDeltaR(float eta1,float phi1,float eta2,float phi2);
        bool Sample_isUsed_forTraining();
		bool Sample_isUsed_forGammaConv();
		double comp_MT_met_lep(TLorentzVector, double, double);
		double Compute_METLD(double);
		int Read_Predefined_Categories();
		void Apply_FakeRate_Weight(TString);
		void Apply_FlipRate_Weight();
		double Get_Distance(double, double, double, double, double, double);
		bool is_GammaConv_Event();
		bool is_Electron_Matched_To_Gamma(int);
		int Ele_To_Lepton_Matching(int);
		int Ele_To_Gamma_Matching(int);

		bool TwoLeptonSelection_THQ2lSS_TrainingSelection(int evt);
		bool ThreeLeptonSelection_THQ3l_TrainingSelection(int evt);

		void Apply_ScaleFactors(int);

		//NEW
		// TString Convert_Number_To_TString(double, int/*=10*/)
		// double Convert_TString_To_Number(TString)
		// bool Check_File_Existence(const TString&)

		ScaleFactors* sf;

        TChain *fChain;   //!pointer to the analyzed TTree or TChain

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
        std::vector<Jet>      *vJetLoose;
		std::vector<Truth>    *vTruth;
		// std::vector<Truth>    *vGenJet;

		std::vector<Lepton>   vLeptonLoose;
		std::vector<Lepton>   vLeptonFakeable;
		std::vector<Lepton>   vLeptonTight;

        //New categories
        std::vector<Jet>	  vLooseBTagJets; //Loose-CSV jets (+pT/eta cuts)
		std::vector<Jet>      vLightJets; //Non-loose CSV jets
		std::vector<Jet>      vLightJets_FwdPtCut; //Non-loose CSV jets (with pT>40 cut on forward jets, SR)

		//Nof objects
        float nLooseBJets;
        float nMediumBJets;
		float nLightJets;
		float nLightJets_Fwd40;
		float nTightLep;
		float nFakeableLep;

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
		Char_t is_tHq_ttWctrl;
		Char_t is_tHq_ttWctrl_SR;
		Char_t is_tHq_ttWctrl_Fake;
        Char_t is_tHq_ttWctrl_Flip;
        Char_t is_tHq_ttWctrl_GammaConv;
		Char_t is_tHq_ttZctrl;
		Char_t is_tHq_ttZctrl_SR;
		Char_t is_tHq_ttZctrl_Fake;
		Char_t is_tHq_WZctrl;
		Char_t is_tHq_WZctrl_SR;
		Char_t is_tHq_WZctrl_Fake;

		//-- Predefined ttH 2017 categories (implemented at NTP level)
        Char_t is_ttH_2lSS;
        Char_t is_ttH_2lSS_SR;
        Char_t is_ttH_2lSS_SR_Data;
		Char_t is_ttH_2lSS_Fake;
        Char_t is_ttH_2lSS_Flip;
        Char_t is_ttH_2lSS_Flip_Data;
        Char_t is_ttH_3l;
        Char_t is_ttH_3l_SR;
        Char_t is_ttH_3l_SR_Data;
		Char_t is_ttH_3l_Fake;
		Char_t is_ttH_ttWctrl;
		Char_t is_ttH_ttWctrl_SR;
        Char_t is_ttH_ttWctrl_SR_Data;
		Char_t is_ttH_ttWctrl_Fake;
        Char_t is_ttH_ttWctrl_Flip;
        Char_t is_ttH_ttWctrl_Flip_Data;
		Char_t is_ttH_ttZctrl;
		Char_t is_ttH_ttZctrl_SR;
        Char_t is_ttH_ttZctrl_SR_Data;
		Char_t is_ttH_ttZctrl_Fake;
		Char_t is_ttH_WZctrl;
		Char_t is_ttH_WZctrl_SR;
        Char_t is_ttH_WZctrl_SR_Data;
		Char_t is_ttH_WZctrl_Fake;


		Char_t is_trigger;
        Char_t is_trigger_ttH;

        Char_t is_hasJetTransitionRegion;

		//3l : 0 uuu, 1 uue, 2 eeu, 3 eee
		//2l : 0 uu, 1 eu+ue, 2 ee
		float channel;

        TTree* tOutput;
        Float_t event_id;
		Float_t event_run;


		//NB -- input vars (tHq2016, ttH2017) are declared in SignalExtractionMVA.h !
        //Additionnal variables, for control plots
        float lep1Pt, lep2Pt, lep3Pt, inv_mll, hardestBjetPt, hardestBjetEta, fwdJetPt;
		float MET, metLD;




        // BTagCalibrationX       calib;
        // BTagCalibrationXReader reader;

        // all weights
        Float_t weight;
        Float_t weightfake;
        Float_t weightflip;
        Float_t mc_weight;                                                                           // weight MC
        Float_t weight_PV;                                                                           // PU reweighting from PV distribution
        Float_t weight_trigger;
        Float_t weight_QF_ee;                                                                        // weight charge flip (only to 2lss where l=e)
        Float_t weight_FR_2lss, weight_FR_3l;                                                        // weight fake rate (to 2lss and 3l)
        // Float_t weight_scale_muF0p5, weight_scale_muF2, weight_scale_muR0p5, weight_scale_muR2;
        std::vector<Float_t> weights_pdf;
        std::vector<std::string> ids_pdf;

        Float_t weight_csv_up, weight_csv_down;



		//--- Lots of variables, used for MEM computation
        Int_t mc_ttZhypAllowed;
        Int_t catJets;

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
        Float_t multilepton_mHT;
        Double_t multilepton_mETcov00;
        Double_t multilepton_mETcov01;
        Double_t multilepton_mETcov10;
        Double_t multilepton_mETcov11;


    private:
        TFile * outputfile;

        TFile *f_CSVwgt_HF, *f_CSVwgt_LF, *f_BTag_eff, *f_QFwgt, *f_FRwgt;

        TString _sampleName;
        TString _process; // Needed for Combine
        TString _outputFileName;

        bool _isdata;
        float _xsec;
        float _lumi;
        int   _nowe; // number of weighted events
        int   _nmax; // max number of events to process

        TFile * _file_PVreweighting;
        TH1F* _h_PV;

        std::ofstream fout_MC;

        std::ofstream fout_RECO;

        std::string del;
        std::string trig;
};

#endif
