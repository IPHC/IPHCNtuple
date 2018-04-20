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
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>
#include "Event.h"
#include "Muon.h"
#include "Tau.h"
#include "Electron.h"
#include "Jet.h"
#include "Truth.h"
#include <TObject.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>

// #include "HistoManager.h"
#include "Lepton.h"
#include "BTagCalibrationXStandalone.h"

class tHqMultileponAnalysis
{

    public:

        tHqMultileponAnalysis(); //Default constructor
        tHqMultileponAnalysis(TString inputFileName, TString sampleName, TString treeName, TString outputFileName, bool isdata, bool doSystCombine, float xsec, float lumi, int nowe, int nmax); //Constructor
        ~tHqMultileponAnalysis(); //Destructor

        // void Debug_Selection(int evt);

        //-- 3l selections
        void ThreeLeptonSelection_THQ3l_SR(int evt);
        void ThreeLeptonSelection_THQ3l_Training(int evt);
        void ThreeLeptonSelection_THQ3l_Z_CR(int evt);

        //-- 2l selections
        void TwoLeptonSelection_THQ2l_SR(int evt);
        void TwoLeptonSelection_THQ2l_Training(int evt);
        void TwoLeptonSelection_THQ2l_EMU_OS_CR(int evt);

		//-- Custom samples : fakes, gamma-conv, charge flips
		void ThreeLeptonSelection_ApplicationFakes(int evt);
		void TwoLeptonSelection_ApplicationFakes(int evt);
		void TwoLeptonSelection_ApplicationChargeFlip(int evt);

        //-- Gamma conv
        void ThreeLeptonSelection_GammaConv(int evt);
        void TwoLeptonSelection_GammaConv(int evt);

		void InitFiles();
		void InitEvent();
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

		//NEW
		// TString Convert_Number_To_TString(double, int/*=10*/)
		// double Convert_TString_To_Number(TString)
		// bool Check_File_Existence(const TString&)

        TChain *fChain;   //!pointer to the analyzed TTree or TChain

        std::vector<Event>    *vEvent;
        std::vector<Electron> *vElectron;
        std::vector<Muon>     *vMuon;
        std::vector<Tau>      *vTau;
        std::vector<Jet>      *vJet;
        std::vector<Truth>    *vTruth;

        std::vector<Lepton>   vLeptons;
        std::vector<Lepton>   vSelectedLeptons;
        std::vector<Muon>	  vSelectedMuons;
        std::vector<Electron> vSelectedElectrons;
        std::vector<Tau>      vSelectedTaus;

        std::vector<Muon>	  vFakeMuons;     // inverted MVA cut
        std::vector<Electron> vFakeElectrons; // inverted MVA cut
        std::vector<Lepton>   vFakeLeptons;   // inverted MVA cut

        std::vector<Jet>	  vSelectedJets; //All jets with pT>25

        //New categories
        std::vector<Jet>	  vLooseBTagJets; //Loose-CSV jets (+pT/eta cuts)
		std::vector<Jet>      vLightJets; //Non-loose CSV jets
		std::vector<Jet>      vLightJets_FwdPtCut; //Non-loose CSV jets (with pT>40 cut on forward jets, SR)
        std::vector<Lepton>   vFakeableLeptons;
        std::vector<Lepton>   vLooseLeptons;
        std::vector<Lepton>   vTightLeptons;

        float nLooseBJets;
        float nMediumBJets;
		float nLightJets;
		float nLightJets_Fwd40;

		float nTightLep;
		float nFakeableLep;


        //-- Categories are encoded as Char_t (not booleans, bc vector<bool> is 'broken' in C++)
        Char_t is_3l_THQ_SR;    // Category : training events tHQ3l analysis
        Char_t is_3l_THQ_Training; //Category : training events tHQ3l analysis
        Char_t is_3l_Z_CR;
        Char_t is_2l_THQ_SR;
        Char_t is_2l_THQ_Training;
        Char_t is_2l_EMU_CR;
        Char_t is_3l_AppFakes_SR;
        Char_t is_3l_GammaConv_SR;
	    Char_t is_2l_AppFakes_SR;
        Char_t is_2l_GammaConv_SR;
        Char_t is_2l_QFlip_SR;

        bool is_trigger;

		//3l : 0 uuu, 1 uue, 2 eeu, 3 eee
		//2l : 0 uu, 1 eu+ue, 2 ee
		float channel;

        TTree* tOutput;
        Float_t event_id;
		Float_t event_run;


        //-- NB : input variables declared in SignalExtractionMVA.h
        //Additionnal variables, for control plots
        float lep1Pt, lep2Pt, lep3Pt, inv_mll, hardestBjetPt, hardestBjetEta, fwdJetPt;
		float MET;




        //--- Lots of variables, used for MEM computation
        BTagCalibrationX       calib;
        BTagCalibrationXReader reader;

        // all weights
        Float_t weight;
        Float_t weightfake;
        Float_t weightflip;
        Float_t mc_weight;                                                                           // weight MC
        Float_t weight_PV;                                                                           // PU reweighting from PV distribution
        Float_t weight_trigger;
        Float_t weight_QF_ee;                                                                        // weight charge flip (only to 2lss where l=e)
        Float_t weight_FR_2lss, weight_FR_3l;                                                        // weight fake rate (to 2lss and 3l)
        Float_t weight_scale_muF0p5, weight_scale_muF2, weight_scale_muR0p5, weight_scale_muR2;
        std::vector<Float_t> weights_pdf;
        std::vector<std::string> ids_pdf;

        Float_t weight_csv_up, weight_csv_down;

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
